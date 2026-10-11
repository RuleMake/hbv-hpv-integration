"""Check the release and optionally rebuild figures in an isolated directory.

Default checks use only the Python standard library and never import analysis sources.
Rebuilding needs requirements-figures.txt and the font described in README.md.
"""
from __future__ import annotations
import argparse
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def inside(rel):
    p = (ROOT / rel).resolve()
    if not p.is_relative_to(ROOT):
        raise ValueError(f'Path outside the release: {rel}')
    return p

def check_release():
    manifest = json.loads((ROOT / 'release_manifest.json').read_text(encoding='utf-8'))
    files = manifest['files']
    for item in files:
        p = inside(item['path'])
        if not p.is_file() or p.stat().st_size != item['bytes'] or sha(p) != item['sha256']:
            raise ValueError(f'Missing or changed release file: {item["path"]}')
    checked_python = 0
    for item in files:
        if item['path'].endswith('.py'):
            ast.parse(inside(item['path']).read_text(encoding='utf-8-sig'), filename=item['path'])
            checked_python += 1

    table_specs = json.loads((ROOT / 'tables/TABLES.json').read_text(encoding='utf-8'))
    assert len(table_specs) == 7
    for item in table_specs:
        path = inside(item['path'])
        rows = read_csv(path)
        assert sha(path) == item['sha256'] and len(rows) == item['rows'], item['id']
        assert list(rows[0]) == item['columns'], item['id']

    claims = read_csv(ROOT / 'provenance/numeric_claim_source_map.csv')
    cache = {}
    for row in claims:
        p = inside(row['source_path'])
        assert sha(p) == row['source_sha256'], row['claim']
        if p not in cache:
            cache[p] = json.loads(p.read_text(encoding='utf-8')) if p.suffix == '.json' else read_csv(p)
        if p.suffix == '.json':
            value = cache[p][row['field']]
        else:
            line = int(row['source_location'].split()[-1])
            value = cache[p][line - 2][row['field']]
        assert str(value) == row['value'], f'Claim/source mismatch: {row["claim"]}'

    register = read_csv(ROOT / 'provenance/current_result_register.csv')
    for row in register:
        target = read_csv(ROOT / 'tables' / row['table'])[int(row['table_csv_line']) - 2]
        assert target[row['p_field']] == row['p_value'], row['result_id']
        if row['adjusted_p_field']:
            assert target[row['adjusted_p_field']] == row['adjusted_p_value'], row['result_id']
        assert sha(inside(row['source_path'])) == row['source_sha256'], row['result_id']

    figures = json.loads((ROOT / 'figures/FIGURES.json').read_text(encoding='utf-8'))['figures']
    assert {x['id'] for x in figures} == {f'F{i}' for i in range(1,7)} | {f'S{i}' for i in range(1,12)}
    for fig in figures:
        folder = inside(fig['dir'])
        assert sha(folder / fig['reference']) == fig['reference_sha256'], fig['id']
        for row in read_csv(folder / 'INPUTS.csv'):
            assert sha(folder / row['path']) == row['sha256'], f'{fig["id"]}: {row["path"]}'
        for cmd in fig['commands']:
            assert cmd[0] == 'python' and (folder / cmd[1]).is_file(), fig['id']
        for rel in fig['tables']:
            assert (folder / rel).is_file(), f'{fig["id"]}: {rel}'
        generated, retained = set(fig['generated_tables']), set(fig['input_tables'])
        assert not generated & retained and generated | retained == set(fig['tables']), fig['id']
    return {'status':'PASS', 'manifest_files':len(files), 'python_sources_parsed':checked_python,
            'tables':7, 'mapped_claims':len(claims), 'registered_result_rows':len(register),
            'figures':len(figures), 'analysis_sources_imported':0}, figures

# Block writes outside the isolated build. Set MPLCONFIGDIR and disable bytecode so
# imported libraries also keep their caches inside this directory.
GUARD = r'''
import os, pathlib, runpy, sys
base=pathlib.Path(sys.argv[1]).resolve()
script=pathlib.Path(sys.argv[2]).resolve()
def check_write(path):
    if isinstance(path,(str,bytes,os.PathLike)):
        p=pathlib.Path(os.fsdecode(path)).resolve()
        if not p.is_relative_to(base):
            raise PermissionError('Build attempted a write outside its isolated directory: '+str(p))
def audit(event,args):
    if event=='open':
        path,mode,flags=args
        if (isinstance(mode,str) and any(c in mode for c in 'wax+')) or (isinstance(flags,int) and flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)):
            check_write(path)
    elif event in ('os.remove','os.rmdir','os.mkdir','os.chmod','os.utime'):
        check_write(args[0])
    elif event in ('os.rename','os.link','os.symlink'):
        check_write(args[0]); check_write(args[1])
sys.addaudithook(audit)
sys.path.insert(0,str(script.parent))
sys.argv=[str(script)]
runpy.run_path(str(script),run_name='__main__')
'''

def images_equal(left, right):
    import fitz
    a, b = fitz.open(left), fitz.open(right)
    try:
        if len(a) != len(b):
            return False
        for p, q in zip(a,b):
            x = p.get_pixmap(matrix=fitz.Matrix(1.5,1.5), alpha=False)
            y = q.get_pixmap(matrix=fitz.Matrix(1.5,1.5), alpha=False)
            if (x.width,x.height,x.n,x.samples)!=(y.width,y.height,y.n,y.samples):
                return False
        return True
    finally:
        a.close(); b.close()

def rebuild_figure(fig, out):
    import pandas as pd
    from pandas.testing import assert_frame_equal
    started = time.monotonic()
    source = inside(fig['dir'])
    with tempfile.TemporaryDirectory(prefix=fig['id']+'_', dir=out) as temp:
        work = Path(temp) / 'figure'
        shutil.copytree(source, work, ignore=shutil.ignore_patterns('__pycache__','rendered','review'))
        # A matching copied file is not evidence of regeneration. Remove every
        # declared output table before running its preparation/build command.
        for rel in fig['generated_tables']:
            (work / rel).unlink()
        for rel in [fig['output'], fig.get('pdf_output')]:
            if rel:
                (work / rel).unlink(missing_ok=True)
        if fig['id'] == 'S11':
            (work / 'saved_tests.json').unlink(missing_ok=True)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLCONFIGDIR=str(Path(temp)/'matplotlib_cache'))
        for command in fig['commands']:
            run = subprocess.run([sys.executable, '-B', '-c', GUARD, temp, str(work/command[1])],
                                 cwd=work, env=env, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=600)
            if run.returncode:
                raise RuntimeError(f'{fig["id"]} {command[1]} failed:\n{run.stdout}\n{run.stderr}')
        for rel in fig['tables']:
            expected, obtained = source/rel, work/rel
            if sha(expected) != sha(obtained):
                assert_frame_equal(pd.read_csv(expected), pd.read_csv(obtained), check_exact=True, check_dtype=False)
        for rel in fig['input_tables']:
            assert sha(source/rel) == sha(work/rel), f'{fig["id"]}: retained input changed: {rel}'
        if not images_equal(source/fig['reference'], work/fig['output']):
            raise ValueError(f'{fig["id"]}: SVG pixels differ; check the recorded fonts and package versions')
        pdf_checked = bool(fig.get('pdf_reference'))
        if pdf_checked and not images_equal(source/fig['pdf_reference'], work/fig['pdf_output']):
            raise ValueError(f'{fig["id"]}: PDF pixels differ; check the recorded fonts and package versions')
        dest = out/fig['id']
        dest.mkdir(exist_ok=True)
        shutil.copy2(work/fig['output'],dest/Path(fig['output']).name)
        if pdf_checked:
            shutil.copy2(work/fig['pdf_output'],dest/Path(fig['pdf_output']).name)
        return {'id':fig['id'],'status':'PASS','tables':len(fig['tables']),
                'generated_tables':len(fig['generated_tables']), 'retained_tables':len(fig['input_tables']),
                'svg_pixel_match':True,'pdf_pixel_match':True if pdf_checked else None,
                'write_guard':'PASS','seconds':round(time.monotonic()-started,2)}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rebuild',action='store_true',help='rebuild all 17 figures from saved inputs')
    parser.add_argument('--figures',nargs='+',help='rebuild selected IDs, for example F3 S11')
    parser.add_argument('--output-dir',type=Path,default=Path('verification_output'))
    args=parser.parse_args()
    args.output_dir=args.output_dir.resolve()
    args.output_dir.mkdir(parents=True,exist_ok=True)
    report={'release':'2026.10.10','python':sys.version,'checks':None,'rebuilds':[]}
    try:
        report['checks'], figures=check_release()
        print(json.dumps(report['checks']),flush=True)
        if args.rebuild or args.figures:
            if args.figures:
                unknown=set(args.figures)-{f['id'] for f in figures}
                if unknown:
                    raise ValueError(f'Unknown figure IDs: {sorted(unknown)}')
                figures=[f for f in figures if f['id'] in args.figures]
            for fig in figures:
                result=rebuild_figure(fig,args.output_dir)
                report['rebuilds'].append(result)
                print(json.dumps(result),flush=True)
        report['status']='PASS'
    except Exception as exc:
        report['status']='FAIL'
        report['error']=f'{type(exc).__name__}: {exc}'
        print(report['error'],file=sys.stderr)
    (args.output_dir/'verification_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return 0 if report['status']=='PASS' else 1

if __name__=='__main__':
    raise SystemExit(main())
