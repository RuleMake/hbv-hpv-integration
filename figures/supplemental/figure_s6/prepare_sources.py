"""Prepare presentation values from saved results; no statistical analyses are run."""
from pathlib import Path
import csv, hashlib, json, math, re
ROOT=Path(__file__).resolve().parent
(ROOT/'review').mkdir(exist_ok=True)
(ROOT/'plotted_data').mkdir(exist_ok=True)
NUMBER=6
KEY='s5'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
src=ROOT/f'source_data/S{NUMBER}_source_data.csv'
with src.open(encoding='utf-8-sig',newline='') as h:
 reader=csv.DictReader(h);original=list(reader);fields=reader.fieldnames
rows=[dict(r) for r in original]
for i,r in enumerate(rows):
 r['source_row']=str(i+2)
 if NUMBER in (5,6):
  r['estimate']=r['fold'];r['lower']=r['ci_lo'];r['upper']=r['ci_hi']
  r['interval_type']='95% spatial-null reference range' if r['panel'] in ('a','c') and NUMBER==5 or r['panel']=='a' and NUMBER==6 else ('95% bootstrap CI' if NUMBER==6 and r['panel']=='b' else 'not applicable')
 if NUMBER==5 and r['panel']=='c':
  t=next(x for x in original if x['panel']=='c' and x['virus']==r['virus'] and x['subset']=='ALL (joint reference)')
  r['coordinate_share_pct']=str(100*float(r['n_coords'])/float(t['n_coords']))
  r['cluster_share_pct']=str(100*float(r['n_clusters'])/float(t['n_clusters']))
  r['total_coordinates']=t['n_coords'];r['total_clusters']=t['n_clusters']
 if NUMBER==7:
  for f in ('excess','mde','excess_20k','mde_excess_20k'):
   r[f+'_pp']=str(100*float(r[f])) if r[f] else ''
  if r['panel']=='a':
   r['lower_pp']=str(100*(float(r['excess'])-1.95996398454*float(r['null_sd'])))
   r['upper_pp']=str(100*(float(r['excess'])+1.95996398454*float(r['null_sd'])))
   r['interval_type']='95% normal-approximation spatial-null reference range'
  if r['panel']=='c':
   r['quantity_type']='arithmetic contrast' if r['rung'] in ('arithmetic sum of the three singles','breakage-density step') else 'separately evaluated matching scheme'
   r['display_measured']='False' if r['quantity_type']=='arithmetic contrast' else 'True'
   r['classification_basis']='The sum of individual reductions and the density increment are arithmetic comparisons, not separately fitted matching schemes.'
  if r['panel']=='b' and not r['numeral']:
   for f in ('excess_pp','mde_pp','excess_20k_pp','mde_excess_20k_pp'):r[f]=''
output=ROOT/f'plotted_data/S{NUMBER}_plotted_values.csv'
allfields=list(dict.fromkeys(k for row in rows for k in row))
with output.open('w',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=allfields);w.writeheader();w.writerows(rows)
assert all(all(r[k]==s[k] for k in fields) for r,s in zip(rows,original))
extra=[]
if NUMBER==6:
 c=[x for x in original if x['panel']=='c']; flag=next(x for x in c if x['arm']=='what the flags identify')
 m=re.search(r'spanning ([\d,]+) of ([\d,]+) caseids',flag['note'])
 ncases,totalcases=[int(x.replace(',','')) for x in m.groups()]
 extra=[dict(label='Annotated calls',count=int(flag['n_events']),total=int(float(flag['n_coords'])),type='annotation'),dict(label='Affected cases',count=ncases,total=totalcases,type='cases')]
 for cl in ('DEL','DUP','INV','TRA'):
  t=next(x for x in c if x['svclass']==cl); count=int(float(t['observed']));total=count+int(t['n_events'])
  extra.append(dict(label=cl,count=count,total=total,type='events removed'))
 extra.append(dict(label='All classes',count=sum(t['count'] for t in extra[2:]),total=sum(t['total'] for t in extra[2:]),type='events removed'))
 for t in extra:t['percent']=100*t['count']/t['total']
 ep=ROOT/'plotted_data/S6_exclusion_accounting.csv'
 with ep.open('w',encoding='utf-8',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(extra[0]));w.writeheader();w.writerows(extra)
maps=[]
for p in sorted((ROOT/'plotted_data').glob('*.csv')):
 maps.append(dict(original_source_relative_to_package=src.relative_to(ROOT).as_posix(),local_source=src.relative_to(ROOT).as_posix(),source_sha256=sha(src),plotted_file=p.relative_to(ROOT).as_posix(),plotted_sha256=sha(p),transformation=('Derived exclusion accounting: counts and denominators are extracted from saved numeric columns and the original case-count note; percent = 100 * count / total. All original source fields are retained separately in S6_plotted_values.csv.' if p.name=='S6_exclusion_accounting.csv' else 'Original columns preserved verbatim; deterministic presentation aliases and percentage conversion only.')))
with (ROOT/'SOURCE_MAP.csv').open('w',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(maps[0]));w.writeheader();w.writerows(maps)
(ROOT/'review/PREPARATION_QA.json').write_text(json.dumps(dict(source_rows=len(original),original_columns_preserved=True,source_sha256=sha(src),new_statistical_analyses=False),indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(figure=NUMBER,rows=len(rows),original_columns_preserved=True)))
