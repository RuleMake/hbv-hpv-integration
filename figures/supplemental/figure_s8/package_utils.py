"""Validate bundled saved data and record reproducible local export checks."""
from pathlib import Path
import hashlib,json
import pandas as pd
import fitz

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def finalize(root,stem,qa,checks):
    root=Path(root)
    (root/'review').mkdir(exist_ok=True)
    rows=[]
    for item in json.loads((root/'source_data/SOURCE_HASHES.json').read_text(encoding='utf-8')):
        raw=root/item['source'];plot=root/item['plotted']
        assert sha(raw)==sha(plot)==item['sha256'],raw.name
        rows.append(dict(local_source=item['source'],
            original_source_relative_to_package=item['source'],
            source_sha256=sha(raw),plotted_file=item['plotted'],
            plotted_sha256=sha(plot),operation='Byte-identical saved source table; selection and display transformations are explicit in the build script. No statistical re-estimation.'))
    pd.DataFrame(rows).to_csv(root/'SOURCE_MAP.csv',index=False)
    doc=fitz.open(root/'rendered'/f'{stem}.pdf')
    assert len(doc)==1
    doc[0].get_pixmap(dpi=220,alpha=False).save(root/'review'/f'{stem}_pdf_render.png')
    checks.update(bundled_source_hashes_equal=True,pdf_pages=len(doc),width_mm=round(doc[0].rect.width*25.4/72,3),
        min_font_pt=qa['min_font_pt'],text_overlap_count=len(qa['text_overlaps']),
        text_outside_canvas=qa['text_outside_canvas'],label_warnings=qa['label_warnings'])
    doc.close()
    (root/'review/NUMERICAL_CHECKS.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(checks,ensure_ascii=False))
