"""Calculate baseline-relative enrichment changes for the display tables."""
from pathlib import Path
import hashlib,json
import pandas as pd
ROOT=Path(__file__).resolve().parent
src=ROOT/'plotted_data/Figure_3B_plotted_values.csv'
df=pd.read_csv(src,float_precision='round_trip').sort_values('excl_kb').copy()
baseline=df.loc[df.excl_kb.eq(0)].iloc[0]
for region,field in [('boundary','edge_fold'),('interior','interior_fold')]:
    df[region+'_baseline_fold']=float(baseline[field])
    df[region+'_relative_change_pct']=100*(df[field]/baseline[field]-1)
df['display_unit']='percent change in fold relative to the same region at 0 kb'
target=ROOT/'plotted_data/Figure_3E_display_values.csv'
df.to_csv(target,index=False,float_format='%.17g')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
(ROOT/'DISPLAY_PROVENANCE.json').write_text(json.dumps(dict(
    input='plotted_data/'+src.name,input_sha256=sha(src),output='plotted_data/'+target.name,output_sha256=sha(target),
    formula='100 * (fold_at_radius / fold_at_0_kb - 1)',new_tests=False),indent=2)+'\n',encoding='utf-8')
print(df[['excl_kb','boundary_relative_change_pct','interior_relative_change_pct']].to_string(index=False))
