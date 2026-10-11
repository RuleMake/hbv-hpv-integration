"""Prepare presentation values from saved results; no statistical analyses are run."""
from pathlib import Path
import csv, hashlib, json, math, re
ROOT=Path(__file__).resolve().parent
(ROOT/'review').mkdir(exist_ok=True)
(ROOT/'plotted_data').mkdir(exist_ok=True)
NUMBER=7
KEY='s10'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
src=ROOT/f'source_data/S{NUMBER}_source_data.csv'
with src.open(encoding='utf-8-sig',newline='') as h:
 reader=csv.DictReader(h);original=list(reader);fields=reader.fieldnames
with (ROOT/'source_data/conditioning_ladder_saved.csv').open(encoding='utf-8-sig',newline='') as h:
 saved_formula_rows={(q['rung'],q['virus']):q for q in csv.DictReader(h)}
reconstructed=[]
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
   r['classification_basis']='The sum of individual reductions and the density increment are arithmetic comparisons rather than separately fitted matching schemes.'
  if r['panel']=='b' and r['numeral']:
   saved=saved_formula_rows[(r['rung'],r['virus'])]
   for key in ('excess','z','axis_observed','axis_null_mean','n_clusters'):
    assert abs(float(saved[key])-float(r[key]))<1e-12
   effect=abs(float(r['excess'])); zabs=abs(float(r['z']))
   assert effect>0 and zabs>.005 and not r['mde']
   mde=280*effect/zabs
   lower=280*(effect-.000005)/(zabs+.005)
   upper=280*(effect+.000005)/(zabs-.005)
   flag=r['clears_mde'].lower()=='true'
   assert (100*effect>=mde)==flag
   assert (zabs-.005>=2.8) if flag else (zabs+.005<2.8)
   r['approx_mde_reconstructed_pp']=str(mde)
   r['approx_mde_rounding_lower_pp']=str(lower)
   r['approx_mde_rounding_upper_pp']=str(upper)
   r['approx_mde_signed_position_pp']=str(mde if r['virus']=='HBV' else -mde)
   r['approx_mde_basis']='100 * 2.80 * abs(saved rounded excess / saved rounded z); original MDE field remains empty'
   r['approx_mde_classification_stable']='True'
   reconstructed.append({'virus':r['virus'],'rung':r['rung'],'mde_pp':mde,'lower_pp':lower,'upper_pp':upper,'classification_stable':True})
  if r['panel']=='c':
   ordering=['1_replication_timing','2_gene_density','3_fragile_site','4_stratum_18way','5_joint_142way','6_joint_plus_density','arithmetic sum of the three singles','breakage-density step']
   r['panel_C_display_order']=str(ordering.index(r['rung']))
   r['panel_C_position']=str([0,1,2,3,4,5,7.5,8.5][ordering.index(r['rung'])])
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
for p in [output]:
 maps.append(dict(original_source_relative_to_package=src.relative_to(ROOT).as_posix(),local_source=src.relative_to(ROOT).as_posix(),source_sha256=sha(src),plotted_file=p.relative_to(ROOT).as_posix(),plotted_sha256=sha(p),transformation=('Original columns preserved verbatim; deterministic presentation aliases and percentage conversion only.' + (' S7 distinguishes arithmetic comparisons from evaluated matching schemes.' if NUMBER==7 else ''))))
legend_source=ROOT/'legend.md'
maps.append(dict(original_source_relative_to_package=legend_source.relative_to(ROOT).as_posix(),local_source=legend_source.relative_to(ROOT).as_posix(),source_sha256=sha(legend_source),plotted_file=output.relative_to(ROOT).as_posix(),plotted_sha256=sha(output),transformation='Classification provenance only: the baseline caption identifies the sum of single factors and density increment as arithmetic quantities. Numeric source values and the legacy measured column are preserved.'))
for filename,note in [
 ('conditioning_ladder_saved.csv','Saved per-virus effect and z; original per-virus MDE is empty. Approximate display threshold reconstructed as 100 * 2.80 * abs(excess/z), with rounding bounds.'),
 ('method_notes.md','Formulas, saved numerical precision, detectable-effect markers and rounding checks used in the display.')]:
 p=ROOT/'source_data'/filename
 maps.append(dict(original_source_relative_to_package=p.relative_to(ROOT).as_posix(),local_source=p.relative_to(ROOT).as_posix(),source_sha256=sha(p),plotted_file=output.relative_to(ROOT).as_posix(),plotted_sha256=sha(output),transformation=note))
preserved_source=ROOT/'source_data/conditioning_ladder_saved.csv'
preserved_display=ROOT/'plotted_data/S7_conditioning_ladder_saved_values.csv'
preserved_display.write_bytes(preserved_source.read_bytes())
maps.append(dict(original_source_relative_to_package=preserved_source.relative_to(ROOT).as_posix(),local_source=preserved_source.relative_to(ROOT).as_posix(),source_sha256=sha(preserved_source),plotted_file=preserved_display.relative_to(ROOT).as_posix(),plotted_sha256=sha(preserved_display),transformation='Complete 21-row saved conditioning-ladder table preserved byte-for-byte as provenance; not an additional plotted analysis. Fourteen per-virus rows supply the saved rounded excess and z used for the approximate display ticks.'))
with (ROOT/'SOURCE_MAP.csv').open('w',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(maps[0]));w.writeheader();w.writerows(maps)
(ROOT/'review/PREPARATION_QA.json').write_text(json.dumps(dict(source_rows=len(original),original_columns_preserved=True,source_sha256=sha(src),new_statistical_analyses=False),indent=2)+'\n',encoding='utf-8')
assert len(reconstructed)==14
(ROOT/'review/MDE_RECONSTRUCTION_QA.json').write_text(json.dumps({'n_reconstructed':14,'original_mde_fields_unchanged_and_empty':True,'all_nonzero_z':True,'all_classifications_stable_under_rounding':True,'max_rounding_error_pp':max(max(q['mde_pp']-q['lower_pp'],q['upper_pp']-q['mde_pp']) for q in reconstructed),'rows':reconstructed},indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(figure=NUMBER,rows=len(rows),original_columns_preserved=True,reconstructed_thresholds=14)))
