"""Prepare Figure 5 display tables from saved results.

Source columns and numerical precision are retained; additional columns
provide display labels and unit conversions.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
(HERE/'review').mkdir(exist_ok=True)
(HERE/'plotted_data').mkdir(exist_ok=True)
SRC=HERE/'source_data'
OUT=HERE/'plotted_data'
OUT.mkdir(exist_ok=True)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(name):
    return pd.read_csv(SRC/name,float_precision='round_trip')

def save(frame,name):
    frame.to_csv(OUT/name,index=False,float_format='%.17g')

a=read('Figure_5A.csv')
assert len(a)==3 and set(a.track)=={'consensus','helas3','hepg2'}
a['estimate']=a.excess
a['lower']=a.excess_reference_low
a['upper']=a.excess_reference_high
for col in ['estimate','lower','upper']:
    a[col+'_pp']=a[col]*100.0
a['estimate_unit']='composition-score difference'
a['display_unit']='percentage points'
a['track_label']=a.track.map({'consensus':'Consensus','helas3':'HeLa-S3','hepg2':'HepG2'})
save(a,'Figure_5A_plotted_values.csv')

b=read('Figure_5B.csv')
assert len(b)==6
b['estimate']=b.fold
b['lower']=b.fold_null_low
b['upper']=b.fold_null_high
b['estimator_label']=b.estimator.map({
    'cluster_event_incidence':'Cluster event incidence',
    'calibrated_cluster_sign':'Calibrated cluster sign',
})
assert b.estimator_label.notna().all()
b['arm_label']=b.arm.map({'joint_geometry':'Joint geometry','gene_status':'Gene status','gene_length':'Gene length'})
assert b.arm_label.notna().all()
b['estimate_unit']='observed / matched-null fold'
b['display_unit']='fold'
save(b,'Figure_5B_plotted_values.csv')

c=read('Figure_5C.csv')
assert len(c)==3
assert c.analysis.tolist()==['HBV_pooled_donor_mean','HPV_pooled_donor_mean','HPV_common_support_histology_standardized']
c['mean_detected_pct']=c.mean_detected*100.0
c['mean_comparison_pct']=c.mean_comparison*100.0
c['estimate']=c.mean_ratio
c['lower']=c.ratio_ci_lo
c['upper']=c.ratio_ci_hi
c['mean_display_type']=np.where(c.analysis=='HPV_common_support_histology_standardized','histology-standardized','donor-average')
c['interval_type']='95% donor-bootstrap CI'
c['estimate_unit']='ratio of group means'
c['mean_display_unit']='percent'
c['ratio_display_unit']='ratio'
save(c,'Figure_5C_plotted_values.csv')

for name,df in [('A',a),('B',b),('C',c)]:
    assert (df.lower<=df.estimate).all() and (df.estimate<=df.upper).all(),name
    # Round-trip checking verifies all preserved original columns, including NaN.
    src=read(f'Figure_5{name}.csv')
    out=pd.read_csv(OUT/f'Figure_5{name}_plotted_values.csv',float_precision='round_trip')
    pd.testing.assert_frame_equal(src,out[src.columns],check_dtype=False,check_exact=True)

maps=[]
for panel,operation in [
    ('A','Preserve three replication-timing contrasts and saved spatial reference ranges; add estimate/lower/upper aliases and ×100 percentage-point display columns.'),
    ('B','Preserve six gene-matching estimates and saved spatial reference ranges; add readable estimator/arm labels and fold aliases.'),
    ('C','Preserve three donor mean ratios and donor-bootstrap confidence intervals; add group-mean percentages and donor-average versus histology-standardized labels.'),
]:
    source=f'Figure_5{panel}.csv'
    plotted=f'Figure_5{panel}_plotted_values.csv'
    maps.append({'panel':panel,'local_source':'source_data/'+source,
                 'original_source_relative_to_package':'source_data/'+source,
                 'source_sha256':sha(SRC/source),'plotted_file':'plotted_data/'+plotted,
                 'plotted_sha256':sha(OUT/plotted),'operation':operation})
pd.DataFrame(maps).to_csv(HERE/'SOURCE_MAP.csv',index=False)
qa={'status':'PASS','rows':{'A':len(a),'B':len(b),'C':len(c)},
    'all_original_columns_exact_after_round_trip':True,
    'new_statistical_tests':False,
    'A_units':'Original score difference; displayed in percentage points by multiplying estimate and both saved endpoints by 100.',
    'B_units':'Unchanged observed / matched-null fold.',
    'C_units':'Means shown as percent; ratios and donor-bootstrap confidence endpoints unchanged.',
    'C_mean_types':c.mean_display_type.tolist(),
    'interval_distinction':'A/B retain spatial-reference semantics; C is a 95% donor-bootstrap CI.'}
(OUT/'preparation_QA.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
print(json.dumps(qa,indent=2))
