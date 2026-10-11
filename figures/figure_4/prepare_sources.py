"""Prepare Figure 4 display tables from saved results.

Source columns are retained. Additional columns supply display labels and
percentage-point conversions for the reported estimates and reference ranges.
"""
from pathlib import Path
import hashlib
import json
import pandas as pd
import numpy as np

HERE = Path(__file__).resolve().parent
(HERE/'review').mkdir(exist_ok=True)
(HERE/'plotted_data').mkdir(exist_ok=True)
SRC = HERE / 'source_data'
OUT = HERE / 'plotted_data'
OUT.mkdir(exist_ok=True)

def read(name):
    return pd.read_csv(SRC / name, float_precision='round_trip')

def save(frame, name):
    frame.to_csv(OUT / name, index=False, float_format='%.17g')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

SPATIAL = '95% spatial reference range; reversed, centred null quantiles; not donor-sampling CI'
POINT_UNIT = 'proportion difference'
DISPLAY_UNIT = 'percentage points'

# A: original 24-cell enrichment atlas. Preserve the original estimators and
# multiplicity columns; the label is for display only.
a = read('Figure_4A.csv')
assert len(a) == 24
a['estimator_label'] = a.estimator.map({
    'bc_uniqSV_blockNull': 'Call count',
    'f_clusterSignTest': 'Cluster sign',
})
assert a.estimator_label.notna().all()
a['estimate'] = a.fold
a['lower'] = a.fold_null_interval_low
a['upper'] = a.fold_null_interval_high
a['interval_type'] = '95% spatial reference range; observed / reversed null quantiles; not donor-sampling CI'
a['estimate_unit'] = 'observed / matched-null fold'
save(a, 'Figure_4A_plotted_values.csv')

# B detail: class-fraction displacements within the same paired-call,
# stratum-matched, full-rung configuration as the saved total axis.
full = read('class_composition.csv')
classes = ['DEL', 'DUP', 'INV', 'TRA']
details = full[
    full.virus.isin(['HBV', 'HPV']) & (full.arm == 'stratum') &
    (full.rung == 'full') & (full.variant == 'paired_sample') &
    (full.window_bp == 1000000) &
    full.statistic.isin(['prop_' + x for x in classes])
].copy()
assert len(details) == 8
details['svclass'] = details.statistic.str.removeprefix('prop_')
details['estimate'] = details.excess
details['lower'] = details.excess - (details.null_reference_q975 - details.null_mean)
details['upper'] = details.excess - (details.null_reference_q025 - details.null_mean)
details['source_interval_type'] = details.interval_type
details['interval_type'] = SPATIAL

total = read('Figure_4B.csv')
assert len(total) == 2
assert set(total.virus) == {'HBV', 'HPV'}
assert (total.arm == 'stratum').all() and (total.rung == 'full').all()
assert (total.variant == 'paired_sample').all()
total['estimate'] = total.excess
total['lower'] = total.excess - (total.null_reference_q975 - total.null_mean)
total['upper'] = total.excess - (total.null_reference_q025 - total.null_mean)
total['source_interval_type'] = total.interval_type
total['interval_type'] = SPATIAL

# C: preserve the four existing contrasts and their original reference rule.
c = read('Figure_4C.csv')
assert len(c) == 4
c['estimate'] = c.observed
c['lower'] = c.observed - c.null_reference_q975
c['upper'] = c.observed - c.null_reference_q025
c['source_interval_type'] = c.interval_type
c['interval_type'] = '95% spatial reference range; observed contrast minus reversed null quantiles; not donor-sampling CI'

for df in [details, total, c]:
    for field in ['estimate', 'lower', 'upper']:
        df[field + '_pp'] = df[field] * 100.0
    df['estimate_unit'] = POINT_UNIT
    df['display_unit'] = DISPLAY_UNIT
    assert (df.lower <= df.estimate).all() and (df.estimate <= df.upper).all()

checks = []
for virus in ['HBV', 'HPV']:
    sub = details[details.virus == virus].set_index('svclass').loc[classes]
    assert len(sub) == 4 and sub.index.is_unique
    reconstructed = sub.loc['TRA','estimate'] + sub.loc['INV','estimate'] - sub.loc['DEL','estimate'] - sub.loc['DUP','estimate']
    summed = float(sub.estimate.sum())
    target = total.loc[total.virus == virus, 'estimate'].item()
    assert np.isclose(summed, 0, atol=1e-12, rtol=0)
    assert np.isclose(reconstructed, target, atol=1e-12, rtol=0)
    # This exact linear identity holds for observed values and null means.
    # Marginal class ranges are NOT summed to obtain the total-axis interval.
    checks.append({'virus':virus, 'four_class_displacement_sum':summed,
                   'reconstructed_axis_displacement':float(reconstructed),
                   'saved_axis_displacement':float(target),
                   'difference':float(reconstructed-target), 'pass':True})

details = details.sort_values(['virus','svclass'], key=lambda col: col.map({'HBV':0,'HPV':1}) if col.name=='virus' else col.map({v:i for i,v in enumerate(classes)}))
save(details, 'Figure_4B_detail_plotted_values.csv')
save(total, 'Figure_4B_total_plotted_values.csv')
save(c, 'Figure_4C_plotted_values.csv')

origins = [
    ('Figure_4A.csv', 'source_data/Figure_4A.csv', 'Figure_4A_plotted_values.csv', 'Original 24 estimates/ranges; readable estimator labels and aliases only'),
    ('class_composition.csv', '../../../class_composition/axis/outputs/paired_sample/class_composition.csv', 'Figure_4B_detail_plotted_values.csv', 'HBV/HPV, paired_sample, stratum, full, 1 Mb, prop_DEL/DUP/INV/TRA; saved displacements and reversed centred ranges; multiply by 100 for pp'),
    ('Figure_4B.csv', 'source_data/Figure_4B.csv', 'Figure_4B_total_plotted_values.csv', 'Two saved axis displacements; reversed centred ranges; multiply by 100 for pp'),
    ('Figure_4C.csv', 'source_data/Figure_4C.csv', 'Figure_4C_plotted_values.csv', 'Four saved contrasts; observed minus reversed null quantiles; multiply by 100 for pp'),
]
rows=[]
for source,origin,plotted,operation in origins:
    rows.append({'local_source':'source_data/'+source,'original_source_relative_to_package':'source_data/'+source,
                 'source_sha256':sha(SRC/source),'plotted_file':'plotted_data/'+plotted,
                 'plotted_sha256':sha(OUT/plotted),'operation':operation})
pd.DataFrame(rows).to_csv(HERE/'SOURCE_MAP.csv',index=False)
qa={'status':'PASS','new_statistical_tests':False,'rows':{'A':len(a),'B_detail':len(details),'B_total':len(total),'C':len(c)},
    'reconstruction_checks':checks,
    'display_conversion':'estimate_pp/lower_pp/upper_pp = original proportion difference × 100',
    'interval_note':'Saved spatial reference quantiles; no donor confidence intervals and no sum of marginal class intervals.'}
(OUT/'preparation_QA.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
print(json.dumps(qa,indent=2))
