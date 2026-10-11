"""Prepare Figure 6 tables and reproduce the displayed regressions.

Panel C uses unweighted ordinary least squares on sample means and two-sided
t intervals with n-2 degrees of freedom. Earlier normal-approximation summaries
are retained in columns with the legacy_ prefix for comparison.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy import stats

HERE=Path(__file__).resolve().parent
(HERE/'review').mkdir(exist_ok=True)
(HERE/'plotted_data').mkdir(exist_ok=True)
SRC=HERE/'source_data'
OUT=HERE/'plotted_data'
OUT.mkdir(exist_ok=True)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(name):return pd.read_csv(SRC/name,float_precision='round_trip')
def save(frame,name):frame.to_csv(OUT/name,index=False,float_format='%.17g')

a=read('Figure_6A.csv')
assert len(a)==2 and set(a.subset)=={'PRODUCTIVE','SILENT'}
a['detection_label']=a.subset.map({'PRODUCTIVE':'Fusion transcript detected','SILENT':'No fusion transcript detected'})
a['estimate']=a.fold
a['lower']=a.fold_ci_lo
a['upper']=a.fold_ci_hi
a['interval_type']='95% spatial reference range; not donor-sampling CI'
a['estimate_unit']='observed / matched-null fold'
save(a,'Figure_6A_plotted_values.csv')

b=read('Figure_6B.csv')
assert len(b)==2 and b.plot_label.tolist()==['All amplicons','Circular class']
b['display_label']=b.plot_label
b['observed_label']='Observed integration clusters'
b['reference_label']='Matched-null mean integration clusters'
b['estimate']=b.fold
b['lower']=b.fold_ci_lo
b['upper']=b.fold_ci_hi
b['reference_upper_available']=b.fold_ci_hi.notna()
b['interval_type']='95% spatial reference range; not donor-sampling CI'
b['estimate_unit']='observed / matched-null fold'
assert b.reference_upper_available.tolist()==[True,False]
assert pd.isna(b.loc[b.amp_class=='Circular','upper']).all()
save(b,'Figure_6B_plotted_values.csv')

c=read('Figure_6C.csv')
assert len(c)==63 and set(c.arm)=={'amplicon','placebo'}
assert c.groupby('arm').size().to_dict()=={'amplicon':33,'placebo':30}
assert not c.duplicated(['arm','sample']).any()
inside_samples=set(c.loc[c.arm=='amplicon','sample'])
outside_samples=set(c.loc[c.arm=='placebo','sample'])
assert outside_samples.issubset(inside_samples)
nmax=float(c.n_integrations.max())
c['marker_area']=8.0+46.0*c.n_integrations/nmax
c['marker_area_unit']='points squared'
c['marker_area_rule']='8 + 46 * n_integrations / common maximum across both groups'
c['marker_area_common_max_integrations']=nmax
c['arm_label']=c.arm.map({'amplicon':'Inside amplicons','placebo':'Outside amplicons'})
save(c,'Figure_6C_plotted_values.csv')

legacy=read('Figure_6C_fits.csv')
assert set(legacy.arm)==set(c.arm) and len(legacy)==2
fit_rows=[]
line_rows=[]
expected_display={
    'amplicon':[1.02753696,.34996055,1.70511337,.00417316],
    'placebo':[.02882561,-.42007972,.47773095,.89629300],
}
for legacy_row in legacy.to_dict('records'):
    arm=legacy_row['arm']
    g=c[c.arm==arm]
    n=len(g)
    fit=stats.linregress(g.log_cn.to_numpy(),g.log_support.to_numpy())
    df=n-2
    crit=float(stats.t.ppf(.975,df))
    lower=float(fit.slope-crit*fit.stderr)
    upper=float(fit.slope+crit*fit.stderr)
    current=np.array([fit.slope,lower,upper,fit.pvalue])
    assert np.allclose(current,expected_display[arm],atol=5e-8,rtol=0),(arm,current)
    row={'legacy_'+key:value for key,value in legacy_row.items()}
    row.update(arm=arm,arm_label='Inside amplicons' if arm=='amplicon' else 'Outside amplicons',
               slope=float(fit.slope),intercept=float(fit.intercept),stderr=float(fit.stderr),
               df=df,t_critical=crit,ci_lo=lower,ci_hi=upper,p=float(fit.pvalue),
               n_units=n,n_integrations=int(g.n_integrations.sum()),
               estimator='Unweighted OLS on sample means',unit='sample mean',
               interval_type='95% OLS slope t interval with n_units - 2 degrees of freedom',
               fit_source='Figure_6C.csv: log_cn and log_support, grouped by arm',
               algorithm='scipy.stats.linregress; scipy.stats.t.ppf(0.975, n_units-2)',
               existing_analysis_reproduced=True,new_statistical_test=False,
               x_min=float(g.log_cn.min()),x_max=float(g.log_cn.max()))
    fit_rows.append(row)
    for i,x in enumerate(np.linspace(g.log_cn.min(),g.log_cn.max(),50)):
        line_rows.append({'arm':arm,'point_index':i,'x':float(x),'y':float(fit.intercept+fit.slope*x),
                          'x_variable':'log_cn','y_variable':'fitted log_support',
                          'line_type':'OLS fitted mean line; not a confidence band',
                          'x_extent':'observed within-group minimum to maximum'})
fits=pd.DataFrame(fit_rows)
lines=pd.DataFrame(line_rows)
save(fits,'Figure_6C_fits_plotted_values.csv')
save(lines,'Figure_6C_line_data.csv')
assert fits.set_index('arm').n_integrations.to_dict()=={'amplicon':87,'placebo':191}
assert lines.groupby('arm').size().to_dict()=={'amplicon':50,'placebo':50}

# Original A/B/C columns must survive CSV round trips exactly, including NaN.
for letter in ['A','B','C']:
    original=read(f'Figure_6{letter}.csv')
    current=pd.read_csv(OUT/f'Figure_6{letter}_plotted_values.csv',float_precision='round_trip')
    pd.testing.assert_frame_equal(original,current[original.columns],check_exact=True,check_dtype=False)
current=pd.read_csv(OUT/'Figure_6C_fits_plotted_values.csv',float_precision='round_trip')
pd.testing.assert_frame_equal(legacy,current[['legacy_'+col for col in legacy]].rename(columns=lambda col:col.removeprefix('legacy_')),check_exact=True,check_dtype=False)

mapping=[
    ('A','Figure_6A.csv','Figure_6A_plotted_values.csv','Original values, detection labels and saved spatial-range aliases only.'),
    ('B','Figure_6B.csv','Figure_6B_plotted_values.csv','Original values and descriptive labels; missing circular-class upper endpoint remains NaN.'),
    ('C points','Figure_6C.csv','Figure_6C_plotted_values.csv','Original 63 sample means; original shared marker-area rule.'),
    ('C fit legacy','Figure_6C_fits.csv','Figure_6C_fits_plotted_values.csv','All original fit columns preserved with legacy_ prefix; they do not define the current t intervals.'),
    ('C fit current','Figure_6C.csv','Figure_6C_fits_plotted_values.csv','Reproduce existing OLS slope/intercept/p and t confidence intervals from displayed sample means.'),
    ('C line','Figure_6C.csv','Figure_6C_line_data.csv','50 points on each reproduced fitted line, restricted to observed within-group x range; no confidence bands.'),
]
maps=[]
for panel,source,plotted,operation in mapping:
    maps.append({'panel':panel,'local_source':'source_data/'+source,
                 'original_source_relative_to_package':'source_data/'+source,
                 'source_sha256':sha(SRC/source),'plotted_file':'plotted_data/'+plotted,
                 'plotted_sha256':sha(OUT/plotted),'operation':operation})
pd.DataFrame(maps).to_csv(HERE/'SOURCE_MAP.csv',index=False)
qa={'status':'PASS','new_statistical_tests':False,'existing_OLS_display_reproduced':True,
    'source_original_columns_exact_after_round_trip':True,
    'C_legacy_columns_preserved_with_prefix':True,
    'rows':{'A':len(a),'B':len(b),'C_points':len(c),'C_fits':len(fits),'C_line':len(lines)},
    'B_reference_upper_available':b.reference_upper_available.tolist(),
    'C_marker_area_common_max_integrations':nmax,
    'C_unique_samples_across_both_groups':int(c['sample'].nunique()),
    'C_outside_samples_are_subset_of_inside_samples':True,
    'C_sample_dependence_note':'63 displayed sample means from 33 unique samples; separate unweighted OLS fits, not 63 independent samples or a direct between-slope test.',
    'C_fits':fits[['arm','slope','intercept','stderr','df','t_critical','ci_lo','ci_hi','p','n_units','n_integrations']].to_dict('records'),
    'C_line_x_extent':'Observed minimum and maximum within each group; 50 points; no extrapolation or confidence bands.',
    'software':{'scipy':__import__('scipy').__version__,'pandas':pd.__version__,'numpy':np.__version__}}
(OUT/'preparation_QA.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
print(json.dumps(qa,indent=2))
