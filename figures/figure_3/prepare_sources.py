"""Prepare Figure 3 composition tables from saved results.

Null composition is reconstructed from the saved six-decimal balance
displacement and exact overlap counts, rather than a ratio of expected counts.
The enrichment tables are copied unchanged; the displacement table retains
its source columns and adds percentage-point conversions.
"""
from pathlib import Path
import hashlib
import json
import shutil
import pandas as pd
import numpy as np

HERE=Path(__file__).resolve().parent
(HERE/'review').mkdir(exist_ok=True)
(HERE/'plotted_data').mkdir(exist_ok=True)
SRC=HERE/'source_data'
OUT=HERE/'plotted_data'
OUT.mkdir(exist_ok=True)

def read(path):return pd.read_csv(path,float_precision='round_trip')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(frame,name):frame.to_csv(OUT/name,index=False,float_format='%.17g')

# Copy the saved enrichment tables and retain source columns in the
# displacement table. A and B are byte-for-byte copies.
for panel in 'AB':
    name=f'Figure_3{panel}_plotted_values.csv'
    shutil.copyfile(SRC/'display_source'/name,OUT/name)
    assert sha(SRC/'display_source'/name)==sha(OUT/name)
a=read(OUT/'Figure_3A_plotted_values.csv')
b=read(OUT/'Figure_3B_plotted_values.csv')
c0=read(SRC/'display_source/Figure_3C_plotted_values.csv')
assert [len(a),len(b),len(c0)]==[2,4,12]
c=c0.copy()
for column in ['displacement','lower','upper','mde_displacement']:
    c[column+'_pp']=c[column]*100.0
c['displacement_display_unit']='percentage points of the boundary-minus-interior balance'
save(c,'Figure_3C_plotted_values.csv')
c_back=read(OUT/'Figure_3C_plotted_values.csv')
pd.testing.assert_frame_equal(c0,c_back[c0.columns],check_exact=True,check_dtype=False)

geometry=read(SRC/'integrase_positive_geometry.csv')
source=geometry[(geometry.statistic=='edge_vs_interior')&geometry.virus.isin(['HPV','HBV','HTLV_all'])]
assert len(source)==12
for r in c.itertuples():
    original=source[(source.virus==r.virus)&(source.excl_kb==r.excl_kb)].iloc[0]
    for field in ['displacement','n_edge','n_interior','n_coords','n_clusters','edge_fold','interior_fold','null_sd','mde_displacement','z','emp_p','n_perm']:
        assert getattr(r,field)==original[field],(r.virus,r.excl_kb,field)
    lower=original.displacement-(original.null_ci_hi-original.expected)
    upper=original.displacement-(original.null_ci_lo-original.expected)
    assert np.isclose(r.lower,lower,atol=1e-15,rtol=0)
    assert np.isclose(r.upper,upper,atol=1e-15,rtol=0)

composition=[]
checks=[]
for virus in ['HPV','HBV','HTLV_all']:
    row=c[(c.virus==virus)&(c.excl_kb==0)].iloc[0]
    src=source[(source.virus==virus)&(source.excl_kb==0)].iloc[0]
    n_edge,n_interior=int(row.n_edge),int(row.n_interior)
    total=n_edge+n_interior
    observed_boundary=n_edge/total
    null_boundary=observed_boundary-row.displacement/2.0
    recovered_balance=2*null_boundary-1
    assert round(recovered_balance,4)==round(float(src.expected),4),(virus,recovered_balance,src.expected)
    assert 0<=null_boundary<=1 and 0<=observed_boundary<=1
    result={'virus':virus,'label':row.label,'n_coords':int(row.n_coords),'n_clusters':int(row.n_clusters),
            'n_edge':n_edge,'n_interior':n_interior,'n_total':total,
            'observed_boundary_share':observed_boundary,'null_boundary_share':null_boundary,
            'observed_interior_share':1-observed_boundary,'null_interior_share':1-null_boundary,
            'observed_boundary_pct':100*observed_boundary,'null_boundary_pct':100*null_boundary,
            'observed_interior_pct':100*(1-observed_boundary),'null_interior_pct':100*(1-null_boundary),
            'stored_expected_balance':float(src.expected),
            'precision_note':'Null shares reconstructed from the six-decimal saved displacement and exact overlap counts; not a ratio of expected counts.'}
    composition.append(result)
    checks.append({'virus':virus,'recovered_null_balance':float(recovered_balance),
                   'stored_expected_balance':float(src.expected),'equal_after_rounding_to_4_decimals':True,
                   'balance_displacement_reconstructed':float(2*(observed_boundary-null_boundary)),
                   'stored_displacement':float(row.displacement)})
composition=pd.DataFrame(composition)
save(composition,'Figure_3_composition_plotted_values.csv')

# Record each source and display table. A source may contribute to
# several display tables.
records=[]
for local in sorted(SRC.rglob('*.csv')):
    rel=local.relative_to(SRC).as_posix()
    if rel.startswith('display_source/'):
        targets=['plotted_data/'+local.name]
        operation='Preserve preceding displayed columns exactly; C additionally receives ×100 display aliases.'
    else:
        if local.name=='Figure_3A_enrichment.csv':targets=['plotted_data/Figure_3A_plotted_values.csv']
        elif local.name=='Figure_3B.csv':targets=['plotted_data/Figure_3B_plotted_values.csv']
        elif local.name=='integrase_positive_geometry.csv':
            targets=['plotted_data/Figure_3C_plotted_values.csv','plotted_data/Figure_3_composition_plotted_values.csv']
        else:targets=['plotted_data/Figure_3C_plotted_values.csv']
        operation='Preserved result source or historical displayed-reference table; no new statistics.'
        if local.name=='integrase_positive_geometry.csv':operation='Verify saved displacement/ranges; baseline composition recovered from overlap counts and half the saved balance displacement.'
    for target in targets:
        records.append({'local_source':'source_data/'+rel,'original_source_relative_to_package':'source_data/'+rel,
                        'source_sha256':sha(local),'plotted_file':target,'plotted_sha256':sha(HERE/target),
                        'operation':operation})
pd.DataFrame(records).to_csv(HERE/'SOURCE_MAP.csv',index=False)
qa={'status':'PASS','new_statistical_tests':False,'new_distance_bands':False,
    'rows':{'A':len(a),'B':len(b),'C':len(c),'composition':len(composition)},
    'A_B_byte_identical_to_previous_display':True,'C_original_columns_exact_after_round_trip':True,
    'C_complete_original_ranges_preserved':True,'composition_reference_checks':checks,
    'composition_unit':'coordinate-interval overlap proportion; not donor or coordinate positivity',
    'composition_null':'Mean boundary share reconstructed from observed share minus half the saved balance displacement; expected-count ratios are not used.',
    'precision':'Saved displacement has six decimals; recovered null balance agrees with saved expected at four decimals.',
    'source_files':len(list(SRC.rglob('*.csv'))),'source_map_rows':len(records)}
(OUT/'preparation_QA.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
print(json.dumps(qa,indent=2))
