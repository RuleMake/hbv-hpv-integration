"""Prepare Supplemental Figure S1 catalogue counts and counting units.

The script selects and labels saved counts for display. Source tables retain
their original values; the display tables identify the counts used in the figure.
"""
from pathlib import Path
import hashlib
import json
import pandas as pd

HERE=Path(__file__).resolve().parent
(HERE/'review').mkdir(exist_ok=True)
(HERE/'plotted_data').mkdir(exist_ok=True)
SRC=HERE/'source_data'
OUT=HERE/'plotted_data'
OUT.mkdir(exist_ok=True)
COUNTS=['HPV','HBV','EBV','total']

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(name):return pd.read_csv(SRC/name,float_precision='round_trip')
def save(frame,name):frame.to_csv(OUT/name,index=False,float_format='%.17g')

legacy=read('catalogue_flow_input.csv')
dup=read('Figure_1B_duplication_reference.csv')
cohort=read('cohort_units.csv')
main_units=read('Figure_1B_units.csv')
dup_reconciliation=json.loads((SRC/'canonical_DUP_reconciliation.json').read_text(encoding='utf-8'))
checks=[]
def check(name,condition):
    assert bool(condition),name
    checks.append({'check':name,'passed':True})

# A: five catalogue stages, with the counting unit made explicit.
a0=legacy[legacy.panel=='a'].reset_index(drop=True)
check('five_catalogue_stages',len(a0)==5)
a=a0[COUNTS].astype('int64').copy()
a.insert(0,'stage_id',range(1,len(a)+1))
a.insert(1,'stage_label',['Published records','GRCh37-mapped records','Unique coordinates','Null-eligible coordinates','After peak masking'])
a.insert(2,'unit',['records','records','coordinates','coordinates','coordinates'])
check('A_totals_equal_sum_of_viral_counts',(a[['HPV','HBV','EBV']].sum(axis=1)==a.total).all())
save(a,'S1A_catalogue_flow.csv')

# B: retain all three source sets; the mapped record set is a processing
# stage shown in A rather than a separate analysis set.
b0=legacy[legacy.panel=='b'].reset_index(drop=True)
check('three_historical_sets_retained',len(b0)==3)
b=b0[COUNTS].astype('int64').copy()
b.insert(0,'set_label',['After peak masking','Null-eligible coordinates','GRCh37-mapped records'])
b['unit']=['coordinates','coordinates','records']
b['shown_in_panel_B']=[True,True,False]
b['use']=['Class atlas and peak-masked corrected-enrichment analyses',
          'Counting and spatial-reference comparisons before peak masking',
          'Mapped catalogue stage retained for provenance; already displayed in panel A']
check('B_sets_reproduce_A_stages',all(b.loc[i,COUNTS].tolist()==a.loc[j,COUNTS].tolist() for i,j in [(0,4),(1,3),(2,1)]))
check('B_totals_equal_sum_of_viral_counts',(b[['HPV','HBV','EBV']].sum(axis=1)==b.total).all())
check('B_only_two_analysis_sets_shown',b.shown_in_panel_B.sum()==2)
save(b,'S1B_analysis_inputs.csv')

# C: list paired calls, supporting mixed records and the large-duplication
# catalogue with their respective units.
selections=[
    ('paired_all','All paired SV calls','normalized_sample_events','all','Paired SV calls'),
    ('paired_excluded','Integration-excluded paired SV calls','normalized_sample_events','whole_positive_donor_exclusion','Paired SV calls'),
    ('support_all','All supporting SV records','legacy_exact','all','Interval and endpoint records'),
    ('support_excluded','Integration-excluded supporting SV records','legacy_exact','legacy_aliquot_exclusion','Interval and endpoint records'),
]
c_rows=[]
for row_id,label,dataset,exclusion,unit_label in selections:
    chosen=cohort[(cohort.dataset==dataset)&(cohort.exclusion_rule==exclusion)&(cohort.svclass=='ALL')]
    check(row_id+'_source_selection_unique',len(chosen)==1)
    r=chosen.iloc[0]
    c_rows.append({'row_id':row_id,'label':label,'donors':int(r.n_donors),'aliquots':int(r.n_aliquots),
                   'records':int(r.n_units),'record_unit':unit_label,'source_file':'cohort_units.csv',
                   'source_selection':f'dataset={dataset}; exclusion_rule={exclusion}; svclass=ALL'})
chosen=dup[dup.unit=='large_tandem_duplication']
check('large_duplication_source_selection_unique',len(chosen)==1)
r=chosen.iloc[0]
c_rows.append({'row_id':'large_duplication_catalogue','label':'Large tandem duplications',
               'donors':int(r.n_donors),'aliquots':int(r.n_aliquots),'records':int(r.n_records),
               'record_unit':'Large tandem-duplication records','source_file':'Figure_1B_duplication_reference.csv',
               'source_selection':'unit=large_tandem_duplication; complete saved catalogue including chromosome Y'})
c=pd.DataFrame(c_rows)
check('C_five_distinct_reference_representations',len(c)==5 and c.row_id.is_unique)
check('C_aliquots_at_least_donors',(c.aliquots>=c.donors).all())
reference=main_units[(main_units.dataset=='normalized_sample_events')&(main_units.exclusion_rule=='all')&(main_units.svclass=='ALL')]
check('main_figure_paired_total_unique',len(reference)==1)
rr=reference.iloc[0]; cc=c[c.row_id=='paired_all'].iloc[0]
check('C_primary_totals_match_accepted_main_figure',[cc.donors,cc.aliquots,cc.records]==[rr.n_donors,rr.n_aliquots,rr.n_units])
check('C_primary_class_counts_sum_to_all',int(main_units.loc[main_units.svclass!='ALL','n_units'].sum())==int(rr.n_units))
old_count=int(legacy.loc[(legacy.panel=='c')&(legacy.stage=='large tandem duplications'),'HPV'].item())
new_count=int(c.loc[c.row_id=='large_duplication_catalogue','aliquots'].item())
check('legacy_duplication_aliquot_value_not_used_for_display',new_count!=old_count)
catalogue_records=int(c.loc[c.row_id=='large_duplication_catalogue','records'].item())
check('large_duplication_records_match_saved_reconciliation',catalogue_records==dup_reconciliation['canonical_events'])
check('large_duplication_full_minus_strict_equals_saved_exclusion',catalogue_records-dup_reconciliation['strict_chr1_22_X_donor_point_pairs']==dup_reconciliation['canonical_excluded_from_strict_set'])
check('large_duplication_excluded_chromosome_is_Y',dup_reconciliation['excluded_chromosomes']==['chrY'])
save(c,'S1C_reference_units.csv')

mapping=[
    ('catalogue_flow_input.csv','S1A_catalogue_flow.csv','Select original panel a counts; five stages and explicit record/coordinate units.'),
    ('catalogue_flow_input.csv','S1B_analysis_inputs.csv','Retain all three original panel b sets; mapped records are not drawn as an analysis input.'),
    ('cohort_units.csv','S1C_reference_units.csv','Select four ALL-class rows using the explicit dataset/exclusion rules recorded in the display table.'),
    ('Figure_1B_duplication_reference.csv','S1C_reference_units.csv','Use donor and aliquot counts from the complete large-duplication catalogue.'),
    ('Figure_1B_units.csv','S1C_reference_units.csv','Independent consistency check against the accepted Figure 1 paired-call counts.'),
    ('canonical_DUP_reconciliation.json','S1C_reference_units.csv','Caption provenance: complete canonical large-DUP set versus chr1-22/X-only donor-pair set; excluded records are on chromosome Y.'),
]
records=[]
for source,plotted,operation in mapping:
    records.append({'local_source':'source_data/'+source,'original_source_relative_to_package':'source_data/'+source,
                    'source_sha256':sha(SRC/source),'plotted_file':'plotted_data/'+plotted,
                    'plotted_sha256':sha(OUT/plotted),'operation':operation})
pd.DataFrame(records).to_csv(HERE/'SOURCE_MAP.csv',index=False)
qa={'status':'PASS','new_statistical_tests':False,'all_displayed_counts_selected_from_sources':True,
    'rows':{'A':len(a),'B':len(b),'C':len(c)},'B_displayed_rows':int(b.shown_in_panel_B.sum()),
    'legacy_large_duplication_aliquot_value_preserved_only_in_raw_source':old_count,
    'audited_large_duplication_aliquots_displayed':new_count,
    'unit_note':'Paired-call and interval/end-point supporting representations retain their own donor/aliquot universes and exclusion rules.',
    'large_duplication_scope':'Complete saved large-duplication catalogue including chromosome Y.',
    'large_duplication_scope_provenance':{key:dup_reconciliation[key] for key in ['canonical_events','strict_chr1_22_X_donor_point_pairs','canonical_excluded_from_strict_set','excluded_chromosomes']},
    'passed':len(checks),'failed':0,'checks':checks,'reference_units':c.to_dict('records')}
(OUT/'preparation_QA.json').write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:qa[k] for k in ['status','rows','B_displayed_rows','passed','failed','reference_units']},indent=2))
