"""Reconstruct the structural-variant catalogue independently from BEDPE files and metadata.

Source inputs are read only; the script writes its own JSON check summary."""
from pathlib import Path
import hashlib, json, csv, gzip, sys
from collections import defaultdict, Counter
HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'event_ledger';DATA=SRC/'data';OUT=SRC/'outputs'
sys.path.insert(0,str(SRC/'runtime_deps'))
import pandas as pd
import numpy as np
checks=[];facts={}
def check(name,ok,detail=None):checks.append(dict(check=name,passed=bool(ok),detail=detail))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def table(p,sep=','):return pd.read_csv(p,sep=sep,keep_default_na=False,dtype=str)
manifest=table(DATA/'input_manifest.csv')
bad=[]
for r in manifest.itertuples():
    p=SRC/r.copy_path
    if p.stat().st_size!=int(r.bytes) or sha(p)!=r.sha256:bad.append(r.copy_path)
check('all independent input hashes',not bad,{'n_files':len(manifest),'bad':bad})
canon=table(DATA/'input_pooled_canonical_events.tsv','\t')
sheet=table(DATA/'pcawg_sample_sheet.tsv','\t')
cohort=table(DATA/'cohort.csv')
inv=table(DATA/'raw_source_inventory.csv')
atoms=defaultdict(set)
for row in canon.itertuples():
    for a in row.aliquot_metadata.split('|'):
        if a.strip():atoms[a.strip()].add(row.donor)
official=defaultdict(set)
for r in sheet.itertuples():
    if r.aliquot_id and r.icgc_donor_id:official[r.aliquot_id].add(r.icgc_donor_id)
check('atomic mappings unambiguous and agree',all(len(v)==1 and (not official.get(a) or official[a]==v) for a,v in atoms.items()))
donor={a:next(iter(v)) for a,v in atoms.items()}
strings=set(canon.aliquot_metadata)
chosen=inv[inv.aliquot.isin(atoms)].drop_duplicates('aliquot',keep='last')
check('all atomic aliquots have chosen source',set(chosen.aliquot)==set(atoms))
facts['selection']={'canonical_records':len(canon),'metadata_strings':len(strings),'atomic_aliquots':len(atoms),'donors':len(set(donor.values())),'unexpanded_aliquots':len(set(chosen.aliquot)&strings),'omitted_aliquots':sorted(set(chosen.aliquot)-strings)}

fields=['sv_id','svclass','chrom1','start1','end1','strand1','chrom2','start2','end2','strand2']
allrows=[];contents={};duplicates=0;conflicts=[]
for r in chosen.itertuples():
    with gzip.open(SRC/r.staged_relative,'rt',encoding='utf-8') as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        for lineno,row in enumerate(reader,1):
            key=(r.aliquot,row['sv_id'])
            content=tuple(sorted(row.items()))
            if key in contents:
                duplicates+=1
                if contents[key]!=content:conflicts.append(key)
            contents[key]=content
            allrows.append([r.aliquot,donor[r.aliquot],lineno]+[row[f] for f in fields])
raw=pd.DataFrame(allrows,columns=['aliquot','donor','source_data_row']+fields)
del allrows,contents
for c in ['start1','end1','start2','end2']:raw[c]=raw[c].astype('int64')
raw['event_id']=raw.aliquot+'::'+raw.sv_id
check('no duplicate or conflicting source IDs',duplicates==0 and not conflicts,{'duplicate_rows':duplicates,'conflicts':conflicts})
check('all source IDs nonblank',raw.sv_id.str.strip().ne('').all())
check('all source donors nonblank',raw.donor.str.strip().ne('').all())
facts['raw_rows']=len(raw)
saved=pd.read_parquet(OUT/'source_rows_all.parquet',columns=list(raw.columns)).set_index('event_id').sort_index()
manual=raw.set_index('event_id').sort_index()
check('every raw paired identity and metadata field restored',saved.astype(str).equals(manual.astype(str)))

allowed={str(i) for i in range(1,23)}|{'X'}
cn1=raw.chrom1.str.removeprefix('chr');cn2=raw.chrom2.str.removeprefix('chr')
cls=raw.svclass.str.upper().replace({'H2HINV':'INV','T2TINV':'INV'})
size=(raw.start2-raw.start1).abs()
allowedpair=cn1.isin(allowed)&cn2.isin(allowed)
validintervals=(raw.start1>=0)&(raw.start2>=0)&(raw.end1>raw.start1)&(raw.end2>raw.start2)
proper=((cls=='TRA')&(cn1!=cn2)) | (cls.isin(['DEL','DUP','INV'])&(cn1==cn2)&(size>100000)&(size<=7500000))
raw['analysis_class']=cls
raw['eligible']=allowedpair&validintervals&proper
eligible=raw[raw.eligible].copy()
stored=pd.read_parquet(OUT/'normalized_events_all.parquet',columns=['event_id','analysis_eligible','analysis_class','integration_positive_donor','legacy_positive_aliquot']).set_index('event_id')
check('all raw eligibility flags and class labels reconstructed',np.array_equal(stored.loc[raw.event_id,'analysis_eligible'],raw.eligible)&np.array_equal(stored.loc[raw.event_id,'analysis_class'],raw.analysis_class))
analysis=pd.read_parquet(OUT/'normalized_analysis_sample_events.parquet',columns=['event_id'])
check('analysis exact event-ID set',set(analysis.event_id)==set(eligible.event_id))
facts['normalized_counts']=eligible.analysis_class.value_counts().sort_index().to_dict()
facts['normalized_total']=len(eligible)
rawdup=raw[cls.eq('DUP') & cn1.eq(cn2) & size.gt(100000) & size.le(7500000)]
rawkeys=set();strictkeys=set()
for r in rawdup.itertuples():
    key=(r.donor,*sorted([(r.chrom1.removeprefix('chr'),r.start1,r.strand1),(r.chrom2.removeprefix('chr'),r.start2,r.strand2)]))
    rawkeys.add(key)
    if r.eligible:strictkeys.add(key)
canonicalkeys={(r.donor,*sorted([(r.chrom.removeprefix('chr'),int(r.left_pos0),r.left_strand),(r.chrom.removeprefix('chr'),int(r.right_pos0),r.right_strand)])) for r in canon.itertuples()}
check('canonical and raw donor point-pair sets identical',rawkeys==canonicalkeys)
check('canonical strict exclusion eight chrY pairs',len(canonicalkeys-strictkeys)==8 and all(k[1][0]=='Y' and k[2][0]=='Y' for k in canonicalkeys-strictkeys))
check('raw BEDPE intervals all one bp',raw.end1.sub(raw.start1).eq(1).all() and raw.end2.sub(raw.start2).eq(1).all())
facts['canonical_reconciliation']={'raw_large_DUP_sample_calls':len(rawdup),'raw_unique_donor_point_pairs':len(rawkeys),'canonical_point_pairs':len(canonicalkeys),'strict_point_pairs':len(strictkeys),'excluded_chrY_pairs':len(canonicalkeys-strictkeys)}
facts['chromosome_restriction']={'one_allowed_TRA':int(((cls=='TRA') & (cn1.isin(allowed)^cn2.isin(allowed))).sum()),'invalid_intervals':int((~validintervals).sum()),'min_valid_interval_span':int(size[raw.eligible & cls.ne('TRA')].min()),'max_valid_interval_span':int(size[raw.eligible & cls.ne('TRA')].max()),'equal_100000':int((size.eq(100000)&cn1.eq(cn2)).sum()),'equal_7500000':int((size.eq(7500000)&cn1.eq(cn2)).sum())}

# Verify both raw endpoints of every normalized call, including excluded calls.
ends=pd.read_parquet(OUT/'normalized_breakends_all.parquet',columns=['event_id','side','breakend_id','chrom','start','end','strand','aliquot','donor'])
check('exactly two separately identified ends per raw call',len(ends)==2*len(raw) and ends.breakend_id.is_unique and ends.groupby('event_id').size().eq(2).all())
for side in [1,2]:
    p=ends[ends.side.eq(side)].set_index('event_id').loc[raw.event_id]
    good=True
    for key in ['chrom','start','end','strand']:
        good=good and np.array_equal(p[key].astype(str),raw[f'{key}{side}'].astype(str))
    good=good and np.array_equal(p.donor,raw.donor) and np.array_equal(p.aliquot,raw.aliquot)
    check('all endpoint '+str(side)+' fields preserved',good)
eligible_ends=pd.read_parquet(OUT/'normalized_analysis_breakends.parquet',columns=['event_id','side'])
check('analysis breakends exact 2x paired calls',len(eligible_ends)==2*len(eligible) and set(eligible_ends.event_id)==set(eligible.event_id))

# Reconstruct the earlier mixed representation from direct tuples independently of production keys.
legacy_all=[]
for r in raw.itertuples():
    c1=r.chrom1.removeprefix('chr');c2=r.chrom2.removeprefix('chr')
    typ=r.analysis_class
    out=[]
    if typ=='TRA' and c1!=c2:
        out=[(c1,r.start1,r.start1),(c2,r.start2,r.start2)]
    elif typ in ['DEL','DUP','INV'] and c1==c2 and 100000<abs(r.start2-r.start1)<=7500000:
        out=[(c1,min(r.start1,r.start2),max(r.start1,r.start2))]
    for c,lo,hi in out:
        if c in allowed:legacy_all.append((r.aliquot,typ,'chr'+c,lo,hi))
expanded=set(legacy_all);old={r for r in expanded if r[0] in strings}
def keyset(path):
    df=pd.read_parquet(path,columns=['aliquot','svclass','chrom','left_pos0','right_pos0'])
    return set(df.itertuples(index=False,name=None))
check('legacy every tuple exact against archived cache',old==keyset(DATA/'legacy_sv_class_ledger.parquet'))
check('expanded every tuple exact',expanded==keyset(OUT/'legacy_expanded_selection.parquet'))
facts['legacy']={'original_records':len(old),'expanded_records':len(expanded),'added_records':len(expanded-old),'expanded_precollapse_rows':len(legacy_all),'collapsed_rows':len(legacy_all)-len(expanded)}

# Exact donor signatures retain raw SV subclass and attached strand with each endpoint.
sigcounts=Counter();sigaliquots=defaultdict(set)
for r in eligible.itertuples():
    e1=(r.chrom1.removeprefix('chr'),r.start1,r.end1,r.strand1)
    e2=(r.chrom2.removeprefix('chr'),r.start2,r.end2,r.strand2)
    key=(r.donor,r.svclass.upper(),*sorted([e1,e2]))
    sigcounts[key]+=1;sigaliquots[key].add(r.aliquot)
facts['donor_signatures']={'n_exact_groups':len(sigcounts),'excess_calls':len(eligible)-len(sigcounts),'multi_aliquot_groups':sum(len(v)>1 for v in sigaliquots.values()),'repeated_groups':sum(n>1 for n in sigcounts.values())}
sg=pd.read_parquet(OUT/'donor_event_exactsignature_groups.parquet')
check('donor exact signature aggregate counts',len(sg)==len(sigcounts) and int(sg.n_sample_events.sum())==len(eligible) and int(sg.n_aliquots.gt(1).sum())==facts['donor_signatures']['multi_aliquot_groups'])

posdon=set(cohort.donor);posaliq={p.replace('\\','/').rsplit('/',1)[-1].split('.')[0] for p in cohort.sv_path if p}
check('positive-donor and aliquot masks all calls',np.array_equal(stored.loc[raw.event_id,'integration_positive_donor'],raw.donor.isin(posdon)) and np.array_equal(stored.loc[raw.event_id,'legacy_positive_aliquot'],raw.aliquot.isin(posaliq)))
extra=eligible[eligible.donor.isin(posdon)&~eligible.aliquot.isin(posaliq)]
facts['positive_exclusion']={'cohort_rows':len(cohort),'positive_donors':len(posdon),'positive_path_aliquots':len(posaliq),'extra_positive_aliquot_calls':len(extra),'eligible_negative_calls':int((~eligible.donor.isin(posdon)).sum()),'eligible_negative_aliquots':eligible.loc[~eligible.donor.isin(posdon),'aliquot'].nunique(),'eligible_negative_donors':eligible.loc[~eligible.donor.isin(posdon),'donor'].nunique()}

summary=pd.read_csv(OUT/'counts_by_unit_class_and_exclusion.csv')
for rule in ['all','legacy_aliquot_exclusion','whole_positive_donor_exclusion']:
    filt=eligible if rule=='all' else eligible[~eligible.aliquot.isin(posaliq)] if rule=='legacy_aliquot_exclusion' else eligible[~eligible.donor.isin(posdon)]
    for typ in ['ALL','DEL','DUP','INV','TRA']:
        sub=filt if typ=='ALL' else filt[filt.analysis_class.eq(typ)]
        expected=summary[(summary.dataset=='normalized_sample_events')&(summary.exclusion_rule==rule)&(summary.svclass==typ)].iloc[0]
        check('published count '+rule+' '+typ, len(sub)==expected.n_units and 2*len(sub)==expected.n_breakend_observations and sub.aliquot.nunique()==expected.n_aliquots and sub.donor.nunique()==expected.n_donors)

result={'scope':'Independent paired-call and denominator audit; no class-axis computed','passed':all(r['passed'] for r in checks),'n_checks':len(checks),'checks':checks,'facts':facts,'producer_code_hashes':{p.name:sha(p) for p in [SRC/'ledger_lib.py',SRC/'run_audit.py',SRC/'protocol.md']},'review_input_output_hashes':{p.name:sha(p) for p in [OUT/'audit_summary.json',OUT/'normalized_analysis_sample_events.parquet',OUT/'counts_by_unit_class_and_exclusion.csv']}}
(HERE/'ledger_independent_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'passed':result['passed'],'n_checks':len(checks),'failures':[r for r in checks if not r['passed']],'facts':facts},ensure_ascii=False,indent=2))
