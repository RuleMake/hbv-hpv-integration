"""Verify donor proportions and their statistical summaries from saved inputs.

The calculation is independent of the production functions and writes its own
check summary beside this script."""
from pathlib import Path
import sys, hashlib, json, bisect
from collections import Counter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / 'donor_duplication_statistics'
DATA, OUT = SRC / 'data', SRC / 'outputs'
checks, facts = [], {}
def check(name, ok, detail=None):
    checks.append(dict(check=name, passed=bool(ok), detail=detail))
def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
manifest = pd.read_csv(SRC/'input_manifest.csv')
for r in manifest.itertuples():
    p=SRC/r.relative_path
    check('input '+r.relative_path, p.stat().st_size==r.bytes and digest(p)==r.sha256)
run=json.loads((OUT/'statistical_run.json').read_text())
for key,p in [('protocol_sha256',SRC/'retrospective_correction_protocol.md'),('script_sha256',SRC/'correct_donor_statistics.py'),('input_manifest_sha256',SRC/'input_manifest.csv'),('donor_scores_sha256',OUT/'donor_scores.csv')]:
    check('run '+key, digest(p)==run[key])

ev=pd.read_csv(DATA/'canonical_dup.tsv',sep='\t')
co=pd.read_csv(OUT/'coordinate_eligibility.csv')
ds=pd.read_csv(OUT/'donor_scores.csv',float_precision='round_trip')
es=pd.read_csv(OUT/'event_scores.csv',float_precision='round_trip')
ss=pd.read_csv(OUT/'corrected_statistics.csv',float_precision='round_trip').set_index('analysis')
h=pd.read_excel(DATA/'histology.xlsx')
integ=pd.read_csv(DATA/'integration_cohort.csv')
check('canonical global ID unique',ev.unique_td_id.is_unique)
check('canonical donor coordinate orientation key unique',not ev.duplicated(['donor','chrom','left_pos0','right_pos0','left_strand','right_strand']).any())
check('no ID shared across donors',ev.groupby('unique_td_id').donor.nunique().max()==1)
check('canonical >100kb and DUP orientation',((ev.right_pos0-ev.left_pos0>100000)&ev.left_strand.eq('-')&ev.right_strand.eq('+')).all())
facts['canonical']={'rows':len(ev),'donors':ev.donor.nunique(),'missing_donors':int(ev.donor.isna().sum()),'min_span':int((ev.right_pos0-ev.left_pos0).min()),'nonfinite_coordinates':int((~np.isfinite(ev[['left_pos0','right_pos0']])).sum().sum())}
check('no missing canonical event fields',ev[['unique_td_id','donor','chrom','left_pos0','right_pos0']].notna().all().all())
check('unique scored donor per virus',not ds.duplicated(['donor','virus_coordinate_set']).any())
check('unique scored event per virus',not es.duplicated(['unique_td_id','virus_coordinate_set']).any())

labels=h.groupby('icgc_donor_id').histology_abbreviation.agg(lambda z:sorted(set(z.dropna())))
resolved={d:v[0] for d,v in labels.items() if len(v)==1}
conflicts={d:v for d,v in labels.items() if len(v)>1}
facts['histology_conflicts']=conflicts
check('nine histology conflicts',len(conflicts)==9)
check('all conflicting donors excluded',not ds.loc[ds.donor.isin(conflicts),'included_corrected'].any())

histsets={'HBV':{'Liver-HCC'},'HPV':{'Cervix-SCC','Cervix-AdenoCA','Head-SCC'}}
arrays={}
for virus in ['HBV','HPV']:
    use=co.loc[co.joint_eligible & co.virus_class.eq(virus)]
    refs={str(c):sorted(set(g.pos0.astype(int))) for c,g in use.groupby('chrom')}
    truth=[]
    for r in ev.itertuples():
        chrom=str(r.chrom)
        if not chrom.startswith('chr'):chrom='chr'+chrom
        ref=refs.get(chrom,[])
        dist=float('inf')
        for p in (int(r.left_pos0),int(r.right_pos0)):
            k=bisect.bisect_left(ref,p)
            if k<len(ref):dist=min(dist,abs(p-ref[k]))
            if k:dist=min(dist,abs(p-ref[k-1]))
        truth.append(dist)
    manual=ev[['unique_td_id','donor']].copy()
    manual['distance']=truth
    manual['hit']=manual.distance.le(100000)
    exported=es.loc[es.virus_coordinate_set.eq(virus)].set_index('unique_td_id').loc[manual.unique_td_id]
    check(virus+' all event endpoint distances',np.allclose(exported.nearest_endpoint_distance_bp,truth,rtol=0,atol=0))
    check(virus+' all event OR hits',np.array_equal(exported.near,manual.hit))
    check(virus+' event donors intact',np.array_equal(exported.donor,manual.donor))
    agg=manual.groupby('donor').hit.agg(['sum','count'])
    agg['r']=agg['sum']/agg['count']
    written=ds.loc[ds.virus_coordinate_set.eq(virus)].set_index('donor')
    joined=written.join(agg)
    check(virus+' donor numerators exact',np.array_equal(joined.near_DUP.fillna(-1),joined['sum'].fillna(-1)))
    check(virus+' donor denominators exact',np.array_equal(joined.total_DUP.fillna(-1),joined['count'].fillna(-1)))
    check(virus+' donor fractions exact to export precision',np.allclose(joined.fraction,joined.r,rtol=0,atol=5e-12,equal_nan=True))
    check(virus+' missing denominator is not zero fraction',joined.loc[joined['count'].isna(),'fraction'].isna().all())
    check(virus+' all true zeros retained',joined.loc[joined['sum'].eq(0),'true_zero_fraction'].all())
    positive=set(integ.loc[integ.virus_class.eq(virus),'donor'])
    selected=agg.loc[[d for d in agg.index if resolved.get(d) in histsets[virus]]].copy()
    selected['hist']=[resolved[d] for d in selected.index]
    selected['positive']=selected.index.isin(positive)
    check(virus+' corrected inclusion reconstructed',set(selected.index)==set(written.index[written.included_corrected]))
    check(virus+' integration grouping reconstructed',np.array_equal(written.integration_detected_for_virus,written.index.isin(positive)))
    a=selected.loc[selected.positive,'r'].to_numpy();b=selected.loc[~selected.positive,'r'].to_numpy()
    arrays[virus+'_pooled_donor_mean']=[(1.,a,b)]
    facts[virus]={'coordinates':len(use),'unique_coordinate_positions':len(use[['chrom','pos0']].drop_duplicates()),'n_detected':len(a),'n_comparison':len(b),'detected_near_events':int(selected.loc[selected.positive,'sum'].sum()),'detected_total_events':int(selected.loc[selected.positive,'count'].sum()),'comparison_near_events':int(selected.loc[~selected.positive,'sum'].sum()),'comparison_total_events':int(selected.loc[~selected.positive,'count'].sum()),'n_scored_roster':len(agg),'n_roster_no_denominator':int(joined['count'].isna().sum()),'detected_any_event_fraction':float((a>0).mean()),'comparison_any_event_fraction':float((b>0).mean())}
    if virus=='HPV':
        support=[];good=[]
        for hist in ['Cervix-SCC','Cervix-AdenoCA','Head-SCC']:
            x=selected[selected['hist'].eq(hist)]
            aa=x.loc[x.positive,'r'].to_numpy();bb=x.loc[~x.positive,'r'].to_numpy()
            support.append(dict(histology=hist,n_detected=len(aa),n_comparison=len(bb)))
            if min(len(aa),len(bb))>=2:good.append((aa,bb))
        total=sum(len(a)+len(b) for a,b in good)
        arrays['HPV_common_support_histology_standardized']=[((len(a)+len(b))/total,a,b) for a,b in good]
        facts['HPV_support']=support

# Independent sum/sum-of-squares statistic, not producer functions.
def effect(cells):
    ma=mb=var=0.
    for w,a,b in cells:
        ma+=w*np.sum(a)/a.size;mb+=w*np.sum(b)/b.size
        va=(np.dot(a,a)-np.sum(a)**2/a.size)/(a.size-1)
        vb=(np.dot(b,b)-np.sum(b)**2/b.size)/(b.size-1)
        var+=w*w*(va/a.size+vb/b.size)
    d=ma-mb
    return np.array([ma,mb,d,ma/mb,d/np.sqrt(var)])
mc={}
for i,(name,cells) in enumerate(arrays.items()):
    row=ss.loc[name];point=effect(cells)
    check(name+' n',sum(len(a) for w,a,b in cells)==row.n_detected and sum(len(b) for w,a,b in cells)==row.n_comparison)
    check(name+' means difference ratio Welch independently reconstructed',np.allclose(point,row[['mean_detected','mean_comparison','mean_difference','mean_ratio','welch_statistic']].to_numpy(float),rtol=0,atol=5e-11),point.tolist())
    # Same documented permutation seed; independently recalculate moments for every allocation.
    rng=np.random.default_rng(np.random.SeedSequence([20260908,i,2]));ct=cd=0
    for iteration in range(99999):
        perm=[]
        for w,a,b in cells:
            q=rng.permutation(np.r_[a,b]);perm.append((w,q[:len(a)],q[len(a):]))
        stat=effect(perm)
        ct+=abs(stat[4])>=abs(point[4])-1e-14
        cd+=abs(stat[2])>=abs(point[2])-1e-14
    pp=np.array([(ct+1)/100000,(cd+1)/100000])
    reported=row[['p_two_sided_studentized_permutation','p_two_sided_raw_difference_permutation']].to_numpy(float)
    check(name+' both permutation P within one Monte Carlo count',np.max(abs(pp-reported))<=0.00001000001,{'exact_count_fraction_P':pp.tolist(),'reported_from_serialized_fractions':reported.tolist()})
    if name=='HPV_common_support_histology_standardized':
        # Recheck the boundary-sensitive arm against exact serialized input.
        # Results use 17 significant digits and round-trip parsing to preserve precision.
        hp=ds.loc[ds.included_corrected & ds.virus_coordinate_set.eq('HPV')]
        cell_csv=[]
        for hist in ['Cervix-SCC','Head-SCC']:
            g=hp[hp.histology_corrected.eq(hist)]
            aa=g.loc[g.integration_detected_for_virus,'fraction'].to_numpy()
            bb=g.loc[~g.integration_detected_for_virus,'fraction'].to_numpy()
            cell_csv.append(((len(aa)+len(bb))/58,aa,bb))
        obs_csv=effect(cell_csv)
        gen=np.random.default_rng(np.random.SeedSequence([20260908,i,2]));nt=nd=0
        for _ in range(99999):
            perm=[]
            for w,a,b in cell_csv:
                q=gen.permutation(np.r_[a,b]);perm.append((w,q[:len(a)],q[len(a):]))
            val=effect(perm)
            nt+=abs(val[4])>=abs(obs_csv[4])-1e-14
            nd+=abs(val[2])>=abs(obs_csv[2])-1e-14
        serialized_p=np.array([(nt+1)/100000,(nd+1)/100000])
        check(name+' exact serialized-input permutation P',np.array_equal(serialized_p,reported),serialized_p.tolist())
    # Independently scheduled donor bootstrap uses a different seed. A replication check,
    # not a replacement estimate or additional inferential test.
    rng=np.random.default_rng(85400+i);nb=25000;ma=np.zeros(nb);mb=np.zeros(nb)
    for w,a,b in cells:
        ma+=w*np.mean(rng.choice(a,(nb,len(a)),replace=True),axis=1)
        mb+=w*np.mean(rng.choice(b,(nb,len(b)),replace=True),axis=1)
    ci=np.r_[np.quantile(ma-mb,[.025,.975]),np.quantile(ma/mb,[.025,.975])]
    old=row[['diff_ci_lo','diff_ci_hi','ratio_ci_lo','ratio_ci_hi']].to_numpy(float)
    mc[name]={'seed':85400+i,'replicates':nb,'difference_ci':ci[:2].tolist(),'ratio_ci':ci[2:].tolist(),'max_difference_from_producer':float(np.max(abs(ci-old))),'zero_comparison_means':int((mb==0).sum())}
    print(name, 'point and permutations checked', flush=True)

check('Holm HBV',abs(ss.loc['HBV_pooled_donor_mean','p_holm_two_pooled_tests']-2*ss.loc['HBV_pooled_donor_mean','p_two_sided_studentized_permutation'])<1e-12)
check('Holm HPV',ss.loc['HPV_pooled_donor_mean','p_holm_two_pooled_tests']==ss.loc['HPV_pooled_donor_mean','p_two_sided_studentized_permutation'])
facts['independent_bootstrap_sanity']=mc
result={'scope':'DONOR_STATISTICS independent output and statistical audit; no producer functions imported','passed':all(c['passed'] for c in checks),'n_checks':len(checks),'checks':checks,'facts':facts,'inspected_output_hashes':{p.name:digest(p) for p in OUT.iterdir() if p.is_file()}}
(HERE/'donor_statistics_independent_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'passed':result['passed'],'n_checks':len(checks),'failures':[c for c in checks if not c['passed']],'facts':facts},ensure_ascii=False,indent=2))
