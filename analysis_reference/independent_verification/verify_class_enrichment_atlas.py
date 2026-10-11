"""Verify class-enrichment results with an independently derived counting algorithm.

Intervals in query-start space provide the event counts. The production estimator
is not imported. Run with Python -B; check summaries are written beside this script."""
from pathlib import Path
import sys,hashlib,json,time
from functools import lru_cache
HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; WS=ROOT.parents[1]
sys.path.insert(0,str(WS/'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
sys.path.insert(0,str(ROOT/'shared'))
import numpy as np
import pandas as pd
import shared_analysis_context as C
OUT=ROOT/'atlas/outputs'; checks=[]
def check(name,ok,detail=None):
    checks.append(dict(check=name,passed=bool(ok),detail=detail))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return pd.read_csv(path,float_precision='round_trip')
table=read(OUT/'atlas_all_variants_full_precision.csv')
diagnostic=read(OUT/'event_geometry_and_union_diagnostics.csv')
for r in table.itertuples():
    tag=f'{r.variant} {r.svclass} {r.virus} {r.estimator}'
    a=np.load(OUT/f'null_{r.variant}_{r.svclass}_{r.virus}.npz')
    key='block' if r.estimator.startswith('bc') else 'sign'
    null=a[key+'_null'];obs=int(a[key+'_observed'])
    hi=(np.count_nonzero(null>=obs)+1)/(len(null)+1)
    lo=(np.count_nonzero(null<=obs)+1)/(len(null)+1)
    qlo,qhi=np.percentile(null,[2.5,97.5]);mu=null.mean();sd=null.std(ddof=1)
    wanted={'observed':obs,'expected':mu,'null_sd':sd,'fold':obs/mu,
      'fold_null_interval_low':obs/qhi,'fold_null_interval_high':obs/qlo,
      'null_quantile_025':qlo,'null_quantile_975':qhi,'z':(obs-mu)/sd,
      'emp_p':hi,'emp_p_lower':lo,'emp_p_two_sided_equal_tail':min(1,2*min(lo,hi)),
      'emp_p_bonferroni_24':min(1,24*hi),
      'emp_p_two_sided_bonferroni_24':min(1,48*min(lo,hi)),
      'bonferroni_24_p_lt_threshold':hi<.05/24,
      'two_sided_bonferroni_24_p_lt_threshold':2*min(lo,hi)<.05/24}
    wrong={k:[float(v),float(getattr(r,k))] for k,v in wanted.items() if v!=getattr(r,k)}
    check(tag+' summary',not wrong,wrong)
    clu=read(OUT/f'clusters_{r.variant}_{r.svclass}_{r.virus}.csv')
    cobs=int(clu.obs.sum()) if key=='block' else int((clu.obs>clu.null_mean).sum())
    check(tag+' cluster aggregate',cobs==obs and len(null)==20000)
for virus,g in diagnostic.groupby('virus'):
    check(virus+' same atlas integer stream',g.rng_stream_sha256.nunique()==1 and g.rng_integers_calls.nunique()==1)

ctx=C.load_context(); feat=ctx['feat']; mask=ctx['mask'];nul=ctx['nul']
# Reconstruct candidate pools from feature tables, without nuller.pool/draw.
clean=feat.loc[[(str(r.chrom),int(r.bin_start)) not in mask for r in feat.itertuples()]]
fine={(str(c),int(s)):g.bin_start.to_numpy(np.int64) for (c,s),g in clean.groupby(['chrom','stratum_ref'])}
basefeat=pd.read_parquet(C.F.FEATURES)
basefeat['chrom']=basefeat.chrom.astype(str).map(lambda x:x if x.startswith('chr') else 'chr'+x)
basefeat=basefeat[(basefeat.stratum>=0)&basefeat.chrom.ne('chrY')]
base={(str(c),int(s)):g.bin_start.to_numpy(np.int64) for (c,s),g in basefeat.groupby(['chrom','stratum'])}
def pool(r):
    cand=np.unique(fine.get((r.chrom,int(r.stratum)),np.array([],dtype=np.int64))+int(r.within_bin_offset))
    cand=cand[cand!=r.pos0]
    if len(cand)<5:
        cand=np.unique(base[(r.chrom,int(r.stratum)//nul.divisor)]+int(r.within_bin_offset))
        cand=cand[cand!=r.pos0]
    return cand
selected=[]
for variant,cls,virus in [('paired_sample','TRA','HBV'),('donor_signature','DEL','HPV')]:
    started=time.time();tag=f'{variant} {cls} {virus}'
    geo,_=C.get_variant_geometry(variant,cls)
    coords=read(OUT/f'coordinates_{virus}.csv')
    gg={str(c):g for c,g in geo.groupby('chrom')}
    @lru_cache(maxsize=96)
    def starts_ends(chrom,length):
        g=gg[chrom];left=g.left_pos0.to_numpy(np.int64);right=g.right_pos0.to_numpy(np.int64)
        joined=right-left<=length
        # A query [x,x+L] hits an endpoint iff x is in [endpoint-L,endpoint].
        # Merge the two such intervals for spans <=L to count that event once.
        starts=np.concatenate((left-length,(right-length)[~joined]))
        ends=np.concatenate((np.where(joined,right,left),right[~joined]))
        return np.sort(starts),np.sort(ends)
    def counts(chrom,length,low):
        starts,ends=starts_ends(chrom,int(length))
        return np.searchsorted(starts,low,side='right')-np.searchsorted(ends,low,side='left')
    rng=np.random.default_rng(7);bn=np.zeros(20000,np.int64);sn=bn.copy();bo=so=0;rows=[];union=set()
    for cid,g in coords.groupby('cluster',sort=True):
        chrom=str(g.chrom.iloc[0]); p=g.pos0.to_numpy(np.int64)
        if chrom not in gg:continue
        j=rng.integers(0,len(p),size=20000);delta=np.empty(20000,np.int64)
        for k,r in enumerate(g.itertuples()):
            choose=j==k;n=int(choose.sum())
            if n:
                candidates=pool(r)
                delta[choose]=candidates[rng.integers(0,len(candidates),size=n)]-int(r.pos0)
        low=int(p.min())-100000;high=int(p.max())+100000
        observed=int(counts(chrom,high-low,np.array([low]))[0]);null=counts(chrom,high-low,low+delta)
        bo+=observed;bn+=null;so+=int(observed>null.mean());sn+=(null>(null.sum()-null)/19999)
        rows.append((int(cid),observed,float(null.mean())))
        ggeo=gg[chrom]
        hit=ggeo.left_pos0.between(low,high)|ggeo.right_pos0.between(low,high)
        union.update(ggeo.loc[hit,'event_id'])
    production=np.load(OUT/f'null_{variant}_{cls}_{virus}.npz')
    for key,array in [('block_null',bn),('sign_null',sn),('block_observed',bo),('sign_observed',so)]:
        check(tag+' independent full '+key,np.array_equal(array,production[key]))
    pcl=read(OUT/f'clusters_{variant}_{cls}_{virus}.csv')
    arr=np.asarray(rows)
    check(tag+' every cluster observed/null mean',np.array_equal(arr,pcl[['cluster','obs','null_mean']].to_numpy()))
    d=diagnostic[(diagnostic.variant==variant)&(diagnostic.svclass==cls)&(diagnostic.virus==virus)].iloc[0]
    check(tag+' independent global union',len(union)==d.observed_catalogue_unique_events and bo-len(union)==d.cross_cluster_repeat_incidences)
    selected.append(dict(variant=variant,svclass=cls,virus=virus,elapsed_seconds=time.time()-started,observed=bo,sign_observed=so,global_union=len(union),block_null_sha256=hashlib.sha256(bn.tobytes()).hexdigest(),sign_null_sha256=hashlib.sha256(sn.tobytes()).hexdigest()))
    print(tag,selected[-1],flush=True)
result=dict(n_checks=len(checks),n_passed=sum(c['passed'] for c in checks),all_passed=all(c['passed'] for c in checks),checks=checks,selected_full_recomputations=selected,reviewed_files={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'atlas/atlas_core.py',ROOT/'atlas/run_atlas.py',OUT/'atlas_all_variants_full_precision.csv',OUT/'artifact_manifest.csv']})
(HERE/'atlas_result_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:result[k] for k in ('n_checks','n_passed','all_passed')}),flush=True)
if not result['all_passed']:print([c for c in checks if not c['passed']]);raise SystemExit(1)
