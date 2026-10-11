"""Verify the four-class virus contrast and paired-sample HPV spatial-null draws.

The check does not import production interaction functions. Run with Python -B."""
from pathlib import Path
from functools import lru_cache
import sys,json,hashlib,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;WS=ROOT.parents[1]
sys.path.insert(0,str(WS/'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
sys.path.insert(0,str(ROOT/'shared'))
import numpy as np
import pandas as pd
import shared_analysis_context as C
OUT=ROOT/'axis/outputs';checks=[];full=[]
def check(name,ok,detail=None):checks.append(dict(check=name,passed=bool(ok),detail=detail))
def read(p):return pd.read_csv(p,float_precision='round_trip')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for variant in C.VARIANTS:
    arrays=np.load(OUT/variant/'class_interaction_null_draws.npz');tab=read(OUT/variant/'class_interaction_interaction.csv')
    counts=read(OUT/variant/'class_interaction_virus_counts.csv')
    provenance=read(OUT/variant/'class_interaction_call_provenance.csv')
    for virus,g in provenance.groupby('virus'):
        check(variant+' '+virus+' coordinates/pool-call identity',g.coordinate_sha256.nunique()==1 and g.n_pool_calls.nunique()==1 and g.n_pool_fallback.nunique()==1)
    for r in counts.itertuples():
        tag=f'{variant} {r.estimator} {r.svclass} {r.virus}'
        nl=arrays[f'{r.estimator}__{r.svclass}__{r.virus}__null_counts'];ob=float(arrays[f'{r.estimator}__{r.svclass}__{r.virus}__observed_count'])
        check(tag+' marginal counts and nonzero floor',ob==r.observed and nl.mean()==r.null_mean and ob/nl.mean()==r.fold and np.count_nonzero(nl<1)==r.n_null_below_one and r.n_null_below_one==0)
        if r.virus=='HBV':
            a=np.load(ROOT/f'atlas/outputs/null_{variant}_{r.svclass}_HBV.npz');key='block' if r.estimator=='bc' else 'sign'
            check(tag+' seed [7,0] atlas stream equality',np.array_equal(nl,a[key+'_null']))
    for r in tab.itertuples():
        tag=f'{variant} {r.estimator} {r.statistic}'
        terms=[(r.svclass,1)] if r.statistic.startswith('HPV_vs_HBV') else [('DEL',1),('DUP',1),('TRA',-1),('INV',-1)] if r.statistic=='crossed_all' else [('DUP',1),('TRA',-1),('INV',-1)]
        obs=0.;null=np.zeros(20000)
        for cls,sign in terms:
            h=arrays[f'{r.estimator}__{cls}__HPV__null_counts'];b=arrays[f'{r.estimator}__{cls}__HBV__null_counts']
            oh=float(arrays[f'{r.estimator}__{cls}__HPV__observed_count']);ob=float(arrays[f'{r.estimator}__{cls}__HBV__observed_count'])
            obs+=sign*(np.log(oh/h.mean())-np.log(ob/b.mean()))
            null+=sign*(np.log(np.maximum(h,1)/h.mean())-np.log(np.maximum(b,1)/b.mean()))
        check(tag+' complete log contrast array',np.array_equal(null,arrays[f'{r.estimator}__{r.statistic}__log_null']))
        mu=null.mean();sd=null.std(ddof=1);qlo,qhi=np.percentile(null,[2.5,97.5]);hi=(np.count_nonzero(null>=obs)+1)/20001;lo=(np.count_nonzero(null<=obs)+1)/20001
        vals={'log_ratio':obs,'ratio_of_folds':np.exp(obs),'null_mean':mu,'null_sd':sd,'null_reference_q025':qlo,'null_reference_q975':qhi,'z_legacy_zero_centered':obs/sd,'z_null_mean_centered_diagnostic':(obs-mu)/sd,'emp_p_greater':hi,'emp_p_less':lo,'emp_p_two_sided':min(1,2*min(hi,lo)),'mde_log_ratio_approx_2p80sd':2.8*sd}
        wrong={k:[float(v),float(getattr(r,k))] for k,v in vals.items() if v!=getattr(r,k)}
        check(tag+' all statistics',not wrong,wrong)
ctx=C.load_context();feat=ctx['feat'];mask=ctx['mask'];joint=ctx['nul']
clean=feat.loc[[(r.chrom,int(r.bin_start)) not in mask for r in feat.itertuples()]]
fine={(str(c),int(s)):g.bin_start.to_numpy(np.int64) for (c,s),g in clean.groupby(['chrom','stratum_ref'])}
bf=pd.read_parquet(C.F.FEATURES);bf['chrom']=bf.chrom.astype(str).map(lambda x:x if x.startswith('chr') else 'chr'+x)
bf=bf[(bf.stratum>=0)&bf.chrom.ne('chrY')]
base={(str(c),int(s)):g.bin_start.to_numpy(np.int64) for (c,s),g in bf.groupby(['chrom','stratum'])}
def pool(r):
    cand=np.unique(fine.get((r.chrom,int(r.stratum)),np.array([],np.int64))+int(r.within_bin_offset));cand=cand[cand!=r.pos0]
    if len(cand)<5:
        cand=np.unique(base[(r.chrom,int(r.stratum)//joint.divisor)]+int(r.within_bin_offset));cand=cand[cand!=r.pos0]
    return cand
coords=read(ROOT/'atlas/outputs/coordinates_HPV.csv');saved=np.load(OUT/'paired_sample/class_interaction_null_draws.npz')
for cls in C.CLASSES:
    started=time.time();geo,_=C.get_variant_geometry('paired_sample',cls);gg={str(c):g for c,g in geo.groupby('chrom')}
    @lru_cache(maxsize=96)
    def bounds(chrom,length):
        g=gg[chrom];left=g.left_pos0.to_numpy(np.int64);right=g.right_pos0.to_numpy(np.int64);joined=right-left<=length
        return np.sort(np.r_[left-length,(right-length)[~joined]]),np.sort(np.r_[np.where(joined,right,left),right[~joined]])
    def count(c,length,low):
        start,end=bounds(c,length)
        return np.searchsorted(start,low,side='right')-np.searchsorted(end,low,side='left')
    rng=np.random.default_rng([7,1]);bn=np.zeros(20000,np.int64);sn=bn.copy();bo=so=0;stream=hashlib.sha256()
    for cid,g in coords.groupby('cluster',sort=True):
        chrom=g.chrom.iloc[0];p=g.pos0.to_numpy(np.int64)
        if chrom not in gg:continue
        j=rng.integers(0,len(p),size=20000);delta=np.empty(20000,np.int64)
        for k,r in enumerate(g.itertuples()):
            choose=j==k;n=int(choose.sum())
            if n:
                cand=pool(r);delta[choose]=cand[rng.integers(0,len(cand),size=n)]-int(r.pos0)
        stream.update(delta.tobytes())
        lo=int(p.min())-100000;hi=int(p.max())+100000
        o=int(count(chrom,hi-lo,np.array([lo]))[0]);null=count(chrom,hi-lo,lo+delta)
        bo+=o;bn+=null;so+=int(o>null.mean());sn+=null>(null.sum()-null)/19999
    for est,nl,o in [('bc',bn,bo),('sign',sn,so)]:
        check(f'paired_sample {cls} HPV independent full {est}',np.array_equal(nl,saved[f'{est}__{cls}__HPV__null_counts']) and o==saved[f'{est}__{cls}__HPV__observed_count'])
    full.append(dict(svclass=cls,variant='paired_sample',virus='HPV',n_perm=20000,observed=bo,sign_observed=so,elapsed_seconds=time.time()-started,actual_cluster_displacement_stream_sha256=stream.hexdigest()))
    print(full[-1],flush=True)
check('paired sample all four HPV classes identical actual displacement stream',len({r['actual_cluster_displacement_stream_sha256'] for r in full})==1)
result=dict(n_checks=len(checks),n_passed=sum(c['passed'] for c in checks),all_passed=all(c['passed'] for c in checks),checks=checks,full_recomputations=full,reviewed_files={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'axis/interaction_recompute.py',ROOT/'axis/PROTOCOL.md']})
(HERE/'class_interaction_result_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print({k:result[k] for k in ('n_checks','n_passed','all_passed')},flush=True)
if not result['all_passed']:print([c for c in checks if not c['passed']]);raise SystemExit(1)
