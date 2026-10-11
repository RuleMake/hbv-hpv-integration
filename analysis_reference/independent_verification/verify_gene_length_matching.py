"""Verify gene-length-matched estimates using direct endpoint counts.

The check does not import the production matching or estimation functions.
The --baseline-only option checks reference inputs; the default also requires
completed sensitivity-analysis outputs."""
from pathlib import Path
import sys,json,hashlib,argparse,time,gzip,csv,re
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).resolve().parent;WS=ROOT.parents[1]
sys.path.insert(0,str(WS/'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
import numpy as np
import pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--baseline-only',action='store_true');args=ap.parse_args()
checks=[];began=time.time()
def ck(name,condition,detail=None):
    checks.append(dict(name=name,passed=bool(condition),**({'detail':detail} if detail is not None else {})))
    if not condition:print('FAIL',name,detail,flush=True)
def near(a,b):return np.allclose(a,b,atol=2e-12,rtol=2e-12,equal_nan=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

class EndpointUnion:
    """Direct per-record OR of endpoint/window predicates, batched 128 windows.

    Does not reuse the production sorted-endpoint inclusion-exclusion algorithm.
    One row per representation ID makes the Boolean row sum a distinct-call count.
    """
    def __init__(self,df):
        if df.event_id.duplicated().any():raise ValueError('Expected one row per intrachromosomal DUP representation ID')
        if not (df.right_pos0>=df.left_pos0).all():raise ValueError('Unordered DUP pair')
        self.pairs={c:(g.left_pos0.to_numpy(np.int64),g.right_pos0.to_numpy(np.int64)) for c,g in df.groupby('chrom')}
    def count(self,chrom,starts,L):
        q=np.asarray(starts,dtype=np.int64);left,right=self.pairs[chrom];result=np.empty(q.shape,np.int64)
        for k in range(0,len(q),128):
            lo=q[k:k+128][None,:];hi=lo+L
            hits=((left[:,None]>=lo)&(left[:,None]<=hi))|((right[:,None]>=lo)&(right[:,None]<=hi))
            result[k:k+128]=hits.sum(axis=0,dtype=np.int64)
        return result
    def event_ids(self,chrom,lo,hi,df):
        g=df[df.chrom.eq(chrom)]
        return set(g.loc[g.left_pos0.between(lo,hi)|g.right_pos0.between(lo,hi),'event_id'])

def stats(o,n):
    n=np.asarray(n);mu=float(n.mean());sd=float(n.std(ddof=1));q=np.percentile(n,[2.5,97.5]);u=(np.count_nonzero(n>=o)+1)/(len(n)+1);l=(np.count_nonzero(n<=o)+1)/(len(n)+1)
    return dict(observed=float(o),expected=mu,null_sd=sd,fold=float(o/mu),null_quantile_025=float(q[0]),null_quantile_975=float(q[1]),
        fold_null_low=float(o/q[1]),fold_null_high=float(o/q[0]),p_upper=float(u),p_lower=float(l),p_two_sided=float(min(1,2*min(u,l))),
        z=(o-mu)/sd,n_permutations=len(n))

# Meaningful oracle counterexamples, including both ends, boundary equality and crossing-only.
syn=pd.DataFrame(dict(chrom=['chr1']*5,left_pos0=[0,10,19,-10,2],right_pos0=[20,11,30,40,18],event_id=list('abcde')))
oracle=EndpointUnion(syn)
for L in [0,1,10,20]:
    x=np.arange(-15,46)
    ref=np.array([int(((syn.left_pos0.between(int(z),int(z+L)))|(syn.right_pos0.between(int(z),int(z+L)))).sum()) for z in x])
    ck('independent_endpoint_direct_union_L'+str(L),np.array_equal(oracle.count('chr1',x,L),ref))

variants=['canonical_legacy','class_legacy','donor_signature','paired_sample']
geoms={v:pd.read_parquet(ROOT/f'baseline/outputs/geometry_{v}.parquet') for v in variants}
indexes={v:EndpointUnion(g) for v,g in geoms.items()}
baseline=pd.read_csv(ROOT/'baseline/outputs/all_baselines_full_precision.csv',float_precision='round_trip')
fields={'observed':'observed','expected':'expected','null_sd':'null_sd','fold':'fold','null_quantile_025':'null_quantile_025','null_quantile_975':'null_quantile_975',
        'fold_null_low':'fold_null_interval_low','fold_null_high':'fold_null_interval_high','p_upper':'emp_p','p_lower':'emp_p_lower','p_two_sided':'emp_p_two_sided_equal_tail','z':'z','n_permutations':'n_perm'}
null_by_variant={v:np.load(ROOT/f'baseline/outputs/null_{v}.npz') for v in variants}
for v in variants:
    for subset in ['ALL','HOTSPOT','NON-HOTSPOT']:
        per=pd.read_csv(ROOT/f'baseline/outputs/clusters_{v}_{subset}.csv',float_precision='round_trip')
        obs=np.array([indexes[v].count(str(r.chrom),[int(r.start)],int(r.end-r.start))[0] for r in per.itertuples()])
        ck(f'baseline_{v}_{subset}_every_cluster_observed_endpoint_union',np.array_equal(obs,per.obs.to_numpy()))
        for key,est in [('event_count','bc_uniqSV_blockNull'),('cluster_sign','f_clusterSignTest')]:
            n=null_by_variant[v][subset+'__'+key]
            o=int(obs.sum()) if key=='event_count' else int(np.count_nonzero(obs>per.null_mean.to_numpy()))
            row=baseline[(baseline.variant==v)&(baseline.subset==subset)&(baseline.estimator==est)].iloc[0]
            for metric,value in stats(o,n).items():ck(f'baseline_{v}_{subset}_{key}_{metric}',near(value,row[fields[metric]]))
for subset in ['ALL','HOTSPOT','NON-HOTSPOT']:
    for key in ['event_count','cluster_sign']:
        ck(f'baseline_canonical_donor_exact_20000_{subset}_{key}',np.array_equal(null_by_variant['canonical_legacy'][subset+'__'+key],null_by_variant['donor_signature'][subset+'__'+key]))
for z in null_by_variant.values():z.close()
summary={'baseline_rows_recovered':len(baseline),'baseline_observed_endpoint_windows':sum(len(pd.read_csv(ROOT/f'baseline/outputs/clusters_{v}_{s}.csv')) for v in variants for s in ['ALL','HOTSPOT','NON-HOTSPOT'])}

if not args.baseline_only:
    run=ROOT/'analysis/outputs';meta=json.loads((run/'run_metadata.json').read_text(encoding='utf-8'))
    ck('new_run_complete',meta['completed'])
    freeze=json.loads((ROOT/'analysis/freeze_record.json').read_text(encoding='utf-8'))
    for p,h in freeze['files'].items():ck('frozen_hash_'+p,sha(ROOT/p)==h)
    table=pd.read_csv(run/'estimates.csv',float_precision='round_trip').set_index(['scope','variant','arm','estimator'])
    per_all=pd.read_csv(run/'per_cluster_results.csv',float_precision='round_trip').set_index(['scope','variant','cluster','arm'])
    diffs=pd.read_csv(run/'paired_changes.csv',float_precision='round_trip').set_index(['scope','variant','estimator','before','after'])
    union_saved=pd.read_csv(run/'event_union_diagnostics.csv').set_index(['scope','variant'])
    arms=['joint_geometry','gene_status','gene_length'];N=20000
    recalculated=[]
    for scope in ['protein_coding','all_gene']:
        co=pd.read_csv(ROOT/f'design/common_coordinates_{scope}.csv',float_precision='round_trip')
        pool=np.load(ROOT/f'design/pools_{scope}.npz');digest=hashlib.sha256()
        block={v:np.zeros((3,N),np.int64) for v in ['donor_signature','paired_sample']}
        sign={v:np.zeros((3,N),np.int64) for v in block};observed={v:0 for v in block};observed_sign={v:np.zeros(3,np.int64) for v in block}
        unions={v:set() for v in block};cluster_issues=[]
        for ii,(cid,g) in enumerate(co.groupby('subset_cluster',sort=True)):
            p=g.pos0.to_numpy(np.int64);chrom=str(g.chrom.iloc[0]);K=len(p);start=int(p.min())-100000;L=int(p.max()-p.min())+200000
            rng=np.random.default_rng(np.random.SeedSequence([20260910,int(cid)]));anchor=rng.integers(0,K,size=N);u=rng.random(N)
            d=np.empty((3,N),np.int64)
            for k in range(K):
                selected=np.flatnonzero(anchor==k)
                for ai,arm in enumerate(arms):
                    values=pool[f'c{int(cid)}_k{k}_{arm}'];d[ai,selected]=values[(u[selected]*len(values)).astype(np.int64)]
            digest.update(np.array([cid],np.int64).tobytes());digest.update(d.tobytes())
            unique,inverse=np.unique(d,return_inverse=True)
            for v in block:
                index=indexes[v];o=int(index.count(chrom,[start],L)[0]);observed[v]+=o
                unions[v].update(index.event_ids(chrom,start,start+L,geoms[v]))
                values=index.count(chrom,start+unique,L);nl=values[inverse].reshape(3,N)
                block[v]+=nl
                # Independent exact integer LOO inequality, no production sign function.
                totals=nl.sum(axis=1);so=o*N>totals;sn=nl*(N-1)>totals[:,None]-nl
                observed_sign[v]+=so.astype(np.int64);sign[v]+=sn.astype(np.int64)
                for ai,arm in enumerate(arms):
                    expected=dict(observed=o,expected=float(nl[ai].mean()),null_sd=float(nl[ai].std(ddof=1)),observed_cluster_sign=int(so[ai]),expected_cluster_sign=float(sn[ai].mean()),window_start=start,window_end=start+L,n_coordinates=K)
                    row=per_all.loc[(scope,v,cid,arm)]
                    for metric,value in expected.items():
                        if not near(value,row[metric]):cluster_issues.append([v,int(cid),arm,metric,value,float(row[metric])])
            if ii%700==0:print(scope,'independent full-draw rebuild',ii,time.time()-began,flush=True)
        pool.close()
        ck(scope+'_all_percluster_counts_and_calibrated_sign_moments',not cluster_issues,{'n_issues':len(cluster_issues),'examples':cluster_issues[:5]})
        for v in block:
            original=np.load(run/f'null_{scope}_{v}.npz')
            for key,a in [('block',block[v]),('sign',sign[v]),('observed_block',observed[v]),('observed_sign',observed_sign[v])]:
                ck(f'{scope}_{v}_{key}_complete_rebuilt_arrays_exact',np.array_equal(a,original[key]))
            original.close()
            ck(scope+'_'+v+'_placement_digest',digest.hexdigest()==union_saved.loc[(scope,v),'placement_sha256'])
            ck(scope+'_'+v+'_global_observed_event_union',len(unions[v])==union_saved.loc[(scope,v),'global_unique_events_observed'])
            ck(scope+'_'+v+'_cross_cluster_repeat_incidence',observed[v]-len(unions[v])==union_saved.loc[(scope,v),'cross_cluster_event_repeat_incidences'])
            for ai,arm in enumerate(arms):
                for est,o,nl in [('cluster_event_incidence',observed[v],block[v][ai]),('calibrated_cluster_sign',observed_sign[v][ai],sign[v][ai])]:
                    s=stats(o,nl);row=table.loc[(scope,v,arm,est)]
                    for metric,value in s.items():
                        if metric not in ['null_sd','z']:ck(f'{scope}_{v}_{arm}_{est}_{metric}',near(value,row[metric]))
                    bonf=min(1.,2*s['p_two_sided']) if scope=='protein_coding' and v=='donor_signature' and arm=='gene_length' else np.nan
                    ck(f'{scope}_{v}_{arm}_{est}_primary_Bonf2',near(bonf,row.p_two_sided_bonferroni2))
                    recalculated.append(dict(scope=scope,variant=v,arm=arm,estimator=est,**s))
            for a,b in [(0,1),(1,2),(0,2)]:
                for est,os,nl in [('cluster_event_incidence',np.repeat(observed[v],3),block[v]),('calibrated_cluster_sign',observed_sign[v],sign[v])]:
                    mu=nl.mean(axis=1);fold=os/mu;cov=np.cov(nl[[a,b]],ddof=1)
                    deriv=np.array([os[a]/mu[a]**2,-os[b]/mu[b]**2]);se=np.sqrt(deriv@cov@deriv/N)
                    expected=dict(fold_before=fold[a],fold_after=fold[b],fold_change=fold[b]-fold[a],log_fold_change=np.log(fold[b]/fold[a]),
                        observed_minus_expected_change=(os[b]-mu[b])-(os[a]-mu[a]),expected_change=mu[b]-mu[a],
                        expected_change_MC_SE=np.sqrt((cov[0,0]+cov[1,1]-2*cov[0,1])/N) if est=='cluster_event_incidence' else np.nan,
                        fold_change_MC_SE=se if est=='cluster_event_incidence' else np.nan)
                    row=diffs.loc[(scope,v,est,arms[a],arms[b])]
                    for metric,value in expected.items():ck(f'{scope}_{v}_{est}_{arms[a]}_to_{arms[b]}_{metric}',near(value,row[metric]))
            np.savez_compressed(OUT/f'independent_null_{scope}_{v}.npz',block=block[v],sign=sign[v],observed_block=np.array(observed[v]),observed_sign=observed_sign[v])
    pd.DataFrame(recalculated).to_csv(OUT/'independent_estimates.csv',index=False,float_format='%.17g')
    # Independent reconstruction of original masked joint pools and their coarse fallback.
    # Read the fixed feature and mask inputs independently of production scoring routines.
    data01=WS/'analysis_inputs/structural_variant_catalogue/donor_duplication_statistics/data'
    base=pd.read_parquet(data01/'genome_features_50kb.parquet')
    base['chrom']=base.chrom.astype(str).map(lambda c:c if c.startswith('chr') else 'chr'+c)
    base=base[(base.stratum>=0)&base.chrom.ne('chrY')]
    base_bins={(str(c),int(s)):np.sort(g.bin_start.to_numpy(np.int64)) for (c,s),g in base.groupby(['chrom','stratum'])}
    masked=set()
    with gzip.open(data01/'gistic_peaks.txt.gz','rt') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            if not next(iter(row.values())).startswith('Amplification Peak'):continue
            m=re.match(r'chr([0-9XY]+):(\d+)-(\d+)',row['Wide Peak Limits'].strip())
            if m:
                c='chr'+m[1];start=max(0,(int(m[2])-100000)//50000)*50000;end=((int(m[3])+100000)//50000)*50000
                masked.update((c,b) for b in range(start,end+50000,50000))
    features=pd.read_parquet(WS/'analysis_inputs/class_composition/data/joint_features.parquet')
    keep=np.array([(str(r.chrom),int(r.bin_start)) not in masked for r in features.itertuples()])
    joint={(str(c),int(s)):np.sort(g.bin_start.to_numpy(np.int64)) for (c,s),g in features[keep].groupby(['chrom','stratum_ref'])}
    def old_pool(chrom,stratum,offset,pos):
        p=np.unique(joint.get((chrom,int(stratum)),np.empty(0,np.int64))+int(offset));p=p[p!=pos]
        if len(p)>=5:return p,False
        p=np.unique(base_bins[(chrom,int(stratum)//12)]+int(offset));return p[p!=pos],True
    for scope in ['protein_coding','all_gene']:
        co=pd.read_csv(ROOT/f'design/common_coordinates_{scope}.csv',float_precision='round_trip')
        rng=np.random.default_rng(7);count={v:np.zeros(N,np.int64) for v in ['donor_signature','paired_sample']}
        signed={v:np.zeros(N,np.int64) for v in count};oo={v:0 for v in count};ss={v:0 for v in count};pool_calls=0;fallbacks=0;issues=[]
        saved_per={v:pd.read_csv(run/f'old_null_common_clusters_{scope}_{v}.csv',float_precision='round_trip') for v in count}
        for newcid,(cid,g) in enumerate(co.groupby('subset_cluster',sort=True)):
            p=g.pos0.to_numpy(np.int64);chrom=str(g.chrom.iloc[0]);start=int(p.min())-100000;L=int(p.max()-p.min())+200000
            j=rng.integers(0,len(p),size=N);d=np.empty(N,np.int64)
            for k,r in enumerate(g.itertuples()):
                sel=j==k;nk=int(sel.sum())
                if nk:
                    pool,fb=old_pool(chrom,r.stratum,r.within_bin_offset,int(p[k]));pool_calls+=1;fallbacks+=int(fb)
                    d[sel]=pool[rng.integers(0,len(pool),size=nk)]-int(p[k])
            uq,inv=np.unique(d,return_inverse=True)
            for v in count:
                nl=indexes[v].count(chrom,start+uq,L)[inv];o=int(indexes[v].count(chrom,[start],L)[0]);total=int(nl.sum())
                count[v]+=nl;signed[v]+=(nl*(N-1)>total-nl).astype(np.int64);oo[v]+=o;ss[v]+=int(o*N>total)
                r=saved_per[v].iloc[newcid]
                if not near([o,nl.mean(),len(p),start,start+L],[r.obs,r.null_mean,r.n_coords,r.start,r.end]):issues.append([v,int(cid)])
            if newcid%700==0:print(scope,'independent original-null rebuild',newcid,time.time()-began,flush=True)
        ck(scope+'_original_null_all_percluster_records',not issues,{'failures':issues[:5]})
        for v in count:
            z=np.load(run/f'old_null_common_{scope}_{v}.npz')
            for key,val in [('block',count[v]),('sign',signed[v]),('observed_block',oo[v]),('observed_sign',ss[v])]:
                ck(f'old_common_{scope}_{v}_{key}_independent_all_draws_exact',np.array_equal(val,z[key]))
            z.close()
            np.savez_compressed(OUT/f'independent_old_common_{scope}_{v}.npz',block=count[v],sign=signed[v],observed_block=np.array(oo[v]),observed_sign=np.array(ss[v]))
        summary[scope+'_historical_pool_calls']=pool_calls;summary[scope+'_historical_fallback_calls']=fallbacks
    common=pd.read_csv(run/'original_null_common_support.csv',float_precision='round_trip')
    for r in common.itertuples():
        z=np.load(run/f'old_null_common_{r.scope}_{r.variant}.npz');kind='block' if r.estimator=='cluster_event_incidence' else 'sign'
        for metric,value in stats(int(z['observed_'+kind]),z[kind]).items():
            if metric not in ['null_sd','z']:ck(f'old_common_{r.scope}_{r.variant}_{kind}_{metric}',near(value,getattr(r,metric)))
        z.close()
    summary['new_full_rebuild']='All 4 scope/representation cells, 3 arms, 20000 draws; direct endpoint union + independent RNG/LOO'
    summary['common_historical_null']='All four complete original-null streams rebuilt from frozen features/GISTIC with independently implemented fallback and direct endpoint Boolean sets'

suffix='baseline' if args.baseline_only else 'numeric_results'
result=dict(n_checks=len(checks),passed=sum(c['passed'] for c in checks),failed=sum(not c['passed'] for c in checks),summary=summary,checks=checks,elapsed_seconds=time.time()-began,
  independent_script_sha256=sha(__file__),input_hashes={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'baseline/outputs/all_baselines_full_precision.csv']+([] if args.baseline_only else [ROOT/'analysis/outputs/estimates.csv',ROOT/'analysis/outputs/paired_changes.csv',ROOT/'analysis/outputs/run_metadata.json',ROOT/'analysis/freeze_record.json'])})
(OUT/(suffix+'_checks.json')).write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({k:result[k] for k in ['n_checks','passed','failed','summary','elapsed_seconds']},ensure_ascii=False))
