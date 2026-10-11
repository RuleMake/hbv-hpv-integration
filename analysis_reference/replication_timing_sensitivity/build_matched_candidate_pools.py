"""Build shared feature partitions, cluster-placement pools and balance summaries.

Only coordinates, genomic features and replication-timing tracks enter this step.
Structural-variant outcomes are scored after the design has been fixed. The script
refuses to change an existing fixed design. Run with Python -B."""
from pathlib import Path
import sys,json,hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]; WS=ROOT.parents[1]
sys.path.insert(0,str(ROOT/'baseline'))
import baseline_api as B
import numpy as np
import pandas as pd
OUT=ROOT/'design'; TRACKS=('consensus','helas3','hepg2'); VIRUSES=('HBV','HPV')
BIN=50000; MIN_POOL=5

def key(virus,cid,k,track):return f'{virus}_c{int(cid)}_k{int(k)}_{track}'
def write_json(name,v):(OUT/name).write_text(json.dumps(v,indent=2,ensure_ascii=False),encoding='utf-8')
def indices(chrom,points,lookup):
    b=np.asarray(points,dtype=np.int64)//BIN; src=lookup[chrom]
    valid=(b>=0)&(b<len(src)); result=np.full(b.shape,-1,np.int64)
    result[valid]=src[b[valid]];return result

def features():
    f=pd.read_parquet(ROOT/'baseline/outputs/fixed_nonRT_features.parquet')
    r=pd.read_parquet(ROOT/'tracks/outputs/rt_comparison_50kb.parquet')
    f=f.merge(r,on=['chrom','bin_start'],validate='one_to_one',suffixes=('','_rt'))
    assert len(f)==60716
    assert np.array_equal(f.rt_consensus,f.rt_consensus_legacy)
    checks=[]; lookup={}
    for t,col in zip(TRACKS,('rt_consensus_legacy','rt_helas3_z','rt_hepg2_z')):
        f[t+'_rt']=f[col]
        f[t+'_rt_tier']=pd.qcut(f[col],3,labels=False,duplicates='drop').astype(int)
        f[t+'_grad']=np.nan
        for chrom,g in f.groupby('chrom',sort=False):
            ix=g.sort_values('bin_start').index.to_numpy(); v=f.loc[ix,col].to_numpy()
            assert np.all(np.diff(f.loc[ix,'bin_start'])==BIN)
            f.loc[ix[1:-1],t+'_grad']=np.abs((v[2:]-v[:-2])/2)
        finite=f[t+'_grad'].notna(); f[t+'_grad_tier']=2
        f.loc[finite,t+'_grad_tier']=pd.qcut(f.loc[finite,t+'_grad'],2,labels=False).astype(int)
        pre=f[t+'_rt_tier']*6+f.original_gd_tier*2+f.in_cfs.astype(int)
        f[t+'_joint']=(pre*3+f[t+'_grad_tier'])*4+f.original_map_tier
        checks.append(dict(track=t,n_bins=len(f),n_rt_tiers=f[t+'_rt_tier'].nunique(),
            rt_cut_q1=float(f[col].quantile(1/3)),rt_cut_q2=float(f[col].quantile(2/3)),
            gradient_median=float(f.loc[finite,t+'_grad'].median()),missing_gradient=int((~finite).sum()),
            n_joint_strata=f[t+'_joint'].nunique()))
    assert np.array_equal(f.consensus_rt_tier,f.original_rt_tier)
    assert np.array_equal(f.consensus_grad_tier,f.original_grad_tier)
    both=f[['consensus_joint','original_joint_stratum']].drop_duplicates()
    assert both.consensus_joint.is_unique and both.original_joint_stratum.is_unique
    for chrom,g in f.groupby('chrom',sort=False):
        a=np.full(int(g.bin_start.max())//BIN+1,-1,np.int64)
        good=g[~g.gistic_masked]; a[good.bin_start.to_numpy(np.int64)//BIN]=good.index.to_numpy()
        lookup[str(chrom)]=a
    f.to_parquet(OUT/'features.parquet',index=False)
    pd.DataFrame(checks).to_csv(OUT/'partition_parameters.csv',index=False,float_format='%.17g')
    cross=[]
    for a in TRACKS:
        for b in TRACKS:
            cross.append(dict(track_a=a,track_b=b,pearson_rt=f[a+'_rt'].corr(f[b+'_rt']),
                spearman_rt=f[a+'_rt'].rank(method='average').corr(f[b+'_rt'].rank(method='average')),
                rt_tier_agreement=float((f[a+'_rt_tier']==f[b+'_rt_tier']).mean()),
                grad_tier_agreement=float((f[a+'_grad_tier']==f[b+'_grad_tier']).mean())))
    pd.DataFrame(cross).to_csv(OUT/'track_similarity.csv',index=False,float_format='%.17g')
    return f,lookup

def summarize_distribution(values,ow,nw):
    values=np.asarray(values,float); valid=np.isfinite(values)
    missing_o=1-ow[valid].sum()/ow.sum(); missing_n=1-nw[valid].sum()/nw.sum()
    v=values[valid]; o=ow[valid]/ow[valid].sum(); n=nw[valid]/nw[valid].sum()
    om=float(v@o);nm=float(v@n); ov=float(((v-om)**2)@o);nv=float(((v-nm)**2)@n)
    den=np.sqrt((ov+nv)/2); order=np.argsort(v,kind='stable'); vv=v[order]
    # Evaluate only at final index of a tied value, so KS never splits a tie.
    last=np.r_[vv[1:]!=vv[:-1],True]
    ks=float(np.max(np.abs((np.cumsum(o[order])-np.cumsum(n[order]))[last])))
    smd=(om-nm)/den if den else (0. if om==nm else float(np.copysign(np.inf,om-nm)))
    return dict(observed_mean=om,null_mean=nm,smd=smd,degenerate_unequal_means=bool(den==0 and om!=nm),ks=ks,
        observed_missing_fraction=float(missing_o),null_missing_fraction=float(missing_n))

def diagnostics(est,virus,co,pools,f,lookup):
    names=['gene_density','in_cfs','gpos100','mapab','cen_dist','telo_dist']
    names += [t+x for t in TRACKS for x in ('_rt','_grad','_rt_tier','_grad_tier')]
    obs={w:np.zeros(len(f)) for w in ('coordinate','cluster')}
    nul={(t,w):np.zeros(len(f)) for t in TRACKS for w in obs}
    for cid,g in co.groupby('original_cluster',sort=True):
        chrom=str(g.chrom.iloc[0]);p=g.pos0.to_numpy(np.int64);K=len(p)
        ix=indices(chrom,p,lookup);assert np.all(ix>=0)
        for w,scale in [('coordinate',1.),('cluster',1./K)]:np.add.at(obs[w],ix,scale)
        anchors=range(K) if est=='atlas' else range(1)
        for k in anchors:
            for t in TRACKS:
                delta=pools[key(virus,cid,k,t)];ci=indices(chrom,delta[:,None]+p[None,:],lookup)
                freq=np.bincount(ci.ravel(),minlength=len(f))/len(anchors)/len(delta)
                nul[t,'coordinate']+=freq;nul[t,'cluster']+=freq/K
    rows=[]
    for (t,w),mass in nul.items():
        assert np.isclose(mass.sum(),obs[w].sum(),atol=1e-8,rtol=0)
        for name in names:rows.append(dict(estimand=est,virus=virus,track=t,weighting=w,covariate=name,
            **summarize_distribution(f[name].to_numpy(float),obs[w],mass)))
    pd.DataFrame(rows).to_csv(OUT/f'balance_{est}_{virus}.csv',index=False,float_format='%.17g')
    np.savez_compressed(OUT/f'balance_bin_weights_{est}_{virus}.npz',
        **{f'observed_{w}':m for w,m in obs.items()},**{f'{t}_{w}':m for (t,w),m in nul.items()})

def main():
    if (ROOT/'analysis/freeze_record.json').exists():raise RuntimeError('Design already frozen')
    OUT.mkdir(exist_ok=True);f,lookup=features()
    sizes=pd.read_csv(WS/'analysis_inputs/gene_length_matching/annotation/outputs/chrom_sizes.csv').set_index('chrom').chrom_size_bp.to_dict()
    candidate={t:{(str(c),int(s)):g.bin_start.to_numpy(np.int64) for (c,s),g in f[~f.gistic_masked].groupby(['chrom',t+'_joint'])} for t in TRACKS}
    summaries=[];allcoverage=[]
    for virus in VIRUSES:
        co=B.load_coordinates(virus);pools={};coverage=[]
        for cid,g in co.groupby('original_cluster',sort=True):
            chrom=str(g.chrom.iloc[0]); p=g.pos0.to_numpy(np.int64); oi=indices(chrom,p,lookup)
            observed_ok=bool(np.all(oi>=0))
            for k,r in enumerate(g.itertuples()):
                row=dict(virus=virus,cluster=int(cid),anchor_index=k,n_members=len(g),
                    coordinate_id=r.tian_coordinate_id,chrom=chrom,observed_all_members_eligible=observed_ok)
                for t in TRACKS:
                    if oi[k]>=0:
                        bins=candidate[t].get((chrom,int(f.loc[oi[k],t+'_joint'])),np.array([],np.int64))
                        delta=np.unique(bins+int(r.within_bin_offset)-int(r.pos0));delta=delta[delta!=0]
                    else:delta=np.array([],np.int64)
                    row[t+'_anchor_pool']=len(delta)
                    shifted=delta[:,None]+p[None,:]
                    good=np.all((indices(chrom,shifted,lookup)>=0)&(shifted>=0)&(shifted<sizes[chrom]),axis=1)
                    delta=delta[good] if observed_ok else np.array([],np.int64)
                    pools[key(virus,cid,k,t)]=delta;row[t+'_geometry_pool']=len(delta)
                coverage.append(row)
        cv=pd.DataFrame(coverage); cv['all_three_ge5']=cv[[t+'_geometry_pool' for t in TRACKS]].min(axis=1)>=MIN_POOL
        for est in ('composition','atlas'):
            anchor_cv=cv[cv.anchor_index==0] if est=='composition' else cv
            support=anchor_cv.groupby('cluster').all_three_ge5.all();ids=support.index[support]
            sub=co[co.original_cluster.isin(ids)].copy();sub.to_csv(OUT/f'common_coordinates_{est}_{virus}.csv',index=False)
            cv['in_'+est+'_support']=cv.cluster.isin(ids)
            summaries.append(dict(estimand=est,virus=virus,full_coordinates=len(co),full_clusters=co.original_cluster.nunique(),
                common_coordinates=len(sub),common_clusters=len(ids),dropped_clusters=int((~support).sum()),
                min_candidates=MIN_POOL,no_fallback=True,SV_outcomes_accessed=False))
            diagnostics(est,virus,sub,pools,f,lookup)
            print(est,virus,'support',len(ids),'clusters',len(sub),'points',flush=True)
        np.savez_compressed(OUT/f'pools_{virus}.npz',**pools)
        allcoverage.append(cv)
    pd.concat(allcoverage).to_csv(OUT/'anchor_coverage.csv',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'support_summary.csv',index=False)
    write_json('design_gate.json',dict(outcomes_read=False,consensus_partition_exact=True,full_grid_bins=len(f),
        gistic_masked_bins=int(f.gistic_masked.sum()),tracks=TRACKS,minimum_candidates=MIN_POOL))

if __name__=='__main__':main()
