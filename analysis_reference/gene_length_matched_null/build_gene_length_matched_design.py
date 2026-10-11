"""Construct genomic-context and gene-length matching pools before scoring outcomes.

The design uses coordinate and annotation tables without reading structural-
variant geometry, event counts or spatial-null outcomes."""
from pathlib import Path
import sys,json,hashlib,argparse,time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
WS=ROOT.parents[1]
sys.path.insert(0,str(ROOT/'baseline'))
sys.path.insert(0,str(ROOT/'annotation'))
import load_inputs_and_original_null as B
import numpy as np
import pandas as pd
from matching_core import length_candidate_masks

OUT=ROOT/'design'
SCOPES=('protein_coding','all_gene')
ARMS=('joint_geometry','gene_status','gene_length')
CALIPER=.20
MIN_POOL=5
FEATURES=('rt_consensus','gene_density','in_cfs','gpos100','grad','mapab','cen_dist','telo_dist')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def savejson(name,x):(OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def load():
    ctx=B.load_context();co=B.load_coordinates()
    feat=ctx['feat'].sort_values(['chrom','bin_start']).reset_index(drop=True)
    lookup={}
    for chrom,g in feat.groupby('chrom',sort=False):
        ix=np.full(int(g.bin_start.max())//50000+1,-1,dtype=np.int64)
        for i,r in g.iterrows():
            if (str(chrom),int(r.bin_start)) not in ctx['mask']:
                ix[int(r.bin_start)//50000]=i
        lookup[str(chrom)]=ix
    ctx['reference_sizes']=pd.read_csv(ROOT/'annotation/outputs/chrom_sizes.csv').set_index('chrom').chrom_size_bp.to_dict()
    return ctx,co,feat,lookup

def indices(chrom,points,lookup):
    b=np.asarray(points,dtype=np.int64)//50000
    src=lookup[str(chrom)];out=np.full(b.shape,-1,dtype=np.int64)
    valid=(b>=0)&(b<len(src));out[valid]=src[b[valid]]
    return out

def geometry(ctx,co,lookup):
    pools={};rows=[]
    for cid,g in co.groupby('subset_cluster',sort=True):
        p=g.pos0.to_numpy(np.int64);chrom=str(g.chrom.iloc[0])
        if np.any(indices(chrom,p,lookup)<0):raise ValueError('Observed point outside eligible bins')
        for k,r in enumerate(g.itertuples()):
            bins=ctx['nul'].bins.get((chrom,int(r.stratum)),np.array([],dtype=np.int64))
            delta=np.unique(bins+int(r.within_bin_offset)-int(r.pos0));delta=delta[delta!=0]
            n_joint=len(delta)
            shifted=delta[:,None]+p[None,:]
            ok=np.all((indices(chrom,shifted,lookup)>=0)&(shifted>=0)&(shifted<ctx['reference_sizes'][chrom]),axis=1)
            delta=delta[ok]
            key=f'c{int(cid)}_k{k}'
            pools[key]=delta
            rows.append(dict(cluster=int(cid),anchor_index=k,tian_coordinate_id=r.tian_coordinate_id,
                             chromosome=chrom,n_members=len(p),n_joint_no_fallback=n_joint,
                             n_geometry=len(delta),geometry_pool_ge5=len(delta)>=MIN_POOL))
    np.savez_compressed(OUT/'geometry_pools.npz',**pools)
    pd.DataFrame(rows).to_csv(OUT/'geometry_anchor_coverage.csv',index=False)
    return pools

def annotation_adapter(engine,chrom,points,scope):
    # Gene annotation is provided by the independently tested annotation module.
    flat=np.asarray(points).reshape(-1)
    frame=engine.annotate_arrays(np.full(len(flat),chrom),flat,scope=scope,compact=True)
    if np.any(frame['n_overlapping_genes']<0):raise ValueError('unsupported reference point in design')
    return (frame['n_overlapping_genes'].astype(np.int64).reshape(np.asarray(points).shape),
            frame['median_log10_span'].reshape(np.asarray(points).shape))

def weighted_distribution(values,weights):
    d={}
    for v,w in zip(np.asarray(values).reshape(-1),np.asarray(weights).reshape(-1)):
        if np.isfinite(v):d[float(v)]=d.get(float(v),0.)+float(w)
    return d

def merge_dist(a,b):
    for k,v in b.items():a[k]=a.get(k,0.)+v

def distribution_summary(obs,ow,null_dict):
    obs=np.asarray(obs,float);ow=np.asarray(ow,float);valid=np.isfinite(obs)
    ow=ow[valid];obs=obs[valid]
    if not len(obs) or not null_dict:return dict(smd=None,ks=None,observed_mean=None,null_mean=None)
    x=np.array(sorted(null_dict));w=np.array([null_dict[k] for k in x]);w=w/w.sum();ow=ow/ow.sum()
    om=float(np.sum(obs*ow));nm=float(np.sum(x*w))
    ov=float(np.sum((obs-om)**2*ow));nv=float(np.sum((x-nm)**2*w));den=np.sqrt((ov+nv)/2)
    sv=np.argsort(obs);ox=obs[sv];oc=np.cumsum(ow[sv]);nc=np.cumsum(w)
    grid=np.unique(np.concatenate((ox,x)))
    oi=np.searchsorted(ox,grid,side='right')-1;ni=np.searchsorted(x,grid,side='right')-1
    oprob=np.where(oi<0,0,oc[np.maximum(0,oi)]);nprob=np.where(ni<0,0,nc[np.maximum(0,ni)])
    def quantile(vals,cumulative,q):return float(vals[min(np.searchsorted(cumulative,q),len(vals)-1)])
    return dict(smd=float((om-nm)/den) if den>0 else (0. if om==nm else None),
                ks=float(np.max(np.abs(oprob-nprob))),observed_mean=om,null_mean=nm,
                observed_q10=quantile(ox,oc,.1),observed_q50=quantile(ox,oc,.5),observed_q90=quantile(ox,oc,.9),
                null_q10=quantile(x,nc,.1),null_q50=quantile(x,nc,.5),null_q90=quantile(x,nc,.9))

def build(engine,scope,ctx,co,feat,lookup,basepools):
    pools={};rows=[];observed=[];start=time.time()
    # Generate diagnostics after the final common-support set has been defined.
    for number,(cid,g) in enumerate(co.groupby('subset_cluster',sort=True)):
        p=g.pos0.to_numpy(np.int64);chrom=str(g.chrom.iloc[0]);oc,ol=annotation_adapter(engine,chrom,p,scope)
        for r,count,length in zip(g.itertuples(),oc,ol):
            observed.append(dict(tian_coordinate_id=r.tian_coordinate_id,cluster=int(cid),chrom=chrom,
                                 pos0=int(r.pos0),gene_count=int(count),median_log10_span=float(length)))
        for k in range(len(p)):
            basekey=f'c{int(cid)}_k{k}';delta=basepools[basekey]
            if len(delta):
                cc,cl=annotation_adapter(engine,chrom,delta[:,None]+p[None,:],scope)
                status,length=length_candidate_masks(oc,ol,cc,cl,CALIPER)
            else:status=length=np.zeros(0,dtype=bool)
            for arm,d in zip(ARMS,(delta,delta[status],delta[length])):
                pools[basekey+'_'+arm]=d
            rows.append(dict(scope=scope,cluster=int(cid),anchor_index=k,n_members=len(p),
                             n_joint_geometry=len(delta),n_gene_status=int(status.sum()),n_gene_length=int(length.sum())))
        if number%400==0:print(f'{scope} annotation {number}/2345 clusters {time.time()-start:.1f}s',flush=True)
    coverage=pd.DataFrame(rows);minimum=coverage.groupby('cluster')[['n_joint_geometry','n_gene_status','n_gene_length']].min()
    support=minimum.index[minimum.n_gene_length>=MIN_POOL].to_numpy(np.int64)
    coverage['in_common_support']=coverage.cluster.isin(support)
    coverage.to_csv(OUT/f'anchor_coverage_{scope}.csv',index=False)
    pd.DataFrame(observed).to_csv(OUT/f'observed_annotations_{scope}.csv',index=False)
    np.savez_compressed(OUT/f'pools_{scope}.npz',**pools)
    sub=co[co.subset_cluster.isin(support)].copy()
    sub.to_csv(OUT/f'common_coordinates_{scope}.csv',index=False)
    savejson(f'support_{scope}.json',dict(scope=scope,full_coordinates=len(co),full_clusters=co.subset_cluster.nunique(),
        geometry_clusters_ge5=int((minimum.n_joint_geometry>=MIN_POOL).sum()),
        status_clusters_ge5=int((minimum.n_gene_status>=MIN_POOL).sum()),
        length_clusters_ge5=len(support),common_coordinates=len(sub),common_clusters=support.tolist(),
        caliper_log10=CALIPER,min_candidates_per_every_original_anchor=MIN_POOL,
        dropped_clusters=minimum.index[minimum.n_gene_length<MIN_POOL].astype(int).tolist(),
        no_fallback=True,SV_outcome_accessed=False))
    diagnostics(engine,scope,sub,feat,lookup,pools)

def diagnostics(engine,scope,sub,feat,lookup,pools):
    # Exactly weight every anchor equally, every eligible candidate within anchor equally,
    # then each member point. Thus large candidate lists do not upweight their anchors.
    names=('median_log10_span','genic','multi_overlap')+FEATURES
    store={(arm,weight,name):{} for arm in ARMS for weight in ('coordinate','cluster') for name in names}
    moments=[];obsvals={name:[] for name in names};obsweights={'coordinate':[],'cluster':[]}
    total_null_mass={(arm,weight,name):0. for arm in ARMS for weight in ('coordinate','cluster') for name in names}
    reuse=[]
    for cid,g in sub.groupby('subset_cluster',sort=True):
        p=g.pos0.to_numpy(np.int64);chrom=str(g.chrom.iloc[0]);K=len(p)
        oc,ol=annotation_adapter(engine,chrom,p,scope);oi=indices(chrom,p,lookup)
        observed=[ol,(oc>0).astype(float),(oc>1).astype(float)]+[feat.loc[oi,f].to_numpy(float) for f in FEATURES]
        for name,v in zip(names,observed):obsvals[name].extend(v.tolist())
        obsweights['coordinate'].extend([1.]*K);obsweights['cluster'].extend([1./K]*K)
        for arm in ARMS:
            per_moment={name:np.zeros((K,3)) for name in names}
            union=[]
            for k in range(K):
                delta=pools[f'c{int(cid)}_k{k}_{arm}'];N=len(delta);union.extend(delta.tolist())
                points=delta[:,None]+p[None,:];cc,cl=annotation_adapter(engine,chrom,points,scope)
                ii=indices(chrom,points,lookup)
                vals=[cl,(cc>0).astype(float),(cc>1).astype(float)]+[feat[f].to_numpy(float)[ii] for f in FEATURES]
                for name,v in zip(names,vals):
                    valid=np.isfinite(v)
                    per_moment[name][:,0]+=np.where(valid,v,0).sum(axis=0)/(K*N)
                    per_moment[name][:,1]+=np.where(valid,v*v,0).sum(axis=0)/(K*N)
                    per_moment[name][:,2]+=valid.sum(axis=0)/(K*N)
                    values,counts=np.unique(v[valid],return_counts=True)
                    for weight,scale in [('coordinate',1.),('cluster',1./K)]:
                        mass=counts*scale/(K*N)
                        key=(arm,weight,name)
                        merge_dist(store[key],dict(zip(values.tolist(),mass.tolist())))
                        total_null_mass[key]+=K*scale/K
            unique=np.unique(union)
            reuse.append(dict(scope=scope,cluster=int(cid),arm=arm,distinct_displacements=len(unique),
                              anchor_candidate_entries=len(union),n_anchors=K,
                              same_placement_reused_across_anchors=len(union)-len(unique)))
            for j,r in enumerate(g.itertuples()):
                for name in names:
                    m=per_moment[name][j]
                    moments.append(dict(scope=scope,cluster=int(cid),coordinate_id=r.tian_coordinate_id,arm=arm,
                                        covariate=name,expected_sum=float(m[0]),expected_sum_squared=float(m[1]),finite_probability=float(m[2])))
    rows=[]
    for (arm,weight,name),dist in store.items():
        obs=np.array(obsvals[name]);ow=np.array(obsweights[weight]);nullmass=total_null_mass[(arm,weight,name)]
        finite_mass=sum(dist.values());observed_missing=float(np.sum(ow[~np.isfinite(obs)])/ow.sum())
        r=dict(scope=scope,arm=arm,weighting=weight,covariate=name,
               observed_missing_fraction=observed_missing,
               null_missing_fraction=1-finite_mass/nullmass if nullmass else None,
               **distribution_summary(obs,ow,dist))
        rows.append(r)
    pd.DataFrame(rows).to_csv(OUT/f'balance_{scope}.csv',index=False,float_format='%.17g')
    pd.DataFrame(moments).to_csv(OUT/f'point_expected_covariates_{scope}.csv',index=False,float_format='%.17g')
    pd.DataFrame(reuse).to_csv(OUT/f'candidate_reuse_{scope}.csv',index=False)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--geometry-only',action='store_true');ap.add_argument('--scope',choices=SCOPES)
    args=ap.parse_args();OUT.mkdir(exist_ok=True)
    if (ROOT/'analysis/freeze_record.json').exists():raise RuntimeError('Design is frozen; do not overwrite')
    ctx,co,feat,lookup=load()
    if args.geometry_only:
        geometry(ctx,co,lookup);print('Outcome-free geometry pools ready');return
    from gene_span_annotation import GeneSpanAnnotation
    engine=GeneSpanAnnotation.from_csv()
    archive=np.load(OUT/'geometry_pools.npz');basepools={k:archive[k] for k in archive.files};archive.close()
    for scope in ([args.scope] if args.scope else SCOPES):build(engine,scope,ctx,co,feat,lookup,basepools)
    savejson('design_manifest.json',dict(script_sha256=sha(__file__),core_sha256=sha(ROOT/'analysis/matching_core.py'),
        coordinate_sha256=sha(ROOT/'baseline/outputs/coordinates_NON-HOTSPOT.csv'),
        feature_sha256=sha(WS/'analysis_inputs/class_composition/data/joint_features.parquet'),
        gene_span_sha256=sha(ROOT/'annotation/outputs/gencode_v19_gene_spans.csv'),outcome_blind_candidate_construction=True,
        no_SV_files_read=True,retrospective_prior_results_already_known=True))

if __name__=='__main__':main()
