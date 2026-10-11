"""Compare three replication-timing references using a fixed common-support design.

Aggregate null streams are saved. Per-cluster draws can be reconstructed from the
candidate pools and recorded seed. Input hashes are checked before results are written."""
from pathlib import Path
import sys,json,hashlib,time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'baseline'))
import baseline_api as B
import numpy as np
import pandas as pd
from build_matched_candidate_pools import TRACKS,VIRUSES,key
from null_core import joint_displacements,axis_of,signed_summary,interaction,summarize_counts
OUT=ROOT/'outputs';N=20000;SEED=20260911

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def csv(rows,name):pd.DataFrame(rows).to_csv(OUT/name,index=False,float_format='%.17g')
def common(est,virus):return pd.read_csv(ROOT/f'design/common_coordinates_{est}_{virus}.csv')
def load_pools(virus):
    with np.load(ROOT/f'design/pools_{virus}.npz') as z:return {k:z[k] for k in z.files}
def draws(est,virus,cid,K,pools):
    seed=[SEED,('composition','atlas').index(est),VIRUSES.index(virus),int(cid)]
    anchors=range(1) if est=='composition' else range(K)
    candidates=[[pools[key(virus,cid,k,t)] for t in TRACKS] for k in anchors]
    delta,anchor,u=joint_displacements(candidates,N,seed)
    h=hashlib.sha256()
    for arr in (delta,anchor,u):h.update(str(arr.dtype).encode());h.update(str(arr.shape).encode());h.update(arr.tobytes())
    return delta,h.hexdigest()

def composition_candidate_scores(chrom,positions,displacements,maps):
    cum,nb,pg=maps;shifted=np.asarray(displacements)[:,None]+np.asarray(positions)[None,:]
    counts=B.E.counts_at(cum,nb,np.full(shifted.size,chrom),shifted.ravel())
    return B.E.props(counts,pg).reshape(len(displacements),len(positions),4).mean(axis=1)

def atlas_candidate_scores(chrom,positions,displacements,index):
    lo=int(np.min(positions))-100000;L=int(np.max(positions))-int(np.min(positions))+200000
    if chrom not in index:return np.zeros(len(displacements),np.int64)
    return index[chrom].count(L,lo+np.asarray(displacements,np.int64))

def run_composition(variant,streams,metadata):
    maps=B.endpoint_maps(B.load_endpoints(variant));saved={};rows=[];detail=[];candidate_archive={}
    observed={};nulls={}
    for virus in VIRUSES:
        co=common('composition',virus);pools=load_pools(virus);C=co.original_cluster.nunique()
        ob=np.zeros(4);nl=np.zeros((3,N,4))
        for j,(cid,g) in enumerate(co.groupby('original_cluster',sort=True)):
            chrom=str(g.chrom.iloc[0]);p=g.pos0.to_numpy(np.int64)
            oc=composition_candidate_scores(chrom,p,np.array([0]),maps)[0];ob+=oc/C
            delta,digest=draws('composition',virus,cid,len(p),pools)
            streams.append(dict(estimand='composition',variant=variant,virus=virus,cluster=int(cid),draw_sha256=digest))
            for t,track in enumerate(TRACKS):
                d=pools[key(virus,cid,0,track)];pc=composition_candidate_scores(chrom,p,d,maps)
                nl[t]+=pc[np.searchsorted(d,delta[t])]/C
                candidate_archive[f'{virus}_c{cid}_{track}']=pc
                detail.append(dict(variant=variant,virus=virus,cluster=int(cid),n_members=len(p),track=track,
                    observed_axis=float(axis_of(oc)),candidate_expected_axis=float(axis_of(pc).mean()),
                    **{f'observed_{c}':float(oc[i]) for i,c in enumerate(B.CLASSES)}))
        observed[virus]=float(axis_of(ob));nulls[virus]=np.stack([axis_of(a) for a in nl])
        saved[virus+'_observed_composition']=ob;saved[virus+'_null_composition']=nl;saved[virus+'_null_axis']=nulls[virus]
        for ti,t in enumerate(TRACKS):rows.append(dict(variant=variant,track=t,contrast=virus,
            n_clusters=C,n_coordinates=len(co),**signed_summary(observed[virus],nulls[virus][ti])))
        metadata.append(dict(variant=variant,estimand='composition',virus=virus,n_coordinates=len(co),n_clusters=C,
            n_reference_breakends=sum(int(c[:,-1].sum()) for c in maps[0].values()),global_prior=maps[2].tolist()))
        print('composition complete',variant,virus,flush=True)
    for ti,t in enumerate(TRACKS):
        oi,ni=interaction(observed['HBV'],observed['HPV'],nulls['HBV'][ti],nulls['HPV'][ti])
        saved[t+'_interaction_null']=ni;saved[t+'_interaction_observed']=np.array(oi)
        r=dict(variant=variant,track=t,contrast='HBV-minus-HPV',**signed_summary(oi,ni))
        r['p_bonferroni_3']=min(1.,3*r['p_two_sided']) if variant=='paired_sample' else np.nan
        rows.append(r)
    np.savez_compressed(OUT/f'composition_null_{variant}.npz',**saved)
    np.savez_compressed(OUT/f'composition_candidate_compositions_{variant}.npz',**candidate_archive)
    csv(detail,f'composition_clusters_{variant}.csv')
    return rows

def run_atlas(variant,streams,metadata):
    idx={c:B.make_index(B.load_geometry(variant,c))[0] for c in B.CLASSES}
    rows=[];detail=[]
    for virus in VIRUSES:
        co=common('atlas',virus);pools=load_pools(virus);ob=np.zeros(4,np.int64)
        total=np.zeros((3,4,N),np.int64);signobs=np.zeros((3,4),np.int64);signnull=np.zeros_like(total)
        candidates={};clusterids=[];clusterobserved=[]
        for j,(cid,g) in enumerate(co.groupby('original_cluster',sort=True)):
            chrom=str(g.chrom.iloc[0]);p=g.pos0.to_numpy(np.int64)
            oc=np.array([atlas_candidate_scores(chrom,p,np.array([0]),idx[c])[0] for c in B.CLASSES]);ob+=oc
            clusterids.append(int(cid));clusterobserved.append(oc)
            delta,digest=draws('atlas',virus,cid,len(p),pools)
            streams.append(dict(estimand='atlas',variant=variant,virus=virus,cluster=int(cid),draw_sha256=digest))
            for ti,t in enumerate(TRACKS):
                d=np.unique(np.concatenate([pools[key(virus,cid,k,t)] for k in range(len(p))]))
                counts=np.stack([atlas_candidate_scores(chrom,p,d,idx[c]) for c in B.CLASSES],axis=1)
                candidates[f'c{cid}_{t}_delta']=d;candidates[f'c{cid}_{t}_counts']=counts
                nl=counts[np.searchsorted(d,delta[ti])].T
                total[ti]+=nl;so=oc>nl.mean(axis=1);sn=nl*N>nl.sum(axis=1)[:,None]
                signobs[ti]+=so;signnull[ti]+=sn.astype(np.int64)
                for ci,c in enumerate(B.CLASSES):detail.append(dict(variant=variant,virus=virus,cluster=int(cid),
                    n_members=len(p),track=t,svclass=c,observed=int(oc[ci]),null_mean=float(nl[ci].mean()),
                    null_sd=float(nl[ci].std(ddof=1)),observed_sign=int(so[ci]),expected_sign=float(sn[ci].mean())))
            if j%800==0:print('atlas',variant,virus,j,flush=True)
        saved=dict(observed=ob,null_counts=total,sign_observed=signobs,sign_null=signnull,
            original_cluster=np.array(clusterids),observed_by_cluster=np.array(clusterobserved))
        np.savez_compressed(OUT/f'atlas_null_{variant}_{virus}.npz',**saved)
        np.savez_compressed(OUT/f'atlas_candidate_counts_{variant}_{virus}.npz',**candidates)
        for ti,t in enumerate(TRACKS):
            for ci,c in enumerate(B.CLASSES):
                for estimator,o,n in [('event_incidence',ob[ci],total[ti,ci]),('calibrated_sign',signobs[ti,ci],signnull[ti,ci])]:
                    r=dict(variant=variant,virus=virus,track=t,svclass=c,estimator=estimator,
                        n_clusters=len(clusterids),n_coordinates=len(co),**summarize_counts(o,n))
                    r['p_bonferroni_48']=min(1.,48*r['p_two_sided']) if variant=='paired_sample' else np.nan
                    rows.append(r)
        metadata.append(dict(variant=variant,estimand='atlas',virus=virus,n_coordinates=len(co),n_clusters=len(clusterids)))
    csv(detail,f'atlas_clusters_{variant}.csv');return rows

def differences(composition,atlas):
    rows=[]
    for variant in B.VARIANTS:
        with np.load(OUT/f'composition_null_{variant}.npz') as z:
            for contrast in (*VIRUSES,'HBV-minus-HPV'):
                group=composition[(composition.variant==variant)&(composition.contrast==contrast)].set_index('track')
                for ti,t in enumerate(TRACKS[1:],start=1):
                    if contrast in VIRUSES:diff=z[contrast+'_null_axis'][ti]-z[contrast+'_null_axis'][0]
                    else:diff=(z['HBV_null_axis'][ti]-z['HBV_null_axis'][0])-(z['HPV_null_axis'][ti]-z['HPV_null_axis'][0])
                    rows.append(dict(estimand='composition',variant=variant,contrast=contrast,track=t,reference='consensus',
                        metric='signed_excess',difference=group.loc[t,'excess']-group.loc['consensus','excess'],
                        monte_carlo_se=float(diff.std(ddof=1)/np.sqrt(N)),uncertainty_type='MC precision of null means; not effect-change CI'))
        for virus in VIRUSES:
            with np.load(OUT/f'atlas_null_{variant}_{virus}.npz') as z:
                for ci,c in enumerate(B.CLASSES):
                    for ti,t in enumerate(TRACKS[1:],start=1):
                        for estimator in ('event_incidence','calibrated_sign'):
                            g=atlas[(atlas.variant==variant)&(atlas.virus==virus)&(atlas.svclass==c)&(atlas.estimator==estimator)].set_index('track')
                            mc=np.nan
                            if estimator=='event_incidence':
                                n0=z['null_counts'][0,ci].astype(float);nt=z['null_counts'][ti,ci].astype(float);o=float(z['observed'][ci])
                                mc=float((-o/nt.mean()**2*nt+o/n0.mean()**2*n0).std(ddof=1)/np.sqrt(N))
                            rows.append(dict(estimand='atlas',variant=variant,contrast=virus,svclass=c,estimator=estimator,track=t,
                                reference='consensus',metric='fold',difference=g.loc[t,'fold']-g.loc['consensus','fold'],monte_carlo_se=mc,
                                uncertainty_type='delta-method MC precision' if estimator=='event_incidence' else 'not estimated: MC-calibrated sign is discontinuous'))
    csv(rows,'track_differences.csv')

def main():
    freeze=json.loads((ROOT/'analysis/freeze_record.json').read_text(encoding='utf-8'))
    for rel,expected in freeze['sha256'].items():
        if sha(ROOT/rel)!=expected:raise RuntimeError('Frozen input changed: '+rel)
    if (OUT/'analysis_complete.json').exists():raise RuntimeError('Completed analysis cannot be overwritten')
    OUT.mkdir(exist_ok=True);start=time.time();er=[];ar=[];streams=[];metadata=[]
    for variant in B.VARIANTS:
        er+=run_composition(variant,streams,metadata);ar+=run_atlas(variant,streams,metadata)
    e=pd.DataFrame(er);a=pd.DataFrame(ar);csv(er,'composition_summary.csv');csv(ar,'atlas_summary.csv')
    csv(streams,'draw_streams.csv');csv(metadata,'reference_metadata.csv');differences(e,a)
    s=pd.DataFrame(streams)
    assert (s.groupby(['estimand','virus','cluster']).draw_sha256.nunique()==1).all()
    (OUT/'analysis_complete.json').write_text(json.dumps(dict(n_permutations=N,seed=SEED,tracks=TRACKS,
        class_order=B.CLASSES,seconds=time.time()-start,freeze_sha256=sha(ROOT/'analysis/freeze_record.json'),
        shared_streams_exact=True,output_rows_composition=len(e),output_rows_atlas=len(a)),indent=2),encoding='utf-8')

if __name__=='__main__':main()
