"""Verify class-composition estimates using direct endpoint searches and independent pools.

The checks compare saved statistics and three complete 20,000-draw streams without
importing the estimator under review. Run with Python -B; results are written to JSON."""
from pathlib import Path
import sys,json,hashlib,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;WS=ROOT.parents[1]
sys.path.insert(0,str(WS/'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
sys.path.insert(0,str(ROOT/'shared'))
import numpy as np
import pandas as pd
import shared_analysis_context as C
OUT=ROOT/'axis/outputs';checks=[];full=[];CLASSES=('DEL','DUP','INV','TRA')
def read(p):return pd.read_csv(p,float_precision='round_trip')
def check(name,ok,detail=None):checks.append(dict(check=name,passed=bool(ok),detail=detail))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def labels(chrom,pos):
    out=[];c0=None;p0=None;k=-1
    for c,p in zip(chrom,pos):
        if c!=c0 or p-p0>100000:k+=1;c0=c;p0=p
        out.append(k)
    return np.array(out,np.int64)
ctx=C.load_context();joint=ctx['nul'];feature=ctx['feat'];mask=ctx['mask']
clean=feature.loc[[(str(r.chrom),int(r.bin_start)) not in mask for r in feature.itertuples()]]
fine={(str(c),int(s)):g.bin_start.to_numpy(np.int64) for (c,s),g in clean.groupby(['chrom','stratum_ref'])}
bf=pd.read_parquet(C.F.FEATURES);bf['chrom']=bf.chrom.astype(str).map(lambda x:x if x.startswith('chr') else 'chr'+x)
bf=bf[(bf.stratum>=0)&bf.chrom.ne('chrY')]
base={(str(c),int(s)):g.bin_start.to_numpy(np.int64) for (c,s),g in bf.groupby(['chrom','stratum'])}
def jointpool(c,s,off,p):
    cand=np.unique(fine.get((c,s),np.array([],np.int64))+off);cand=cand[cand!=p]
    if len(cand)<5:
        cand=np.unique(base[(c,s//joint.divisor)]+off);cand=cand[cand!=p]
    return cand
ind=read(HERE/'independent_composition_observed.csv')
for variant in C.VARIANTS:
    where=OUT/variant;arrays=np.load(where/'null_draws.npz')
    tables=[read(where/n) for n in ('class_composition.csv','axis.csv','strict_sensitivity.csv')]
    bp=C.load_breakpoints(variant,exclude_positive=True)
    pg=bp.svclass.value_counts().reindex(CLASSES).to_numpy(float)/len(bp)
    comp=read(where/'background_composition.csv')
    check(variant+' global endpoint counts/priors',np.array_equal(comp.n_breakpoints,bp.svclass.value_counts().reindex(CLASSES)) and np.array_equal(comp.global_fraction,pg))
    for r in ind[ind.variant==variant].itertuples():
        oc=arrays[f'stratum__{r.virus}__observed_class']
        check(variant+' '+r.virus+' direct endpoint observed',np.array_equal(oc,[getattr(r,'observed_prop_'+c) for c in CLASSES]))
        d=tables[0][(tables[0].virus==r.virus)&(tables[0].arm=='stratum')&(tables[0].statistic=='breakpoint_density')].iloc[0]
        check(variant+' '+r.virus+' direct endpoint density',np.isclose(d.observed,r.mean_endpoint_density,rtol=0,atol=1e-12))
    for table in tables:
        for r in table.itertuples():
            arm='density_strict' if r.rung=='strict' else r.arm
            tag=f'{variant} {r.virus} {arm} {r.statistic}'
            if r.virus=='HBV-vs-HPV':
                nl=arrays[f'{arm}__interaction_null']
                h=arrays[f'{arm}__HBV__axis_null'];p=arrays[f'{arm}__HPV__axis_null']
                reconstructed=(h-h.mean())-(p-p.mean())
                check(tag+' centered independent streams',np.array_equal(nl,reconstructed))
            elif r.statistic=='class_axis_TRAINV_minus_DELDUP':nl=arrays[f'{arm}__{r.virus}__axis_null']
            elif r.statistic=='breakpoint_density':nl=arrays[f'{arm}__{r.virus}__density_null']
            else:nl=arrays[f'{arm}__{r.virus}__class_null'][:,CLASSES.index(r.statistic[5:])]
            hi=(np.count_nonzero(nl>=r.observed)+1)/20001;lo=(np.count_nonzero(nl<=r.observed)+1)/20001
            qlo,qhi=np.percentile(nl,[2.5,97.5]);mu=nl.mean();sd=nl.std(ddof=1)
            vals={'null_mean':mu,'null_sd':sd,'null_reference_q025':qlo,'null_reference_q975':qhi,'excess':r.observed-mu,'z':(r.observed-mu)/sd,'emp_p_greater':hi,'emp_p_less':lo,'emp_p_two_sided':min(1,2*min(hi,lo)),'mde_excess_approx_2p80sd':2.8*sd}
            wrong={k:[float(v),float(getattr(r,k))] for k,v in vals.items() if v!=getattr(r,k)}
            check(tag+' saved summary',not wrong,wrong)
    # Direct endpoint searches at all feature-bin centers independently define density.
    positions={(c,cl):np.sort(g.pos0.to_numpy(np.int64)) for (c,cl),g in bp.groupby(['chrom','svclass'])}
    def direct(c,b):
        b=np.asarray(b,np.int64);out=np.zeros((b.size,4),np.int64)
        for k,cl in enumerate(CLASSES):
            p=positions.get((c,cl),np.array([],np.int64))
            out[:,k]=np.searchsorted(p,(b.ravel()+21)*50000,side='left')-np.searchsorted(p,(b.ravel()-20)*50000,side='left')
        return out
    member=read(where/'density_bin_membership.csv')
    dd=np.zeros(len(clean),np.int64)
    for c,g in clean.groupby('chrom'):
        ii=np.flatnonzero(clean.chrom.to_numpy()==c)
        dd[ii]=direct(c,g.bin_start.to_numpy(np.int64)//50000).sum(axis=1)
    order=np.argsort(dd,kind='stable');ranks=np.empty(len(dd),np.int64);ranks[order]=np.arange(1,len(dd)+1)
    tiers=pd.qcut(ranks,10,labels=False)
    check(variant+' full density bin counts and stable tie ranks',np.array_equal(member.breakpoint_density,dd) and np.array_equal(member.tier,tiers) and np.array_equal(member.bin_start,clean.bin_start) and np.array_equal(member.chrom,clean.chrom))
    tiermap={(str(c),int(b)):int(t) for c,b,t in zip(clean.chrom,clean.bin_start,tiers)}
    db={}
    for c,s,b,t in zip(clean.chrom,clean.stratum_ref,clean.bin_start,tiers):db.setdefault((str(c),int(s),int(t)),[]).append(int(b))
    db={k:np.array(v,np.int64) for k,v in db.items()}
    for virus in C.VIRUSES:
        co=ctx['co_joint'];co=co[co.virus_class==virus].sort_values(['chrom','pos0']).reset_index(drop=True)
        lab=labels(co.chrom,co.pos0);lead=np.r_[0,np.flatnonzero(np.diff(lab))+1]
        sizes=[]
        for i in lead:
            r=co.iloc[i];off=int(r.pos0)%50000;t=tiermap[(r.chrom,int(r.pos0)//50000*50000)]
            cand=np.unique(db.get((r.chrom,int(r.stratum),t),np.array([],np.int64))+off);cand=cand[cand!=r.pos0]
            sizes.append(len(cand))
        sizes=np.array(sizes);s=read(where/f'{virus}_strict_pool_audit.csv')
        check(variant+' '+virus+' strict support',np.array_equal(s.strict_pool_size,sizes) and np.array_equal(s.strict,sizes>=5))
    selections=[]
    if variant=='paired_sample':selections=[('HBV','stratum'),('HPV','density')]
    if variant=='donor_signature':selections=[('EBV','density_strict')]
    for virus,arm in selections:
        started=time.time();co=ctx['co_joint'];co=co[co.virus_class==virus].sort_values(['chrom','pos0']).reset_index(drop=True)
        lab=labels(co.chrom,co.pos0)
        if arm=='density_strict':
            lead=np.r_[0,np.flatnonzero(np.diff(lab))+1];keep=[]
            for i in lead:
                r=co.iloc[i];t=tiermap[(r.chrom,int(r.pos0)//50000*50000)]
                cand=np.unique(db.get((r.chrom,int(r.stratum),t),np.array([],np.int64))+int(r.pos0)%50000)
                keep.append(np.count_nonzero(cand!=r.pos0)>=5)
            co=co[np.array(keep)[lab]].reset_index(drop=True);lab=labels(co.chrom,co.pos0)
        lead=np.r_[0,np.flatnonzero(np.diff(lab))+1];K=len(lead);pos=co.pos0.to_numpy(np.int64);chrom=co.chrom.to_numpy(str)
        weights=1/(K*np.bincount(lab)[lab]);delta=pos-pos[lead][lab];pools=[]
        for i in lead:
            r=co.iloc[i];off=int(r.pos0)%50000;s=int(r.stratum)
            cand=np.array([],np.int64)
            if arm!='stratum':
                t=tiermap[(r.chrom,int(r.pos0)//50000*50000)]
                cand=np.unique(db.get((r.chrom,s,t),np.array([],np.int64))+off);cand=cand[cand!=r.pos0]
            if len(cand)<5:cand=jointpool(r.chrom,s,off,int(r.pos0))
            pools.append(cand)
        size=np.array([len(p) for p in pools]);start=np.r_[0,np.cumsum(size)[:-1]];flat=np.concatenate(pools)
        # Maps built by independent binary searches, not production cumulative counts.
        maps={};indices={}
        for c,g in bp.groupby('chrom'):
            n=int(g.pos0.max())//50000+2;cnt=direct(str(c),np.arange(n+21))
            maps[str(c)]=((cnt+30*pg)/(cnt.sum(axis=1)[:,None]+30),cnt.sum(axis=1))
        for c in np.unique(chrom):indices[c]=np.flatnonzero(chrom==c)
        rng=np.random.default_rng([7,C.VIRUSES.index(virus)]);null=np.empty((20000,4));den=np.empty(20000)
        for done in range(0,20000,250):
            draws=start[None,:]+np.floor(rng.random((250,K))*size[None,:]).astype(np.int64)
            new=np.maximum(flat[draws][:,lab]+delta[None,:],0)
            pp=np.empty((250,len(co),4));density=np.empty((250,len(co)),np.int64)
            for c,ii in indices.items():
                if c not in maps:pp[:,ii,:]=pg;density[:,ii]=0;continue
                pm,dm=maps[c];bins=np.minimum(new[:,ii]//50000,len(pm)-1)
                pp[:,ii,:]=pm[bins];density[:,ii]=dm[bins]
            null[done:done+250]=(pp*weights[None,:,None]).sum(axis=1)
            den[done:done+250]=(density*weights[None,:]).sum(axis=1)
        pnull=arrays[f'{arm}__{virus}__class_null'];pdens=arrays[f'{arm}__{virus}__density_null']
        check(f'{variant} {virus} {arm} independent 20000 full class null',np.array_equal(null,pnull),float(np.max(abs(null-pnull))))
        check(f'{variant} {virus} {arm} independent 20000 full density null',np.array_equal(den,pdens),float(np.max(abs(den-pdens))))
        full.append(dict(variant=variant,virus=virus,arm=arm,n_coords=len(co),n_clusters=K,elapsed_seconds=time.time()-started,class_null_sha256=hashlib.sha256(null.tobytes()).hexdigest()))
        print(full[-1],flush=True)
# Concentration is the top-10 signed excess divided by the total signed excess.
conc=read(ROOT/'summary/top10_cluster_excess.csv')
for r in conc.itertuples():
    cl=read(ROOT/f'atlas/outputs/clusters_{r.variant}_{r.svclass}_{r.virus}.csv');ex=cl.obs-cl.null_mean
    vals=[ex.sum(),ex.nlargest(10).sum(),ex.nlargest(10).sum()/ex.sum() if ex.sum()>0 else np.nan]
    check(f'{r.variant} {r.virus} {r.svclass} top10 signed excess',np.allclose(vals,[r.total_excess,r.top10_excess,r.top10_share_of_excess],rtol=0,atol=1e-11,equal_nan=True))
result=dict(n_checks=len(checks),n_passed=sum(c['passed'] for c in checks),all_passed=all(c['passed'] for c in checks),checks=checks,full_recomputations=full,reviewed_files={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'axis/axis_recompute.py',ROOT/'axis/PROTOCOL.md',ROOT/'summary/top10_cluster_excess.csv']})
(HERE/'axis_result_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print({k:result[k] for k in ('n_checks','n_passed','all_passed')},flush=True)
if not result['all_passed']:print([c for c in checks if not c['passed']]);raise SystemExit(1)
