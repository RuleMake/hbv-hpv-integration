"""Estimate local class composition in 1-Mb windows under four breakend representations.

Shared adapters supply the input tables and spatial-null routines. Run with Python -B."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
# Make the reviewed shared adapter's dependency directory available before pandas
# caches pyarrow availability and version.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "structural_variant_catalogue/event_ledger/runtime_deps"))
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
SHARED = ROOT.parent / "shared"
BIN, WINDOW, HALF = 50_000, 1_000_000, 20
CLUSTER_RADIUS, SHRINK, MIN_POOL = 100_000, 30.0, 5
CLASSES = ("DEL", "DUP", "INV", "TRA")
VIRUSES = ("HBV", "HPV", "EBV")
VARIANTS = ("legacy_exact", "legacy_expanded", "paired_sample", "donor_signature")
SEED, N_PERM, CHUNK = 7, 20_000, 250


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(table, path):
    table.to_csv(path, index=False, float_format="%.17g", na_rep="")


def breakpoint_fingerprint(bp):
    """Fingerprint the endpoint MULTISET, including identical repeated endpoints."""
    canonical = bp[["chrom", "pos0", "svclass"]].sort_values(["chrom", "pos0", "svclass"], kind="stable")
    payload = canonical.to_csv(index=False, lineterminator="\n").encode()
    return hashlib.sha256(payload).hexdigest()


def propensity_map(bp):
    if bp.empty or not set(bp.svclass).issubset(CLASSES):
        raise ValueError("Empty or unknown-class breakpoint input")
    if (bp.pos0 < 0).any():
        raise ValueError("Negative breakpoint positions")
    p_global = bp.svclass.value_counts().reindex(CLASSES, fill_value=0).to_numpy(float)/len(bp)
    cum, nb = {}, {}
    for chrom, group in bp.groupby("chrom"):
        n = int(group.pos0.max()) // BIN + 2
        a = np.zeros((4, n), dtype=np.int64)
        for k, cls in enumerate(CLASSES):
            bins = group.loc[group.svclass == cls, "pos0"].to_numpy(np.int64)//BIN
            np.add.at(a[k], bins, 1)
        cum[chrom] = np.concatenate([np.zeros((4,1),np.int64),np.cumsum(a,axis=1)],axis=1)
        nb[chrom] = n
    return cum, nb, p_global


def counts_at(cum, nb, chrom, pos, half=HALF):
    out = np.zeros((len(pos),4),dtype=np.int64)
    bins = np.asarray(pos,dtype=np.int64)//BIN
    for c in np.unique(chrom):
        mask = chrom == c
        cc, n = cum.get(str(c)), nb.get(str(c))
        if cc is None:
            continue
        lo = np.clip(bins[mask]-half,0,n)
        hi = np.clip(bins[mask]+half+1,0,n)
        out[mask] = (cc[:,hi]-cc[:,lo]).T
    return out


def props(counts, p_global):
    return (counts + SHRINK*p_global)/(counts.sum(axis=1,keepdims=True)+SHRINK)


def lookup_maps(cum, nb, p_global):
    # Positions beyond n+HALF have no breakends in the original clipped window.
    maps = {}
    for chrom, n in nb.items():
        b = np.arange(n+HALF+1,dtype=np.int64)
        counts = counts_at(cum,nb,np.repeat(chrom,len(b)),b*BIN)
        maps[chrom] = (props(counts,p_global), counts.sum(axis=1))
    return maps


class DensityNuller:
    def __init__(self, joint, cum, nb):
        self.joint = joint
        idx = joint.lookup.index
        chrom = idx.get_level_values(0).to_numpy(str)
        start = idx.get_level_values(1).to_numpy(np.int64)
        dens = counts_at(cum,nb,chrom,start+BIN//2).sum(axis=1)
        tier = pd.qcut(pd.Series(dens).rank(method="first"),10,labels=False).to_numpy(np.int64)
        self.tier = pd.Series(tier,index=idx)
        key = pd.DataFrame({"chrom":chrom,"stratum":joint.lookup.to_numpy(np.int64),"tier":tier,"bin_start":start})
        self.bins = {(c,s,t):g.bin_start.to_numpy(np.int64) for (c,s,t),g in key.groupby(["chrom","stratum","tier"])}
        self.n_calls = self.n_fallback = 0
        self.membership = key.assign(breakpoint_density=dens)

    def pool(self, chrom, stratum, offset, pos0, tier=None):
        self.n_calls += 1
        bins = self.bins.get((str(chrom),int(stratum),int(tier)))
        if bins is not None and bins.size:
            cand = np.unique(bins+int(offset))
            cand = cand[cand != int(pos0)]
            if cand.size >= MIN_POOL:
                return cand
        self.n_fallback += 1
        return self.joint.pool(chrom,stratum,offset,pos0)


def leader_indices(labels):
    uniq, inv = np.unique(labels,return_inverse=True)
    leader = np.zeros(len(uniq),dtype=np.int64)
    leader[inv[::-1]] = np.arange(len(inv))[::-1]
    return inv, leader


def strict_mask(sub, labels, dens, tier):
    inv, lead = leader_indices(labels)
    ok = np.zeros(len(lead),dtype=bool)
    pool_sizes = np.zeros(len(lead),dtype=int)
    for k,i in enumerate(lead):
        r = sub.iloc[i]
        if not np.isfinite(tier[i]):
            continue
        b = dens.bins.get((str(r.chrom),int(r.stratum),int(tier[i])))
        if b is not None:
            cand = np.unique(b+int(r.pos0)%BIN)
            cand = cand[cand != int(r.pos0)]
            pool_sizes[k] = len(cand)
            ok[k] = len(cand) >= MIN_POOL
    return ok[inv], ok, pool_sizes


def score(sub, labels, nuller, cum, nb, pg, rng, n_perm, density_tier=None, verify_lookup=False):
    chrom, pos = sub.chrom.to_numpy(str), sub.pos0.to_numpy(np.int64)
    strat = sub.stratum.to_numpy(np.int64)
    inv, lead = leader_indices(labels)
    n_clusters = len(lead)
    size = np.bincount(inv,minlength=n_clusters)
    delta = pos-pos[lead][inv]
    weight = 1.0/(n_clusters*size[inv])
    cnt = counts_at(cum,nb,chrom,pos)
    pr = props(cnt,pg)
    obs = (pr*weight[:,None]).sum(axis=0)
    obs_density = float((cnt.sum(axis=1)*weight).sum())
    pools=[]
    for k,i in enumerate(lead):
        kw = {} if density_tier is None else {"tier":density_tier[i]}
        pools.append(nuller.pool(str(chrom[i]),int(strat[i]),int(pos[i])%BIN,int(pos[i]),**kw))
    sizes=np.array([len(p) for p in pools],dtype=np.int64)
    if (sizes==0).any():
        raise ValueError("Empty candidate pool")
    starts=np.concatenate([[0],np.cumsum(sizes)[:-1]])
    flat=np.concatenate(pools).astype(np.int64)
    maps=lookup_maps(cum,nb,pg)
    coordinate_chrom=chrom[lead][inv]
    chrom_indices={c:np.flatnonzero(coordinate_chrom==c) for c in np.unique(coordinate_chrom)}
    null=np.empty((n_perm,4),dtype=float)
    null_density=np.empty(n_perm,dtype=float)
    max_lookup_error=0.0
    for done in range(0,n_perm,CHUNK):
        count=min(CHUNK,n_perm-done)
        idx=starts[None,:]+(rng.random((count,n_clusters))*sizes[None,:]).astype(np.int64)
        newpos=flat[idx][:,inv]+delta[None,:]
        np.maximum(newpos,0,out=newpos)
        pp=np.empty((count,len(pos),4),dtype=float)
        density=np.empty((count,len(pos)),dtype=np.int64)
        for chrom_name,ix in chrom_indices.items():
            if chrom_name not in maps:
                pp[:,ix,:]=pg
                density[:,ix]=0
            else:
                prob_map,den_map=maps[chrom_name]
                b=np.clip(newpos[:,ix]//BIN,0,len(prob_map)-1)
                pp[:,ix,:]=prob_map[b]
                density[:,ix]=den_map[b]
        if verify_lookup and done==0:
            slow=counts_at(cum,nb,np.tile(coordinate_chrom,count),newpos.ravel())
            slow_prob=props(slow,pg).reshape(count,len(pos),4)
            if not np.array_equal(pp,slow_prob) or not np.array_equal(density,slow.sum(axis=1).reshape(count,len(pos))):
                raise AssertionError("Fast bin lookup differs from original counts_at arithmetic")
            max_lookup_error=float(np.max(np.abs(pp-slow_prob)))
        null[done:done+count]=(pp*weight[None,:,None]).sum(axis=1)
        null_density[done:done+count]=(density*weight[None,:]).sum(axis=1)
    detail=sub[["tian_coordinate_id","chrom","pos0"]].copy() if "tian_coordinate_id" in sub else sub[["chrom","pos0"]].copy()
    detail["cluster"]=labels
    detail["coordinate_weight"]=weight
    detail["observed_axis"]=pr[:,2:4].sum(axis=1)-pr[:,0:2].sum(axis=1)
    detail["observed_breakpoint_density"]=cnt.sum(axis=1)
    return obs,obs_density,null,null_density,detail,max_lookup_error


def summarize(obs, null, name, signed=False):
    obs=float(obs)
    exp,sd=float(null.mean()),float(null.std(ddof=1))
    lo,hi=np.percentile(null,[2.5,97.5])
    n=len(null)
    n_hi,n_lo=int(np.sum(null>=obs)),int(np.sum(null<=obs))
    p_hi,p_lo=(n_hi+1)/(n+1),(n_lo+1)/(n+1)
    return {"statistic":name,"observed":obs,"null_mean":exp,"null_reference_q025":float(lo),"null_reference_q975":float(hi),"interval_type":"randomization_reference_not_sampling_CI","excess":obs-exp,"fold":np.nan if signed or exp<=0 else obs/exp,"null_sd":sd,"z":(obs-exp)/sd if sd>0 else np.nan,"n_null_ge_observed":n_hi,"n_null_le_observed":n_lo,"emp_p_greater":p_hi,"emp_p_less":p_lo,"emp_p_two_sided":min(1.,2*min(p_hi,p_lo)),"mde_excess_approx_2p80sd":2.8*sd,"mde_fold_approx":np.nan if signed or exp<=0 else (exp+2.8*sd)/exp,"clears_legacy_mde":abs(obs-exp)>=2.8*sd,"mde_interpretation":"approximate_sensitivity_not_effect_upper_bound"}


def axis_values(obs,null):
    return float(obs[[3,2]].sum()-obs[[0,1]].sum()), null[:,[3,2]].sum(axis=1)-null[:,[0,1]].sum(axis=1)


def compare_legacy(actual_class,actual_axis,actual_strict):
    results=[]
    references=(("class",actual_class,"legacy_composition_class.csv"),("axis",actual_axis,"legacy_composition_axis.csv"),("strict",actual_strict,"legacy_strict_matching_strict.csv"))
    mapping={"null_ci_lo":"null_reference_q025","null_ci_hi":"null_reference_q975","mde_excess":"mde_excess_approx_2p80sd","mde_fold":"mde_fold_approx","clears_mde":"clears_legacy_mde"}
    rounding={"observed":5,"null_mean":5,"null_ci_lo":5,"null_ci_hi":5,"excess":5,"fold":4,"z":2,"mde_excess":5,"mde_fold":4,"null_stream_corr":4,"null_fallback_pct":2}
    for category,table,filename in references:
        ref=pd.read_csv(ROOT/"provenance"/filename,float_precision="round_trip")
        ref=ref[ref.window_bp==WINDOW]
        keys=["virus","arm","statistic"]+(["rung"] if category=="strict" else [])
        for r in ref.to_dict("records"):
            match=table
            for key in keys:
                match=match[match[key]==r[key]]
            if len(match)!=1:
                raise ValueError(f"Reference row is not unique {category}:{r}")
            got=match.iloc[0]
            for metric,expected in r.items():
                if metric in keys or metric not in got.index and metric not in mapping or pd.isna(expected):
                    continue
                actual=got[mapping.get(metric,metric)]
                if metric == "mde_excess":
                    # Stored results were corrected from K=2.486 *after* the
                    # original cell had been rounded, rather than rerun.
                    actual=round(round(2.486*float(got.null_sd),5)*(2.80/2.486),5)
                elif metric == "mde_fold":
                    exp,sd=float(got.null_mean),float(got.null_sd)
                    old=round((exp+2.486*sd)/exp,4)
                    actual=round(1+(old-1)*(2.80/2.486),5)
                elif metric in rounding:
                    actual=round(float(actual),rounding[metric])
                if isinstance(expected,(float,int,np.number)) and not isinstance(expected,(bool,np.bool_)):
                    passed=bool(np.isclose(actual,expected,rtol=0,atol=1e-13))
                else:
                    passed=bool(actual==expected)
                results.append({"table":category,"virus":r["virus"],"arm":r["arm"],"rung":r.get("rung",""),"statistic":r["statistic"],"metric":metric,"legacy":expected,"recomputed_at_legacy_precision":actual,"passed":passed})
    return pd.DataFrame(results)


def run_variant(common,context,variant,n_perm=N_PERM):
    out=OUT/variant
    out.mkdir(parents=True,exist_ok=True)
    bp=common.load_breakpoints(variant,exclude_positive=True)
    cum,nb,pg=propensity_map(bp)
    joint=context["nul"]
    co=context["co_joint"]
    F=context["F"]
    dens=DensityNuller(joint,cum,nb)
    write_csv(dens.membership,out/"density_bin_membership.csv")
    composition=pd.DataFrame({"svclass":CLASSES,"n_breakpoints":bp.svclass.value_counts().reindex(CLASSES,fill_value=0).to_numpy(),"global_fraction":pg})
    write_csv(composition,out/"background_composition.csv")
    row,axrow,strictrow,draws,null_arrays,checks=[],[],[],{}, {}, []
    t0=time.time()
    for vi,vc in enumerate(VIRUSES):
        sub=co[co.virus_class==vc].sort_values(["chrom","pos0"]).reset_index(drop=True)
        labels=F.leader_clusters(sub.chrom.to_numpy(str),sub.pos0.to_numpy(np.int64),CLUSTER_RADIUS)
        key=pd.MultiIndex.from_arrays([sub.chrom.to_numpy(str),(sub.pos0.to_numpy(np.int64)//BIN)*BIN])
        tier=dens.tier.reindex(key).to_numpy()
        if np.isnan(tier).any():
            raise ValueError("Missing density tier")
        keep_coord,keep_clu,pool_sizes=strict_mask(sub,labels,dens,tier)
        pool_audit=pd.DataFrame({"cluster":np.arange(len(keep_clu)),"strict_pool_size":pool_sizes,"strict":keep_clu})
        write_csv(pool_audit,out/f"{vc}_strict_pool_audit.csv")
        for arm,nuller,dt in (("stratum",joint,None),("density",dens,tier),("density_strict",dens,tier[keep_coord])):
            s=sub if arm!="density_strict" else sub[keep_coord].reset_index(drop=True)
            lb=labels if arm!="density_strict" else F.leader_clusters(s.chrom.to_numpy(str),s.pos0.to_numpy(np.int64),CLUSTER_RADIUS)
            nuller.n_calls=nuller.n_fallback=0
            ob,od,nl,nd,detail,err=score(s,lb,nuller,cum,nb,pg,np.random.default_rng([SEED,vi]),n_perm,density_tier=dt,verify_lookup=True)
            fallback=100*nuller.n_fallback/max(nuller.n_calls,1)
            base={"variant":variant,"virus":vc,"arm":"density" if arm=="density_strict" else arm,"window_bp":WINDOW,"n_coords":len(s),"n_clusters":len(np.unique(lb)),"n_perm":n_perm,"null_fallback_pct":fallback,"rung":"strict" if arm=="density_strict" else "full"}
            oax,nax=axis_values(ob,nl)
            summary={**base,**summarize(oax,nax,"class_axis_TRAINV_minus_DELDUP",True)}
            draws[(arm,vc)]=(oax-nax.mean(),nax-nax.mean())
            null_arrays[f"{arm}__{vc}__class_null"]=nl
            null_arrays[f"{arm}__{vc}__density_null"]=nd
            null_arrays[f"{arm}__{vc}__observed_class"]=ob
            null_arrays[f"{arm}__{vc}__axis_null"]=nax
            checks.append({"virus":vc,"arm":arm,"fast_lookup_max_absolute_error":err,"lookup_identical":err==0,"weights_sum":detail.coordinate_weight.sum()})
            write_csv(detail,out/f"{vc}_{arm}_observed_coordinates.csv")
            if arm!="density_strict":
                axrow.append(summary)
                for k,cls in enumerate(CLASSES):
                    row.append({**base,**summarize(ob[k],nl[:,k],"prop_"+cls)})
                row.append({**base,**summarize(od,nd,"breakpoint_density")})
            if arm!="stratum":
                strictrow.append(summary)
            print(f"{variant} {vc} {arm}: ncl={base['n_clusters']} axis_excess={summary['excess']:.8f} p2={summary['emp_p_two_sided']:.7g} fallback={fallback:.4f}% ({time.time()-t0:.0f}s)",flush=True)
    for arm in ("stratum","density","density_strict"):
        (eh,nh),(ep,npv)=draws[(arm,"HBV")],draws[(arm,"HPV")]
        summary={"variant":variant,"virus":"HBV-vs-HPV","arm":"density" if arm=="density_strict" else arm,"rung":"strict" if arm=="density_strict" else "full","window_bp":WINDOW,"n_coords":np.nan,"n_clusters":np.nan,"n_perm":n_perm,"null_fallback_pct":np.nan,"null_stream_corr":float(np.corrcoef(nh,npv)[0,1]),**summarize(eh-ep,nh-npv,"interaction_axis_HBV_minus_HPV",True)}
        null_arrays[f"{arm}__interaction_null"]=nh-npv
        if arm!="density_strict":
            axrow.append(summary)
        if arm!="stratum":
            strictrow.append(summary)
    tables=[pd.DataFrame(row),pd.DataFrame(axrow),pd.DataFrame(strictrow)]
    for table,name in zip(tables,("class_composition.csv","axis.csv","strict_sensitivity.csv")):
        write_csv(table,out/name)
    write_csv(pd.DataFrame(checks),out/"lookup_arithmetic_checks.csv")
    np.savez_compressed(out/"null_draws.npz",**null_arrays)
    atomic_aliquots={a for value in bp.aliquot.dropna().astype(str) for a in value.split('|') if a} if "aliquot" in bp else set()
    metadata={"variant":variant,"window_bp":WINDOW,"bin_size":BIN,"grid_window_bins":41,"shrink":SHRINK,"seed":SEED,"n_perm":n_perm,"n_background_breakpoints":len(bp),"n_background_donors":int(bp.donor.nunique()) if "donor" in bp else None,"n_background_aliquots":len(atomic_aliquots),"aliquot_count_definition":"atomic union of pipe-delimited contributor aliquots","endpoint_multiset_sha256":breakpoint_fingerprint(bp),"protocol_sha256":sha(ROOT/"PROTOCOL.md"),"script_sha256":sha(__file__),"shared_adapter_sha256":sha(SHARED/"shared_analysis_context.py"),"versions":{"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__}}
    (out/"run.json").write_text(json.dumps(metadata,indent=2),encoding="utf-8")
    if variant=="legacy_exact":
        audit=compare_legacy(*tables)
        write_csv(audit,out/"legacy_reproduction.csv")
        passed=bool(audit.passed.all())
        (out/"legacy_gate.json").write_text(json.dumps({"passed":passed,"n_checks":len(audit),"n_passed":int(audit.passed.sum())},indent=2))
        if not passed:
            print(audit[~audit.passed].to_string(index=False),flush=True)
            raise AssertionError("Class-composition reference check failed")
    return metadata


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--variant",choices=VARIANTS+('all',),default="legacy_exact")
    args=parser.parse_args()
    sys.path.insert(0,str(SHARED))
    import shared_analysis_context as common
    context=common.load_context()
    variants=VARIANTS if args.variant=="all" else (args.variant,)
    if variants[0]!="legacy_exact":
        gate=OUT/"legacy_exact/legacy_gate.json"
        if not gate.exists() or not json.loads(gate.read_text())["passed"]:
            raise ValueError("Run and reproduce legacy_exact first")
    for variant in variants:
        run_variant(common,context,variant)


if __name__=="__main__":
    main()
