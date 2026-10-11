"""Estimate class-specific enrichment in prespecified 100-kb windows."""
from pathlib import Path
import sys
import argparse
import json
import time
import hashlib

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT.parent / "shared"))
import shared_analysis_context as common
import numpy as np
import pandas as pd
from atlas_core import make_index, full_summary, global_observed_union, compare_published, TrackedRNG

OUT = ROOT / "outputs"
CLASS_ORDER, VIRUS_ORDER = ("DUP", "DEL", "INV", "TRA"), ("HBV", "HPV", "EBV")
WINDOW, N_PERM, SEED = 100_000, 20_000, 7
WORKSPACE = ROOT.parents[2]
PUBLISHED = WORKSPACE / "reference_data/results/sv_class_locus_level.csv"


def output(name):
    path = (OUT / name).resolve()
    if not path.is_relative_to(ROOT):
        raise RuntimeError("refusing output outside atlas ownership")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def csv(frame, name):
    frame.to_csv(output(name), index=False, float_format="%.17g")


def select_coordinates(co, nul, mask, virus):
    # Direct safe port of N.eligible, rerun lookup to preserve original semantics even if co enriched.
    sub = co[co.virus_class == virus].sort_values(["chrom", "pos0"]).reset_index(drop=True)
    in_mask = np.array([(str(c), int(p)//50_000*50_000) in mask for c,p in zip(sub.chrom,sub.pos0)])
    sub = sub[~in_mask].reset_index(drop=True)
    key = pd.MultiIndex.from_arrays([sub.chrom.astype(str).to_numpy(), (sub.pos0.to_numpy(np.int64)//50_000)*50_000])
    sub = sub.assign(stratum=nul.lookup.reindex(key).to_numpy())
    sub = sub[sub.stratum.notna()].reset_index(drop=True)
    sub["stratum"] = sub.stratum.astype(int)
    sub["within_bin_offset"] = sub.pos0 - sub.pos0//50_000*50_000
    return sub


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", default=",".join(common.VARIANTS))
    args = parser.parse_args()
    variants = args.variants.split(",")
    if any(v not in common.VARIANTS for v in variants):
        raise ValueError("unknown variant")
    if variants[0] != "legacy_exact" and not output("legacy_replication_checks.csv").exists():
        raise RuntimeError("Exact reference check must pass first")
    context = common.load_context()
    F, nul, co, mask = [context[k] for k in ["F", "nul", "co", "mask"]]
    coords = {v:select_coordinates(co,nul,mask,v) for v in VIRUS_ORDER}
    labels = {v:F.leader_clusters(d.chrom.to_numpy(str),d.pos0.to_numpy(np.int64),WINDOW) for v,d in coords.items()}
    for virus, sub in coords.items():
        csv(sub.assign(cluster=labels[virus]), f"coordinates_{virus}.csv")
    rows, rounded_rows, diagnostics = [], [], []
    start = time.time()
    for variant in variants:
        vr, vround = [], []
        for cls in CLASS_ORDER:
            geometry,n_unique = common.get_variant_geometry(variant,cls)
            idx, actual_unique = make_index(geometry,F.ChromEvents)
            if actual_unique != n_unique:
                raise RuntimeError("adapter unique event count disagrees with geometry")
            for virus in VIRUS_ORDER:
                sub = coords[virus]
                rng = TrackedRNG(SEED)
                obs,null,per,sign_obs,sign_null = F.est_block(idx,nul,rng,sub,WINDOW,N_PERM,labels[virus],unit="uniqSV")
                cell = f"{variant}_{cls}_{virus}"
                np.savez_compressed(output(f"null_{cell}.npz"),block_null=null,sign_null=sign_null,
                    block_observed=np.int64(obs),sign_observed=np.int64(sign_obs))
                csv(per,f"clusters_{cell}.csv")
                union = global_observed_union(idx,per)
                diagnostics.append(dict(variant=variant,svclass=cls,virus=virus,n_unique=n_unique,
                    n_geometry_rows=len(geometry),n_coords=len(sub),n_clusters=len(per),
                    observed_cluster_event_incidences=obs,observed_catalogue_unique_events=union,
                    cross_cluster_repeat_incidences=obs-union,
                    rng_stream_sha256=rng.hexdigest(),rng_integers_calls=rng.n_calls))
                for est,o,nl in [("bc_uniqSV_blockNull",obs,null),("f_clusterSignTest",sign_obs,sign_null)]:
                    meta=dict(variant=variant,svclass=cls,n_events=int(n_unique),virus=virus,estimator=est,n_coords=len(sub),
                              n_clusters=len(per),seed=SEED,window_bp=WINDOW)
                    res=full_summary(o,nl)
                    res["nominal_upper_tail_p_lt_005"]=res["emp_p"]<.05
                    res["bonferroni_24_p_lt_threshold"]=res["emp_p"]<.05/24
                    res["emp_p_bonferroni_24"]=min(1.,24*res["emp_p"])
                    res["emp_p_two_sided_bonferroni_24"]=min(1.,24*res["emp_p_two_sided_equal_tail"])
                    res["two_sided_bonferroni_24_p_lt_threshold"]=res["emp_p_two_sided_equal_tail"]<.05/24
                    vr.append(dict(**meta,**res))
                    vround.append(dict(**meta,**F.summarize(o,nl)))
                print(f"[{cell}] n={n_unique} block fold={vr[-2]['fold']:.8f} p={vr[-2]['emp_p']:.8g}; sign fold={vr[-1]['fold']:.8f} p={vr[-1]['emp_p']:.8g}; elapsed={time.time()-start:.1f}s",flush=True)
        csv(pd.DataFrame(vr),f"atlas_{variant}_full_precision.csv")
        csv(pd.DataFrame(vround),f"atlas_{variant}_legacy_rounding.csv")
        if variant == "legacy_exact":
            published=pd.read_csv(PUBLISHED,float_precision="round_trip")
            checks=compare_published(pd.DataFrame(vround)[published.columns],published)
            csv(checks,"legacy_replication_checks.csv")
            if not checks.passed.all():
                raise RuntimeError("Class-atlas reference mismatch; alternative representations were not evaluated")
            print(f"Exact reference: {len(checks)} per-field comparisons passed",flush=True)
        rows.extend(vr); rounded_rows.extend(vround)
        csv(pd.DataFrame(rows),"atlas_all_variants_full_precision.csv")
        csv(pd.DataFrame(rounded_rows),"atlas_all_variants_legacy_rounding.csv")
        csv(pd.DataFrame(diagnostics),"event_geometry_and_union_diagnostics.csv")
        for virus in VIRUS_ORDER:
            streams = {r["rng_stream_sha256"] for r in diagnostics if r["virus"] == virus}
            if len(streams) != 1:
                raise RuntimeError(f"class/variant common random numbers mismatch for {virus}")
    final=dict(variants=variants,n_cells=len(rows)//2,n_estimator_rows=len(rows),n_perm=N_PERM,seed=SEED,
        window_bp=WINDOW,observed_only_global_union_is_diagnostic=True,
        block_estimand="sum over clusters of within-cluster unique event counts",
        significance_tail="one-sided upper",interval_type="observed divided by null percentiles",
        elapsed_seconds=time.time()-start,published_sha256=hashlib.sha256(PUBLISHED.read_bytes()).hexdigest())
    output("run_summary.json").write_text(json.dumps(final,indent=2),encoding="utf-8")


if __name__ == "__main__":
    main()
