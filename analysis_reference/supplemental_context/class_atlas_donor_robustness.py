"""Computes sample- and study-removal sensitivity for the class atlas.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F                    # noqa: E402
import joint_correction as J                   # noqa: E402
import amplification_and_histology_control as A    # noqa: E402
import discriminating_tests as N                # noqa: E402
import sv_class_lolo as L                           # noqa: E402
from class_composition_donor_robustness import sample_map     # noqa: E402

PROJ = Path(r"source_project")
OUT = PROJ / "results" / "class_atlas_donor_robustness.csv"
OUT_G = PROJ / "results" / "class_atlas_donor_robustness_groups.csv"
WINDOW = 100_000
CELLS = [("DUP", "HPV"), ("DEL", "HPV"), ("TRA", "HBV"), ("INV", "HBV"), ("TRA", "HPV")]
TOP_N = 15


def loo_folds(o_tot, per_obs, per_null, groups):
    """Exact leave-one-GROUP-out folds. per_null is a list of (n_perm,) arrays, one per cluster."""
    n_perm = len(per_null[0])
    null_stack = np.vstack(per_null)                       # (n_clusters, n_perm)
    n_tot = null_stack.sum(axis=0)                         # (n_perm,)
    full = o_tot / n_tot.mean()
    rows = []
    for g in pd.unique(groups):
        m = groups == g
        k = int(m.sum())
        if k == 0 or k >= len(per_obs):
            continue
        o = o_tot - per_obs[m].sum()
        nn = (n_tot - null_stack[m].sum(axis=0)).mean()
        f = o / nn if nn else np.nan
        rows.append(dict(group=g, n_clusters=k,
                         obs_dropped=int(per_obs[m].sum()),
                         pct_of_excess=100 * float(per_obs[m].sum() - null_stack[m].sum(axis=0).mean())
                         / max(o_tot - n_tot.mean(), 1e-9),
                         fold_loo=float(f),
                         shrinkage_pct=100 * (1 - (f - 1) / (full - 1)) if full != 1 else np.nan,
                         still_above_1=bool(f > 1)))
    return full, pd.DataFrame(rows), n_perm


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    nul = J.JointNuller(feat, F.Nuller(), mask)
    co = F.load_coordinates()
    canon = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
    led = N.build_ledger(set(canon["aliquot_metadata"].astype(str)))
    smap = sample_map()

    rows, grows = [], []
    for cls, vc in CELLS:
        idx, n_ev = N.to_index(led[led.svclass == cls])
        sub = N.eligible(co, nul, mask, vc).merge(smap, on=["chrom", "pos0"], how="left")
        labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                   sub["pos0"].to_numpy(np.int64), WINDOW)
        rng = np.random.default_rng(a.seed)
        o_tot, _n_tot, per_obs, per_null, meta = L.est_block_uniqsv_percluster(
            idx, nul, rng, sub, WINDOW, a.n_perm, labels)

        uniq, inv = np.unique(labels, return_inverse=True)
        tmp = pd.DataFrame({"cluster": inv, "sample": sub["sample"].to_numpy(),
                            "study": sub["data_source"].to_numpy()})
        dom_s = tmp.groupby("cluster")["sample"].agg(
            lambda s: s.value_counts().idxmax() if s.notna().any() else "UNKNOWN").to_numpy()
        dom_t = tmp.groupby("cluster")["study"].agg(
            lambda s: s.value_counts().idxmax() if s.notna().any() else "UNKNOWN").to_numpy()
        mixed = int(tmp.groupby("cluster")["sample"].nunique().gt(1).sum())
        print(f"\n[{cls} {vc}] {len(sub)} coords / {len(per_obs)} clusters, "
              f"{sub['sample'].nunique()} samples, {sub['data_source'].nunique()} studies; "
              f"{mixed} clusters span >1 sample", flush=True)

        for level, g in (("sample", dom_s), ("study", dom_t)):
            full, df, n_perm = loo_folds(o_tot, per_obs, per_null, g)
            if df.empty:
                continue
            df = df.sort_values("fold_loo")            # worst = smallest remaining fold
            worst = df.iloc[0]
            rows.append(dict(
                svclass=cls, virus=vc, level=level, n_groups=len(df),
                n_clusters=len(per_obs), n_perm=n_perm, fold_full=round(full, 4),
                worst_group=str(worst["group"]), worst_group_clusters=int(worst["n_clusters"]),
                fold_loo_min=round(float(df["fold_loo"].min()), 4),
                fold_loo_max=round(float(df["fold_loo"].max()), 4),
                worst_shrinkage_pct=round(float(worst["shrinkage_pct"]), 2),
                worst_pct_of_excess=round(float(worst["pct_of_excess"]), 2),
                n_still_above_1=int(df["still_above_1"].sum()),
                all_still_above_1=bool(df["still_above_1"].all())))
            df.insert(0, "level", level); df.insert(0, "virus", vc); df.insert(0, "svclass", cls)
            grows.append(df.head(TOP_N))
            print(f"  [{level:6s}] full {full:.4f}; worst drop {worst['group']} "
                  f"({int(worst['n_clusters'])} clusters, {worst['pct_of_excess']:.0f}% of excess) "
                  f"-> {worst['fold_loo']:.4f}; still >1 in "
                  f"{int(df['still_above_1'].sum())}/{len(df)}", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    pd.concat(grows, ignore_index=True).to_csv(OUT_G, index=False)
    pd.set_option("display.width", 250)
    print("\n" + d.to_string(index=False))
    print(f"\nwrote {OUT}\nwrote {OUT_G}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
