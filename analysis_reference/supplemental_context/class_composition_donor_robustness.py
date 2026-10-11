"""Computes sample- and study-removal sensitivity for the composition axis.
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
import locus_class_propensity as E              # noqa: E402

PROJ = Path(r"source_project")
SRC = Path(r"source_project")
OUT = PROJ / "results" / "class_composition_donor_robustness.csv"
OUT_G = PROJ / "results" / "class_composition_donor_robustness_groups.csv"
TOP_N = 25


def sample_map() -> pd.DataFrame:
    up = pd.read_csv(SRC / "data/processed/tian2022_viral_integrations.csv", low_memory=False)
    up = up[up["lifted"] & up["pos"].notna()].copy()
    up["chrom"] = "chr" + up["chrom"].astype(str).str.replace("chr", "", regex=False)
    up["pos0"] = up["pos"].astype(np.int64)
    return up[["chrom", "pos0", "sample", "data_source"]].drop_duplicates(["chrom", "pos0"])


def group_loo(obs_c: np.ndarray, null_c: np.ndarray, groups: np.ndarray):
    """Exact leave-one-GROUP-out for a weighted-sum statistic.

    Returns a frame with, per group, the clusters it owns, its contribution, and the excess that
    remains after dropping it (renormalised by n/(n-k) so the statistic stays a mean over clusters).
    """
    n = len(obs_c)
    s_obs = obs_c.sum()
    exp_k = null_c.mean(axis=0)
    s_exp = exp_k.sum()
    full = s_obs - s_exp
    rows = []
    for g in pd.unique(groups):
        m = groups == g
        k = int(m.sum())
        if k == 0 or k >= n:
            continue
        r = n / (n - k)
        loo = (s_obs - obs_c[m].sum() - (s_exp - exp_k[m].sum())) * r
        rows.append(dict(group=g, n_clusters=k,
                         contribution=float(obs_c[m].sum() - exp_k[m].sum()),
                         pct_of_full=100 * float(obs_c[m].sum() - exp_k[m].sum()) / full,
                         loo_excess=float(loo),
                         shrinkage_pct=100 * (1 - abs(loo) / abs(full)),
                         same_sign=bool(np.sign(loo) == np.sign(full))))
    return full, pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=5000)
    ap.add_argument("--window", type=int, default=1_000_000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    cum, nb, p_global = E.propensity_map()
    half = a.window // E.BIN
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), J.WINDOW)
    joint = J.JointNuller(feat, F.Nuller(), mask)
    dens = E.DensityNuller(joint, cum, nb, half)
    co = F.load_coordinates()
    smap = sample_map()

    rows, grows = [], []
    for vi, vc in enumerate(("HBV", "HPV")):
        _sub0, sub, tier, labels = E.prepare(co, vc, joint, dens, mask)
        sub = sub.merge(smap, on=["chrom", "pos0"], how="left")
        miss = int(sub["sample"].isna().sum())
        uniq, inv = np.unique(labels, return_inverse=True)
        n_clu = len(uniq)
        # cluster -> dominant sample / study
        tmp = pd.DataFrame({"cluster": inv, "sample": sub["sample"].to_numpy(),
                            "study": sub["data_source"].to_numpy()})
        dom_s = tmp.groupby("cluster")["sample"].agg(
            lambda s: s.value_counts().idxmax() if s.notna().any() else "UNKNOWN")
        dom_t = tmp.groupby("cluster")["study"].agg(
            lambda s: s.value_counts().idxmax() if s.notna().any() else "UNKNOWN")
        mixed = int(tmp.groupby("cluster")["sample"].nunique().gt(1).sum())
        print(f"\n[{vc}] {len(sub)} coords / {n_clu} clusters; "
              f"{sub['sample'].nunique()} samples, {sub['data_source'].nunique()} studies; "
              f"{mixed} clusters span >1 sample; {miss} coords unjoined")

        for arm, nuller, dt in (("stratum", joint, None), ("density", dens, tier)):
            rng = np.random.default_rng([a.seed, vi])
            *_, per_cluster = E.score(sub, labels, nuller, cum, nb, half, p_global, rng,
                                      a.n_perm, density_tier=dt, keep_per_cluster=True)
            obs_c, null_c = per_cluster
            null_c = null_c.astype(np.float64)

            for level, g in (("sample", dom_s.to_numpy()), ("study", dom_t.to_numpy())):
                full, df = group_loo(obs_c, null_c, g)
                if df.empty:
                    continue
                df = df.sort_values("shrinkage_pct", ascending=False)
                worst = df.iloc[0]
                rows.append(dict(
                    virus=vc, arm=arm, level=level, n_groups=len(df), n_clusters=n_clu,
                    full_excess=round(full, 5),
                    worst_group=str(worst["group"]), worst_group_clusters=int(worst["n_clusters"]),
                    worst_loo_excess=round(float(worst["loo_excess"]), 5),
                    worst_shrinkage_pct=round(float(worst["shrinkage_pct"]), 2),
                    worst_pct_of_full=round(float(worst["pct_of_full"]), 2),
                    n_same_sign=int(df["same_sign"].sum()),
                    all_same_sign=bool(df["same_sign"].all())))
                df.insert(0, "arm", arm); df.insert(0, "level", level); df.insert(0, "virus", vc)
                grows.append(df.head(TOP_N))
                print(f"  [{arm} | {level:6s}] full {full:+.5f}; worst drop "
                      f"{worst['group']} ({int(worst['n_clusters'])} clusters) -> "
                      f"{worst['loo_excess']:+.5f} ({worst['shrinkage_pct']:.1f}% shrink); "
                      f"same sign {int(df['same_sign'].sum())}/{len(df)}", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    pd.concat(grows, ignore_index=True).to_csv(OUT_G, index=False)
    pd.set_option("display.width", 250)
    print("\n" + d.to_string(index=False))
    print(f"\nwrote {OUT}\nwrote {OUT_G}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
