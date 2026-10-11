"""Estimate host structural-variant class composition around integration loci.

The propensity map excludes donors with called viral integrations. Within a window,
the class fraction is (count_C + a*p_C) / (n + a), with a = 30; sparse windows are
therefore shrunk toward the genome-wide class composition. Values are evaluated
on the 50-kb grid used for matched placement.

Coordinate clusters receive equal weight. The joint-stratum null preserves each
cluster's internal offsets, with an additional analysis matched on local breakpoint-
density decile. Outputs include class-specific estimates, signed class-composition
contrasts, the HBV-HPV interaction and breakpoint-density diagnostics."""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F              # noqa: E402
import joint_correction as J             # noqa: E402
import amplification_and_histology_control as A   # noqa: E402

PROJ = Path(r"source_project")
SRC = Path(r"upstream_project")
LEDGER = PROJ / "results" / "cache" / "sv_class_ledger.parquet"
COHORT = SRC / "data" / "processed" / "cohort.csv"
OUT = PROJ / "results" / "locus_class_propensity.csv"
OUT_AX = PROJ / "results" / "locus_class_propensity_axes.csv"

BIN = F.BIN_SIZE                 # 50 kb, the grid the null moves on
CLUSTER_RADIUS = 100_000
CLASSES = ["DEL", "DUP", "INV", "TRA"]
SHRINK = 30.0
DENSITY_TIERS = 10
MIN_POOL = 5
# One signed class axis, oriented a priori from the class atlas (sv_class_locus_level.csv):
# HBV's classes minus HPV's. HBV coordinates should sit above their null on it, HPV below, EBV on.
# Scoring two mirrored axes would be redundant -- they are exact negations of each other.
AXIS = (("TRA", "INV"), ("DEL", "DUP"))
WINDOWS = (200_000, 500_000, 1_000_000, 2_000_000)   # 1 Mb is primary; the rest are a scale profile


# ---------------------------------------------------------------- virus-free propensity map
def propensity_map():
    """Per chromosome, cumulative virus-negative breakpoint counts by class on the 50-kb grid.

    Returns (cum, n_bins, p_global) where cum[chrom] has shape (4, n_bins+1), so counts in a
    bin range [i, j) are cum[:, j] - cum[:, i] and any window size can be read off the same map.
    """
    cohort = pd.read_csv(COHORT)
    drop = {os.path.basename(str(p)).split(".")[0] for p in cohort["sv_path"]}

    d = pd.read_parquet(LEDGER)
    n_all, al_all = len(d), set(d["aliquot"].unique())
    d = d[~d["aliquot"].isin(drop)]
    print(f"[map]    virus-negative {len(d):,}/{n_all:,} events "
          f"({d['aliquot'].nunique()} aliquots); {n_all - len(d):,} events from the "
          f"{len(drop & al_all)} integration-positive donors removed")

    tra = d["svclass"] == "TRA"
    bp = pd.concat([
        d[["chrom", "svclass"]].assign(pos=d["left_pos0"].to_numpy(np.int64)),
        d.loc[~tra, ["chrom", "svclass"]].assign(pos=d.loc[~tra, "right_pos0"].to_numpy(np.int64)),
    ], ignore_index=True)

    p_global = (bp["svclass"].value_counts().reindex(CLASSES).to_numpy(float)
                / len(bp))
    print(f"[map]    {len(bp):,} breakpoints; genome composition "
          + ", ".join(f"{c} {p:.3f}" for c, p in zip(CLASSES, p_global)))

    cum, nb = {}, {}
    for c, g in bp.groupby("chrom"):
        n = int(g["pos"].max()) // BIN + 2
        a = np.zeros((len(CLASSES), n), dtype=np.int64)
        for k, cls in enumerate(CLASSES):
            gg = g.loc[g["svclass"] == cls, "pos"].to_numpy(np.int64) // BIN
            if gg.size:
                np.add.at(a[k], gg, 1)
        cum[c] = np.concatenate([np.zeros((len(CLASSES), 1), np.int64),
                                 np.cumsum(a, axis=1)], axis=1)
        nb[c] = n
    return cum, nb, p_global


def counts_at(cum, nb, half, chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """(n_positions, 4) virus-negative breakpoint counts in +/-window of each position."""
    out = np.zeros((len(pos), len(CLASSES)), dtype=np.int64)
    b = pos // BIN
    for c in np.unique(chrom):
        m = chrom == c
        cc, n = cum.get(str(c)), nb.get(str(c))
        if cc is None:
            continue
        lo = np.clip(b[m] - half, 0, n)
        hi = np.clip(b[m] + half + 1, 0, n)
        out[m] = (cc[:, hi] - cc[:, lo]).T
    return out


def props(cnt: np.ndarray, p_global: np.ndarray) -> np.ndarray:
    tot = cnt.sum(axis=1, keepdims=True)
    return (cnt + SHRINK * p_global) / (tot + SHRINK)


def set_shrink(v: float) -> None:
    """Set shrinkage toward the genome-wide class composition.
    
    Comparisons at low shrinkage assess how strongly the regularization of sparse
    windows affects the estimated association with local breakpoint density."""
    global SHRINK
    SHRINK = float(v)


# ---------------------------------------------------------------- density-matched null
class DensityNuller:
    """Joint-stratum pool further restricted to the same local-density decile."""

    def __init__(self, joint, cum, nb, half):
        self.joint = joint
        idx = joint.lookup.index
        chrom = idx.get_level_values(0).to_numpy(str)
        start = idx.get_level_values(1).to_numpy(np.int64)
        dens = counts_at(cum, nb, half, chrom, start + BIN // 2).sum(axis=1)
        tier = pd.qcut(pd.Series(dens).rank(method="first"), DENSITY_TIERS,
                       labels=False).to_numpy(np.int64)
        self.tier = pd.Series(tier, index=idx)
        strat = joint.lookup.to_numpy(np.int64)
        key = pd.DataFrame({"chrom": chrom, "stratum": strat, "tier": tier,
                            "bin_start": start})
        self.bins = {(c, s, t): g["bin_start"].to_numpy(np.int64)
                     for (c, s, t), g in key.groupby(["chrom", "stratum", "tier"])}
        self.n_calls = self.n_fallback = 0
        print(f"[dens]   {DENSITY_TIERS} density deciles x joint strata -> "
              f"{len(self.bins):,} pools")

    def pool(self, chrom, stratum, offset, pos0, tier=None):
        self.n_calls += 1
        b = self.bins.get((str(chrom), int(stratum), int(tier)))
        if b is not None and b.size:
            cand = np.unique(b + int(offset))
            cand = cand[cand != int(pos0)]
            if cand.size >= MIN_POOL:
                return cand
        self.n_fallback += 1
        return self.joint.pool(chrom, stratum, offset, pos0)


# ---------------------------------------------------------------- estimator
def score(sub: pd.DataFrame, labels: np.ndarray, nuller, cum, nb, half,
          p_global, rng, n_perm: int, density_tier=None, chunk: int = 250,
          keep_per_cluster: bool = False):
    """Observed and null cluster-mean propensity vectors.

    Mean-over-clusters-of-mean-over-members collapses exactly to one weighted sum over
    coordinates with weight 1/(n_clusters * cluster_size), which is what makes the permutation
    affordable: no per-permutation aggregation is needed, only a dot product.

    Returns (obs (4,), obs_density, null (n_perm,4), null_density (n_perm,), per_cluster) where
    per_cluster is None unless `keep_per_cluster`, in which case it carries the observed and null
    class-axis contribution of every cluster for leave-one-cluster-out.
    """
    chrom = sub["chrom"].to_numpy(str)
    pos = sub["pos0"].to_numpy(np.int64)
    strat = sub["stratum"].to_numpy(np.int64)

    uniq, inv = np.unique(labels, return_inverse=True)
    n_clu = len(uniq)
    size = np.bincount(inv, minlength=n_clu)
    # cluster leader = its first coordinate (sub is sorted by chrom, pos); members keep their
    # offset to it, so the null moves each cluster rigidly rather than scattering its members
    lead_i = np.zeros(n_clu, dtype=np.int64)
    lead_i[inv[::-1]] = np.arange(len(inv))[::-1]
    delta = pos - pos[lead_i][inv]
    w = 1.0 / (n_clu * size[inv])                       # exact cluster-mean-of-means weights

    cnt = counts_at(cum, nb, half, chrom, pos)
    pr = props(cnt, p_global)
    obs = (pr * w[:, None]).sum(axis=0)
    obs_dens = float((cnt.sum(axis=1) * w).sum())

    # per-leader null pools, flattened so every permutation is a single vectorised draw
    pools = []
    for k in range(n_clu):
        i = lead_i[k]
        kw = {} if density_tier is None else {"tier": density_tier[i]}
        pools.append(nuller.pool(str(chrom[i]), int(strat[i]), int(pos[i]) % BIN,
                                 int(pos[i]), **kw))
    sizes = np.array([p.size for p in pools], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    flat = np.concatenate(pools).astype(np.int64)
    lead_chrom_by_coord = chrom[lead_i][inv]

    own_i = [CLASSES.index(c) for c in AXIS[0]]
    oth_i = [CLASSES.index(c) for c in AXIS[1]]
    null = np.empty((n_perm, len(CLASSES)))
    null_dens = np.empty(n_perm)
    per_cluster_null = np.empty((n_perm, n_clu), np.float32) if keep_per_cluster else None

    done = 0
    while done < n_perm:
        c = min(chunk, n_perm - done)
        idx = starts[None, :] + (rng.random((c, n_clu)) * sizes[None, :]).astype(np.int64)
        newpos = flat[idx][:, inv] + delta[None, :]
        np.maximum(newpos, 0, out=newpos)
        cc = counts_at(cum, nb, half, np.tile(lead_chrom_by_coord, c), newpos.ravel())
        cc = cc.reshape(c, len(pos), len(CLASSES))
        pp = props(cc.reshape(-1, len(CLASSES)), p_global).reshape(c, len(pos), len(CLASSES))
        null[done:done + c] = (pp * w[None, :, None]).sum(axis=1)
        null_dens[done:done + c] = (cc.sum(axis=2) * w[None, :]).sum(axis=1)
        if keep_per_cluster:
            ax = (pp[:, :, own_i].sum(axis=2) - pp[:, :, oth_i].sum(axis=2)) * w[None, :]
            # bincount per permutation, not np.add.at on a 2-D fancy index: same result, ~20x faster
            per_cluster_null[done:done + c] = np.stack(
                [np.bincount(inv, weights=ax[i], minlength=n_clu) for i in range(c)]
            ).astype(np.float32)
        done += c
        if done % 5000 < chunk:
            print(f"          perm {done}/{n_perm}", flush=True)

    per_cluster = None
    if keep_per_cluster:
        ax_obs = (pr[:, own_i].sum(axis=1) - pr[:, oth_i].sum(axis=1)) * w
        agg = np.zeros(n_clu)
        np.add.at(agg, inv, ax_obs)
        per_cluster = (agg, per_cluster_null)
    return obs, obs_dens, null, null_dens, per_cluster


def summarize(obs: float, null: np.ndarray, name: str, *, signed: bool = False) -> dict:
    """`signed` statistics (the class axis, the interaction) can straddle zero, so a fold is
    undefined for them and the MDE is reported as an absolute displacement instead."""
    obs = float(obs)
    exp, sd = float(null.mean()), float(null.std(ddof=1))
    lo, hi = np.percentile(null, [2.5, 97.5])
    n = len(null)
    p_hi = (int(np.sum(null >= obs)) + 1) / (n + 1)
    p_lo = (int(np.sum(null <= obs)) + 1) / (n + 1)
    mde_abs = 2.80 * sd            # two-sided alpha 0.05, 80% power (1.96 + 0.84); see mde_constant_correction.py
    return {
        "statistic": name,
        "observed": round(obs, 5), "null_mean": round(exp, 5),
        "null_ci_lo": round(float(lo), 5), "null_ci_hi": round(float(hi), 5),
        "excess": round(obs - exp, 5),
        "fold": np.nan if signed or exp <= 0 else round(obs / exp, 4),
        "z": round((obs - exp) / sd, 2) if sd > 0 else np.nan,
        "emp_p_greater": p_hi, "emp_p_less": p_lo,
        "emp_p_two_sided": min(1.0, 2 * min(p_hi, p_lo)),
        "mde_excess": round(mde_abs, 5),
        "mde_fold": np.nan if signed or exp <= 0 else round((exp + mde_abs) / exp, 4),
        "clears_mde": bool(abs(obs - exp) >= mde_abs),
    }


def prepare(co, vc, joint, dens, mask):
    """Coordinates for one virus: amplicon-masked, joint-stratum-labelled, leader-clustered."""
    sub0 = co[co["virus_class"] == vc].sort_values(["chrom", "pos0"]).reset_index(drop=True)
    in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                        for c, p in zip(sub0["chrom"], sub0["pos0"])])
    sub = sub0[~in_mask].reset_index(drop=True)
    key = pd.MultiIndex.from_arrays([sub["chrom"].astype(str).to_numpy(),
                                     (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
    sub = sub.assign(stratum=joint.lookup.reindex(key).to_numpy())
    keep = sub["stratum"].notna().to_numpy()
    sub = sub[keep].reset_index(drop=True)
    sub["stratum"] = sub["stratum"].astype(int)
    tier = dens.tier.reindex(pd.MultiIndex.from_arrays(
        [sub["chrom"].astype(str).to_numpy(),
         (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])).to_numpy()
    labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                               sub["pos0"].to_numpy(np.int64), CLUSTER_RADIUS)
    return sub0, sub, tier, labels


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--windows", type=int, nargs="*", default=list(WINDOWS))
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--shrink", type=float, default=SHRINK,
                    help="shrinkage alpha; set near 0 to test the shrinkage-bias sensitivity")
    a = ap.parse_args()
    t0 = time.time()
    set_shrink(a.shrink)
    # Write non-default shrinkage results separately to preserve the primary outputs.
    global OUT, OUT_AX
    if abs(SHRINK - 30.0) > 1e-9:
        tag = f"_shrink{SHRINK:g}"
        OUT = OUT.with_name(OUT.stem + tag + ".csv")
        OUT_AX = OUT_AX.with_name(OUT_AX.stem + tag + ".csv")
    print(f"[cfg]    shrinkage alpha = {SHRINK}  ->  {OUT.name}")

    cum, nb, p_global = propensity_map()
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), J.WINDOW)
    co = F.load_coordinates()
    joint = J.JointNuller(feat, F.Nuller(), mask)

    own_i = [CLASSES.index(c) for c in AXIS[0]]
    oth_i = [CLASSES.index(c) for c in AXIS[1]]
    rows, axrows = [], []

    for window in a.windows:
        half = window // BIN
        dens = DensityNuller(joint, cum, nb, half)
        axis_draws = {}                       # (arm, virus) -> (obs_excess, null axis draws)
        for vi, vc in enumerate(("HBV", "HPV", "EBV")):
            sub0, sub, tier, labels = prepare(co, vc, joint, dens, mask)
            n_clu = int(len(np.unique(labels)))
            print(f"\n[{window//1000} kb | {vc}] {len(sub0)} -> {len(sub)} coordinates, "
                  f"{n_clu} clusters")

            for arm, nuller, dt in (("stratum", joint, None), ("density", dens, tier)):
                # Seed per virus, not globally. The interaction test differences one virus's null
                # draws from the other's, and a shared stream would correlate them: pool bins are
                # position-sorted, so the same uniform draw lands at a similarly-ranked genomic
                # position in both viruses' pools. Correlated nulls shrink Var(difference) and
                # make the interaction p anti-conservative. Independent streams remove the path.
                rng = np.random.default_rng([a.seed, vi])
                nuller.n_calls = nuller.n_fallback = 0
                obs, obs_d, null, null_d, _ = score(sub, labels, nuller, cum, nb, half,
                                                    p_global, rng, a.n_perm, density_tier=dt)
                base = dict(virus=vc, arm=arm, window_bp=window, n_coords=len(sub),
                            n_clusters=n_clu, n_perm=a.n_perm,
                            null_fallback_pct=round(100 * nuller.n_fallback
                                                    / max(nuller.n_calls, 1), 2))
                for ci, cls in enumerate(CLASSES):
                    rows.append({**base, **summarize(obs[ci], null[:, ci], f"prop_{cls}")})
                rows.append({**base, **summarize(obs_d, null_d, "breakpoint_density")})

                o_ax = obs[own_i].sum() - obs[oth_i].sum()
                n_ax = null[:, own_i].sum(axis=1) - null[:, oth_i].sum(axis=1)
                axrows.append({**base, **summarize(o_ax, n_ax, "class_axis_TRAINV_minus_DELDUP",
                                                   signed=True)})
                axis_draws[(arm, vc)] = (o_ax - n_ax.mean(), n_ax - n_ax.mean())
                print(f"  [{arm}] " + "  ".join(
                    f"{c} {obs[i]:.4f}/{null[:, i].mean():.4f}"
                    for i, c in enumerate(CLASSES))
                    + f"   axis {o_ax - n_ax.mean():+.5f}"
                    + f"   fallback {base['null_fallback_pct']}%")

        # crossed interaction: the two viruses' displacements on the same axis, opposite signs.
        # HBV and HPV are disjoint coordinate sets scored against independent null draws, so
        # pairing draw t of one with draw t of the other is a valid joint null.
        for arm in ("stratum", "density"):
            (eh, nh), (ep, npv) = axis_draws[(arm, "HBV")], axis_draws[(arm, "HPV")]
            r = float(np.corrcoef(nh, npv)[0, 1])   # must be ~0, else the pairing is not a null
            axrows.append({
                "virus": "HBV-vs-HPV", "arm": arm, "window_bp": window,
                "n_coords": np.nan, "n_clusters": np.nan, "n_perm": a.n_perm,
                "null_fallback_pct": np.nan, "null_stream_corr": round(r, 4),
                **summarize(eh - ep, nh - npv, "interaction_axis_HBV_minus_HPV", signed=True)})
            print(f"  [{arm}] interaction {eh - ep:+.5f}   null-stream corr {r:+.4f}")

    df, dfa = pd.DataFrame(rows), pd.DataFrame(axrows)
    df.to_csv(OUT, index=False)
    dfa.to_csv(OUT_AX, index=False)
    cols = ["window_bp", "virus", "arm", "statistic", "observed", "null_mean", "fold",
            "excess", "z", "emp_p_greater", "emp_p_two_sided", "clears_mde"]
    print("\n" + df[df.statistic != "breakpoint_density"][cols].to_string(index=False))
    print("\n" + df[df.statistic == "breakpoint_density"][cols].to_string(index=False))
    print("\n" + dfa[cols].to_string(index=False))
    print(f"\nwrote {OUT}\nwrote {OUT_AX}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
