"""Assess covariate balance in the joint-matched spatial null.

Observed integration coordinates are compared with uniformly sampled eligible bins
and with actual rigid-block draws from the matched null. Standardized mean differences
summarize replication timing, timing gradient, gene density, fragile-site status and
mappability. Centromere and telomere distances provide additional diagnostics for
features not explicitly matched."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402
import joint_correction as J  # noqa: E402
import amplification_and_histology_control as A  # noqa: E402

OUT = J.OUT
BIN = F.BIN_SIZE
WINDOW = 100_000
SEED = 7

# (column, label, matched-on?)
COVS = [("rt_consensus", "replication timing (OK-seq consensus)", True),
        ("gene_density", "gene density (genes per 50-kb bin)", True),
        ("in_cfs", "fragile-site status (in common fragile site)", True),
        ("mapab", "mappability (CRG Align100)", True),
        ("grad", "replication-timing gradient", True),
        ("cen_dist", "centromere distance [NOT matched: negative control]", False),
        ("telo_dist", "telomere distance [NOT matched: negative control]", False)]


def smd(obs: np.ndarray, ref: np.ndarray) -> tuple[float, float, float]:
    o = obs[np.isfinite(obs)]
    r = ref[np.isfinite(ref)]
    if len(o) < 2 or len(r) < 2:
        return np.nan, np.nan, np.nan
    sd = np.sqrt((o.var(ddof=1) + r.var(ddof=1)) / 2.0)
    d = (o.mean() - r.mean()) / sd if sd > 0 else np.nan
    return float(o.mean()), float(r.mean()), float(d)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rep", type=int, default=200)
    ap.add_argument("--viruses", default="HBV,HPV,EBV")
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), 100_000)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)
    co = F.load_coordinates()

    # eligible = the bins the null can draw from at all: in the feature table, not amplicon-masked
    fe = feat.copy()
    fe["_masked"] = [(str(c), int(b) // BIN * BIN) in mask
                     for c, b in zip(fe["chrom"], fe["bin_start"])]
    elig = fe[~fe["_masked"]].reset_index(drop=True)
    print(f"[pool]   {len(feat):,} feature bins, {len(elig):,} eligible after the amplicon mask "
          f"({100*len(elig)/len(feat):.1f}%)", flush=True)
    key = pd.MultiIndex.from_arrays([feat["chrom"].astype(str).to_numpy(),
                                     feat["bin_start"].to_numpy(np.int64)])
    FT = feat.set_index(key)

    rows = []
    for virus in a.viruses.split(","):
        sub = co[co["virus_class"] == virus].sort_values(
            ["chrom", "pos0"]).reset_index(drop=True)
        in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                            for c, p in zip(sub["chrom"], sub["pos0"])])
        sub = sub[~in_mask].reset_index(drop=True)
        k = pd.MultiIndex.from_arrays([sub["chrom"].astype(str).to_numpy(),
                                       (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
        sub = sub.assign(stratum=nul.lookup.reindex(k).to_numpy())
        sub = sub[sub["stratum"].notna()].reset_index(drop=True)
        sub["stratum"] = sub["stratum"].astype(int)
        sub["within_bin_offset"] = sub["pos0"] - (sub["pos0"] // BIN) * BIN
        labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                   sub["pos0"].to_numpy(np.int64), WINDOW)
        n_clu = int(len(np.unique(labels)))

        # observed covariates, one row per coordinate
        okey = pd.MultiIndex.from_arrays([sub["chrom"].astype(str).to_numpy(),
                                          (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
        obs_f = FT.reindex(okey)

        # matched draws: the actual joint null, rigid blocks, n_rep replicates
        nul.n_calls = nul.n_fallback = 0
        rng = np.random.default_rng(SEED)
        df = sub.copy()
        df["_cluster"] = labels
        chroms, bins_ = [], []
        for _, g in df.groupby("_cluster", sort=True):
            c = str(g["chrom"].iloc[0])
            p = g["pos0"].to_numpy(np.int64)
            st = g["stratum"].to_numpy()
            off = g["within_bin_offset"].to_numpy()
            j = rng.integers(0, len(p), size=a.n_rep)
            delta = np.empty(a.n_rep, dtype=np.int64)
            for kk in range(len(p)):
                m = j == kk
                nk = int(m.sum())
                if nk:
                    delta[m] = nul.draw(rng, c, st[kk], off[kk], p[kk], nk) - p[kk]
            for pi in p:
                q = (pi + delta) // BIN * BIN
                chroms.append(np.full(a.n_rep, c))
                bins_.append(q)
        mkey = pd.MultiIndex.from_arrays([np.concatenate(chroms),
                                          np.concatenate(bins_).astype(np.int64)])
        mat_f = FT.reindex(mkey)
        fb = 100 * nul.n_fallback / max(nul.n_calls, 1)
        print(f"\n[{virus}]  {len(sub):,} coordinates, {n_clu:,} clusters, "
              f"{len(mat_f):,} matched draws ({a.n_rep} per coordinate), "
              f"pool fallback {fb:.2f}%")
        print(f"         {'covariate':<46}{'obs':>11}{'unmatched':>11}{'SMD':>8}"
              f"{'matched':>11}{'SMD':>8}  balance")
        for col, lab, matched_on in COVS:
            o = obs_f[col].to_numpy(float) if col != "in_cfs" else obs_f[col].to_numpy(float)
            u = elig[col].to_numpy(float)
            m = mat_f[col].to_numpy(float)
            mo, mu, d_u = smd(o, u)
            _, mm, d_m = smd(o, m)
            good = (abs(d_m) < 0.10) if np.isfinite(d_m) else False
            improved = (abs(d_m) < abs(d_u)) if np.isfinite(d_m) and np.isfinite(d_u) else False
            verdict = ("balanced" if good else "IMBALANCED") + ("" if matched_on else " (not matched)")
            print(f"         {lab:<46}{mo:>11.4g}{mu:>11.4g}{d_u:>+8.3f}"
                  f"{mm:>11.4g}{d_m:>+8.3f}  {verdict}")
            rows.append(dict(virus=virus, covariate=col, label=lab, matched_on=matched_on,
                             n_coords=len(sub), n_clusters=n_clu, n_draws=len(mat_f),
                             pct_pool_fallback=round(fb, 2),
                             mean_observed=round(mo, 6), mean_unmatched=round(mu, 6),
                             smd_unmatched=round(d_u, 4), mean_matched=round(mm, 6),
                             smd_matched=round(d_m, 4),
                             balanced_matched=bool(good), improved_by_matching=bool(improved)))

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "covariate_balance.csv", index=False)
    m = d[d["matched_on"]]
    print(f"\n=== SUMMARY ===")
    print(f"matched-on covariates: {int((m['smd_matched'].abs() < 0.10).sum())}/{len(m)} "
          f"balanced at |SMD| < 0.10 after matching; "
          f"{int((m['smd_unmatched'].abs() >= 0.10).sum())}/{len(m)} were imbalanced before it; "
          f"worst matched |SMD| {m['smd_matched'].abs().max():.3f}")
    nc = d[~d["matched_on"]]
    print(f"negative controls (not matched on): worst matched |SMD| "
          f"{nc['smd_matched'].abs().max():.3f}")
    print(f"\nwrote {OUT/'covariate_balance.csv'}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
