"""Combine genomic-covariate matching with recurrent-amplification masking.

The matching scheme combines 18 context strata with replication-timing-gradient
and mappability tiers. Pools below MIN_POOL = 5 fall back to the corresponding
18-stratum pool; fallback use is recorded. Amplification peaks are excluded from
both integration coordinates and null candidate positions."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402
import amplification_and_histology_control as A  # noqa: E402
from minimal_bigwig_reader import BigWig  # noqa: E402

SRC = Path(r"upstream_project")
OUT = Path(r"source_project\results")
BIN = F.BIN_SIZE
WINDOW = 100_000
MIN_POOL = 5
GRAD_TIERS, MAP_TIERS = 2, 3


def build_joint_features():
    feat = pd.read_parquet(F.FEATURES)
    feat = feat.assign(chrom=feat["chrom"].astype(str).map(
        lambda c: c if c.startswith("chr") else f"chr{c}"))
    feat = feat[(feat["stratum"] >= 0) & (feat["chrom"] != "chrY")]
    feat = feat.sort_values(["chrom", "bin_start"]).reset_index(drop=True)

    grads = []
    for _, g in feat.groupby("chrom", sort=False):
        v = g["rt_consensus"].to_numpy(float)
        d = np.full(len(v), np.nan)
        if len(v) > 2:
            d[1:-1] = (v[2:] - v[:-2]) / 2.0
        grads.append(pd.Series(np.abs(d), index=g.index))
    feat["grad"] = pd.concat(grads).reindex(feat.index)
    ok = feat["grad"].notna()
    gt = np.full(len(feat), GRAD_TIERS, dtype=np.int64)
    gt[ok.to_numpy()] = pd.qcut(feat.loc[ok, "grad"], GRAD_TIERS, labels=False,
                                duplicates="drop").to_numpy()

    bw = BigWig(SRC / "data/raw/annotation/wgEncodeCrgMapabilityAlign100mer.bigWig")
    v = bw.validate()
    if not v["passed"]:
        raise SystemExit("bigWig reader failed self-validation")
    binned, _, _ = bw.binned_means(BIN)
    bw.close()
    mp = pd.concat([pd.DataFrame(dict(chrom=(n if n.startswith("chr") else f"chr{n}"),
                                      bin_start=s, mapab=m))
                    for n, (s, m) in binned.items()], ignore_index=True).dropna()
    feat = feat.merge(mp, on=["chrom", "bin_start"], how="left")
    ok2 = feat["mapab"].notna()
    mt = np.full(len(feat), MAP_TIERS, dtype=np.int64)
    mt[ok2.to_numpy()] = pd.qcut(feat.loc[ok2, "mapab"], MAP_TIERS, labels=False,
                                 duplicates="drop").to_numpy()

    feat["stratum_ref"] = (feat["stratum"].astype(np.int64) * (GRAD_TIERS + 1) + gt) \
        * (MAP_TIERS + 1) + mt
    print(f"[joint]  bigWig self-validation rel_diff {v['rel_diff']:.2e}; "
          f"{feat['stratum_ref'].nunique()} joint strata")
    return feat


class JointNuller(F.Nuller):
    def __init__(self, feat, base, mask=None):
        self.base = base
        self.divisor = (GRAD_TIERS + 1) * (MAP_TIERS + 1)
        self.n_calls = self.n_fallback = 0
        f = feat
        if mask:
            keep = ~pd.MultiIndex.from_arrays(
                [f["chrom"].astype(str), f["bin_start"].astype(np.int64)]).isin(mask)
            removed = int((~keep).sum())
            f = f[keep]
            print(f"[joint]  amplification-peak mask removes {removed:,}/{len(feat):,} bins "
                  f"({100*removed/len(feat):.1f}%)")
        self.bins = {(str(c), int(s)): g["bin_start"].to_numpy(np.int64)
                     for (c, s), g in f.groupby(["chrom", "stratum_ref"])}
        self.chrom_size = {str(c): int(g["bin_start"].max()) + BIN
                           for c, g in feat.groupby("chrom")}
        self.lookup = pd.Series(
            f["stratum_ref"].to_numpy(np.int64),
            index=pd.MultiIndex.from_arrays([f["chrom"].to_numpy(str),
                                             f["bin_start"].to_numpy(np.int64)]))

    def pool(self, chrom, stratum, offset, pos0):
        self.n_calls += 1
        b = self.bins.get((str(chrom), int(stratum)))
        if b is not None and b.size:
            cand = np.unique(b + int(offset))
            cand = cand[cand != int(pos0)]
            if cand.size >= MIN_POOL:
                return cand
        self.n_fallback += 1
        return self.base.pool(chrom, int(stratum) // self.divisor, offset, pos0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    feat = build_joint_features()
    peaks = A.amp_peaks()
    mask = A.masked_bins(peaks, WINDOW)
    idx, _ = F.load_events()
    co = F.load_coordinates()
    base = F.Nuller()
    nul_joint = JointNuller(feat, base, mask)

    rows, t0 = [], time.time()
    for vc in ("HBV", "HPV", "EBV"):
        sub0 = co[co["virus_class"] == vc].sort_values(
            ["chrom", "pos0"]).reset_index(drop=True)
        in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                            for c, p in zip(sub0["chrom"], sub0["pos0"])])
        sub = sub0[~in_mask].reset_index(drop=True)
        k = pd.MultiIndex.from_arrays([sub["chrom"].astype(str).to_numpy(),
                                       (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
        sub = sub.assign(stratum=nul_joint.lookup.reindex(k).to_numpy())
        sub = sub[sub["stratum"].notna()].reset_index(drop=True)
        sub["stratum"] = sub["stratum"].astype(int)

        labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                   sub["pos0"].to_numpy(np.int64), WINDOW)
        n_clu = int(len(np.unique(labels)))
        print(f"\n[{vc}] {len(sub0)} -> {len(sub)} coordinates after amp-peak exclusion and "
              f"joint-stratum lookup; {n_clu} clusters")

        nul_joint.n_calls = nul_joint.n_fallback = 0
        rng = np.random.default_rng(a.seed)
        o, nl, _, s_obs, s_null = F.est_block(idx, nul_joint, rng, sub, WINDOW,
                                              a.n_perm, labels, unit="uniqSV")
        fb = round(100 * nul_joint.n_fallback / max(nul_joint.n_calls, 1), 2)
        for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                         ("f_clusterSignTest", F.summarize(s_obs, s_null))):
            rows.append(dict(virus=vc, estimator=est, n_coords=len(sub), n_clusters=n_clu,
                             pct_pool_fallback=fb, **res))
        print(f"  bc={rows[-2]['fold']} p={rows[-2]['emp_p']:.5g} | "
              f"sign={rows[-1]['fold']} p={rows[-1]['emp_p']:.5g} | fallback {fb}% "
              f"({time.time()-t0:.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "fig3_joint_correction.csv", index=False)
    pd.set_option("display.width", 240)
    print("\nJOINT (multiplicity + amp-peaks + RT gradient + mappability)")
    print(df[["virus", "estimator", "n_coords", "n_clusters", "observed", "expected",
              "fold", "fold_ci_lo", "fold_ci_hi", "emp_p", "pct_pool_fallback"]]
          .to_string(index=False))
    print("\nreference, 18 strata only: HBV bc 1.0580 / sign 1.0873 ; "
          "HPV bc 1.1157 / sign 1.1203 ; EBV bc 1.0120 / sign 1.1307 (ns)")
    print(f"\nwrote {OUT / 'fig3_joint_correction.csv'}")


if __name__ == "__main__":
    main()
