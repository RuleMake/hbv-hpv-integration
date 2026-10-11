# Replication-timing gradient and mappability matching for spatial-null pools.
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import pandas as pd
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WS = ROOT.parents[1]
DATA01 = WS / 'analysis_inputs/structural_variant_catalogue/donor_duplication_statistics/data'
sys.path.insert(0, str(WS / 'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
import spatial_null_core as F
from bigwig_reader import BigWig
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

    bw = BigWig(ROOT / "data/wgEncodeCrgMapabilityAlign100mer.bigWig")
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