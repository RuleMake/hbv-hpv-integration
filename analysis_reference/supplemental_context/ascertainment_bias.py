"""Examines integration read support, local copy number and detection-related covariates.
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

SRC_ROOT = Path(r"source_project")
sys.path.insert(0, str(SRC_ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import mannwhitneyu, spearmanr  # noqa: E402

from src.common.config import (  # noqa: E402
    BIN_SIZE, chrom_norm, HG19_CHROM_SIZES, LARGE_TD_MIN, LARGE_TD_MAX,
)
from src.common.genome_features import load_features  # noqa: E402

OUT = Path(r"source_project\results")
CLCA = SRC_ROOT / "data" / "interim" / "clca"
WIN = 1_000_000
N_PERM = 2000
SEED = 7


def parse_intervals(s):
    out = []
    if not isinstance(s, str):
        return out
    for tok in s.split(","):
        m = re.match(r"^(chr[\dXYM]+):(\d+)-(\d+)$", tok.strip())
        if m:
            out.append((chrom_norm(m.group(1)), int(m.group(2)), int(m.group(3))))
    return out


def main() -> None:
    feat = load_features()
    pools = {}
    for c, g in feat[feat["stratum"] >= 0].groupby("chrom"):
        pools[c] = {int(x): sg["bin_start"].to_numpy() for x, sg in g.groupby("stratum")}

    sv = pd.read_csv(CLCA / "sv_api.csv")
    sv["c1"] = sv["chra"].map(chrom_norm)
    sv["c2"] = sv["chrb"].map(chrom_norm)
    sv = sv[sv["c1"] == sv["c2"]].copy()
    sv["span"] = (sv["posb"] - sv["posa"]).abs()
    sv["sv_class"] = ["DUP" if (a == "-" and b == "+") else
                      "DEL" if (a == "+" and b == "-") else
                      "INV" if a == b else "OTHER"
                      for a, b in zip(sv["strand_a"], sv["strand_b"])]
    svL = sv[(sv["span"] >= LARGE_TD_MIN) & (sv["span"] <= LARGE_TD_MAX)]

    ints = pd.read_csv(CLCA / "integrations.csv")
    ints["chrom"] = ints["chrom"].astype(str).map(
        lambda c: chrom_norm(f"chr{c}" if not str(c).startswith("chr") else c))
    ints["pos"] = ints["pos"].astype(int)
    ints["bin_start"] = (ints["pos"] // BIN_SIZE * BIN_SIZE).astype(int)
    ints = ints.merge(feat[["chrom", "bin_start", "stratum"]], on=["chrom", "bin_start"], how="left")
    ints = ints[ints["stratum"] >= 0].copy()
    ints["support"] = ints["split_reads"].fillna(0) + ints["pair_reads"].fillna(0)

    print("=" * 78)
    print("ASCERTAINMENT-BIAS CONTROL")
    print("=" * 78)
    print(f"integrations with stratum: {len(ints)} in {ints['sample'].nunique()} samples")

    # ---------------- TEST 1: amplicon membership
    amp = pd.read_csv(CLCA / "amplicons.csv")
    amp_by_sample: dict[str, list] = {}
    amp_copy: dict[str, list] = {}
    for _, r in amp.iterrows():
        iv = parse_intervals(r["intervals"])
        amp_by_sample.setdefault(r["sample"], []).extend(iv)
        for x in iv:
            amp_copy.setdefault(r["sample"], []).append((x, r["avg_copy"]))

    in_amp, amp_cn = [], []
    for _, r in ints.iterrows():
        hit, cn = False, np.nan
        for (c, s, e), v in amp_copy.get(r["sample"], []):
            if c == r["chrom"] and s <= r["pos"] <= e:
                hit, cn = True, v
                break
        in_amp.append(hit)
        amp_cn.append(cn)
    ints["in_amplicon"] = in_amp
    ints["amp_avg_copy"] = amp_cn

    a = ints.loc[ints.in_amplicon, "support"]
    b = ints.loc[~ints.in_amplicon, "support"]
    print("\n--- TEST 1: does amplicon membership predict read support? ---")
    print(f"  inside amplicon : n={len(a):4d}  median support={a.median():.1f}  mean={a.mean():.1f}")
    print(f"  outside amplicon: n={len(b):4d}  median support={b.median():.1f}  mean={b.mean():.1f}")
    if len(a) > 5 and len(b) > 5:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        print(f"  Mann-Whitney p={p:.4g}  -> "
              f"{'SUPPORT DEPENDS ON COPY NUMBER (ascertainment risk)' if p < 0.05 else 'no significant dependence'}")
    if ints["amp_avg_copy"].notna().sum() > 10:
        sub = ints.dropna(subset=["amp_avg_copy"])
        rho, p = spearmanr(sub["amp_avg_copy"], sub["support"])
        print(f"  support vs amplicon avg_copy: Spearman rho={rho:.3f}, p={p:.4g} (n={len(sub)})")

    # ---------------- TEST 2: local SV density vs support
    bymap_all = {}
    for smp, g in sv.groupby("caseid"):
        m = {}
        for _, t in g.iterrows():
            m.setdefault(t["c1"], []).extend([int(t["posa"]), int(t["posb"])])
        bymap_all[smp] = {c: np.array(sorted(v)) for c, v in m.items()}

    dens = []
    for _, r in ints.iterrows():
        arr = bymap_all.get(r["sample"], {}).get(r["chrom"])
        if arr is None or arr.size == 0:
            dens.append(0)
        else:
            dens.append(int(np.searchsorted(arr, r["pos"] + WIN, "right")
                            - np.searchsorted(arr, r["pos"] - WIN, "left")))
    ints["local_sv_density"] = dens
    rho, p = spearmanr(ints["local_sv_density"], ints["support"])
    print("\n--- TEST 2: does local rearrangement density predict read support? ---")
    print(f"  Spearman rho={rho:.3f}, p={p:.4g} (n={len(ints)})")
    print(f"  -> {'DENSITY PREDICTS SUPPORT (ascertainment risk)' if p < 0.05 and rho > 0 else 'no positive dependence'}")

    # Enrichment stratified by junction-read support.
    ints["support_tertile"] = pd.qcut(ints["support"], 3,
                                      labels=["low", "mid", "high"], duplicates="drop")
    bymapL = {}
    for smp, g in svL.groupby("caseid"):
        d = {}
        for cls, gg in g.groupby("sv_class"):
            m = {}
            for _, t in gg.iterrows():
                m.setdefault(t["c1"], []).extend([int(t["posa"]), int(t["posb"])])
            d[cls] = {c: np.array(sorted(v)) for c, v in m.items()}
        bymapL[smp] = d

    def n_near(arr, pos, win=WIN):
        if arr is None or arr.size == 0:
            return np.zeros(np.shape(pos), dtype=int)
        return (np.searchsorted(arr, np.asarray(pos) + win, "right")
                - np.searchsorted(arr, np.asarray(pos) - win, "left"))

    print("\n--- Enrichment by read-support tertile, +/-1 Mb ---")
    print("  If ascertainment drives the signal, fold should RISE with support.")
    rows = []
    for tert in ["low", "mid", "high"]:
        sub = ints[ints.support_tertile == tert]
        sup_range = f"{sub['support'].min():.0f}-{sub['support'].max():.0f}"
        print(f"\n  tertile={tert} (n={len(sub)} sites, support {sup_range})")
        for cls in ("DUP", "DEL", "INV"):
            rng = np.random.default_rng(SEED)
            obs = 0
            null = np.zeros(N_PERM, dtype=np.int64)
            for _, r in sub.iterrows():
                arr = bymapL.get(r["sample"], {}).get(cls, {}).get(r["chrom"])
                if arr is None or arr.size == 0:
                    continue
                obs += int(n_near(arr, np.array([r["pos"]]))[0])
                pool = pools.get(r["chrom"], {}).get(int(r["stratum"]))
                if pool is None or pool.size == 0:
                    mpos = rng.integers(0, HG19_CHROM_SIZES[r["chrom"]], size=N_PERM)
                else:
                    mpos = (pool[rng.integers(0, pool.size, size=N_PERM)]
                            + rng.integers(0, BIN_SIZE, size=N_PERM))
                null += n_near(arr, mpos)
            exp = float(null.mean())
            p_ = (int(np.sum(null >= obs)) + 1) / (N_PERM + 1)
            fold = obs / exp if exp > 0 else np.nan
            rows.append({"tertile": tert, "sv_class": cls, "n_sites": len(sub),
                         "observed": obs, "exp": round(exp, 2),
                         "fold": round(fold, 3), "emp_p": round(p_, 5)})
            print(f"    {cls}: obs={obs:5d} exp={exp:8.2f} fold={fold:7.3f} p={p_:.5f}")

    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "ascertainment_bias.csv", index=False)
    ints.to_csv(OUT / "ascertainment_integration_annot.csv", index=False)

    print("\n" + "=" * 78)
    print("Fold enrichment across support tertiles")
    print("=" * 78)
    for cls in ("DUP", "DEL", "INV"):
        k = df[df.sv_class == cls].set_index("tertile")["fold"]
        lo, mid, hi = (k.get(t, np.nan) for t in ("low", "mid", "high"))
        trend = "RISES with support (ascertainment risk)" if hi > lo * 1.5 else \
                "flat / no upward trend (ascertainment NOT driving it)"
        print(f"  {cls}: low={lo:6.2f}  mid={mid:6.2f}  high={hi:6.2f}   -> {trend}")
    print(f"\nWrote {OUT / 'ascertainment_bias.csv'}")


if __name__ == "__main__":
    main()
