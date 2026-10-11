"""Summarizes Hi-C contact comparisons with donor-level aggregation.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SRC = Path(r"source_project")
PAIRS = SRC / "results" / "supplementary" / "exp1_2_hic_pair_oe.csv"
OUT = Path(r"source_project\results")
SEED = 7
NBOOT = 20000
STRATA = ((100_000, "<=100kb"), (1_000_000, "<=1Mb"), (5_000_000, "<=5Mb"))

# carried forward from the original run; needs the 22 GB .hic file to recompute
ORIGINAL_PERM_P = {"<=1Mb": 0.8580709645177411, "<=5Mb": 0.7581209395302348}
ORIGINAL_NULL_MEDIAN = {"<=1Mb": 0.990784924240707, "<=5Mb": 1.1191959144132813}
# Reference values used to check the pair-level calculation.
GATE = {"<=1Mb": (16, 0.8307), "<=5Mb": (47, 0.9957)}


def mde_fold(x: np.ndarray) -> float:
    """80% power, alpha 0.05, two-sided, on the log2 scale -> returned as a fold."""
    if len(x) < 2:
        return np.nan
    return float(2 ** ((1.96 + 0.84) * x.std(ddof=1) / np.sqrt(len(x))))


def main() -> None:
    argparse.ArgumentParser().parse_args()
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    d = pd.read_csv(PAIRS)
    d = d[np.isfinite(d["oe"])].copy()
    d["log2oe"] = np.log2(d["oe"].replace(0, np.nan))
    print(f"[in]      {len(d)} pairs with a defined O/E, {d.donor.nunique()} donors, "
          f"line {d.line.iloc[0]}")
    print(f"[in]      pairs per donor: {d.donor.value_counts().to_dict()}")
    print(f"[in]      {int((d.oe == 0).sum())} pair(s) with O/E exactly 0 "
          f"(dropped from log2-scale statistics only)")

    rows = []
    for thr, lab in STRATA:
        s = d[d.dist_bp <= thr]
        sl = s[np.isfinite(s.log2oe)]
        n, nd = len(s), int(s.donor.nunique())

        # Check pair-level summaries against the reference values.
        med = float(s.oe.median()) if n else np.nan
        gate = GATE.get(lab)
        drift = abs(med - gate[1]) if gate else np.nan
        ok = gate is not None and n == gate[0] and drift < 5e-4
        if gate:
            print(f"[check]    {lab}: n={n} (report {gate[0]}), median O/E {med:.4f} "
                  f"(report {gate[1]:.4f}, drift {drift:.5f}) -> {'PASS' if ok else 'FAIL'}")
        rows.append(dict(analysis="pair_level", stratum=lab, n_pairs=n, n_donors=nd,
                         median_oe=round(med, 4) if n else np.nan,
                         frac_oe_gt1=round(float((s.oe > 1).mean()), 4) if n else np.nan,
                         wilcoxon_p=round(float(stats.wilcoxon(sl.log2oe).pvalue), 5)
                         if len(sl) >= 3 else np.nan,
                         mde_fold_80pct=round(mde_fold(sl.log2oe.to_numpy(float)), 4)
                         if len(sl) >= 2 else np.nan,
                         original_perm_p=ORIGINAL_PERM_P.get(lab),
                         original_null_median_oe=(round(ORIGINAL_NULL_MEDIAN[lab], 4)
                                                  if lab in ORIGINAL_NULL_MEDIAN else np.nan),
                         reproduction_gate=("PASS" if ok else ("FAIL" if gate else "n/a"))))

        # Aggregate contacts within donors.
        if nd < 2 or len(sl) < 2:
            continue
        per = sl.groupby("donor").log2oe.median()
        k = len(per)
        n_up = int((per > 0).sum())
        sign_p = float(stats.binomtest(n_up, k, 0.5).pvalue)
        draws = np.empty(NBOOT)
        vals = per.to_numpy(float)
        for i in range(NBOOT):
            draws[i] = np.median(vals[rng.integers(0, k, k)])
        lo, hi = np.percentile(draws, [2.5, 97.5])
        rows.append(dict(analysis="donor_level", stratum=lab, n_donors=k,
                         n_pairs=int(len(sl)),
                         median_of_donor_medians_oe=round(float(2 ** np.median(vals)), 4),
                         boot_ci_lo_oe=round(float(2 ** lo), 4),
                         boot_ci_hi_oe=round(float(2 ** hi), 4),
                         n_donors_above_1=n_up,
                         donor_sign_test_p=round(sign_p, 5),
                         mde_fold_80pct=round(mde_fold(vals), 4)))
        print(f"[donor]   {lab}: {k} donors, median-of-medians O/E "
              f"{2**np.median(vals):.4f} [{2**lo:.4f},{2**hi:.4f}], "
              f"{n_up}/{k} above 1 (sign P={sign_p:.4f}), MDE {mde_fold(vals):.4f}")

    # Compare saved table values with the reference summaries.
    pw = SRC / "results" / "supplementary" / "exp1_2_power.csv"
    if pw.exists():
        p = pd.read_csv(pw)
        for _, r in p.iterrows():
            lab = str(r["window"]).strip()
            g = GATE.get(lab)
            if g and (int(r["n"]) != g[0] or abs(float(r["obs_median_oe"]) - g[1]) > 5e-3):
                rows.append(dict(analysis="stale_erratum", stratum=lab,
                                 stale_n=int(r["n"]), stale_median_oe=float(r["obs_median_oe"]),
                                 corrected_n=g[0], corrected_median_oe=g[1],
                                 note="exp1_2_power.csv still ships pre-erratum values; "
                                      "report md was corrected 2026-07-30, the CSV was not"))
                print(f"[stale]   exp1_2_power.csv {lab}: ships n={int(r['n'])}, "
                      f"median {float(r['obs_median_oe'])} -- corrected values are "
                      f"n={g[0]}, median {g[1]}")

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "hic_contact_donor_aware.csv", index=False)
    d.to_csv(OUT / "hic_pair_oe_source.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    for an in out.analysis.unique():
        s = out[out.analysis == an].dropna(axis=1, how="all")
        print(f"\n=== {an} ===")
        print(s.drop(columns=["analysis"]).to_string(index=False))
    print(f"\n=== per donor (<=5 Mb) ===")
    s = d[(d.dist_bp <= 5e6) & np.isfinite(d.log2oe)]
    print(s.groupby("donor").agg(n_pairs=("oe", "size"), median_oe=("oe", "median"),
                                 frac_gt1=("oe", lambda x: (x > 1).mean())).round(4).to_string())
    print(f"\nwrote {OUT/'hic_contact_donor_aware.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
