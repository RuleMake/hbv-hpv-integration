"""Recalculates Hi-C resolution summaries while retaining zero-contact pairs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SRC = Path(r"source_project")
PAIR_OE = SRC / "results" / "supplementary" / "exp1_2_hic_pair_oe.csv"
OUT = Path(r"source_project") / "results" / "hic_power_regenerated.csv"


def power_at(log2oe, delta, rng, n_sim=2000, alpha=0.05):
    centred = log2oe - np.median(log2oe)
    hits = 0
    for _ in range(n_sim):
        s = rng.choice(centred, size=len(log2oe), replace=True) + delta
        try:
            p = stats.wilcoxon(s)[1]
        except ValueError:
            p = 1.0
        hits += p < alpha
    return hits / n_sim


def mde_for(l2, rng):
    for delta in np.arange(0.05, 1.51, 0.05):
        if power_at(l2, delta, rng) >= 0.80:
            return float(delta)
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()

    d = pd.read_csv(PAIR_OE)
    n_all, n_na, n_zero = len(d), int(d["oe"].isna().sum()), int((d["oe"] == 0).sum())
    print(f"[input]  {n_all} pairs: {n_na} missing O/E, {n_zero} with O/E exactly 0, "
          f"{int((d['oe'] > 0).sum())} positive")
    print(f"         the {n_zero} zero-contact pairs all sit at window "
          f"{sorted(d.loc[d['oe'] == 0, 'window_bp'].unique())}")

    rows = []
    for keep_zero in (True, False):
        base = d[d["oe"].notna()] if keep_zero else d[d["oe"].notna() & (d["oe"] > 0)]
        filt = "keep_zero_contact_pairs" if keep_zero else "drop_zero_contact_pairs_log2_required"
        for W, lab in ((1_000_000, "<=1Mb"), (5_000_000, "<=5Mb"), (None, "all_windows")):
            sub = base if W is None else base[base["window_bp"] <= W]
            if len(sub) < 3:
                continue
            oe = sub["oe"].to_numpy(float)
            med_oe = float(np.median(oe))
            donor_med = float(np.median(sub.groupby("donor")["oe"].median().to_numpy()))
            rng = np.random.default_rng(a.seed)
            boot = [np.median(rng.choice(oe, size=len(oe), replace=True)) for _ in range(3000)]
            lo, hi = np.percentile(boot, [2.5, 97.5])
            # the MDE is a log2-scale quantity, so it can only be computed on positive pairs
            pos = oe[oe > 0]
            l2 = np.log2(pos)
            mde = mde_for(l2, np.random.default_rng(a.seed)) if len(l2) >= 3 else None
            rows.append(dict(
                filter=filt, window=lab, n=len(sub), n_donors=int(sub["donor"].nunique()),
                median_oe=round(med_oe, 4), ci95_lo=round(float(lo), 4),
                ci95_hi=round(float(hi), 4),
                donor_aware_median_oe=round(donor_med, 4),
                n_positive_for_mde=len(pos),
                mde_log2_80pct=None if mde is None else round(mde, 2),
                mde_fold_80pct=None if mde is None else round(2 ** mde, 2)))
            print(f"  [{filt[:18]:18s} {lab:11s}] n={len(sub):>3}  median O/E {med_oe:.4f} "
                  f"[{lo:.3f}, {hi:.3f}]  donor-aware {donor_med:.4f}  "
                  f"MDE {'-' if mde is None else f'{2**mde:.2f}x'}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    pd.set_option("display.width", 250)
    print("\n" + df.to_string(index=False))

    k = df[(df["filter"].str.startswith("keep")) & (df["window"] == "all_windows")].iloc[0]
    p = df[(df["filter"].str.startswith("drop")) & (df["window"] == "all_windows")].iloc[0]
    print(f"\nIncluding zero-contact pairs: n = {k['n']}, median O/E = {k['median_oe']}. "
          f"Positive pairs used for the log2 analysis: n = {p['n']}, median O/E = {p['median_oe']}.")
    print(f"Donor-level estimates are closer ({k['donor_aware_median_oe']} vs "
          f"{p['donor_aware_median_oe']}) after summarizing contact values, including zeros, within donors.")
    print(f"\nSaved {OUT}")


if __name__ == "__main__":
    main()
