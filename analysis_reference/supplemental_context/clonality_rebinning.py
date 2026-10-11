"""Rebins integration clonality by distance to the nearest structural rearrangement.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))

import timing_clonality_order as TC  # noqa: E402

OUT = TC.OUT
BANDS = ["JUNCTION", "SEAM", "FAR", "NONE"]

PUBLISHED_GATE = {  # band -> (n, n_subclonal)
    "JUNCTION": (9, 1),
    "SEAM": (62, 9),
    "FAR": (28, 15),
    "NONE": (19, 12),
}


def published_band(d: pd.DataFrame, w: int = 100) -> pd.Series:
    """The shipped rule: the OUTERMOST occupied shell."""
    j0 = d[f"n_sv_{w}k_j0"] > 0
    j10 = d[f"n_sv_{w}k_j10000"] > 0
    j50 = d[f"n_sv_{w}k_j50000"] > 0
    return pd.Series(np.select([j50, j10 & ~j50, j0 & ~j10], ["FAR", "SEAM", "JUNCTION"],
                               default="NONE"), index=d.index)


def nearest_band(d: pd.DataFrame, w: int = 100) -> pd.Series:
    """The sensitivity rule: the shell holding the NEAREST rearrangement."""
    n0 = d[f"n_sv_{w}k_j0"]
    n10 = d[f"n_sv_{w}k_j10000"]
    n50 = d[f"n_sv_{w}k_j50000"]
    in10 = (n0 - n10) > 0            # something within 10 kb
    in10_50 = (n10 - n50) > 0        # nothing within 10 kb, something within 50
    in50_100 = n50 > 0
    return pd.Series(np.select([in10, in10_50 & ~in10, in50_100 & ~in10_50 & ~in10],
                               ["JUNCTION", "SEAM", "FAR"], default="NONE"), index=d.index)


def profile(d: pd.DataFrame, band: pd.Series, rule: str) -> tuple[list[dict], float, float]:
    rows = []
    for b in BANDS:
        s = d[band == b]
        n = len(s)
        nsub = int((s.clonality == "Subclonal").sum())
        rows.append(dict(rule=rule, band=b, n=n, n_donors=int(s.donor.nunique()),
                         n_subclonal=nsub,
                         pct_subclonal=round(100.0 * nsub / n, 1) if n else np.nan,
                         median_ccf=round(float(s.ccf.median()), 4) if n else np.nan))
    tab = np.array([[r["n_subclonal"], r["n"] - r["n_subclonal"]] for r in rows if r["n"] > 0])
    chi2, p, dof, exp = stats.chi2_contingency(tab)
    return rows, float(chi2), float(p), float(exp.min())


def main() -> int:
    t0 = time.time()
    d = pd.read_csv(TC.CIS)
    d["clonality"] = d["clonality"].astype(str)
    d = d[d["clonality"].isin(["Clonal", "Subclonal"])].reset_index(drop=True)
    print(f"[in]      {len(d)} breakpoints, {d.donor.nunique()} donors")

    pub = published_band(d)
    for b, (n, nsub) in PUBLISHED_GATE.items():
        got_n = int((pub == b).sum())
        got_s = int(((pub == b) & (d.clonality == "Subclonal")).sum())
        if (got_n, got_s) != (n, nsub):
            print(f"[CHECK FAILED] published band {b}: got n={got_n} sub={got_s}, expected {n}/{nsub}")
            return 1
    print("[check]    Reference counts match in every distance band")

    rows_pub, chi_pub, p_pub, exp_pub = profile(d, pub, "published_outermost_shell")
    near = nearest_band(d)
    rows_near, chi_near, p_near, exp_near = profile(d, near, "sensitivity_nearest_rearrangement")

    # the rule-invariant contrast, which is what the paragraph actually rests on
    prox = d[f"n_sv_100k_j0"] > 0
    tab = np.array([[int((prox & (d.clonality == "Subclonal")).sum()),
                     int((prox & (d.clonality == "Clonal")).sum())],
                    [int((~prox & (d.clonality == "Subclonal")).sum()),
                     int((~prox & (d.clonality == "Clonal")).sum())]])
    orr, p_fisher = stats.fisher_exact(tab)

    rows = rows_pub + rows_near
    for r, (chi, p, e) in ((rows_pub, (chi_pub, p_pub, exp_pub)),
                           (rows_near, (chi_near, p_near, exp_near))):
        for x in r:
            x["chi2_omnibus"] = round(chi, 4)
            x["p_omnibus"] = float(f"{p:.6g}")
            x["min_expected_cell"] = round(e, 3)
    rows.append(dict(rule="rule_invariant_presence_contrast", band="any_SV_within_100kb",
                     n=int(prox.sum()), n_donors=int(d[prox].donor.nunique()),
                     n_subclonal=int(tab[0, 0]),
                     pct_subclonal=round(100.0 * tab[0, 0] / prox.sum(), 1),
                     median_ccf=np.nan, chi2_omnibus=np.nan,
                     p_omnibus=float(f"{p_fisher:.6g}"), min_expected_cell=np.nan))
    rows.append(dict(rule="rule_invariant_presence_contrast", band="no_SV_within_100kb",
                     n=int((~prox).sum()), n_donors=int(d[~prox].donor.nunique()),
                     n_subclonal=int(tab[1, 0]),
                     pct_subclonal=round(100.0 * tab[1, 0] / (~prox).sum(), 1),
                     median_ccf=np.nan, chi2_omnibus=np.nan,
                     p_omnibus=float(f"{p_fisher:.6g}"), min_expected_cell=np.nan))

    # Presence comparison using the shared coordinate-clustering rule.
    # Methods states that clusters within one patient are not independent, and the band
    # contrasts are tested by within-donor permutation. The presence contrast is the one the
    # paragraph keeps, so it is tested the same way here rather than pooled only.
    don = d["donor"].to_numpy()
    subc = (d["clonality"] == "Subclonal").to_numpy()
    both = [k for k in np.unique(don) if len(np.unique(prox.to_numpy()[don == k])) == 2]
    m = np.isin(don, both)
    tr = np.array([[int((prox & subc & m).sum()), int((prox & ~subc & m).sum())],
                   [int((~prox & subc & m).sum()), int((~prox & ~subc & m).sum())]])
    or_r, p_r = stats.fisher_exact(tr)

    rng = np.random.default_rng(7)
    idx = {k: np.where(don == k)[0] for k in np.unique(don)}
    pv = prox.to_numpy()

    def diff(lab):
        a = subc[lab].mean() if lab.sum() else np.nan
        b = subc[~lab].mean() if (~lab).sum() else np.nan
        return a - b

    obs_d = diff(pv)
    n_perm, hits = 20000, 0
    for _ in range(n_perm):
        pl = pv.copy()
        for k, ix in idx.items():
            pl[ix] = rng.permutation(pv[ix])
        v = diff(pl)
        if np.isfinite(v) and abs(v) >= abs(obs_d):
            hits += 1
    p_wd = (hits + 1) / (n_perm + 1)

    rows.append(dict(rule="presence_contrast_within_donor_permutation",
                     band="any_vs_no_SV_within_100kb", n=int(len(d)),
                     n_donors=int(len(both)), n_subclonal=int((prox & subc).sum()),
                     pct_subclonal=round(100.0 * obs_d, 1), median_ccf=np.nan,
                     chi2_omnibus=np.nan, p_omnibus=float(f"{p_wd:.6g}"),
                     min_expected_cell=np.nan))
    rows.append(dict(rule="presence_contrast_donors_in_both_arms",
                     band="any_vs_no_SV_within_100kb", n=int(m.sum()),
                     n_donors=int(len(both)), n_subclonal=int(tr[0, 0]),
                     pct_subclonal=round(100.0 * tr[0, 0] / tr[0].sum(), 1),
                     median_ccf=np.nan, chi2_omnibus=np.nan,
                     p_omnibus=float(f"{p_r:.6g}"), min_expected_cell=np.nan))

    out = pd.DataFrame(rows)
    out["odds_ratio_presence"] = round(float(orr), 4)
    out["odds_ratio_presence_donors_in_both_arms"] = round(float(or_r), 4)
    out["n_donors_single_arm_only"] = int(d.donor.nunique() - len(both))
    out.to_csv(OUT / "clonality_rebinning.csv", index=False)
    print(out.to_string(index=False))
    print(f"\n[published]   chi2 {chi_pub:.4f}, P {p_pub:.3g}, min expected cell {exp_pub:.3f}")
    print(f"[re-binned]   chi2 {chi_near:.4f}, P {p_near:.3g}, min expected cell {exp_near:.3f}")
    print(f"[invariant]   {tab[0,0]}/{int(prox.sum())} vs {tab[1,0]}/{int((~prox).sum())}, "
          f"Fisher P {p_fisher:.4f}, OR {orr:.2f}")
    print(f"[clustered]   within-donor permutation P {p_wd:.4f} on {len(both)} informative donors; "
          f"restricted to those donors {tr[0,0]}/{tr[0].sum()} vs {tr[1,0]}/{tr[1].sum()}, "
          f"Fisher P {p_r:.3f}, OR {1/or_r:.2f}; "
          f"{d.donor.nunique()-len(both)} of {d.donor.nunique()} donors contribute to one arm only")
    print(f"\nwrote {OUT/'clonality_rebinning.csv'} ({time.time()-t0:.0f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
