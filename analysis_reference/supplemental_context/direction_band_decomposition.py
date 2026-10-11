"""Examines clonality distance bands and copy-number/ploidy adjustment.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import timing_clonality_order as TC  # noqa: E402

OUT = TC.OUT
SEED = 7
BANDS = ("JUNCTION", "SEAM", "FAR", "NONE")


def attach_cn(d: pd.DataFrame) -> pd.DataFrame:
    """Attach local consensus copy number using the shared clonality-analysis definitions."""
    pur = pd.read_csv(TC.PURITY, sep="\t").set_index("samplename")
    sites = pd.read_csv(TC.SITES)
    ali = {str(r.project): str(r.aliquot_id) for r in sites.itertuples()}
    cn = pd.read_csv(TC.CNBYGENE, sep="\t", low_memory=False).set_index("Gene Symbol")
    gp = pd.read_parquet(TC.GENES)
    gp = gp[gp["gene_type"] == "protein_coding"].copy()
    gp["chrom"] = gp["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
    gx = {c: (s.sort_values("start")["start"].to_numpy(np.int64),
              s.sort_values("start")["end"].to_numpy(np.int64),
              s.sort_values("start")["gene_name"].to_numpy(str))
          for c, s in gp.groupby("chrom")}
    loc, rel, pur_v = [], [], []
    for r in d.itertuples():
        al = ali.get(str(r.donor))
        g = TC.nearest_gene(gx, str(r.chrom), int(r.pos0))
        v = np.nan
        if al and g and al in cn.columns and g in cn.index:
            try:
                v = float(np.atleast_1d(cn.loc[g, al])[0])
            except Exception:  # noqa: BLE001
                v = np.nan
        plo = float(pur.loc[al, "ploidy"]) if al in pur.index else np.nan
        loc.append(v)
        rel.append(v / plo if np.isfinite(v) and np.isfinite(plo) and plo > 0 else np.nan)
        pur_v.append(float(pur.loc[al, "purity"]) if al in pur.index else np.nan)
    d = d.copy()
    d["local_cn"], d["cn_over_ploidy"], d["purity"] = loc, rel, pur_v
    return d


def assign_band(d: pd.DataFrame, w: int = 100) -> pd.Series:
    j0 = d[f"n_sv_{w}k_j0"] > 0
    j10 = d[f"n_sv_{w}k_j10000"] > 0
    j50 = d[f"n_sv_{w}k_j50000"] > 0
    return pd.Series(np.select([j50, j10 & ~j50, j0 & ~j10], ["FAR", "SEAM", "JUNCTION"],
                               default="NONE"), index=d.index, name="band")


def within_donor_perm(d: pd.DataFrame, stat_fn, n_perm: int, rng) -> tuple[float, float, int]:
    """Permute band labels within donor. Returns (observed, two-sided p, n informative donors)."""
    obs = stat_fn(d["band"].to_numpy(), d)
    groups = [g.index.to_numpy() for _, g in d.groupby("donor", sort=True)]
    informative = sum(1 for ix in groups if d.loc[ix, "band"].nunique() > 1)
    lab = d["band"].to_numpy().copy()
    pos = {ix_val: i for i, ix_val in enumerate(d.index)}
    gidx = [np.array([pos[i] for i in ix]) for ix in groups]
    hits = 0
    for _ in range(n_perm):
        p = lab.copy()
        for gi in gidx:
            p[gi] = rng.permutation(p[gi])
        s = stat_fn(p, d)
        if np.isfinite(s) and abs(s) >= abs(obs) - 1e-12:
            hits += 1
    return float(obs), (hits + 1) / (n_perm + 1), informative


def frac_sub(labels: np.ndarray, d: pd.DataFrame, band: str) -> float:
    m = labels == band
    return float((d["clonality"].to_numpy()[m] == "Subclonal").mean()) if m.sum() else np.nan


def median_cn(labels: np.ndarray, d: pd.DataFrame, band: str) -> float:
    v = d["cn_over_ploidy"].to_numpy(float)[(labels == band) & d["cn_over_ploidy"].notna().to_numpy()]
    return float(np.median(v)) if len(v) else np.nan


def mde(a: np.ndarray, b: np.ndarray) -> float:
    """80% power, alpha 0.05, two-sided -- the project's MDE convention."""
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return float((1.96 + 0.84) * np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    d = pd.read_csv(TC.CIS)
    d["clonality"] = d["clonality"].astype(str)
    d = d[d["clonality"].isin(["Clonal", "Subclonal"])].reset_index(drop=True)
    d = attach_cn(d)
    d["band"] = assign_band(d)
    print(f"[in]      {len(d)} breakpoints, {d.donor.nunique()} donors, "
          f"CN attached to {int(d.local_cn.notna().sum())}")
    print(f"[bands]   {d.band.value_counts().reindex(BANDS).to_dict()}")

    # Check that the distance bands reproduce the saved two-group contrasts.
    prox50 = d["n_sv_100k_j50000"] > 0
    prox10 = d["n_sv_100k_j10000"] > 0
    for name, m in (("j=50kb", prox50), ("j=10kb", prox10)):
        t = np.array([[int((m & (d.clonality == "Subclonal")).sum()), int((m & (d.clonality == "Clonal")).sum())],
                      [int((~m & (d.clonality == "Subclonal")).sum()), int((~m & (d.clonality == "Clonal")).sum())]])
        orr, pp = stats.fisher_exact(t)
        print(f"[check]    Clonality contrast {name}: n_prox={int(m.sum())} OR={orr:.4f} p={pp:.5f}")

    rows = []

    # ---- 1. band profile
    for b in BANDS:
        s = d[d.band == b]
        cnv = s["cn_over_ploidy"].dropna()
        rows.append(dict(analysis="band_profile", band=b, n=len(s), n_donors=s.donor.nunique(),
                         n_subclonal=int((s.clonality == "Subclonal").sum()),
                         frac_subclonal=round(float((s.clonality == "Subclonal").mean()), 4) if len(s) else np.nan,
                         median_ccf=round(float(s.ccf.median()), 4) if len(s) else np.nan,
                         median_local_cn=round(float(s.local_cn.median()), 3) if s.local_cn.notna().any() else np.nan,
                         median_cn_over_ploidy=round(float(cnv.median()), 4) if len(cnv) else np.nan,
                         n_cn=len(cnv)))

    # ---- 2. omnibus across bands
    grp_cn = [d.loc[d.band == b, "cn_over_ploidy"].dropna().to_numpy(float) for b in BANDS]
    grp_ccf = [d.loc[d.band == b, "ccf"].dropna().to_numpy(float) for b in BANDS]
    for label, gs in (("cn_over_ploidy", grp_cn), ("ccf", grp_ccf)):
        gs2 = [g for g in gs if len(g) >= 3]
        if len(gs2) >= 3:
            h = stats.kruskal(*gs2)
            rows.append(dict(analysis="omnibus_kruskal", metric=label, n_groups=len(gs2),
                             statistic=round(float(h.statistic), 4), p=round(float(h.pvalue), 5)))
    ct = pd.crosstab(d.band, d.clonality)
    if ct.shape[0] >= 2:
        chi = stats.chi2_contingency(ct.to_numpy())
        rows.append(dict(analysis="omnibus_clonality", metric="band x clonality",
                         statistic=round(float(chi.statistic), 4), p=round(float(chi.pvalue), 5),
                         n=int(ct.to_numpy().sum())))

    # ---- 3. pairwise: clonality (Fisher + within-donor permutation)
    for x, y in (("SEAM", "FAR"), ("SEAM", "NONE"), ("FAR", "NONE"), ("JUNCTION", "SEAM")):
        s = d[d.band.isin([x, y])]
        if s.band.nunique() < 2 or len(s) < 10:
            continue
        t = pd.crosstab(s.band, s.clonality)
        if t.shape != (2, 2):
            continue
        t = t.reindex(index=[x, y], columns=["Subclonal", "Clonal"]).fillna(0).to_numpy()
        orr, pp = stats.fisher_exact(t)
        obs, permp, ndon = within_donor_perm(
            s.reset_index(drop=True),
            lambda lab, dd, x=x, y=y: frac_sub(lab, dd, x) - frac_sub(lab, dd, y),
            a.n_perm, rng)
        rows.append(dict(analysis="pairwise_clonality", band=x, vs=y,
                         n_x=int((s.band == x).sum()), n_y=int((s.band == y).sum()),
                         frac_subclonal_x=round(float((s[s.band == x].clonality == "Subclonal").mean()), 4),
                         frac_subclonal_y=round(float((s[s.band == y].clonality == "Subclonal").mean()), 4),
                         fisher_or=round(float(orr), 4) if np.isfinite(orr) else np.nan,
                         fisher_p=round(float(pp), 5), diff=round(obs, 4),
                         within_donor_perm_p=round(permp, 5), n_informative_donors=ndon))

    # ---- 4. THE DETECTION-BIAS TEST: is SEAM sitting at high copy number?
    for x, y in (("SEAM", "FAR"), ("SEAM", "NONE"), ("FAR", "NONE")):
        av = d.loc[d.band == x, "cn_over_ploidy"].dropna().to_numpy(float)
        bv = d.loc[d.band == y, "cn_over_ploidy"].dropna().to_numpy(float)
        if len(av) < 3 or len(bv) < 3:
            continue
        u = stats.mannwhitneyu(av, bv, alternative="two-sided")
        s = d[d.band.isin([x, y]) & d.cn_over_ploidy.notna()].reset_index(drop=True)
        obs, permp, ndon = within_donor_perm(
            s, lambda lab, dd, x=x, y=y: median_cn(lab, dd, x) - median_cn(lab, dd, y),
            a.n_perm, rng)
        rows.append(dict(analysis="detection_bias_cn", band=x, vs=y, n_x=len(av), n_y=len(bv),
                         median_cn_over_ploidy_x=round(float(np.median(av)), 4),
                         median_cn_over_ploidy_y=round(float(np.median(bv)), 4),
                         diff=round(obs, 4), mwu_p=round(float(u.pvalue), 5),
                         within_donor_perm_p=round(permp, 5), n_informative_donors=ndon,
                         mde_80pct=round(mde(av, bv), 4)))

    # ---- 5. joint model: does band survive adjustment for copy number?
    s = d[d.cn_over_ploidy.notna() & d.band.isin(["SEAM", "FAR", "NONE"])].copy()
    if s.band.nunique() == 3 and len(s) >= 30:
        import statsmodels.api as sm
        X = pd.DataFrame({
            "is_far": (s.band == "FAR").astype(float).to_numpy(),
            "is_none": (s.band == "NONE").astype(float).to_numpy(),
            "log_cn": np.log2(np.maximum(s.cn_over_ploidy.to_numpy(float), 0.05)),
        })
        y = (s.clonality == "Subclonal").astype(float).to_numpy()
        try:
            m = sm.Logit(y, sm.add_constant(X)).fit(disp=0)
            for term in ("is_far", "is_none", "log_cn"):
                rows.append(dict(analysis="logit_subclonal_ref_SEAM", term=term, n=len(s),
                                 odds_ratio=round(float(np.exp(m.params[term])), 4),
                                 ci_lo=round(float(np.exp(m.conf_int().loc[term, 0])), 4),
                                 ci_hi=round(float(np.exp(m.conf_int().loc[term, 1])), 4),
                                 p=round(float(m.pvalues[term]), 5)))
        except Exception as e:  # noqa: BLE001
            print(f"[warn]    logit failed: {e}")

    # ---- 6. CN-normal stratum, where a coverage-driven bias cannot operate
    n = d[(d.cn_over_ploidy - 1).abs() <= 0.25]
    for b in BANDS:
        s = n[n.band == b]
        if len(s):
            rows.append(dict(analysis="cn_normal_stratum", band=b, n=len(s),
                             n_donors=s.donor.nunique(),
                             frac_subclonal=round(float((s.clonality == "Subclonal").mean()), 4)))

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "direction_band_decomposition.csv", index=False)

    dt = (d.groupby(["donor", "band"])
            .agg(n=("ccf", "size"), n_subclonal=("clonality", lambda s: int((s == "Subclonal").sum())),
                 median_ccf=("ccf", "median"), median_cn_over_ploidy=("cn_over_ploidy", "median"))
            .reset_index().round(4))
    dt.to_csv(OUT / "direction_band_decomposition_donors.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    for an in out.analysis.unique():
        sub = out[out.analysis == an].dropna(axis=1, how="all")
        print(f"\n=== {an} ===")
        print(sub.drop(columns=["analysis"]).to_string(index=False))
    print(f"\n=== per donor ===\n{dt.to_string(index=False)}")
    print(f"\nwrote {OUT/'direction_band_decomposition.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
