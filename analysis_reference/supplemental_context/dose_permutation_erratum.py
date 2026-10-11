"""Computes integration-count associations using the corrected permutation assignment.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

SRC = Path(r"source_project")
PROJ = Path(r"source_project")
BURDEN = PROJ / "results" / "hbv_sample_level_td_burden_donors.csv"
INTS = SRC / "data" / "processed" / "integration_sites_with_sv.csv"
OUT = PROJ / "results" / "dose_permutation_erratum.csv"
SEED = 7
OUTCOMES = ["dup_canonical_100k_7p5M", "dup_10k_50k", "dup_all", "del_all",
            "inv_all", "tra_all", "sv_all"]


def build() -> pd.DataFrame:
    d = pd.read_csv(BURDEN)
    d = d[d["hbv"] == True].reset_index(drop=True)              # noqa: E712
    s = pd.read_csv(INTS, low_memory=False)
    s = s[s["virus"].astype(str).str.upper() == "HBV"].copy()
    s["chrom"] = s["Chromosome"].astype(str).map(
        lambda c: c if c.startswith("chr") else f"chr{c}")
    s["locus"] = s["chrom"] + ":" + (s["chromStart"].to_numpy(np.int64) // 100_000).astype(str)
    d["n_loci"] = d["donor"].map(s.groupby("donor")["locus"].nunique()).fillna(0).astype(int)
    d["dose"] = np.minimum(d["n_loci"], 2)
    return d


def stat_broken(y, prj, lab):
    """The original statistic, reproduced exactly -- integer `lab` used as a positional index."""
    tot = 0.0
    for pp in np.unique(prj):
        k = prj == pp
        if lab[k].sum() in (0, int(k.sum())):
            continue
        r = stats.rankdata(y[k])
        tot += float(r[lab[k]].mean() - r.mean())
    return tot


def stat_boolean(y, prj, lab):
    """What the original meant: within-project rank advantage of exposed donors."""
    tot = 0.0
    for pp in np.unique(prj):
        k = prj == pp
        m = lab[k] > 0
        if m.sum() in (0, int(k.sum())):
            continue
        r = stats.rankdata(y[k])
        tot += float(r[m].mean() - r.mean())
    return tot


def stat_ordinal(y, prj, lab):
    """Sum of within-project Spearman correlations, weighted by project size."""
    tot = 0.0
    for pp in np.unique(prj):
        k = prj == pp
        if len(np.unique(lab[k])) < 2 or len(np.unique(y[k])) < 2:
            continue
        tot += float(k.sum()) * float(stats.spearmanr(lab[k], y[k]).statistic)
    return tot


def perm_p(fn, y, prj, lab, n_perm, rng):
    obs = fn(y, prj, lab)
    null = np.empty(n_perm)
    for i in range(n_perm):
        sh = lab.copy()
        for pp in np.unique(prj):
            k = prj == pp
            sh[k] = rng.permutation(lab[k])
        null[i] = fn(y, prj, sh)
    return obs, null, float((np.abs(null) >= abs(obs)).mean())


def nb2_irr(y, X, lab):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p0 = sm.GLM(y, X, family=sm.families.Poisson()).fit()
        mu = np.maximum(p0.mu, 1e-6)
        alpha = max(float(sm.OLS(((y - mu) ** 2 - y) / mu, mu).fit().params[0]), 1e-6)
        m = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
        return round(float(np.exp(m.params[lab])), 3), float(f"{m.pvalues[lab]:.4g}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    a = ap.parse_args()

    d = build()
    prj = d["project"].to_numpy()
    lab = d["dose"].to_numpy()
    print(f"[data]   {len(d)} HBV+ donors; dose "
          f"{d['dose'].value_counts().sort_index().to_dict()}; "
          f"{d['project'].nunique()} projects")

    rows = []
    for metric in OUTCOMES:
        y = d[metric].to_numpy(float)
        X = pd.DataFrame(dict(dose=lab.astype(float),
                              purity=d["purity"].to_numpy(float),
                              wgd=(d["wgd_status"].astype(str) == "wgd").astype(float),
                              tp53=pd.to_numeric(d.get("tp53_cn", np.nan),
                                                 errors="coerce").fillna(2.0),
                              load=np.log1p(d["hbv_load"].fillna(0).to_numpy(float))))
        for p in d["project"].unique():
            X[p] = (d["project"] == p).astype(float)
        X = X.loc[:, X.columns[~X.apply(lambda c: c.nunique(dropna=False) <= 1)]]
        X["_one"] = 1.0
        irr, glm_p = nb2_irr(y, X, "dose")

        out = dict(metric=metric, n=len(d), irr=irr, glm_p=glm_p)
        for name, fn in (("broken", stat_broken), ("boolean", stat_boolean),
                         ("ordinal", stat_ordinal)):
            rng = np.random.default_rng(SEED)
            npp = 200 if name == "broken" else a.n_perm     # 200 suffices to show a point mass
            obs, null, p = perm_p(fn, y, prj, lab, npp, rng)
            out[f"{name}_obs"] = round(obs, 4)
            out[f"{name}_null_sd"] = round(float(null.std(ddof=1)), 6)
            out[f"{name}_perm_p"] = round(p, 5)
            out[f"{name}_n_perm"] = npp
        rows.append(out)
        print(f"  [{metric:24s}] IRR {irr:5.3f} (GLM p={glm_p:.4g})   "
              f"broken p={out['broken_perm_p']:.3f} (null sd {out['broken_null_sd']:.2e})   "
              f"boolean p={out['boolean_perm_p']:.4f}   ordinal p={out['ordinal_perm_p']:.4f}",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    pd.set_option("display.width", 260)
    print("\n" + df[["metric", "irr", "glm_p", "broken_null_sd", "broken_perm_p",
                     "boolean_perm_p", "ordinal_perm_p"]].to_string(index=False))
    deg = (df["broken_null_sd"] == 0).all()
    print(f"\nbroken statistic has zero null variance in all rows: {deg}"
          f"  -> perm_p = 1.0 was a construction artefact, not a result")
    # Holm across the seven outcomes on the ordinal statistic
    o = df.sort_values("ordinal_perm_p").reset_index(drop=True)
    o["holm"] = [min(1.0, (len(o) - i) * p) for i, p in enumerate(o["ordinal_perm_p"])]
    o["holm"] = o["holm"].cummax()
    print("\nordinal permutation, Holm-corrected across 7 outcomes:")
    print(o[["metric", "irr", "ordinal_perm_p", "holm"]].to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
