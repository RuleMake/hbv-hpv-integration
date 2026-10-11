"""Compares large-duplication burden across HBV detection groups in liver samples.
"""

from __future__ import annotations

import gzip
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

SRC = Path(r"source_project")
OUT = Path(r"source_project\results")
SV_DIRS = [SRC / "data/raw/pcawg_sv/icgc_open", SRC / "data/raw/pcawg_sv/tcga_open"]
HIST = SRC / "data/raw/clinical/pcawg_specimen_histology_August2016_v9.xlsx"
SHEET = SRC / "data/raw/viral_integration/pcawg_sample_sheet.v1.4.2016-09-14.tsv"
ZAP = SRC / "data/raw/viral_integration/Zapatka2020_PCAWG_SupplTables_MOESM3.xlsx"
PURITY = SRC / "data/raw/pcawg_cnv/consensus.20170217.purity.ploidy.txt.gz"
SEED = 7
N_PERM = 20000

SPAN_CLASSES = {
    "dup_1k_10k": (1_000, 10_000),
    "dup_10k_50k": (10_000, 50_000),
    "dup_50k_100k": (50_000, 100_000),
    "dup_canonical_100k_7p5M": (100_000, 7_500_000),
    "dup_gt_7p5M": (7_500_000, np.inf),
    "dup_all": (0, np.inf),
}
OTHER = ["del_all", "inv_all", "tra_all", "sv_all"]


def liver_donors() -> pd.DataFrame:
    h = pd.read_excel(HIST)
    lv = h[h["histology_abbreviation"].astype(str).str.contains("Liver", case=False, na=False)]
    lv = (lv[["icgc_donor_id", "project_code", "histology_abbreviation"]]
          .drop_duplicates("icgc_donor_id")
          .rename(columns={"icgc_donor_id": "donor", "project_code": "project"}))
    print(f"[cohort]  {len(lv)} PCAWG liver donors "
          f"({lv['histology_abbreviation'].value_counts().to_dict()})")
    return lv


def sv_paths() -> dict[str, Path]:
    p = {f.name.split(".")[0]: f for d in SV_DIRS for f in d.rglob("*.bedpe.gz")}
    print(f"[sv]      {len(p)} consensus SV call sets on disk")
    return p


def donor_to_aliquot() -> pd.DataFrame:
    ss = pd.read_csv(SHEET, sep="\t", low_memory=False)
    w = ss[(ss["library_strategy"].astype(str).str.upper() == "WGS")
           & (~ss["dcc_specimen_type"].astype(str).str.contains("Normal", na=False))]
    return w[["icgc_donor_id", "aliquot_id"]].rename(columns={"icgc_donor_id": "donor"})


def count_sv(path: Path) -> dict:
    """Per-donor SV counts. svclass comes from the PCAWG consensus caller; inversions are
    emitted as h2hINV / t2tINV, never as a bare 'INV'."""
    with gzip.open(path, "rt") as fh:
        d = pd.read_csv(fh, sep="\t", low_memory=False)
    row = {k: 0 for k in list(SPAN_CLASSES) + OTHER}
    if d.empty:
        return row
    same = (d["chrom1"].astype(str) == d["chrom2"].astype(str)).to_numpy()
    span = np.abs(d["start2"].to_numpy(np.int64) - d["start1"].to_numpy(np.int64))
    cls = d["svclass"].astype(str).str.upper().to_numpy()
    isdup = cls == "DUP"
    for name, (lo, hi) in SPAN_CLASSES.items():
        row[name] = int((isdup & same & (span > lo) & (span <= hi)).sum())
    row["del_all"] = int(((cls == "DEL") & same).sum())
    row["inv_all"] = int((np.isin(cls, ["H2HINV", "T2TINV", "INV"]) & same).sum())
    row["tra_all"] = int((cls == "TRA").sum())
    row["sv_all"] = int(len(d))
    return row


def zapatka() -> tuple[set[str], pd.DataFrame]:
    cons = pd.read_excel(ZAP, sheet_name="6 WGS Consensus-Tumor")
    tax = cons["NCBI_taxname"].astype(str)
    hbv = set(cons.loc[tax.str.contains("hepatitis B|orthohepadna", case=False, na=False),
                       "icgc_donor_id"].astype(str))
    print(f"[zapatka] {len(hbv)} PCAWG donors HBV-positive by WGS consensus")
    load = pd.read_excel(ZAP, sheet_name="12 Virus load")
    hb = load[load["genus_name"].astype(str).str.contains("hepadna", case=False, na=False)]
    hb = (hb.groupby("icgc_donor_id")["virus_load"].max().rename("hbv_load")
          .reset_index().rename(columns={"icgc_donor_id": "donor"}))
    return hbv, hb


def fit(d: pd.DataFrame, metric: str, offset_col: str | None, rng,
        extra: tuple[str, ...] = ()) -> dict:
    """NB2 GLM: count ~ hbv + project + purity + wgd (+ extra), optional log(total SV) offset."""
    y = d[metric].to_numpy(float)
    X = pd.get_dummies(d["project"], prefix="prj", drop_first=True).astype(float)
    X["purity"] = d["purity"].to_numpy(float)
    X["wgd"] = (d["wgd_status"].astype(str) == "wgd").astype(float)
    for c in extra:
        X[c] = pd.to_numeric(d[c], errors="coerce").fillna(
            pd.to_numeric(d[c], errors="coerce").median()).to_numpy(float)
    X["hbv"] = d["hbv"].astype(float).to_numpy()
    X = sm.add_constant(X, has_constant="add")
    off = np.log(np.maximum(d[offset_col].to_numpy(float), 1.0)) if offset_col else None
    res = dict(metric=metric, offset=offset_col or "none", n=len(d),
               n_hbv_pos=int(d.hbv.sum()), n_hbv_neg=int((~d.hbv).sum()),
               mean_hbv_pos=round(float(y[d.hbv.to_numpy()].mean()), 3),
               mean_hbv_neg=round(float(y[~d.hbv.to_numpy()].mean()), 3),
               median_hbv_pos=float(np.median(y[d.hbv.to_numpy()])),
               median_hbv_neg=float(np.median(y[~d.hbv.to_numpy()])))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:  # estimate alpha from a Poisson fit, then refit NB2
            p0 = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
            mu = np.maximum(p0.mu, 1e-6)
            aux = ((y - mu) ** 2 - y) / mu
            alpha = max(float(sm.OLS(aux, mu).fit().params[0]), 1e-6)
            m = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha), offset=off).fit()
            b, se = float(m.params["hbv"]), float(m.bse["hbv"])
            res.update(alpha=round(alpha, 3), IRR=round(float(np.exp(b)), 4),
                       ci_lo=round(float(np.exp(b - 1.96 * se)), 4),
                       ci_hi=round(float(np.exp(b + 1.96 * se)), 4),
                       wald_p=float(f"{m.pvalues['hbv']:.4g}"))
        except Exception as e:  # noqa: BLE001
            res.update(alpha=np.nan, IRR=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                       wald_p=np.nan, note=str(e)[:60])

    # distribution-free backup: permute the HBV label WITHIN project, statistic = sum of
    # within-project rank-sums for the HBV+ group. Assumes nothing about the mean-variance link.
    lab = d["hbv"].to_numpy()
    prj = d["project"].to_numpy()
    val = y / np.maximum(d[offset_col].to_numpy(float), 1.0) if offset_col else y

    def stat(mask):
        s = 0.0
        for p in np.unique(prj):
            k = prj == p
            if mask[k].sum() in (0, int(k.sum())):
                continue
            s += float(stats.rankdata(val[k])[mask[k]].mean() - stats.rankdata(val[k]).mean())
        return s

    obs = stat(lab)
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        sh = lab.copy()
        for p in np.unique(prj):
            k = prj == p
            sh[k] = rng.permutation(lab[k])
        null[i] = stat(sh)
    res["perm_p_2sided"] = round(float((np.abs(null) >= abs(obs)).mean()), 5)
    return res


def main() -> None:
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    lv, paths, al = liver_donors(), sv_paths(), donor_to_aliquot()
    al = al[al["aliquot_id"].isin(paths)]
    m = lv.merge(al, on="donor", how="inner").drop_duplicates("donor")
    print(f"[cohort]  {len(m)}/{len(lv)} liver donors have a tumour WGS SV call set")

    rows = []
    for i, r in enumerate(m.itertuples(), 1):
        rows.append({"donor": r.donor, "project": r.project, "aliquot": r.aliquot_id,
                     **count_sv(paths[r.aliquot_id])})
        if i % 100 == 0:
            print(f"          parsed {i}/{len(m)}", flush=True)
    d = pd.DataFrame(rows)

    pur = pd.read_csv(PURITY, sep="\t")[["samplename", "purity", "ploidy", "wgd_status"]]
    d = d.merge(pur.rename(columns={"samplename": "aliquot"}), on="aliquot", how="left")
    tp = pd.read_csv(SRC / "data/raw/pcawg_cnv/all_samples.consensus_CN.by_gene.170214.txt.gz",
                     sep="\t", low_memory=False).set_index("Gene Symbol")
    if "TP53" in tp.index:
        row = tp.loc["TP53"]
        row = row.iloc[0] if isinstance(row, pd.DataFrame) else row
        d["tp53_cn"] = d["aliquot"].map(row.to_dict())
        print(f"[tp53]    gene-level consensus CN for {int(d['tp53_cn'].notna().sum())} donors "
              f"(COPY NUMBER only -- R249S point mutation is not on disk)")
    del tp

    hbv, load = zapatka()
    d["hbv"] = d["donor"].isin(hbv)
    d = d.merge(load, on="donor", how="left")
    before = len(d)
    d = d[d["purity"].notna() & d["wgd_status"].notna()].reset_index(drop=True)
    d.to_csv(OUT / "hbv_sample_level_td_burden_donors.csv", index=False)

    print(f"\n[model]   {len(d)}/{before} donors with consensus purity+WGD "
          f"(dropped {before-len(d)} for missing covariates)")
    print(f"[exposure] HBV+ {int(d.hbv.sum())} | HBV- {int((~d.hbv).sum())}")
    print(d.groupby(["project", "hbv"]).size().unstack(fill_value=0).to_string())
    a, b = d.loc[d.hbv, "purity"], d.loc[~d.hbv, "purity"]
    print(f"[purity]   median HBV+ {a.median():.3f} vs HBV- {b.median():.3f} "
          f"(MWU p = {stats.mannwhitneyu(a, b).pvalue:.3g})")
    print(f"[wgd]      HBV+ {100*(d.loc[d.hbv,'wgd_status']=='wgd').mean():.1f}% vs "
          f"HBV- {100*(d.loc[~d.hbv,'wgd_status']=='wgd').mean():.1f}%")

    res = []
    for k in list(SPAN_CLASSES) + OTHER:
        res.append(fit(d, k, None, rng))
    for k in ("dup_canonical_100k_7p5M", "dup_10k_50k", "del_all", "inv_all"):
        res.append(fit(d, k, "sv_all", rng))
    # Assess TP53 copy number as a potential covariate. These inputs provide gene-level
    # consensus copy number, but no per-donor point-mutation calls; this comparison
    # therefore adjusts for TP53 deletion rather than mutation status.
    if "tp53_cn" in d.columns and d["tp53_cn"].notna().sum() >= 100:
        for k in ("dup_canonical_100k_7p5M", "dup_all", "sv_all"):
            r = fit(d, k, None, rng, extra=("tp53_cn",))
            r["metric"] = k + " [+TP53 CN]"
            res.append(r)
    out = pd.DataFrame(res)
    out.to_csv(OUT / "hbv_sample_level_td_burden.csv", index=False)

    hp = d[d.hbv & d["hbv_load"].notna()]
    dose = []
    for k in ("dup_canonical_100k_7p5M", "dup_10k_50k", "del_all", "inv_all", "sv_all"):
        if len(hp) >= 15:
            rho, p = stats.spearmanr(hp["hbv_load"], hp[k])
            rho_p, p_p = stats.spearmanr(hp["hbv_load"], hp["purity"])
            dose.append(dict(metric=k, n=len(hp), spearman_rho=round(float(rho), 4),
                             p=round(float(p), 5),
                             rho_load_vs_purity=round(float(rho_p), 4),
                             p_load_vs_purity=round(float(p_p), 5)))

    pd.set_option("display.width", 260)
    cols = ["metric", "offset", "mean_hbv_pos", "mean_hbv_neg", "median_hbv_pos",
            "median_hbv_neg", "IRR", "ci_lo", "ci_hi", "wald_p", "perm_p_2sided"]
    print("\n=== NB2 GLM, adjusted for project + purity + WGD ===")
    print(out[cols].to_string(index=False))
    if dose:
        print("\n=== dose-response within HBV+ (virus load vs burden) ===")
        print(pd.DataFrame(dose).to_string(index=False))
    print(f"\nwrote {OUT/'hbv_sample_level_td_burden.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
