"""Builds clonality and rearrangement-distance comparisons with sample-level controls.
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402
import joint_correction as J  # noqa: E402

SRC = Path(r"source_project")
OUT = J.OUT
CIS = OUT / "sasa_cis_orientation.csv"          # produced by the Sasa orientation analysis
SITES = SRC / "data/processed/integration_sites_with_sv.csv"
PURITY = SRC / "data/raw/pcawg_cnv/consensus.20170217.purity.ploidy.txt.gz"
CNBYGENE = SRC / "data/raw/pcawg_cnv/all_samples.consensus_CN.by_gene.170214.txt.gz"
GENES = SRC / "data/interim/gencode_v19_genes.parquet"
JUNCS = (10_000, 50_000)
WINDOWS = (100_000, 1_000_000)
CCF_SUB = 0.6          # "clearly subclonal"
SEED = 7


def cn_control(d: pd.DataFrame, prox: pd.Series, jx: int, w: int) -> dict:
    """Assess the contribution of local copy number to the clonality comparison.
    
    Compare copy number near and away from structural variants, adjust the proximity
    association for copy number/ploidy, and repeat it within a near-normal copy-number
    stratum."""
    out: dict[str, float] = {}
    a = d.loc[prox, "cn_over_ploidy"].dropna()
    b = d.loc[~prox, "cn_over_ploidy"].dropna()
    if len(a) >= 5 and len(b) >= 5:
        out["median_cnratio_proximal"] = round(float(a.median()), 3)
        out["median_cnratio_distal"] = round(float(b.median()), 3)
        out["cn_mwu_p"] = round(float(stats.mannwhitneyu(a, b).pvalue), 5)
    s = d.loc[d["cn_over_ploidy"].notna()].copy()
    s["prox"] = (s[f"n_sv_{w//1000}k_j{jx}"] > 0).astype(float)
    s["sub"] = (s["clonality"] == "Subclonal").astype(float)
    if s["sub"].nunique() > 1 and s["prox"].nunique() > 1 and len(s) >= 30:
        import statsmodels.api as sm
        X = sm.add_constant(pd.DataFrame({"prox": s["prox"].to_numpy(),
                                          "log_cn": np.log2(np.maximum(
                                              s["cn_over_ploidy"].to_numpy(float), 0.05))}))
        try:
            m = sm.Logit(s["sub"].to_numpy(), X).fit(disp=0)
            out["adj_OR_prox"] = round(float(np.exp(m.params["prox"])), 4)
            out["adj_p_prox"] = round(float(m.pvalues["prox"]), 5)
            out["adj_OR_logcn"] = round(float(np.exp(m.params["log_cn"])), 4)
            out["adj_p_logcn"] = round(float(m.pvalues["log_cn"]), 5)
        except Exception:  # noqa: BLE001
            pass
    # Restrict to |CN/ploidy - 1| <= 0.25 to reduce copy-number differences.
    n = s[(s["cn_over_ploidy"] - 1).abs() <= 0.25]
    if len(n) >= 20 and n["prox"].nunique() > 1:
        t = pd.crosstab(n["prox"], n["sub"])
        if t.shape == (2, 2):
            orr, pp = stats.fisher_exact(t.to_numpy())
            out["cn_normal_n"] = int(len(n))
            out["cn_normal_OR"] = round(float(orr), 4) if np.isfinite(orr) else np.nan
            out["cn_normal_p"] = round(float(pp), 5)
    return out


def nearest_gene(gx, chrom: str, pos: int) -> str | None:
    g = gx.get(chrom)
    if g is None:
        return None
    st, en, name = g
    d = np.maximum(np.maximum(st - pos, pos - en), 0)
    j = int(np.argmin(d))
    return str(name[j]) if d[j] < 500_000 else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.parse_args()
    t0 = time.time()

    if not CIS.exists():
        raise SystemExit("run analysis/sasa_cis_orientation.py first (it writes the cis table)")
    d = pd.read_csv(CIS)
    d["clonality"] = d["clonality"].astype(str)
    d = d[d["clonality"].isin(["Clonal", "Subclonal"])].copy()
    print(f"[cis]     {len(d)} breakpoints, {d.donor.nunique()} donors, "
          f"clonality {d.clonality.value_counts().to_dict()}")
    print(f"[ccf]     median CCF {d['ccf'].median():.3f}; "
          f"Clonal {d.loc[d.clonality=='Clonal','ccf'].median():.3f} vs "
          f"Subclonal {d.loc[d.clonality=='Subclonal','ccf'].median():.3f}")

    # Attach local consensus copy number to assess potential detection differences:
    # higher copy number can increase junction coverage and affect subclonal detection.
    pur0 = pd.read_csv(PURITY, sep="\t").set_index("samplename")
    sites0 = pd.read_csv(SITES)
    ali0 = {str(r.project): str(r.aliquot_id) for r in sites0.itertuples()}
    cn0 = pd.read_csv(CNBYGENE, sep="\t", low_memory=False).set_index("Gene Symbol")
    gp0 = pd.read_parquet(GENES)
    gp0 = gp0[gp0["gene_type"] == "protein_coding"].copy()
    gp0["chrom"] = gp0["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
    gx0 = {c: (s.sort_values("start")["start"].to_numpy(np.int64),
               s.sort_values("start")["end"].to_numpy(np.int64),
               s.sort_values("start")["gene_name"].to_numpy(str))
           for c, s in gp0.groupby("chrom")}
    loc_cn, rel_cn = [], []
    for r in d.itertuples():
        al = ali0.get(str(r.donor))
        g = nearest_gene(gx0, str(r.chrom), int(r.pos0))
        v = np.nan
        if al and g and al in cn0.columns and g in cn0.index:
            try:
                v = float(np.atleast_1d(cn0.loc[g, al])[0])
            except Exception:  # noqa: BLE001
                v = np.nan
        loc_cn.append(v)
        plo = float(pur0.loc[al, "ploidy"]) if al in pur0.index else np.nan
        rel_cn.append(v / plo if np.isfinite(v) and np.isfinite(plo) and plo > 0 else np.nan)
    d["local_cn"] = loc_cn
    d["cn_over_ploidy"] = rel_cn
    print(f"[cn]      local consensus CN attached to {int(d['local_cn'].notna().sum())}/{len(d)} "
          f"breakpoints; median CN/ploidy {np.nanmedian(d['cn_over_ploidy']):.3f}")

    rows = []
    for jx in JUNCS:
        for w in WINDOWS:
            col = f"n_sv_{w//1000}k_j{jx}"
            if col not in d.columns:
                continue
            prox = d[col] > 0
            # ---- readout 1: clonality composition of SV-proximal integrants
            n_prox = int(prox.sum())
            n_sub_prox = int((prox & (d.clonality == "Subclonal")).sum())
            n_sub_dist = int((~prox & (d.clonality == "Subclonal")).sum())
            n_dist = int((~prox).sum())
            tab = np.array([[n_sub_prox, n_prox - n_sub_prox],
                            [n_sub_dist, n_dist - n_sub_dist]])
            odds, pf = stats.fisher_exact(tab) if tab.min() >= 0 and tab.sum() > 0 else (np.nan,
                                                                                        np.nan)
            # ---- readout 2: CCF contrast
            a_ = d.loc[prox, "ccf"].dropna().to_numpy(float)
            b_ = d.loc[~prox, "ccf"].dropna().to_numpy(float)
            if len(a_) >= 5 and len(b_) >= 5:
                sp = np.sqrt(a_.var(ddof=1) / len(a_) + b_.var(ddof=1) / len(b_))
                mde = (1.96 + 0.84) * sp
                u = stats.mannwhitneyu(a_, b_, alternative="two-sided")
                rows.append(dict(junction_excl_kb=jx // 1000, window_kb=w // 1000,
                                 n_sv_proximal=n_prox, n_sv_distal=n_dist,
                                 frac_subclonal_proximal=round(n_sub_prox / max(n_prox, 1), 4),
                                 frac_subclonal_distal=round(n_sub_dist / max(n_dist, 1), 4),
                                 fisher_or=round(float(odds), 4) if np.isfinite(odds) else np.nan,
                                 fisher_p=round(float(pf), 5) if np.isfinite(pf) else np.nan,
                                 median_ccf_proximal=round(float(np.median(a_)), 4),
                                 median_ccf_distal=round(float(np.median(b_)), 4),
                                 ccf_diff=round(float(np.median(a_) - np.median(b_)), 4),
                                 mde_ccf_diff_80pct=round(float(mde), 4),
                                 mwu_p=round(float(u.pvalue), 5),
                                 **cn_control(d, prox, jx, w)))

    # ---- readout 3: the logical-ordering count
    pur = pd.read_csv(PURITY, sep="\t").set_index("samplename")
    sites = pd.read_csv(SITES)
    ali = {str(r.project): str(r.aliquot_id) for r in sites.itertuples()}
    cn = pd.read_csv(CNBYGENE, sep="\t", low_memory=False).set_index("Gene Symbol")
    gp = pd.read_parquet(GENES)
    gp = gp[gp["gene_type"] == "protein_coding"].copy()
    gp["chrom"] = gp["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
    gx = {c: (s.sort_values("start")["start"].to_numpy(np.int64),
              s.sort_values("start")["end"].to_numpy(np.int64),
              s.sort_values("start")["gene_name"].to_numpy(str))
          for c, s in gp.groupby("chrom")}

    proofs = []
    for r in d.itertuples():
        al = ali.get(str(r.donor))
        if al is None or al not in pur.index or al not in cn.columns:
            continue
        p = float(pur.loc[al, "purity"])
        plo = float(pur.loc[al, "ploidy"])
        gname = nearest_gene(gx, str(r.chrom), int(r.pos0))
        if gname is None or gname not in cn.index:
            continue
        try:
            v = np.atleast_1d(cn.loc[gname, al])[0]
        except Exception:  # noqa: BLE001
            continue
        if not np.isfinite(v):
            continue
        # the cell fraction a gain of this size needs, given purity: a gain of (v - ploidy) copies
        # observed in the bulk requires at least (v - ploidy) / (purity * max_extra) of the cells.
        # Use the conservative form: a gain is at least near-clonal if v >= ploidy + 2 / purity.
        clonal_gain = v >= plo + 2.0 / max(p, 0.05)
        if clonal_gain and np.isfinite(r.ccf) and r.ccf < CCF_SUB:
            proofs.append(dict(donor=r.donor, chrom=r.chrom, pos0=r.pos0, gene=gname,
                               integrant_ccf=round(float(r.ccf), 3),
                               clonality=r.clonality, gene_cn=float(v), ploidy=plo,
                               purity=p, verdict="GAIN PRECEDES INTEGRANT"))
    pf_ = pd.DataFrame(proofs)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "timing_clonality_order.csv", index=False)
    pf_.to_csv(OUT / "timing_clonality_order_pairs.csv", index=False)
    pd.set_option("display.width", 250)
    print("\n=== readouts 1 and 2: clonality and CCF of SV-proximal vs SV-distal integrants ===")
    print(out.to_string(index=False))
    print("\n=== Subclonal integrations within gains classified as near-clonal ===")
    if len(pf_):
        print(pf_.to_string(index=False))
    else:
        print("  none found -- no integrant in this set is both clearly subclonal (CCF < "
              f"{CCF_SUB}) and sitting inside a gain large enough to be near-clonal at its "
              "tumour's purity. Report as 'the design found no unambiguous case', not as evidence "
              "for either direction.")
    print(f"\nwrote {OUT/'timing_clonality_order.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
