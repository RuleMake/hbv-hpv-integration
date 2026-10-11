"""Summarizes published long-read HBV junction classes and paired tissue contrasts.
"""

from __future__ import annotations

import argparse
import itertools
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import joint_correction as J  # noqa: E402

SRC = Path(r"source_project")
WANG = SRC / "data" / "external" / "wang2025" / "Supplemental_Tables.xlsx"
OUT = J.OUT
BURDEN = OUT / "hbv_sample_level_td_burden.csv"
SEED = 7
NBOOT = 20000


def parse_locus(s) -> tuple[str | None, int | None]:
    if pd.isna(s):
        return (None, None)
    m = re.search(r"Chr([0-9XYM]+)\s*:\s*([0-9]+)", str(s), re.I)
    return (m.group(1).upper(), int(m.group(2))) if m else (None, None)


def load() -> pd.DataFrame:
    d = pd.read_excel(WANG, sheet_name="Sup. Table S2", header=2)
    d.columns = [str(c).strip() for c in d.columns]
    d["c1"], d["p1"] = zip(*d["HBV integration locus ID1"].map(parse_locus))
    d["c2"], d["p2"] = zip(*d["HBV integration locus ID2"].map(parse_locus))
    d["sample"] = d["Sample"].astype(str).str.strip()
    d["tissue"] = np.where(d["sample"].str.startswith("T"), "Tumour", "Non-tumour")
    d["patient"] = d["sample"].str.extract(r"(\d+)")[0]
    d["n_loci"] = d[["c1", "c2"]].notna().sum(axis=1)
    d["pattern"] = d["HBV integration pattern"].astype(str).str.strip()
    d["cfs1"] = ~d["CFS_1"].astype(str).str.strip().isin(["NA", "nan", "None", ""])
    d["log2cn"] = pd.to_numeric(d["LOG2_RATIO_MEAN_1"], errors="coerce")
    # canonical unordered locus pair, for deduplication
    def key(r):
        if r.n_loci < 2:
            return None
        a, b = (r.c1, r.p1), (r.c2, r.p2)
        return f"{r['sample']}|" + "|".join(f"{c}:{p}" for c, p in sorted([a, b]))
    d["junction_key"] = d.apply(key, axis=1)
    d["locus_key"] = d["sample"] + "|" + d["c1"].astype(str) + ":" + d["p1"].astype(str)
    return d


def cluster_boot(per_sample: pd.DataFrame, num: str, den: str, rng) -> tuple[float, float, float]:
    """Share with a cluster-bootstrap CI: resample SAMPLES with replacement."""
    n, dd = per_sample[num].to_numpy(float), per_sample[den].to_numpy(float)
    obs = n.sum() / dd.sum() if dd.sum() else np.nan
    k = len(n)
    draws = np.empty(NBOOT)
    for i in range(NBOOT):
        ix = rng.integers(0, k, k)
        s = dd[ix].sum()
        draws[i] = n[ix].sum() / s if s else np.nan
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return float(obs), float(lo), float(hi)


def paired_signflip(pairs: list[tuple[float, float]]) -> tuple[float, float, int]:
    """Exact 2^k sign-flip permutation on within-patient tumour-minus-normal differences."""
    diffs = np.array([t - n for t, n in pairs], float)
    diffs = diffs[np.isfinite(diffs)]
    k = len(diffs)
    if k == 0:
        return (np.nan, np.nan, 0)
    obs = float(diffs.mean())
    hits = sum(1 for signs in itertools.product([-1, 1], repeat=k)
               if abs((diffs * np.array(signs)).mean()) >= abs(obs) - 1e-12)
    return obs, hits / (2 ** k), k


def exact_unpaired(values: np.ndarray, is_t: np.ndarray) -> float:
    """Exact permutation of the tissue label over samples."""
    obs = values[is_t].mean() - values[~is_t].mean()
    idx = np.arange(len(values))
    nt = int(is_t.sum())
    hits = tot = 0
    for combo in itertools.combinations(idx, nt):
        m = np.zeros(len(values), bool)
        m[list(combo)] = True
        if abs(values[m].mean() - values[~m].mean()) >= abs(obs) - 1e-12:
            hits += 1
        tot += 1
    return hits / tot


def mde(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return np.nan
    return float((1.96 + 0.84) * np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    rows: list[dict] = []

    d = load()
    print(f"[in]      {len(d)} molecules, {d['sample'].nunique()} samples "
          f"({d[d.tissue=='Tumour']['sample'].nunique()} T / "
          f"{d[d.tissue=='Non-tumour']['sample'].nunique()} N)")

    # ---- counting unit: molecules -> unique junctions
    two = d[d.n_loci == 2].copy()
    uj = two.drop_duplicates(subset=["junction_key"]).copy()
    ul = d.drop_duplicates(subset=["locus_key"]).copy()
    print(f"[unit]    two-locus molecules {len(two)} -> unique junctions {len(uj)} "
          f"(multiplicity {len(two)/max(len(uj),1):.2f}); "
          f"all molecules {len(d)} -> unique loci {len(ul)} "
          f"(multiplicity {len(d)/max(len(ul),1):.2f})")
    rows.append(dict(analysis="counting_unit", unit="two_locus_molecule", n=len(two)))
    rows.append(dict(analysis="counting_unit", unit="unique_junction", n=len(uj),
                     multiplicity=round(len(two) / max(len(uj), 1), 4)))
    rows.append(dict(analysis="counting_unit", unit="all_molecules", n=len(d)))
    rows.append(dict(analysis="counting_unit", unit="unique_locus", n=len(ul),
                     multiplicity=round(len(d) / max(len(ul), 1), 4)))

    for lab, t in (("per_molecule", two), ("per_unique_junction", uj)):
        t = t.copy()
        t["inter"] = (t.c1 != t.c2)
        t["dist"] = np.where(t.c1 == t.c2, (t.p2 - t.p1).abs(), np.nan)

        # ---- readout 1: class composition, tumour arm, cluster bootstrap over samples
        tt = t[t.tissue == "Tumour"]
        ps = tt.groupby("sample").agg(n_inter=("inter", "sum"), n=("inter", "size")).reset_index()
        obs, lo, hi = cluster_boot(ps, "n_inter", "n", rng)
        dd = tt.dist.dropna()
        rows.append(dict(analysis="class_composition", unit=lab, tissue="Tumour",
                         n_junctions=len(tt), n_samples=len(ps),
                         frac_inter_chromosomal=round(obs, 4),
                         boot_ci_lo=round(lo, 4), boot_ci_hi=round(hi, 4),
                         n_intra=int(len(dd)),
                         median_intra_dist_bp=int(dd.median()) if len(dd) else np.nan,
                         frac_intra_gt_100kb=round(float((dd > 1e5).mean()), 4) if len(dd) else np.nan,
                         frac_intra_gt_1Mb=round(float((dd > 1e6).mean()), 4) if len(dd) else np.nan,
                         frac_intra_gt_10Mb=round(float((dd > 1e7).mean()), 4) if len(dd) else np.nan))
        print(f"[class]   {lab:22s} tumour inter-chromosomal {obs:.4f} [{lo:.4f},{hi:.4f}] "
              f"(n={len(tt)}, {len(ps)} samples)")

        # ---- readout 2: tumour vs non-tumour, exact paired sign-flip
        ps_all = t.groupby(["sample", "patient", "tissue"]).agg(
            n_inter=("inter", "sum"), n=("inter", "size")).reset_index()
        ps_all["frac"] = ps_all.n_inter / ps_all.n
        pairs = []
        for pid, g in ps_all.groupby("patient"):
            a = g[g.tissue == "Tumour"]["frac"]
            b = g[g.tissue == "Non-tumour"]["frac"]
            if len(a) and len(b):
                pairs.append((float(a.iloc[0]), float(b.iloc[0])))
        obs_d, pp, k = paired_signflip(pairs)
        rows.append(dict(analysis="inter_frac_T_vs_N", unit=lab, n_pairs=k,
                         mean_diff=round(obs_d, 4) if np.isfinite(obs_d) else np.nan,
                         exact_signflip_p=round(pp, 5) if np.isfinite(pp) else np.nan,
                         min_attainable_p=round(2 / 2 ** k, 5) if k else np.nan))

    # ---- readout 1b: reference -- translocation share of all SVs in HBV+ PCAWG donors
    if BURDEN.exists():
        b = pd.read_csv(BURDEN).set_index("metric")
        if {"tra_all", "sv_all"} <= set(b.index):
            tra = float(b.loc["tra_all", "mean_hbv_pos"])
            sv = float(b.loc["sv_all", "mean_hbv_pos"])
            rows.append(dict(analysis="reference_sv_composition",
                             source="hbv_sample_level_td_burden.csv (HBV+ donors, n=56)",
                             mean_tra=round(tra, 3), mean_sv=round(sv, 3),
                             tra_share_of_all_sv=round(tra / sv, 4)))
            print(f"[ref]     TRA share of all SVs in HBV+ PCAWG donors: {tra/sv:.4f}")

    # ---- readout 2b: two-locus fraction and pattern composition, T vs N
    ps = d.groupby(["sample", "patient", "tissue"]).agg(
        n_two=("n_loci", lambda s: int((s == 2).sum())), n=("n_loci", "size")).reset_index()
    ps["frac_two"] = ps.n_two / ps.n
    pairs = [(float(g[g.tissue == "Tumour"].frac_two.iloc[0]),
              float(g[g.tissue == "Non-tumour"].frac_two.iloc[0]))
             for _, g in ps.groupby("patient")
             if (g.tissue == "Tumour").any() and (g.tissue == "Non-tumour").any()]
    obs_d, pp, k = paired_signflip(pairs)
    vals = ps.frac_two.to_numpy(float)
    ist = (ps.tissue == "Tumour").to_numpy()
    rows.append(dict(analysis="two_locus_fraction_T_vs_N", n_pairs=k,
                     mean_T=round(float(vals[ist].mean()), 4),
                     mean_N=round(float(vals[~ist].mean()), 4),
                     mean_diff=round(obs_d, 4) if np.isfinite(obs_d) else np.nan,
                     exact_signflip_p=round(pp, 5) if np.isfinite(pp) else np.nan,
                     exact_unpaired_p=round(exact_unpaired(vals, ist), 5),
                     mde_80pct=round(mde(vals[ist], vals[~ist]), 4)))

    for pat in sorted(d.pattern.unique()):
        ps2 = d.assign(hit=(d.pattern == pat)).groupby(
            ["sample", "patient", "tissue"]).agg(h=("hit", "sum"), n=("hit", "size")).reset_index()
        ps2["frac"] = ps2.h / ps2.n
        pr = [(float(g[g.tissue == "Tumour"].frac.iloc[0]),
               float(g[g.tissue == "Non-tumour"].frac.iloc[0]))
              for _, g in ps2.groupby("patient")
              if (g.tissue == "Tumour").any() and (g.tissue == "Non-tumour").any()]
        o, p_, kk = paired_signflip(pr)
        v = ps2.frac.to_numpy(float)
        m = (ps2.tissue == "Tumour").to_numpy()
        rows.append(dict(analysis="pattern_T_vs_N", pattern=pat, n_pairs=kk,
                         mean_T=round(float(v[m].mean()), 4), mean_N=round(float(v[~m].mean()), 4),
                         mean_diff=round(o, 4) if np.isfinite(o) else np.nan,
                         exact_signflip_p=round(p_, 5) if np.isfinite(p_) else np.nan,
                         mde_80pct=round(mde(v[m], v[~m]), 4)))

    # ---- readout 3: common fragile sites, two-locus vs single-locus
    ps3 = d.groupby(["sample", "n_loci"]).agg(h=("cfs1", "sum"), n=("cfs1", "size")).reset_index()
    for k_ in (1, 2):
        s = ps3[ps3.n_loci == k_]
        o, lo, hi = cluster_boot(s, "h", "n", rng)
        rows.append(dict(analysis="cfs_by_n_loci", n_loci=k_, n_molecules=int(s.n.sum()),
                         n_samples=len(s), frac_cfs=round(o, 4),
                         boot_ci_lo=round(lo, 4), boot_ci_hi=round(hi, 4)))
    w = ps3.pivot(index="sample", columns="n_loci", values="h").div(
        ps3.pivot(index="sample", columns="n_loci", values="n")).dropna()
    if len(w) >= 2:
        o, p_, kk = paired_signflip(list(zip(w[2], w[1])))
        rows.append(dict(analysis="cfs_two_minus_one", n_samples=kk,
                         mean_diff=round(o, 4), exact_signflip_p=round(p_, 5),
                         min_attainable_p=round(2 / 2 ** kk, 5)))
        print(f"[cfs]     two-minus-one-locus CFS share {o:+.4f}, exact P={p_:.4f} ({kk} samples)")

    # ---- readout 4: detection-bias control on the authors' CNV column
    cn = d.dropna(subset=["log2cn"])
    a = cn.loc[cn.n_loci == 2, "log2cn"].to_numpy(float)
    b_ = cn.loc[cn.n_loci == 1, "log2cn"].to_numpy(float)
    per = cn.groupby(["sample", "n_loci"]).log2cn.median().unstack().dropna()
    o = p_ = kk = np.nan
    if len(per) >= 2:
        o, p_, kk = paired_signflip(list(zip(per[2], per[1])))
    rows.append(dict(analysis="cnv_detection_bias", n_two=len(a), n_one=len(b_),
                     median_log2cn_two=round(float(np.median(a)), 4),
                     median_log2cn_one=round(float(np.median(b_)), 4),
                     pooled_mwu_p=round(float(stats.mannwhitneyu(a, b_).pvalue), 6),
                     sample_level_mean_diff=round(o, 4) if np.isfinite(o) else np.nan,
                     exact_signflip_p=round(p_, 5) if np.isfinite(p_) else np.nan,
                     n_samples=kk, mde_80pct=round(mde(a, b_), 4)))

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "wang2025_longread_classes.csv", index=False)

    smp = (two.assign(inter=two.c1 != two.c2)
              .groupby(["sample", "tissue"])
              .agg(n_two_locus=("inter", "size"), n_inter=("inter", "sum")).reset_index())
    smp["frac_inter"] = (smp.n_inter / smp.n_two_locus).round(4)
    tot = d.groupby("sample").size().rename("n_molecules")
    smp = smp.merge(tot, on="sample").sort_values(["tissue", "sample"])
    smp.to_csv(OUT / "wang2025_longread_classes_samples.csv", index=False)
    uj.assign(inter=uj.c1 != uj.c2)[
        ["sample", "tissue", "patient", "c1", "p1", "c2", "p2", "inter", "pattern",
         "cfs1", "log2cn"]].to_csv(OUT / "wang2025_longread_classes_junctions.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    for an in out.analysis.unique():
        s = out[out.analysis == an].dropna(axis=1, how="all")
        print(f"\n=== {an} ===")
        print(s.drop(columns=["analysis"]).to_string(index=False))
    print(f"\n=== per sample ===\n{smp.to_string(index=False)}")
    print(f"\nwrote {OUT/'wang2025_longread_classes.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
