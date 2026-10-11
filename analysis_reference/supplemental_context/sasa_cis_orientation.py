"""Analyzes viral orientation, viral breakpoint annotation and same-sample structural variants in the Sasa catalogue.
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
import amplification_and_histology_control as A  # noqa: E402
from minimal_bigwig_reader import BigWig  # noqa: E402

SRC = Path(r"source_project")
OUT = J.OUT
SASA = Path(r"source_project\data\sasa2025\MOESM5.xlsx")
ZAP = SRC / "data/raw/viral_integration/Zapatka2020_PCAWG_SupplTables_MOESM3.xlsx"
CHAIN = SRC / "data/raw/annotation/hg38ToHg19.over.chain.gz"
RFD_BW = SRC / "data/raw/okseq/Petryk2016_HeLa_OkaSeq_RFD_hg19.bw"
SITES = SRC / "data/processed/integration_sites_with_sv.csv"
BIN = F.BIN_SIZE
SEED = 7
JUNCS = (0, 10_000, 50_000)
WINDOWS = (100_000, 1_000_000)
# looping needs an intact origin (URR) and E1 helicase
GENE_GROUPS = {
    "EARLY_disrupted_E1_E2": ("E1", "E2", "E1orE2"),
    "LATE_disrupted_L1_L2_E5": ("L1", "L2", "E5", "E5_L2"),
    "URR_E6_E7": ("URR/LCR", "E6", "E7"),
}


def load_sasa() -> pd.DataFrame:
    from pyliftover import LiftOver
    d = pd.read_excel(SASA, sheet_name="Data S2", skiprows=2)
    d = d[d["ID"].notna()].copy()
    d["pos38"] = pd.to_numeric(d["Human Position"], errors="coerce")
    d["chrom"] = d["Human Chromosome"].astype(str)
    lo = LiftOver(str(CHAIN))
    pos19 = []
    for c, p in zip(d["chrom"], d["pos38"]):
        if not np.isfinite(p):
            pos19.append(np.nan)
            continue
        r = lo.convert_coordinate(c, int(p))
        pos19.append(int(r[0][1]) if r and r[0][0] == c else np.nan)
    d["pos0"] = pos19
    d = d[d["pos0"].notna()].copy()
    d["pos0"] = d["pos0"].astype(np.int64)
    hs = d["Human Strand"].astype(str).map({"+": 1, "-": -1}).fillna(0)
    vs = d["HPV Strand"].astype(str).map({"+": 1, "-": -1}).fillna(0)
    d["T"] = (hs * vs).astype(int)
    d["gene_group"] = "OTHER"
    for g, labs in GENE_GROUPS.items():
        d.loc[d["HPV Gene"].astype(str).isin(labs), "gene_group"] = g
    print(f"[sasa]    {len(d)} breakpoints lifted to hg19 | "
          f"{d.gene_group.value_counts().to_dict()} | T {d['T'].value_counts().to_dict()}")
    return d


def map_pcawg(sasa: pd.DataFrame) -> dict[str, str]:
    z = pd.read_excel(ZAP, sheet_name="17 Integration")
    z = z[z["virus"].astype(str).str.contains("papilloma|HPV", case=False, na=False)].copy()
    z["chrom"] = z["Chromosome"].astype(str)
    z["pos"] = pd.to_numeric(z["chromStart"], errors="coerce")
    s = sasa[sasa["ID"].astype(str).str.startswith("PCAWG")]
    best: dict[str, tuple[str, int]] = {}
    for sid, g in s.groupby("ID"):
        for proj, zg in z.groupby("project"):
            n = sum(int(((zg.chrom == c) & ((zg.pos - p).abs() <= 1000)).any())
                    for c, p in zip(g.chrom, g.pos0))
            if n and (sid not in best or n > best[sid][1]):
                best[str(sid)] = (str(proj), n)
    used, out = set(), {}
    for sid, (proj, n) in sorted(best.items(), key=lambda kv: -kv[1][1]):
        if proj in used:
            print(f"[map]     COLLISION: {sid} also matches {proj}; skipped")
            continue
        used.add(proj)
        out[sid] = proj
        print(f"[map]     {sid} -> {proj}  ({n} coordinates matched within 1 kb)")
    return out


def load_donor_sv(path: str, chrom: str) -> np.ndarray:
    with gzip.open(path, "rt") as fh:
        d = pd.read_csv(fh, sep="\t", low_memory=False)
    if d.empty:
        return np.empty((2, 0), dtype=np.int64)
    m = ((d["chrom1"].astype(str) == d["chrom2"].astype(str))
         & ("chr" + d["chrom1"].astype(str).str.replace("chr", "", regex=False) == chrom))
    if not m.any():
        return np.empty((2, 0), dtype=np.int64)
    s1 = d.loc[m, "start1"].to_numpy(np.int64)
    s2 = d.loc[m, "start2"].to_numpy(np.int64)
    return np.stack([np.minimum(s1, s2), np.maximum(s1, s2)])


def run_locus_strata(sasa: pd.DataFrame, n_perm: int) -> pd.DataFrame:
    """Estimate joint-matched enrichment by the HPV gene disrupted at the viral breakpoint."""
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), 100_000)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)
    idx, _ = F.load_events()
    rows = []
    for tag in ["ALL"] + list(GENE_GROUPS):
        s = sasa if tag == "ALL" else sasa[sasa.gene_group == tag]
        s = s.drop_duplicates(["chrom", "pos0"]).sort_values(["chrom", "pos0"])
        in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                            for c, p in zip(s["chrom"], s["pos0"])])
        s = s[~in_mask].copy()
        k = pd.MultiIndex.from_arrays([s["chrom"].astype(str).to_numpy(),
                                       (s["pos0"].to_numpy(np.int64) // BIN) * BIN])
        s["stratum"] = nul.lookup.reindex(k).to_numpy()
        s = s[s["stratum"].notna()].reset_index(drop=True)
        if len(s) < 25:
            print(f"[strata]  {tag}: only {len(s)} eligible coordinates -- skipped")
            continue
        s["stratum"] = s["stratum"].astype(int)
        s["within_bin_offset"] = s["pos0"] - (s["pos0"] // BIN) * BIN
        labels = F.leader_clusters(s["chrom"].to_numpy(str), s["pos0"].to_numpy(np.int64), 100_000)
        rng = np.random.default_rng(SEED)
        o, nl, _, so, sn = F.est_block(idx, nul, rng, s, 100_000, n_perm, labels, unit="uniqSV")
        for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                         ("f_clusterSignTest", F.summarize(so, sn))):
            rows.append(dict(stratum=tag, estimator=est, n_coords=len(s),
                             n_clusters=int(len(np.unique(labels))), **res))
        print(f"  [{tag:26s}] n={len(s):4d} bc={rows[-2]['fold']:.4f} "
              f"[{rows[-2]['fold_ci_lo']:.4f},{rows[-2]['fold_ci_hi']:.4f}] "
              f"p={rows[-2]['emp_p']:.5g} | sign={rows[-1]['fold']:.4f} "
              f"p={rows[-1]['emp_p']:.5g}", flush=True)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    a = ap.parse_args()
    t0 = time.time()

    sasa = load_sasa()
    # Exclude PCAWG-labelled subjects from the external-coordinate comparison because
    # their own structural variants contribute to the shared event background.
    indep = sasa[~sasa["ID"].astype(str).str.startswith("PCAWG")]
    print("\n=== Locus-level enrichment by disrupted HPV gene ===")
    print(f"[indep]   dropping {len(sasa)-len(indep)} PCAWG-labelled breakpoints "
          f"(their own SVs are in the event pool); {len(indep)} remain")
    gs = run_locus_strata(indep, a.n_perm)
    gs.insert(0, "arm", "PCAWG-labelled samples EXCLUDED")
    gs_all = run_locus_strata(sasa, a.n_perm)
    gs_all.insert(0, "arm", "all Sasa samples (NOT independent)")
    gs = pd.concat([gs, gs_all], ignore_index=True)
    gs.to_csv(OUT / "sasa_viral_gene_stratum.csv", index=False)

    print("\n=== Same-sample orientation dataset ===")
    m = map_pcawg(sasa)
    sites = pd.read_csv(SITES)
    path_of = {}
    for r in sites.itertuples():
        path_of[str(r.project)] = str(r.path)
    bw = BigWig(RFD_BW)
    binned, _, _ = bw.binned_means(50_000)
    bw.close()
    rfd = {c: np.sign(np.nan_to_num(v, nan=0.0)).astype(np.int8) for c, (s, v) in binned.items()}

    rows = []
    sv_cache: dict[tuple[str, str], np.ndarray] = {}
    for sid, proj in m.items():
        p = path_of.get(proj)
        if p is None:
            print(f"[cis]     {sid} -> {proj}: no SV path on disk, skipped")
            continue
        g = sasa[sasa["ID"].astype(str) == sid]
        for r in g.itertuples():
            key = (p, r.chrom)
            if key not in sv_cache:
                sv_cache[key] = load_donor_sv(p, r.chrom)
            iv = sv_cache[key]
            rr = rfd.get(r.chrom)
            R = int(rr[min(int(r.pos0) // 50_000, len(rr) - 1)]) if rr is not None else 0
            rec = dict(sasa_id=sid, donor=proj, chrom=r.chrom, pos0=int(r.pos0),
                       hpv_gene=str(r._asdict().get("_9", "")), gene_group=r.gene_group,
                       T=int(r.T), R=R, TR=int(r.T) * R,
                       clonality=str(getattr(r, "Clonality", "")),
                       ccf=float(getattr(r, "CCF", np.nan)) if pd.notna(
                           getattr(r, "CCF", np.nan)) else np.nan)
            if iv.shape[1]:
                dl = np.abs(iv[0] - int(r.pos0))
                dh = np.abs(iv[1] - int(r.pos0))
                for jx in JUNCS:
                    ok = np.minimum(dl, dh) > jx
                    for w in WINDOWS:
                        rec[f"n_sv_{w//1000}k_j{jx}"] = int((((dl <= w) | (dh <= w)) & ok).sum())
            else:
                for jx in JUNCS:
                    for w in WINDOWS:
                        rec[f"n_sv_{w//1000}k_j{jx}"] = 0
            rows.append(rec)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "sasa_cis_orientation.csv", index=False)
    print(f"[cis]     {len(d)} breakpoints across {d.donor.nunique()} donors with SV calls")
    print(f"[cis]     T x R: {d.TR.value_counts().to_dict()}  (TR=-1 head-on, +1 co-directional "
          f"under the stated convention)")

    res = []
    for jx in JUNCS:
        for w in WINDOWS:
            col = f"n_sv_{w//1000}k_j{jx}"
            ho = d.loc[d.TR == -1, col].to_numpy(float)
            cd = d.loc[d.TR == 1, col].to_numpy(float)
            if len(ho) < 10 or len(cd) < 10:
                continue
            # MDE FIRST: smallest ratio of means resolvable at 80% power, two-sided alpha 0.05
            sp = np.sqrt((ho.var(ddof=1) / len(ho)) + (cd.var(ddof=1) / len(cd)))
            mde = 1.0 + (1.96 + 0.84) * sp / max(cd.mean(), 1e-9)
            u = stats.mannwhitneyu(ho, cd, alternative="two-sided")
            res.append(dict(window_kb=w // 1000, junction_excl_kb=jx // 1000,
                            n_headon=len(ho), n_codir=len(cd),
                            mean_headon=round(float(ho.mean()), 3),
                            mean_codir=round(float(cd.mean()), 3),
                            ratio=round(float(ho.mean() / cd.mean()), 4) if cd.mean() else np.nan,
                            mde_ratio_80pct=round(float(mde), 3),
                            mwu_p=round(float(u.pvalue), 5)))
    r = pd.DataFrame(res)
    r.to_csv(OUT / "sasa_cis_orientation_test.csv", index=False)
    pd.set_option("display.width", 250)
    print("\n=== Head-on and co-directional integrations in the same genome ===")
    print(r[["window_kb", "junction_excl_kb", "n_headon", "n_codir", "mde_ratio_80pct",
             "mean_headon", "mean_codir", "ratio", "mwu_p"]].to_string(index=False))
    print("\n=== Locus-level results by disrupted HPV gene ===")
    print(gs[["arm", "stratum", "estimator", "n_coords", "n_clusters", "fold", "fold_ci_lo",
              "fold_ci_hi", "emp_p"]].to_string(index=False))
    print(f"\nwrote {OUT/'sasa_cis_orientation.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
