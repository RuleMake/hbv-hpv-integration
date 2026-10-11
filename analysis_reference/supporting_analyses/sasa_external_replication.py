"""Evaluate duplication enrichment at independently collected Sasa HPV coordinates.

The source is Supplementary Data 2 from Sasa et al. (2025), Nature Communications
16:1052, DOI 10.1038/s41467-025-56150-z. PCAWG-labelled samples are excluded to avoid
overlap with the structural-variant background. OC and OU samples are retained.

Coordinates are converted from GRCh38 to GRCh37 using pyliftover, retaining only
unique, chromosome-consistent mappings. The analysis uses the same event counting,
matching strata and variant background as the primary spatial comparison."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402

SRC = Path(r"upstream_project")
DATA = Path(r"source_project\data\sasa2025\MOESM5.xlsx")
CHAIN = SRC / "data" / "raw" / "annotation" / "hg38ToHg19.over.chain.gz"
OUT = Path(r"source_project\results")
WINDOW = 100_000
BIN = F.BIN_SIZE


def load_sasa(keep_pcawg: bool = False) -> pd.DataFrame:
    d = pd.read_excel(DATA, sheet_name="Data S2", header=2)
    d = d[d["ID"].notna()].copy()
    n_all = len(d)
    is_pc = d["ID"].astype(str).str.startswith("PCAWG")
    if not keep_pcawg:
        d = d[~is_pc]
    print(f"[sasa]   {n_all} breakpoints / {d['ID'].nunique()} samples retained "
          f"({int(is_pc.sum())} PCAWG-labelled breakpoints "
          f"{'kept' if keep_pcawg else 'DROPPED for independence'})")
    return d


def lift(d: pd.DataFrame) -> pd.DataFrame:
    from pyliftover import LiftOver
    lo = LiftOver(str(CHAIN))
    rows, fail, multi, xchrom = [], 0, 0, 0
    for _, r in d.iterrows():
        c38 = str(r["Human Chromosome"])
        p38 = int(r["Human Position"])
        res = lo.convert_coordinate(c38, p38)
        if not res:
            fail += 1
            continue
        if len(res) > 1:
            multi += 1
            continue
        c19, p19 = res[0][0], int(res[0][1])
        if c19 != c38:
            xchrom += 1
            continue
        rows.append(dict(sample=r["ID"], chrom=c19, pos0=p19,
                         human_strand=r["Human Strand"], hpv_strand=r["HPV Strand"],
                         hpv_gene=r["HPV Gene"], hpv_pos=r["HPV Position"]))
    out = pd.DataFrame(rows)
    print(f"[lift]   {len(out)}/{len(d)} lifted uniquely and chromosome-consistently "
          f"(unmapped {fail}, multi-mapping {multi}, chromosome-switch {xchrom})")
    return out


def attach_strata(d: pd.DataFrame) -> pd.DataFrame:
    feat = pd.read_parquet(F.FEATURES)
    feat = feat.assign(chrom=feat["chrom"].astype(str).map(
        lambda c: c if c.startswith("chr") else f"chr{c}"))
    feat = feat[(feat["stratum"] >= 0) & (feat["chrom"] != "chrY")]
    d = d.copy()
    d["bin_start"] = (d["pos0"] // BIN) * BIN
    d = d.merge(feat[["chrom", "bin_start", "stratum"]], on=["chrom", "bin_start"], how="left")
    before = len(d)
    d = d[d["stratum"].notna()].copy()
    d["within_bin_offset"] = d["pos0"] - d["bin_start"]
    d["stratum"] = d["stratum"].astype(int)
    print(f"[strata] {len(d)}/{before} coordinates land in an eligible 18-stratum bin")
    return d.sort_values(["chrom", "pos0"]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--keep-pcawg", action="store_true")
    a = ap.parse_args()

    idx, _ = F.load_events()
    nul = F.Nuller()
    sub = attach_strata(lift(load_sasa(a.keep_pcawg)))

    labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                               sub["pos0"].to_numpy(np.int64), WINDOW)
    n_clu = int(len(np.unique(labels)))
    print(f"[design] {len(sub)} coordinates -> {n_clu} clusters "
          f"(multiplicity {len(sub)/n_clu:.2f}); {sub['sample'].nunique()} samples")

    rows = []
    t0 = time.time()

    rng = np.random.default_rng(a.seed)
    o, nl = F.est_independent(idx, nul, rng, sub, WINDOW, a.n_perm)
    rows.append(dict(catalogue="Sasa2025_HPV16", estimator="naive_published",
                     n_coords=len(sub), n_clusters=n_clu, **F.summarize(o, nl)))

    lead = pd.Series(labels).drop_duplicates(keep="first").index
    rng = np.random.default_rng(a.seed)
    o, nl = F.est_independent(idx, nul, rng, sub.loc[lead].reset_index(drop=True),
                              WINDOW, a.n_perm)
    rows.append(dict(catalogue="Sasa2025_HPV16", estimator="a_dedup_coords",
                     n_coords=len(lead), n_clusters=n_clu, **F.summarize(o, nl)))

    rng = np.random.default_rng(a.seed)
    o, nl, _, s_obs, s_null = F.est_block(idx, nul, rng, sub, WINDOW, a.n_perm,
                                          labels, unit="uniqSV")
    rows.append(dict(catalogue="Sasa2025_HPV16", estimator="bc_uniqSV_blockNull",
                     n_coords=len(sub), n_clusters=n_clu, **F.summarize(o, nl)))
    rows.append(dict(catalogue="Sasa2025_HPV16", estimator="f_clusterSignTest",
                     n_coords=len(sub), n_clusters=n_clu, **F.summarize(s_obs, s_null)))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "sasa_external_replication.csv", index=False)
    pd.set_option("display.width", 240)
    print()
    print(df[["estimator", "n_coords", "n_clusters", "observed", "expected", "fold",
              "fold_ci_lo", "fold_ci_hi", "emp_p"]].to_string(index=False))
    print(f"\nTian HPV reference (same estimator, same pool, same null): "
          f"bc_uniqSV_blockNull 1.1157 (p < 5e-5), sign test 1.1203 (p < 5e-5)")
    print(f"\nwrote {OUT / 'sasa_external_replication.csv'}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
