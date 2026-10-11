"""Tests expression at integration-associated genes and nearby genomic contexts.
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

SRC = Path(r"source_project")
OUT = J.OUT
BIN = F.BIN_SIZE
WINDOW = 100_000
FPKM_ON = 1.0
SEED = 7
FPKM = SRC / "data/raw/pcawg_expression/tophat_star_fpkm.v2_aliquot_gl.tsv.gz"
GENES = SRC / "data/interim/gencode_v19_genes.parquet"
SITES = SRC / "data/processed/integration_sites_with_sv.csv"
SHEET = SRC / "data/raw/viral_integration/pcawg_sample_sheet.v1.4.2016-09-14.tsv"
HIST = SRC / "data/raw/clinical/pcawg_specimen_histology_August2016_v9.xlsx"
CNBYGENE = SRC / "data/raw/pcawg_cnv/all_samples.consensus_CN.by_gene.170214.txt.gz"


def load_genes() -> dict[str, tuple]:
    g = pd.read_parquet(GENES)
    g = g[g["gene_type"] == "protein_coding"].copy()
    g["chrom"] = g["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
    out = {}
    for c, s in g.groupby("chrom"):
        s = s.sort_values("start")
        out[c] = (s["start"].to_numpy(np.int64), s["end"].to_numpy(np.int64),
                  s["gene_id_base"].to_numpy(str), s["gene_name"].to_numpy(str),
                  s["strand"].to_numpy(str))
    print(f"[genes]   {len(g):,} protein-coding genes on {len(out)} chromosomes")
    return out


def near(gx, pos: int, win: int) -> np.ndarray:
    """indices of genes whose BODY overlaps [pos-win, pos+win]"""
    st, en = gx[0], gx[1]
    j = np.searchsorted(st, pos + win, side="right")
    if j == 0:
        return np.empty(0, dtype=int)
    i = np.flatnonzero(en[:j] >= pos - win)
    return i


def donor_rna_map() -> dict[str, str]:
    with gzip.open(FPKM, "rt") as f:
        cols = f.readline().rstrip("\n").split("\t")
    have = set(cols[1:])
    ss = pd.read_csv(SHEET, sep="\t", low_memory=False)
    r = ss[ss["library_strategy"].astype(str).str.upper() == "RNA-SEQ"]
    m = {}
    for d, al in zip(r["icgc_donor_id"], r["aliquot_id"]):
        if al in have and d not in m:
            m[str(d)] = str(al)
    print(f"[rna]     {len(m)} donors have an RNA aliquot in the FPKM matrix")
    return m


def wgs_aliquot_map() -> dict[str, str]:
    ss = pd.read_csv(SHEET, sep="\t", low_memory=False)
    w = ss[(ss["library_strategy"].astype(str).str.upper() == "WGS")
           & (~ss["dcc_specimen_type"].astype(str).str.contains("Normal", na=False))]
    return {str(d): str(a) for d, a in zip(w["icgc_donor_id"], w["aliquot_id"])}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=2000)
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    sites = pd.read_csv(SITES)
    sites["chrom"] = sites["Chromosome"].astype(str)
    sites["pos0"] = sites["chromStart"].astype(np.int64)
    sites["donor"] = sites["donor"].astype(str)

    h = pd.read_excel(HIST)[["icgc_donor_id", "histology_abbreviation"]].drop_duplicates(
        "icgc_donor_id")
    hist = dict(zip(h["icgc_donor_id"].astype(str), h["histology_abbreviation"].astype(str)))
    rna = donor_rna_map()
    wgs = wgs_aliquot_map()
    sites["hist"] = sites["donor"].map(hist)
    sites["rna"] = sites["donor"].map(rna)
    print(f"[sites]   {len(sites)} PCAWG integration sites, {sites.donor.nunique()} donors; "
          f"{int(sites['rna'].notna().sum())} sites in donors with RNA-seq")
    print(f"[sites]   histology: {sites.loc[sites['rna'].notna(),'hist'].value_counts().to_dict()}")

    # ---- load only the FPKM columns we need: every donor of the same histologies
    hist_needed = set(sites.loc[sites["rna"].notna(), "hist"].dropna())
    peers = {d: al for d, al in rna.items() if hist.get(d) in hist_needed}
    keep = ["feature"] + sorted(set(peers.values()))
    print(f"[rna]     loading {len(keep)-1} aliquots of histologies {sorted(hist_needed)}")
    fp = pd.read_csv(FPKM, sep="\t", usecols=keep)
    fp["gid"] = fp["feature"].astype(str).str.split(".").str[0]
    fp = fp.drop_duplicates("gid").set_index("gid").drop(columns=["feature"])
    al2donor = {al: d for d, al in peers.items()}
    fp.columns = [al2donor[c] for c in fp.columns]
    print(f"[rna]     FPKM matrix {fp.shape[0]:,} genes x {fp.shape[1]} donors")

    cn = pd.read_csv(CNBYGENE, sep="\t", low_memory=False)
    cn_cols = {d: wgs[d] for d in fp.columns if wgs.get(d) in cn.columns}
    cn = cn.set_index("Gene Symbol")
    print(f"[cnv]     gene-level consensus CN available for {len(cn_cols)}/{fp.shape[1]} donors")

    # ---- per-histology percentile-rank matrix: rank of each donor for each gene, among peers
    # of the same histology. Precomputed so a null draw is an O(1) lookup, which is what makes
    # the within-donor matched-gene control affordable.
    arr = fp.to_numpy(dtype=np.float32)
    rank_mat = np.full(arr.shape, np.nan, dtype=np.float32)
    donors = list(fp.columns)
    col_of = {d: i for i, d in enumerate(donors)}
    row_of = {g: i for i, g in enumerate(fp.index)}
    for hh in sorted(hist_needed):
        cols = [col_of[d] for d in donors if hist.get(d) == hh]
        if len(cols) < 20:
            print(f"[rank]    histology {hh}: only {len(cols)} donors -- excluded from Q2")
            continue
        sub_a = arr[:, cols]
        r = stats.rankdata(sub_a, axis=1, method="average", nan_policy="omit")
        rank_mat[:, cols] = ((r - 0.5) / len(cols)).astype(np.float32)
        print(f"[rank]    histology {hh}: {len(cols)} donors")

    gx = load_genes()
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)

    rows = []
    null_rank_store: list[np.ndarray] = []
    virus_of_store: list[str] = []
    obs_on = 0
    null_on = np.zeros(a.n_perm, dtype=np.int64)
    n_scored = 0
    for r in sites.itertuples():
        if not isinstance(r.rna, str) or r.chrom not in gx:
            continue
        b = (int(r.pos0) // BIN) * BIN
        if (r.chrom, b) in mask:
            continue
        st = nul.lookup.get((r.chrom, b))
        if st is None or (isinstance(st, float) and np.isnan(st)):
            continue
        st = int(st)
        off = int(r.pos0) - b
        n_scored += 1

        g = gx[r.chrom]
        idx = near(g, int(r.pos0), WINDOW)
        gid = g[2][idx]
        e = fp.reindex(gid)[r.donor].to_numpy(float) if len(idx) else np.empty(0)
        e = e[~np.isnan(e)]
        on = bool(e.size and np.nanmax(e) >= FPKM_ON)
        obs_on += int(on)

        # ---- null positions, reused by BOTH questions
        draws = nul.draw(rng, r.chrom, st, off, int(r.pos0), a.n_perm)
        dcol = col_of.get(r.donor)
        null_rank = np.full(a.n_perm, np.nan, dtype=np.float32)
        for i, p in enumerate(draws):
            ii = near(g, int(p), WINDOW)
            if not len(ii):
                continue
            rr = [row_of.get(x) for x in g[2][ii]]
            rr = [x for x in rr if x is not None]
            if not rr:
                continue
            ee = arr[rr, dcol] if dcol is not None else np.array([np.nan])
            if np.any(ee >= FPKM_ON):          # Q1
                null_on[i] += 1
            if dcol is not None:               # Q2 within-donor matched-gene control
                dd = np.maximum(np.maximum(g[0][ii] - int(p), int(p) - g[1][ii]), 0)
                null_rank[i] = rank_mat[rr[int(np.argmin(dd[:len(rr)]))], dcol]

        # Q2: nearest gene, percentile rank of this donor among same-histology peers
        if len(idx) and dcol is not None:
            d0 = np.maximum(np.maximum(g[0][idx] - int(r.pos0), int(r.pos0) - g[1][idx]), 0)
            k = idx[int(np.argmin(d0))]
            gname, gidb = g[3][k], g[2][k]
            grow = row_of.get(gidb)
            if grow is not None and np.isfinite(rank_mat[grow, dcol]):
                cnv = np.nan
                if gname in cn.index and wgs.get(r.donor) in cn.columns:
                    try:
                        cnv = float(np.atleast_1d(cn.loc[gname, wgs[r.donor]])[0])
                    except Exception:  # noqa: BLE001
                        cnv = np.nan
                peer = arr[grow, [col_of[d] for d in donors if hist.get(d) == r.hist]]
                rows.append(dict(virus=r.virus_class, donor=r.donor, hist=r.hist,
                                 chrom=r.chrom, pos0=int(r.pos0), gene=gname,
                                 dist_bp=int(np.min(d0)), fpkm=float(arr[grow, dcol]),
                                 n_peers=int(np.isfinite(peer).sum()),
                                 peer_median=float(np.nanmedian(peer)),
                                 rank=float(rank_mat[grow, dcol]), gene_cn=cnv,
                                 site_expressed=on,
                                 null_rank_mean=float(np.nanmean(null_rank))
                                 if np.isfinite(null_rank).any() else np.nan))
                null_rank_store.append(null_rank)
                virus_of_store.append(r.virus_class)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "integration_site_transcription_pairs.csv", index=False)

    res = []
    q1 = F.summarize(obs_on, null_on)
    res.append(dict(question="Q1_site_has_expressed_gene_within_100kb", subset="ALL",
                    n_sites=n_scored, **q1))
    nr = np.vstack(null_rank_store) if null_rank_store else np.zeros((0, a.n_perm))
    vs = np.array(virus_of_store)
    for v in ("HBV", "HPV"):
        s = d[d.virus == v]
        if len(s) < 10:
            continue
        w = stats.wilcoxon(s["rank"] - 0.5, alternative="two-sided")
        row = dict(question="Q2_rank_vs_same_histology_peers", subset=v, n_pairs=len(s),
                   mean_rank=round(float(s["rank"].mean()), 4),
                   median_rank=round(float(s["rank"].median()), 4),
                   frac_rank_gt_0p9=round(float((s["rank"] > 0.9).mean()), 4),
                   wilcoxon_p=round(float(w.pvalue), 5),
                   median_fpkm=round(float(s["fpkm"].median()), 3),
                   median_peer_fpkm=round(float(s["peer_median"].median()), 3),
                   frac_site_expressed=round(float(s["site_expressed"].mean()), 4))
        # Compare the same donors' expression ranks at matched null positions.
        m = vs == v
        if m.sum() >= 10:
            with np.errstate(invalid="ignore"):
                nulldist = np.nanmean(nr[m], axis=0)
            nulldist = nulldist[np.isfinite(nulldist)]
            obs_mean = float(s["rank"].mean())
            if nulldist.size:
                pp = (int(np.sum(np.abs(nulldist - nulldist.mean())
                                 >= abs(obs_mean - nulldist.mean()))) + 1) / (nulldist.size + 1)
                lo, hi = np.percentile(nulldist, [2.5, 97.5])
                row.update(matched_gene_null_mean=round(float(nulldist.mean()), 4),
                           matched_null_ci_lo=round(float(lo), 4),
                           matched_null_ci_hi=round(float(hi), 4),
                           excess_over_matched=round(obs_mean - float(nulldist.mean()), 4),
                           matched_null_p=round(float(pp), 5))
        c = s.dropna(subset=["gene_cn"])
        if len(c) >= 15:
            sl, ic, rr, pp, _ = stats.linregress(c["gene_cn"], c["rank"])
            resid = c["rank"] - (ic + sl * c["gene_cn"])
            base_r = float(0.5 - (ic + sl * c["gene_cn"]).mean())
            wr = stats.wilcoxon(resid - base_r, alternative="two-sided")
            row.update(n_with_cn=len(c), cn_rank_slope=round(float(sl), 5),
                       cn_rank_p=round(float(pp), 5),
                       resid_mean_rank=round(float(0.5 + resid.mean() - base_r), 4),
                       resid_wilcoxon_p=round(float(wr.pvalue), 5))
        res.append(row)

    out = pd.DataFrame(res)
    out.to_csv(OUT / "integration_site_transcription.csv", index=False)
    pd.set_option("display.width", 260)
    print("\n=== Q1: is the integration site transcribed? (vs joint-corrected matched null) ===")
    print(f"  {obs_on}/{n_scored} sites ({100*obs_on/max(n_scored,1):.1f}%) have a protein-coding "
          f"gene at FPKM >= {FPKM_ON} within +/-{WINDOW//1000} kb, in the same tumour")
    print(f"  matched null expectation {null_on.mean():.1f} "
          f"({100*null_on.mean()/max(n_scored,1):.1f}%)  ->  fold {q1['fold']}, "
          f"p = {q1['emp_p']:.4g}")
    print("\n=== Q2: is transcription ELEVATED at the site? ===")
    print(out[[c for c in out.columns if c not in ("question",)]].to_string(index=False))
    print(f"\nwrote {OUT/'integration_site_transcription.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
