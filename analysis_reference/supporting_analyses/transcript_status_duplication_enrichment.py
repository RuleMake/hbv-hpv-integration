"""Compare duplication enrichment between detected fusion-transcript categories.

Fan et al. (2023), Cell Genomics 3:100211, provides paired DNA/RNA annotations for
762 HPV integration coordinates from 89 cervical cancers. The source labels
Silent and Productive denote the authors' detection categories; an undetected
transcript is not evidence of biological inactivity.

After unique, chromosome-consistent liftover, all coordinates and each transcript
category are scored against the same PCAWG variant background and joint-matched
cluster-block null. The comparison is also restricted to patients with both categories.

Independent cluster bootstraps estimate the log ratio of category-specific folds.
This bootstrap does not model patient pairing. Intervals and minimum detectable
ratios describe the precision of the category contrast."""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402
import joint_correction as J  # noqa: E402
import amplification_and_histology_control as A  # noqa: E402

SRC = Path(r"upstream_project")
DATA = SRC / "data" / "external" / "fan2023" / "Supplemental_Tables.xlsx"
CHAIN = SRC / "data" / "raw" / "annotation" / "hg38ToHg19.over.chain.gz"
OUT = J.OUT
BIN = F.BIN_SIZE
WINDOW = 100_000
SEED = 7
NBOOT = 20000


def load_fan() -> pd.DataFrame:
    d = pd.read_excel(DATA, sheet_name="Table S2", header=2)
    d.columns = [str(c).strip() for c in d.columns]
    d = d.rename(columns={"Human chromosome integrated": "chrom38",
                          "Human integration location": "pos38",
                          "Inregration status": "status"})       # sic -- author's spelling
    d = d[d["Sample"].notna()].copy()
    d["pos38"] = pd.to_numeric(d["pos38"], errors="coerce")
    d["chrom38"] = d["chrom38"].astype(str).str.strip()
    d["patient"] = d["Sample"].astype(str).str.rsplit("-", n=1).str[0]
    d = d[d["pos38"].notna() & d["status"].isin(["Silent", "Productive"])]
    d = d[d["chrom38"].str.match(r"^chr[0-9XY]+$") & (d["chrom38"] != "chrY")]
    n0 = len(d)

    # Define unique sites before summarizing transcript categories: the same
    # physical site can appear on several rows, and a site labelled both ways is uninterpretable.
    d = d.drop_duplicates(subset=["Sample", "chrom38", "pos38", "status"])
    amb = d.duplicated(subset=["Sample", "chrom38", "pos38"], keep=False)
    n_amb = int(amb.sum())
    d = d[~amb].copy()
    print(f"[fan]    {n0} rows -> {len(d)} unique per-sample coordinates "
          f"({n_amb} dropped as labelled both Silent and Productive)")
    print(f"[fan]    {d.status.value_counts().to_dict()} | {d.patient.nunique()} patients | "
          f"HPV {d['HPV type integrated'].value_counts().head(4).to_dict()}")
    return d.reset_index(drop=True)


def lift(d: pd.DataFrame) -> pd.DataFrame:
    from pyliftover import LiftOver
    lo = LiftOver(str(CHAIN))
    rows, fail, multi, xchrom = [], 0, 0, 0
    for r in d.itertuples():
        res = lo.convert_coordinate(r.chrom38, int(r.pos38))
        if not res:
            fail += 1
            continue
        if len(res) > 1:
            multi += 1
            continue
        c19, p19 = res[0][0], int(res[0][1])
        if c19 != r.chrom38:
            xchrom += 1
            continue
        rows.append(dict(sample=r.Sample, patient=r.patient, chrom=c19, pos0=p19,
                         chrom38=r.chrom38, pos38=int(r.pos38), reads=r.Reads,
                         status=r.status))
    out = pd.DataFrame(rows)
    print(f"[lift]   {len(out)}/{len(d)} lifted uniquely and chromosome-consistently "
          f"(unmapped {fail}, multi-mapping {multi}, chromosome-switch {xchrom}) | "
          f"{out.status.value_counts().to_dict()}")
    return out


def prep(d: pd.DataFrame, nul, mask) -> pd.DataFrame:
    """Amplicon-peak exclusion + joint-strata lookup, identical to the Furuta/Sasa joint arms."""
    in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                        for c, p in zip(d["chrom"], d["pos0"])])
    s = d[~in_mask].copy()
    k = pd.MultiIndex.from_arrays([s["chrom"].astype(str).to_numpy(),
                                   (s["pos0"].to_numpy(np.int64) // BIN) * BIN])
    s["stratum"] = nul.lookup.reindex(k).to_numpy()
    s = s[s["stratum"].notna()].copy()
    s["stratum"] = s["stratum"].astype(int)
    s["within_bin_offset"] = s["pos0"] - (s["pos0"] // BIN) * BIN
    return s.sort_values(["chrom", "pos0"]).reset_index(drop=True)


def run(idx, nul, sub, tag, n_perm):
    labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                               sub["pos0"].to_numpy(np.int64), WINDOW)
    n_clu = int(len(np.unique(labels)))
    nul.n_calls = nul.n_fallback = 0
    rng = np.random.default_rng(SEED)
    o, nl, per_cluster, s_obs, s_null = F.est_block(idx, nul, rng, sub, WINDOW, n_perm,
                                                    labels, unit="uniqSV")
    fb = round(100 * nul.n_fallback / max(nul.n_calls, 1), 2)
    out = []
    for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                     ("f_clusterSignTest", F.summarize(s_obs, s_null))):
        out.append(dict(subset=tag, estimator=est, n_coords=len(sub), n_clusters=n_clu,
                        n_patients=int(sub["patient"].nunique()), pct_pool_fallback=fb, **res))
    print(f"  [{tag:38s}] n={len(sub):4d} clu={n_clu:4d} pat={sub['patient'].nunique():3d}  "
          f"bc={out[0]['fold']:.4f} [{out[0]['fold_ci_lo']:.4f},{out[0]['fold_ci_hi']:.4f}] "
          f"p={out[0]['emp_p']:.5g} | sign={out[1]['fold']:.4f} p={out[1]['emp_p']:.5g}",
          flush=True)
    return out, per_cluster


def contrast(pc_a: pd.DataFrame, pc_b: pd.DataFrame, tag_a: str, tag_b: str, rng) -> dict:
    """Cluster-robust bootstrap for log(fold_a / fold_b), plus the MDE of that ratio."""
    oa, ea = pc_a["obs"].to_numpy(float), pc_a["null_mean"].to_numpy(float)
    ob, eb = pc_b["obs"].to_numpy(float), pc_b["null_mean"].to_numpy(float)
    fa, fb_ = oa.sum() / ea.sum(), ob.sum() / eb.sum()
    obs_ratio = fa / fb_
    ka, kb = len(oa), len(ob)
    draws = np.empty(NBOOT)
    for i in range(NBOOT):
        ia, ib = rng.integers(0, ka, ka), rng.integers(0, kb, kb)
        sa, sb = ea[ia].sum(), eb[ib].sum()
        draws[i] = ((oa[ia].sum() / sa) / (ob[ib].sum() / sb)) if sa > 0 and sb > 0 else np.nan
    lg = np.log(draws[np.isfinite(draws) & (draws > 0)])
    lo, hi = np.exp(np.percentile(lg, [2.5, 97.5]))
    # two-sided bootstrap p for ratio = 1
    p = 2 * min((lg <= 0).mean(), (lg >= 0).mean())
    return dict(analysis="contrast", subset=f"{tag_a} / {tag_b}",
                fold_a=round(fa, 4), fold_b=round(fb_, 4), ratio=round(obs_ratio, 4),
                ratio_ci_lo=round(float(lo), 4), ratio_ci_hi=round(float(hi), 4),
                boot_p=round(float(min(p, 1.0)), 5),
                mde_ratio_80pct=round(float(np.exp((1.96 + 0.84) * lg.std(ddof=1))), 4),
                n_clusters_a=ka, n_clusters_b=kb)


def detectability_confounder(d: pd.DataFrame) -> list[dict]:
    """The named competing explanation for any Silent/Productive difference.

    'Productive' means a viral-human chimeric transcript was DETECTED. Detection needs copy number
    (more copies -> more transcript) and local transcriptional activity -- and amplified,
    rearranged loci supply both. So a Productive arm with more SV nearby can arise with no causal
    role for transcription at all, purely because rearranged loci make the chimera visible.

    Tested on the authors' own gene-level CNV matrix and their own read support, within patient.
    """
    import gzip
    from scipy import stats

    cnv_p = DATA.parent / "cnv_log2.csv"
    gtf = DATA.parent / "gencode.v33.primary_assembly.annotation.gtf.gz"
    if not (cnv_p.exists() and gtf.exists()):
        print("[cnv]    fan2023 CNV matrix or GENCODE v33 absent -- confounder check skipped")
        return []

    genes: dict[str, list] = {}
    with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9 or f[2] != "gene" or 'gene_type "protein_coding"' not in f[8]:
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if m:
                genes.setdefault(f[0], []).append((int(f[3]), int(f[4]), m.group(1)))
    gx = {c: (np.array([g[0] for g in v]), np.array([g[1] for g in v]),
              np.array([g[2] for g in v], dtype=object)) for c, v in genes.items()}
    print(f"[cnv]    GENCODE v33: {sum(len(v) for v in genes.values())} protein-coding genes")

    cnv = pd.read_csv(cnv_p).set_index("gene")
    vals = []
    for r in d.itertuples():
        g = gx.get(str(r.chrom38))
        v = np.nan
        if g is not None:
            st, en, nm = g
            dist = np.maximum(np.maximum(st - r.pos38, r.pos38 - en), 0)
            j = int(np.argmin(dist))
            if dist[j] < 500_000:
                gn = str(nm[j])
                if gn in cnv.index and r.sample in cnv.columns:
                    v = float(np.atleast_1d(cnv.loc[gn, r.sample])[0])
        vals.append(v)
    d = d.assign(local_log2cnv=vals)
    ok = d.dropna(subset=["local_log2cnv"])
    print(f"[cnv]    local gene-level log2 CNV attached to {len(ok)}/{len(d)} coordinates")

    rows = []
    a_ = ok.loc[ok.status == "Productive", "local_log2cnv"].to_numpy(float)
    b_ = ok.loc[ok.status == "Silent", "local_log2cnv"].to_numpy(float)
    rows.append(dict(analysis="confounder_cnv", subset="pooled", n_productive=len(a_),
                     n_silent=len(b_), median_productive=round(float(np.median(a_)), 4),
                     median_silent=round(float(np.median(b_)), 4),
                     mwu_p=round(float(stats.mannwhitneyu(a_, b_).pvalue), 6)))

    # Restrict to patients with both transcript categories.
    per = (ok.groupby(["patient", "status"]).local_log2cnv.median().unstack().dropna())
    if len(per) >= 6:
        w = stats.wilcoxon(per["Productive"], per["Silent"])
        diff = (per["Productive"] - per["Silent"]).to_numpy(float)
        rows.append(dict(analysis="confounder_cnv", subset="within_patient",
                         n_patients=len(per), mean_diff=round(float(diff.mean()), 4),
                         median_diff=round(float(np.median(diff)), 4),
                         frac_patients_productive_higher=round(float((diff > 0).mean()), 4),
                         wilcoxon_p=round(float(w.pvalue), 6),
                         mde_80pct=round(float((1.96 + 0.84) * diff.std(ddof=1) / np.sqrt(len(diff))), 4)))
        print(f"[cnv]    within-patient productive-minus-silent local log2 CNV "
              f"{diff.mean():+.4f} (Wilcoxon P={w.pvalue:.4g}, {len(per)} patients, "
              f"{100*(diff>0).mean():.0f}% higher)")

    # Assess detectability using the published read-support column.
    rd = d.assign(reads=pd.to_numeric(d.reads, errors="coerce")).dropna(subset=["reads"])
    perr = rd.groupby(["patient", "status"]).reads.median().unstack().dropna()
    if len(perr) >= 6:
        w = stats.wilcoxon(perr["Productive"], perr["Silent"])
        dr = np.log2(perr["Productive"].to_numpy(float) / perr["Silent"].to_numpy(float))
        rows.append(dict(analysis="confounder_reads", subset="within_patient",
                         n_patients=len(perr),
                         median_log2_read_ratio=round(float(np.median(dr)), 4),
                         frac_patients_productive_higher=round(float((dr > 0).mean()), 4),
                         wilcoxon_p=round(float(w.pvalue), 6)))
        print(f"[reads]  within-patient log2(productive/silent) read support "
              f"{np.median(dr):+.4f} (Wilcoxon P={w.pvalue:.4g}, {len(perr)} patients)")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    d = lift(load_fan())
    both = (d.groupby("patient").status.nunique() == 2)
    both_pat = set(both[both].index)
    print(f"[design] {len(both_pat)} patients carry BOTH classes "
          f"({int(d.patient.isin(both_pat).sum())} coordinates) -- the balanced arm")

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    idx, _ = F.load_events()
    nul = J.JointNuller(feat, F.Nuller(), mask)

    print("\n[joint-corrected estimator, +/-100 kb, vs the same PCAWG pool as Tian/Sasa/Furuta]")
    rows, pcs = [], {}
    arms = [("ALL Fan coordinates", d),
            ("SILENT", d[d.status == "Silent"]),
            ("PRODUCTIVE", d[d.status == "Productive"]),
            ("SILENT (both-class patients)", d[(d.status == "Silent") & d.patient.isin(both_pat)]),
            ("PRODUCTIVE (both-class patients)",
             d[(d.status == "Productive") & d.patient.isin(both_pat)])]
    for tag, sel in arms:
        sub = prep(sel, nul, mask)
        if len(sub) < 30:
            print(f"  [{tag:38s}] SKIPPED, only {len(sub)} coordinates after joint lookup")
            continue
        r, pc = run(idx, nul, sub, tag, a.n_perm)
        rows += r
        pcs[tag] = pc

    out = pd.DataFrame(rows)
    out.insert(0, "analysis", "arm")

    con = []
    for x, y in (("PRODUCTIVE", "SILENT"),
                 ("PRODUCTIVE (both-class patients)", "SILENT (both-class patients)")):
        if x in pcs and y in pcs:
            c = contrast(pcs[x], pcs[y], x, y, rng)
            con.append(c)
            print(f"\n[contrast] {x} / {y}: ratio {c['ratio']:.4f} "
                  f"[{c['ratio_ci_lo']:.4f}, {c['ratio_ci_hi']:.4f}] boot p={c['boot_p']:.4f} "
                  f"| MDE {c['mde_ratio_80pct']:.4f}")
    conf = detectability_confounder(d)
    parts = [out] + ([pd.DataFrame(con)] if con else []) + ([pd.DataFrame(conf)] if conf else [])
    res = pd.concat(parts, ignore_index=True)
    res.to_csv(OUT / "fan2023_silent_productive.csv", index=False)
    d.to_csv(OUT / "fan2023_silent_productive_coords.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print("\n=== arms ===")
    print(out[out.estimator == "bc_uniqSV_blockNull"][
        ["subset", "n_coords", "n_clusters", "n_patients", "observed", "expected",
         "fold", "fold_ci_lo", "fold_ci_hi", "emp_p"]].to_string(index=False))
    print("\n=== sign-test estimator ===")
    print(out[out.estimator == "f_clusterSignTest"][
        ["subset", "n_clusters", "fold", "fold_ci_lo", "fold_ci_hi", "emp_p"]].to_string(index=False))
    if con:
        print("\n=== contrasts ===")
        print(pd.DataFrame(con).drop(columns=["analysis"]).to_string(index=False))
    print("\nreference under the SAME joint correction:")
    print("  Tian HPV (primary)   bc 1.0836 [1.0521,1.1166] p<5e-5 | sign 1.0791 p=3.5e-4")
    print("  Sasa2025 HPV16       bc 1.3168 [1.1385,1.5371] p=2e-4 | sign 1.2922 p=0.0103")
    print("  Furuta NT (HBV)      bc 1.0213 [0.9709,1.0736] p=0.201")
    print(f"\nwrote {OUT/'fan2023_silent_productive.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
