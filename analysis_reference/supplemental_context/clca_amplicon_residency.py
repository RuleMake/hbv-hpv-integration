"""Measures same-sample integration residency in CLCA amplicon intervals.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F                    # noqa: E402
import joint_correction as J                   # noqa: E402
import amplification_and_histology_control as A    # noqa: E402
from probe_clca_amplicon import parse_intervals, CLCA   # noqa: E402

OUT = Path(r"source_project\results") / "clca_amplicon_residency.csv"
BIN = F.BIN_SIZE
WINDOW = 100_000
HOTSPOT = ("TERT", "CCNE1", "KMT2B", "CLPTM1L")


def merged_amplicons() -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    """(sample, chrom) -> (sorted merged starts, ends). Overlapping calls are merged so that
    membership is one searchsorted rather than a scan, and so a locus covered by two overlapping
    amplicon calls is not counted twice."""
    amp = pd.read_csv(CLCA / "amplicons.csv", low_memory=False)
    rows = []
    for _, r in amp.iterrows():
        for c, lo, hi in parse_intervals(r["intervals"]):
            rows.append((r["sample"], c, lo, hi))
    seg = pd.DataFrame(rows, columns=["sample", "chrom", "start", "end"])
    out = {}
    for (s, c), g in seg.groupby(["sample", "chrom"]):
        iv = np.sort(g[["start", "end"]].to_numpy(np.int64), axis=0)
        order = np.argsort(g["start"].to_numpy(np.int64))
        st = g["start"].to_numpy(np.int64)[order]
        en = g["end"].to_numpy(np.int64)[order]
        ms, me = [st[0]], [en[0]]
        for a, b in zip(st[1:], en[1:]):
            if a <= me[-1]:
                me[-1] = max(me[-1], b)
            else:
                ms.append(a); me.append(b)
        out[(s, c)] = (np.array(ms, np.int64), np.array(me, np.int64))
    print(f"[amp]   {len(seg):,} intervals -> {sum(len(v[0]) for v in out.values()):,} merged, "
          f"{len({k[0] for k in out})} samples")
    return out


def inside(amps, sample: np.ndarray, chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Boolean: does each position fall in a merged amplicon of its own sample?"""
    hit = np.zeros(len(pos), dtype=bool)
    key = pd.Series(list(zip(sample, chrom)))
    for k, idx in key.groupby(key).groups.items():
        iv = amps.get(k)
        if iv is None:
            continue
        i = np.asarray(idx, dtype=np.int64)
        st, en = iv
        j = np.searchsorted(st, pos[i], side="right") - 1
        ok = j >= 0
        hit[i[ok]] = pos[i][ok] <= en[j[ok]]
    return hit


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    amps = merged_amplicons()
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    nul = J.JointNuller(feat, F.Nuller(), mask)

    d = pd.read_csv(CLCA / "integrations.csv", low_memory=False)
    d["chrom"] = "chr" + d["chrom"].astype(str).str.replace("chr", "", regex=False)
    d["pos0"] = d["pos"].astype(np.int64)
    d = d[d["sample"].isin({k[0] for k in amps})].reset_index(drop=True)
    key = pd.MultiIndex.from_arrays([d["chrom"].to_numpy(str),
                                     (d["pos0"].to_numpy(np.int64) // BIN) * BIN])
    d = d.assign(stratum=nul.lookup.reindex(key).to_numpy())
    d = d[d["stratum"].notna()].reset_index(drop=True)
    d["stratum"] = d["stratum"].astype(int)
    d["within_bin_offset"] = d["pos0"] - (d["pos0"] // BIN) * BIN
    print(f"[ints]  {len(d):,} integrations in amplicon-bearing samples, "
          f"{d['sample'].nunique()} samples, after joint-stratum lookup")

    d["support"] = d["split_reads"].fillna(0) + d["pair_reads"].fillna(0)
    rows = []
    # `support_ge_50` IS NOT A VALID DETECTABILITY CONTROL. It is retained, and reported, only so
    # that the reasoning is on record.
    #
    # Amplicon-resident integrations carry about twice the read support of non-resident ones
    # (median 70 vs 36; within-sample paired 78 vs 37, Wilcoxon P = 0.0013, 17/23 samples higher
    # inside). That is predicted both by the hypothesis -- an integrant on an amplified allele is
    # physically present in more copies -- and by its rival, that high copy number makes an
    # integration easier to detect at all. The obvious move is to restrict to well-supported
    # integrations, where sensitivity should not be limiting, and see whether the enrichment erodes.
    #
    # It does not erode; it roughly doubles. But that is what the RIVAL hypothesis predicts too,
    # because read support is a COLLIDER here:
    #
    #       read support  <--  copy number  -->  inside an amplicon
    #
    # conditioning on support selects for high copy number, and high copy number is what defines an
    # amplicon. The arm is inflated by the very bias it was meant to exclude, so its larger fold is
    # not evidence for anything. Detection bias cannot be excluded with the data available: CLCA
    # gives copy number only inside amplicon calls, so there is no equally-amplified non-amplicon
    # comparison region to match against.
    for arm in ("all", "no_hotspot", "support_ge_50"):
        if arm == "all":
            s = d
        elif arm == "no_hotspot":
            s = d[~d["gene"].astype(str).isin(HOTSPOT)].reset_index(drop=True)
        else:
            s = d[d["support"] >= 50].reset_index(drop=True)
        # leader clusters WITHIN sample
        lab = np.empty(len(s), np.int64)
        base = 0
        for smp, g in s.groupby("sample", sort=True):
            gi = g.index.to_numpy()
            l = F.leader_clusters(g["chrom"].to_numpy(str), g["pos0"].to_numpy(np.int64), WINDOW)
            lab[gi] = l + base
            base += l.max() + 1
        uniq, inv = np.unique(lab, return_inverse=True)
        n_clu = len(uniq)
        lead = np.zeros(n_clu, np.int64)
        lead[inv[::-1]] = np.arange(len(inv))[::-1]
        delta = s["pos0"].to_numpy(np.int64) - s["pos0"].to_numpy(np.int64)[lead][inv]

        smp = s["sample"].to_numpy(str)
        chm = s["chrom"].to_numpy(str)
        hit = inside(amps, smp, chm, s["pos0"].to_numpy(np.int64))
        # cluster is a hit if any member is inside; sample is a hit if any cluster is
        clu_hit = np.zeros(n_clu, bool)
        np.logical_or.at(clu_hit, inv, hit)
        obs_clu = int(clu_hit.sum())
        obs_coord = int(hit.sum())
        obs_samp = int(pd.Series(hit).groupby(smp).any().sum())

        pools = [nul.pool(str(chm[lead[k]]), int(s["stratum"].to_numpy()[lead[k]]),
                          int(s["within_bin_offset"].to_numpy()[lead[k]]),
                          int(s["pos0"].to_numpy()[lead[k]])) for k in range(n_clu)]
        sizes = np.array([p.size for p in pools], np.int64)
        starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
        flat = np.concatenate(pools).astype(np.int64)
        lead_chrom = chm[lead][inv]

        rng = np.random.default_rng(a.seed)
        nc = np.empty(a.n_perm); ncoord = np.empty(a.n_perm); nsamp = np.empty(a.n_perm)
        for t in range(a.n_perm):
            idx = starts + (rng.random(n_clu) * sizes).astype(np.int64)
            p = np.maximum(flat[idx][inv] + delta, 0)
            h = inside(amps, smp, lead_chrom, p)
            ch = np.zeros(n_clu, bool)
            np.logical_or.at(ch, inv, h)
            nc[t] = ch.sum(); ncoord[t] = h.sum()
            nsamp[t] = pd.Series(h).groupby(smp).any().sum()
            if (t + 1) % 5000 == 0:
                print(f"        perm {t+1}/{a.n_perm}", flush=True)

        for name, o, nl in (("clusters_in_amplicon", obs_clu, nc),
                            ("coords_in_amplicon", obs_coord, ncoord),
                            ("samples_with_any", obs_samp, nsamp)):
            e, sd = float(nl.mean()), float(nl.std(ddof=1))
            p_hi = (int(np.sum(nl >= o)) + 1) / (a.n_perm + 1)
            lo, hi = np.percentile(nl, [2.5, 97.5])
            rows.append(dict(
                arm=arm, statistic=name, n_coords=len(s), n_clusters=n_clu,
                n_samples=int(s["sample"].nunique()), n_perm=a.n_perm,
                observed=o, expected=round(e, 2),
                fold=round(o / e, 3) if e > 0 else np.nan,
                fold_ci_lo=round(o / hi, 3) if hi > 0 else np.nan,
                fold_ci_hi=round(o / lo, 3) if lo > 0 else np.nan,
                z=round((o - e) / sd, 2) if sd > 0 else np.nan,
                emp_p=p_hi,
                mde_fold=round((e + 2.80 * sd) / e, 3) if e > 0 else np.nan,
                clears_mde=bool(o >= e + 2.80 * sd)))
        print(f"  [{arm}] clusters {obs_clu} vs {nc.mean():.1f} "
              f"(fold {obs_clu/nc.mean():.2f}), coords {obs_coord} vs {ncoord.mean():.1f}",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nwrote {OUT}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
