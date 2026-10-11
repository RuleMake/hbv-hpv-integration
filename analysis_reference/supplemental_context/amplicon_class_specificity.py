"""Compares integration residency across CLCA amplicon classes.
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
import clca_amplicon_residency as M             # noqa: E402
from probe_clca_amplicon import parse_intervals, CLCA   # noqa: E402

OUT = Path(r"source_project\results") / "amplicon_class_specificity.csv"
CLASSES = ["BFB", "Circular", "Heavily rearranged", "Linear"]


def merged_by_class() -> dict[str, dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]]:
    """Merge amplicon intervals within each class, keyed by sample and chromosome."""
    amp = pd.read_csv(CLCA / "amplicons.csv", low_memory=False)
    amp["cls"] = amp["amp_class"].astype(str).str.replace(
        "Heavily rearranged rearranged", "Heavily rearranged", regex=False)
    out = {}
    for cls in CLASSES:
        rows = []
        for _, r in amp[amp["cls"] == cls].iterrows():
            for c, lo, hi in parse_intervals(r["intervals"]):
                rows.append((r["sample"], c, lo, hi))
        seg = pd.DataFrame(rows, columns=["sample", "chrom", "start", "end"])
        d = {}
        for (s, c), g in seg.groupby(["sample", "chrom"]):
            order = np.argsort(g["start"].to_numpy(np.int64))
            st = g["start"].to_numpy(np.int64)[order]
            en = g["end"].to_numpy(np.int64)[order]
            ms, me = [st[0]], [en[0]]
            for x, y in zip(st[1:], en[1:]):
                if x <= me[-1]:
                    me[-1] = max(me[-1], y)
                else:
                    ms.append(x); me.append(y)
            d[(s, c)] = (np.array(ms, np.int64), np.array(me, np.int64))
        out[cls] = d
        print(f"[amp]   {cls:20s} {len(seg):5,} intervals -> "
              f"{sum(len(v[0]) for v in d.values()):5,} merged, "
              f"{len({k[0] for k in d})} samples")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    byclass = merged_by_class()
    allsamples = {k[0] for d in byclass.values() for k in d}
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), M.WINDOW)
    nul = J.JointNuller(feat, F.Nuller(), mask)

    d = pd.read_csv(CLCA / "integrations.csv", low_memory=False)
    d["chrom"] = "chr" + d["chrom"].astype(str).str.replace("chr", "", regex=False)
    d["pos0"] = d["pos"].astype(np.int64)
    d = d[d["sample"].isin(allsamples)].reset_index(drop=True)
    key = pd.MultiIndex.from_arrays([d["chrom"].to_numpy(str),
                                     (d["pos0"].to_numpy(np.int64) // M.BIN) * M.BIN])
    d = d.assign(stratum=nul.lookup.reindex(key).to_numpy())
    d = d[d["stratum"].notna()].reset_index(drop=True)
    d["stratum"] = d["stratum"].astype(int)
    d["within_bin_offset"] = d["pos0"] - (d["pos0"] // M.BIN) * M.BIN
    print(f"[ints]  {len(d):,} integrations, {d['sample'].nunique()} samples")

    lab = np.empty(len(d), np.int64)
    base = 0
    for _s, g in d.groupby("sample", sort=True):
        l = F.leader_clusters(g["chrom"].to_numpy(str), g["pos0"].to_numpy(np.int64), M.WINDOW)
        lab[g.index.to_numpy()] = l + base
        base += l.max() + 1
    uniq, inv = np.unique(lab, return_inverse=True)
    n_clu = len(uniq)
    lead = np.zeros(n_clu, np.int64)
    lead[inv[::-1]] = np.arange(len(inv))[::-1]
    pos = d["pos0"].to_numpy(np.int64)
    delta = pos - pos[lead][inv]
    smp = d["sample"].to_numpy(str)
    chm = d["chrom"].to_numpy(str)

    obs = {}
    for c in CLASSES:
        h = M.inside(byclass[c], smp, chm, pos)
        ch = np.zeros(n_clu, bool)
        np.logical_or.at(ch, inv, h)
        obs[c] = int(ch.sum())
    print(f"[obs]   {n_clu} clusters; hits by class " + str(obs))

    pools = [nul.pool(str(chm[lead[k]]), int(d["stratum"].to_numpy()[lead[k]]),
                      int(d["within_bin_offset"].to_numpy()[lead[k]]), int(pos[lead[k]]))
             for k in range(n_clu)]
    sizes = np.array([p.size for p in pools], np.int64)
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    flat = np.concatenate(pools).astype(np.int64)
    lead_chrom = chm[lead][inv]

    rng = np.random.default_rng(a.seed)
    null = {c: np.empty(a.n_perm) for c in CLASSES}
    chi = np.empty(a.n_perm)
    obs_vec = np.array([obs[c] for c in CLASSES], float)
    for t in range(a.n_perm):
        idx = starts + (rng.random(n_clu) * sizes).astype(np.int64)
        p = np.maximum(flat[idx][inv] + delta, 0)
        v = []
        for c in CLASSES:
            h = M.inside(byclass[c], smp, lead_chrom, p)
            ch = np.zeros(n_clu, bool)
            np.logical_or.at(ch, inv, h)
            null[c][t] = ch.sum()
            v.append(ch.sum())
        if (t + 1) % 5000 == 0:
            print(f"        perm {t+1}/{a.n_perm}", flush=True)
    # omnibus: chi-square of the class vector against the null's mean class vector
    exp_vec = np.array([null[c].mean() for c in CLASSES])
    safe = np.maximum(exp_vec, 0.5)
    chi_obs = float(np.sum((obs_vec - exp_vec) ** 2 / safe))
    for t in range(a.n_perm):
        v = np.array([null[c][t] for c in CLASSES])
        chi[t] = np.sum((v - exp_vec) ** 2 / safe)
    p_omni = (int(np.sum(chi >= chi_obs)) + 1) / (a.n_perm + 1)

    rows = []
    for c in CLASSES:
        nl = null[c]
        e, sd = float(nl.mean()), float(nl.std(ddof=1))
        lo, hi = np.percentile(nl, [2.5, 97.5])
        rows.append(dict(
            amp_class=c, n_clusters=n_clu, n_perm=a.n_perm,
            observed=obs[c], expected=round(e, 2),
            fold=round(obs[c] / e, 3) if e > 0 else np.nan,
            fold_ci_lo=round(obs[c] / hi, 3) if hi > 0 else np.nan,
            fold_ci_hi=round(obs[c] / lo, 3) if lo > 0 else np.nan,
            z=round((obs[c] - e) / sd, 2) if sd > 0 else np.nan,
            emp_p_two_sided=min(1.0, 2 * min(
                (int(np.sum(nl >= obs[c])) + 1) / (a.n_perm + 1),
                (int(np.sum(nl <= obs[c])) + 1) / (a.n_perm + 1))),
            mde_fold=round((e + 2.80 * sd) / e, 3) if e > 0 else np.nan,
            clears_mde=bool(abs(obs[c] - e) >= 2.80 * sd)))
    rows.append(dict(amp_class="OMNIBUS", n_clusters=n_clu, n_perm=a.n_perm,
                     observed=round(chi_obs, 3), expected=round(float(chi.mean()), 3),
                     fold=np.nan, fold_ci_lo=np.nan, fold_ci_hi=np.nan,
                     z=round((chi_obs - chi.mean()) / chi.std(ddof=1), 2),
                     emp_p_two_sided=p_omni, mde_fold=np.nan, clears_mde=bool(p_omni < 0.05)))

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nwrote {OUT}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
