"""Tests repeat-content adjustment of the local class-composition axis.
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
import locus_class_propensity as E              # noqa: E402

RMSK = E.PROJ / "data" / "external_sv" / "rmsk_hg19.txt.gz"
CACHE = E.PROJ / "results" / "cache" / "rmsk_bin_coverage.parquet"
OUT = E.PROJ / "results" / "repeat_content_axis.csv"

REPEATS = ["SINE", "LINE", "LTR", "DNA", "Satellite", "Simple_repeat"]
VIRUSES = ("HBV", "HPV")
PREDICTED = {"SINE": "HPV > HBV (NAHR produces the DEL/DUP side HPV sits on)"}


def build_coverage() -> pd.DataFrame:
    """Base pairs of each repeat class per 50-kb bin. Cached -- the parse is the slow part."""
    if CACHE.exists():
        d = pd.read_parquet(CACHE)
        print(f"[rmsk]   cache hit: {len(d):,} bins x {len(REPEATS)} repeat classes")
        return d
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    print(f"[rmsk]   parsing {RMSK.name} (5.5M records) ...", flush=True)
    d = pd.read_csv(RMSK, sep="\t", header=None, usecols=[5, 6, 7, 11],
                    names=["chrom", "start", "end", "repclass"], compression="gzip",
                    dtype={5: str, 6: np.int64, 7: np.int64, 11: str})
    d = d[d["repclass"].isin(REPEATS)]
    d = d[~d["chrom"].str.contains("_")]
    # a repeat is short relative to a 50-kb bin, so charge its full length to the bin holding its
    # midpoint; the spill for the rare multi-bin LINE/Satellite is immaterial to a coverage tier
    mid = ((d["start"] + d["end"]) // 2).to_numpy(np.int64)
    d = d.assign(bin_start=(mid // E.BIN) * E.BIN, bp=(d["end"] - d["start"]).to_numpy(np.int64))
    g = (d.groupby(["chrom", "bin_start", "repclass"])["bp"].sum()
         .unstack("repclass").fillna(0.0).reset_index())
    for r in REPEATS:
        if r not in g.columns:
            g[r] = 0.0
    g = g[["chrom", "bin_start"] + REPEATS]
    g.to_parquet(CACHE, index=False)
    print(f"[rmsk]   {len(g):,} bins with repeat coverage -> {CACHE.name}")
    return g


def cumulative(cov: pd.DataFrame, feat: pd.DataFrame):
    """Per-chromosome (n_bins+1, n_feat) cumulative coverage, aligned to the null's own bin grid."""
    grid = feat[["chrom", "bin_start"]].drop_duplicates()
    m = grid.merge(cov, on=["chrom", "bin_start"], how="left").fillna(0.0)
    m = m.sort_values(["chrom", "bin_start"])
    cum, nbin = {}, {}
    for c, g in m.groupby("chrom", sort=False):
        v = g[REPEATS].to_numpy(np.float64)
        cum[str(c)] = np.vstack([np.zeros((1, len(REPEATS))), np.cumsum(v, axis=0)])
        nbin[str(c)] = len(v)
    return cum, nbin


def window_cov(cum, nbin, half, chrom, pos):
    """Mean repeat bp per Mb in +/-half bins around each position."""
    out = np.zeros((len(pos), len(REPEATS)))
    b = pos // E.BIN
    for c in pd.unique(chrom):
        m = chrom == c
        if c not in cum:
            continue
        n = nbin[c]
        lo = np.clip(b[m] - half, 0, n)
        hi = np.clip(b[m] + half + 1, 0, n)
        out[m] = (cum[c][hi] - cum[c][lo]) / np.maximum((hi - lo)[:, None] * E.BIN / 1e6, 1e-9)
    return out


def score(sub, labels, nuller, cum, nbin, half, rng, n_perm, chunk=250):
    """Estimate observed and null cluster-mean repeat coverage using rigid displacement."""
    chrom = sub["chrom"].to_numpy(str)
    pos = sub["pos0"].to_numpy(np.int64)
    strat = sub["stratum"].to_numpy(np.int64)
    uniq, inv = np.unique(labels, return_inverse=True)
    n_clu = len(uniq)
    size = np.bincount(inv, minlength=n_clu)
    lead_i = np.zeros(n_clu, dtype=np.int64)
    lead_i[inv[::-1]] = np.arange(len(inv))[::-1]
    delta = pos - pos[lead_i][inv]
    w = 1.0 / (n_clu * size[inv])

    obs = (window_cov(cum, nbin, half, chrom, pos) * w[:, None]).sum(axis=0)

    pools = [nuller.pool(str(chrom[lead_i[k]]), int(strat[lead_i[k]]),
                         int(pos[lead_i[k]]) % E.BIN, int(pos[lead_i[k]]))
             for k in range(n_clu)]
    sizes = np.array([p.size for p in pools], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(sizes)[:-1]])
    flat = np.concatenate(pools).astype(np.int64)
    lead_chrom = chrom[lead_i][inv]

    null = np.empty((n_perm, len(REPEATS)))
    done = 0
    while done < n_perm:
        c = min(chunk, n_perm - done)
        idx = starts[None, :] + (rng.random((c, n_clu)) * sizes[None, :]).astype(np.int64)
        newpos = flat[idx][:, inv] + delta[None, :]
        np.maximum(newpos, 0, out=newpos)
        cc = window_cov(cum, nbin, half, np.tile(lead_chrom, c), newpos.ravel())
        cc = cc.reshape(c, len(pos), len(REPEATS))
        null[done:done + c] = (cc * w[None, :, None]).sum(axis=1)
        done += c
        if done % 5000 < chunk:
            print(f"          perm {done}/{n_perm}", flush=True)
    return obs, null


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--windows", type=int, nargs="*", default=[200_000, 1_000_000, 2_000_000])
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    cov = build_coverage()
    cum_sv, nb_sv, p_global = E.propensity_map()
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), J.WINDOW)
    co = F.load_coordinates()
    joint = J.JointNuller(feat, F.Nuller(), mask)
    cum, nbin = cumulative(cov, feat)

    rows = []
    for win in a.windows:
        half = win // E.BIN
        dens = E.DensityNuller(joint, cum_sv, nb_sv, half)
        draws = {}
        for vi, vc in enumerate(VIRUSES):
            sub0, sub, tier, labels = E.prepare(co, vc, joint, dens, mask)
            rng = np.random.default_rng([a.seed, vi])
            obs, null = score(sub, labels, joint, cum, nbin, half, rng, a.n_perm)
            draws[vc] = (obs, null)
            for i, r in enumerate(REPEATS):
                rows.append(dict(window_bp=win, virus=vc, repeat_class=r,
                                 n_clusters=int(len(np.unique(labels))), n_perm=a.n_perm,
                                 **{k: v for k, v in
                                    E.summarize(obs[i], null[:, i], r).items()
                                    if k != "statistic"}))
            print(f"  [{vc} {win//1000}kb] " + "  ".join(
                f"{r} {obs[i]/null[:, i].mean():.4f}" for i, r in enumerate(REPEATS)), flush=True)

        (oh, nh), (op, npv) = draws["HBV"], draws["HPV"]
        for i, r in enumerate(REPEATS):
            eh = oh[i] - nh[:, i].mean()
            ep = op[i] - npv[:, i].mean()
            inter = ep - eh                      # HPV minus HBV, the predicted direction for SINE
            inull = (npv[:, i] - npv[:, i].mean()) - (nh[:, i] - nh[:, i].mean())
            sd = float(inull.std(ddof=1))
            p = (int(np.sum(np.abs(inull) >= abs(inter))) + 1) / (len(inull) + 1)
            rows.append(dict(window_bp=win, virus="HPV-minus-HBV", repeat_class=r,
                             n_perm=a.n_perm, excess_HBV=round(float(eh), 3),
                             excess_HPV=round(float(ep), 3),
                             observed=round(float(inter), 3), null_sd=round(sd, 3),
                             z=round(inter / sd, 2) if sd else np.nan,
                             emp_p_two_sided=round(p, 5), mde_excess=round(2.80 * sd, 3),
                             clears_mde=bool(abs(inter) >= 2.80 * sd),
                             prediction=PREDICTED.get(r, "")))
            flag = "  <-- PREDICTED" if r in PREDICTED else ""
            print(f"  [HPV-HBV {win//1000}kb] {r:14s} {inter:+9.2f} bp/Mb  z {inter/sd:+5.2f}  "
                  f"P {p:.4f}  MDE {2.80*sd:.2f}  "
                  f"{'clears' if abs(inter) >= 2.80*sd else 'FAILS'}{flag}", flush=True)
        print(flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT, index=False)
    pd.set_option("display.width", 250)
    print(d[d.virus == "HPV-minus-HBV"].to_string(index=False))
    print(f"\nwrote {OUT}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
