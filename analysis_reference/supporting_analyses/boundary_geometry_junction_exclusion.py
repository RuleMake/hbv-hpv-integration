"""Compare duplication-boundary and interior occupancy after junction exclusion.

For a coordinate p, depth(m) counts intervals containing p with at least m base
pairs on each side. With exclusion radius Jx, the boundary count is
  edge(Jx) = depth(Jx) - depth(100 kb),
and the interior count is depth(100 kb). Their signed contrast is
  A(Jx) = (edge(Jx) - interior) / (edge(Jx) + interior).

The same margin is applied to observed and randomized positions. Jx = 0 checks
the reference statistic; larger radii assess sensitivity to events whose endpoints
lie close to the integration junction. This statistic is orientation-independent
and does not require the OK-seq or gene-orientation tracks."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F  # noqa: E402
import joint_correction as J  # noqa: E402
import amplification_and_histology_control as A  # noqa: E402

OUT = J.OUT
BIN = F.BIN_SIZE
RES = 1_000
EDGE_BP = 100_000                       # margin defining the interval interior
EXCL = (0, 1_000, 10_000, 50_000)       # junction-exclusion radii; 0 = published statistic
MARGINS = tuple(sorted(set(EXCL) | {EDGE_BP}))
WINDOW = 100_000
SEED = 7


def _sweep(left: np.ndarray, right: np.ndarray, n: int) -> np.ndarray:
    d = np.zeros(n + 1, dtype=np.int32)
    lo = np.clip(left // RES, 0, n - 1)
    hi = np.clip(right // RES, 0, n - 1)
    np.add.at(d, lo, 1)
    np.add.at(d, hi + 1, -1)
    return np.cumsum(d)[:n].astype(np.int32)


class ChromDepth:
    """Containment depth at every margin in MARGINS, on one chromosome."""

    __slots__ = ("n", "depth")

    def __init__(self, left: np.ndarray, right: np.ndarray, size: int):
        n = size // RES + 2
        self.n = n
        self.depth = {}
        for m in MARGINS:
            if m == 0:
                self.depth[m] = _sweep(left, right, n)
            else:
                keep = (right - left) > 2 * m
                self.depth[m] = _sweep(left[keep] + m, right[keep] - m, n)


def build_geometry(nul) -> dict[str, ChromDepth]:
    ev = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv",
                     sep="\t").drop_duplicates("unique_td_id")
    sizes = {c: int(max(nul.chrom_size.get(c, 0), int(g["right_pos0"].max()) + RES))
             for c, g in ev.groupby("chrom")}
    geo = {str(c): ChromDepth(g["left_pos0"].to_numpy(np.int64),
                              g["right_pos0"].to_numpy(np.int64), sizes[str(c)])
           for c, g in ev.groupby("chrom")}
    span = (ev["right_pos0"] - ev["left_pos0"]).to_numpy(np.int64)
    print(f"[geom]   {len(ev):,} canonical intervals on {len(geo)} chromosomes; median span "
          f"{np.median(span):,.0f} bp; margins {[m//1000 for m in MARGINS]} kb", flush=True)
    return geo


KEYS = tuple(f"depth_{m}" for m in MARGINS)


def stats_at(gm: ChromDepth, pos: np.ndarray) -> dict[str, np.ndarray]:
    q = pos // RES
    ok = ((q >= 0) & (q < gm.n)).astype(np.int64)
    q = np.clip(q, 0, gm.n - 1)
    return {f"depth_{m}": gm.depth[m][q] * ok for m in MARGINS}


def run_block(geo, nul, rng, sub, n_perm, labels):
    """Rigid-block displacement, identical rule to est_block / td_geometry_orientation."""
    df = sub.copy()
    df["_cluster"] = labels
    obs = {k: 0 for k in KEYS}
    null = {k: np.zeros(n_perm, dtype=np.int64) for k in KEYS}
    for _, g in df.groupby("_cluster", sort=True):
        c = str(g["chrom"].iloc[0])
        gm = geo.get(c)
        if gm is None:
            continue
        p = g["pos0"].to_numpy(np.int64)
        st = g["stratum"].to_numpy()
        off = g["within_bin_offset"].to_numpy()
        j = rng.integers(0, len(p), size=n_perm)
        delta = np.empty(n_perm, dtype=np.int64)
        for k in range(len(p)):
            m = j == k
            nk = int(m.sum())
            if nk:
                delta[m] = nul.draw(rng, c, st[k], off[k], p[k], nk) - p[k]
        for pi in p:
            o = stats_at(gm, np.array([pi], dtype=np.int64))
            nl = stats_at(gm, pi + delta)
            for k in KEYS:
                obs[k] += int(o[k][0])
                null[k] += nl[k]
    return obs, null


def asym(o_num, n_num, o_den, n_den) -> dict:
    """Calculate A = (num - den)/(num + den) relative to its null mean.
    
    The two-sided comparison centers on the simulated null mean, which need not be zero."""
    tot_o = o_num + o_den
    tot_n = n_num + n_den
    if tot_o == 0 or np.all(tot_n == 0):
        return dict(observed=np.nan, expected=np.nan, z=np.nan, emp_p=np.nan)
    a_o = (o_num - o_den) / tot_o
    a_n = (n_num - n_den) / np.maximum(tot_n, 1)
    mu, sd = float(a_n.mean()), float(a_n.std(ddof=1))
    p = (int(np.sum(np.abs(a_n - mu) >= abs(a_o - mu))) + 1) / (len(a_n) + 1)
    lo, hi = np.percentile(a_n, [2.5, 97.5])
    return dict(observed=round(float(a_o), 4), expected=round(mu, 4),
                null_ci_lo=round(float(lo), 4), null_ci_hi=round(float(hi), 4),
                z=round((a_o - mu) / max(sd, 1e-12), 2), emp_p=round(p, 5),
                n_edge=int(o_num), n_interior=int(o_den))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--viruses", default="HBV,HPV,EBV")
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), 100_000)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)
    geo = build_geometry(base)
    co = F.load_coordinates()

    rows = []
    for virus in a.viruses.split(","):
        sub = co[co["virus_class"] == virus].sort_values(
            ["chrom", "pos0"]).reset_index(drop=True)
        in_mask = np.array([(str(c), int(p) // BIN * BIN) in mask
                            for c, p in zip(sub["chrom"], sub["pos0"])])
        sub = sub[~in_mask].reset_index(drop=True)
        k = pd.MultiIndex.from_arrays([sub["chrom"].astype(str).to_numpy(),
                                       (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
        sub = sub.assign(stratum=nul.lookup.reindex(k).to_numpy())
        sub = sub[sub["stratum"].notna()].reset_index(drop=True)
        sub["stratum"] = sub["stratum"].astype(int)
        sub["within_bin_offset"] = sub["pos0"] - (sub["pos0"] // BIN) * BIN

        labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                   sub["pos0"].to_numpy(np.int64), WINDOW)
        n_clu = int(len(np.unique(labels)))
        nul.n_calls = nul.n_fallback = 0
        rng = np.random.default_rng(SEED)
        obs, null = run_block(geo, nul, rng, sub, a.n_perm, labels)
        fb = round(100 * nul.n_fallback / max(nul.n_calls, 1), 2)
        meta = dict(virus=virus, window_kb=WINDOW // 1000, n_coords=len(sub),
                    n_clusters=n_clu, pct_pool_fallback=fb, n_perm=a.n_perm)

        o_int, n_int = obs[f"depth_{EDGE_BP}"], null[f"depth_{EDGE_BP}"]
        for jx in EXCL:
            o_edge = obs[f"depth_{jx}"] - o_int
            n_edge = null[f"depth_{jx}"] - n_int
            r = {**meta, "excl_kb": jx // 1000,
                 "statistic": "edge_vs_interior",
                 "edge_fold": round(o_edge / max(float(n_edge.mean()), 1e-9), 4),
                 "interior_fold": round(o_int / max(float(n_int.mean()), 1e-9), 4),
                 **asym(o_edge, n_edge, o_int, n_int)}
            rows.append(r)
            print(f"  [{virus} excl {jx//1000:>2} kb] edge {o_edge:>7} (fold {r['edge_fold']:.4f}) "
                  f"vs interior {o_int:>7} (fold {r['interior_fold']:.4f}) | "
                  f"A {r['observed']:+.4f} null {r['expected']:+.4f} "
                  f"z {r['z']:+.2f} p {r['emp_p']:.5g}", flush=True)
        # containment folds per margin, for the record
        for m in MARGINS:
            rows.append({**meta, "excl_kb": m // 1000, "statistic": f"depth_margin_{m//1000}kb",
                         **F.summarize(obs[f"depth_{m}"], null[f"depth_{m}"])})

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "boundary_junction_excluded.csv", index=False)
    pd.set_option("display.width", 260)
    print("\n=== edge_vs_interior by junction-exclusion radius (0 kb = published statistic) ===")
    print(d[d["statistic"] == "edge_vs_interior"][
        ["virus", "excl_kb", "n_clusters", "n_edge", "n_interior", "edge_fold",
         "interior_fold", "observed", "expected", "z", "emp_p"]].to_string(index=False))
    print(f"\nwrote {OUT/'boundary_junction_excluded.csv'}  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
