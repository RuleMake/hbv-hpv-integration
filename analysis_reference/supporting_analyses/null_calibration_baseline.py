"""Calibrate counting and placement rules using simulated coordinate sets.

Each observed coordinate cluster is translated rigidly to a matched-stratum position,
preserving multiplicity, cluster size and internal spacing. The analysis measures
rejection at nominal alpha = 0.05 for coordinate-level counting, one coordinate per
cluster, and distinct-event counting with a cluster-block null.

The reference estimator is checked against observed-coordinate summaries before
simulation. Summary, replicate-level and reference-check outputs are written separately."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F                    # noqa: E402

OUT = Path(r"source_project\results")
BIN_SIZE = F.BIN_SIZE
ALPHAS = (0.05, 0.01, 0.001)
# Reference fold estimates: corrected_estimator_baseline.PUBLISHED_FOLD / estimator_reference_unclipped.csv
RUNG_TARGET = {"HBV": 1.1590, "HPV": 1.4094}
RUNG_TOL = 0.01


def bin_lookup() -> tuple[dict[tuple[str, int], int], set[tuple[str, int]]]:
    """(chrom, bin_start) -> stratum, for the bins available to the spatial null sampler,
    plus the set of (chrom, stratum) pools that hold a single bin.

    A single-bin pool is unusable: Nuller.pool builds bins + offset and then excludes the
    coordinate's own position, so one bin leaves nothing to draw from. The real catalogue happens
    to contain no such coordinate -- the rung below asserts it -- but a translated block can land
    in one, and a coordinate the estimator cannot score must not enter the synthetic set."""
    feat = pd.read_parquet(F.FEATURES)
    feat = feat.assign(chrom=feat["chrom"].astype(str).map(
        lambda c: c if c.startswith("chr") else f"chr{c}"))
    feat = feat[(feat["stratum"] >= 0) & (feat["chrom"] != "chrY")]
    look = {(str(c), int(b)): int(s)
            for c, b, s in zip(feat["chrom"], feat["bin_start"], feat["stratum"])}
    sizes = feat.groupby(["chrom", "stratum"]).size()
    thin = {(str(c), int(st)) for (c, st), n in sizes.items() if n < 2}
    print(f"[pools]   {len(sizes):,} chrom x stratum pools, {len(thin)} hold a single bin "
          f"and are not drawable")
    return look, thin


def synth_set(rng, nul, sub, labels, strat: dict, thin: set) -> pd.DataFrame:
    """One true-null coordinate set: every cluster translated as a rigid block.

    No clipping. `corrected_estimator_baseline.est_block` documents why an earlier clip made the observed
    and null sides score windows of different effective length; the same argument applies here, so
    a block that leaves the eligible genome simply loses the coordinates that land outside it.
    """
    out = []
    df = sub.assign(_cluster=labels)
    for _, g in df.groupby("_cluster", sort=True):
        c = str(g["chrom"].iloc[0])
        p = g["pos0"].to_numpy(np.int64)
        st = g["stratum"].to_numpy()
        off = g["within_bin_offset"].to_numpy()
        k = int(rng.integers(0, len(p)))
        delta = int(nul.draw(rng, c, st[k], off[k], p[k], 1)[0]) - int(p[k])
        for q in p + delta:
            b = int(q) // BIN_SIZE * BIN_SIZE
            s = strat.get((c, b))
            if s is None or q < 0 or (c, s) in thin:
                continue
            out.append((c, int(q), b, int(q) - b, s))
    return pd.DataFrame(out, columns=["chrom", "pos0", "bin_start", "within_bin_offset",
                                      "stratum"])


def rate(p: np.ndarray, a: float) -> float:
    return round(float(np.mean(p <= a)), 4)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-rep", type=int, default=150)
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--windows", default="100000,1000000")
    ap.add_argument("--viruses", default="HBV,HPV")
    a = ap.parse_args()

    idx, n_ev = F.load_events()
    nul = F.Nuller()
    co = F.load_coordinates()
    strat, thin_pools = bin_lookup()

    # ---------------------------------------------------------------- reproduction rung
    rung = []
    for vc in ("HBV", "HPV"):
        sub = co[co["virus_class"] == vc].reset_index(drop=True)
        rng = np.random.default_rng(a.seed)
        obs, null = F.est_independent(idx, nul, rng, sub, 100_000, 400)
        fold = obs / float(null.mean())
        ok = abs(fold - RUNG_TARGET[vc]) <= RUNG_TOL
        rung.append(dict(virus=vc, n_coords=len(sub), window_bp=100_000, fold=round(fold, 4),
                         published_fold=RUNG_TARGET[vc], tolerance=RUNG_TOL, reproduces=bool(ok)))
        print(f"[reference] {vc}: published estimator on the real coordinates -> fold {fold:.4f} "
              f"against the reference {RUNG_TARGET[vc]:.4f} -- {'OK' if ok else 'FAIL'}",
              flush=True)
    bad = [(c, s) for c, s in zip(co["chrom"].astype(str), co["stratum"].astype(int))
           if (c, s) in thin_pools]
    if bad:
        raise SystemExit(f"{len(bad)} coordinates have no alternative position in their matched pool; "
                         "the reference estimator cannot score these coordinates")
    print(f"[reference] no real coordinate sits in an undrawable pool ({len(co):,} checked)")
    rd = pd.DataFrame(rung)
    rd.to_csv(OUT / "null_calibration_rung.csv", index=False)
    if not rd.reproduces.all():
        raise SystemExit("Reference enrichment check failed; "
                         "calibration was stopped")

    # ---------------------------------------------------------------- calibration
    rows = []
    for win in [int(w) for w in a.windows.split(",")]:
        for vc in a.viruses.split(","):
            sub = co[co["virus_class"] == vc].reset_index(drop=True)
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), win)
            gen = np.random.default_rng(a.seed + 1000)
            t0 = time.time()
            for rep in range(a.n_rep):
                syn = synth_set(gen, nul, sub, labels, strat, thin_pools)
                if len(syn) == 0:
                    continue
                slab = F.leader_clusters(syn["chrom"].to_numpy(str),
                                         syn["pos0"].to_numpy(np.int64), win)
                thin = syn.loc[pd.Series(slab).drop_duplicates(keep="first").index]
                thin = thin.reset_index(drop=True)

                rng = np.random.default_rng(a.seed + rep)
                o1, n1 = F.est_independent(idx, nul, rng, syn, win, a.n_perm)
                rng = np.random.default_rng(a.seed + rep)
                o2, n2 = F.est_independent(idx, nul, rng, thin, win, a.n_perm)
                rng = np.random.default_rng(a.seed + rep)
                o3, n3 = F.est_block(idx, nul, rng, syn, win, a.n_perm, slab, unit="uniqSV")[:2]

                rows.append(dict(
                    window_bp=win, virus=vc, rep=rep, n_coords=len(syn),
                    n_clusters=int(pd.Series(slab).nunique()),
                    fold_published=o1 / float(n1.mean()),
                    p_published=(int(np.sum(n1 >= o1)) + 1) / (a.n_perm + 1),
                    fold_dedup=o2 / float(n2.mean()),
                    p_dedup=(int(np.sum(n2 >= o2)) + 1) / (a.n_perm + 1),
                    fold_corrected=o3 / float(n3.mean()),
                    p_corrected=(int(np.sum(n3 >= o3)) + 1) / (a.n_perm + 1)))
                if (rep + 1) % 10 == 0:
                    el = time.time() - t0
                    print(f"  [{vc} {win // 1000} kb] {rep + 1}/{a.n_rep} reps, "
                          f"{el:.0f}s elapsed, {el / (rep + 1) * (a.n_rep - rep - 1):.0f}s left",
                          flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "null_calibration_detail.csv", index=False)

    out = []
    for (win, vc), g in d.groupby(["window_bp", "virus"]):
        r = dict(window_bp=win, virus=vc, n_rep=len(g), n_coords=int(g.n_coords.median()),
                 n_events=n_ev)
        for tag in ("published", "dedup", "corrected"):
            r[f"{tag}_median_fold"] = round(float(g[f"fold_{tag}"].median()), 4)
            r[f"{tag}_p95_fold"] = round(float(g[f"fold_{tag}"].quantile(0.95)), 4)
            for al in ALPHAS:
                r[f"{tag}_FPR_{al}"] = rate(g[f"p_{tag}"].to_numpy(), al)
        out.append(r)
    s = pd.DataFrame(out)
    s.to_csv(OUT / "null_calibration.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["window_bp", "virus", "n_rep", "n_coords", "published_FPR_0.05",
            "dedup_FPR_0.05", "corrected_FPR_0.05"]
    print("\n" + s[cols].to_string(index=False))
    print("\nwrote null_calibration.csv / _detail.csv / _rung.csv")


if __name__ == "__main__":
    main()
