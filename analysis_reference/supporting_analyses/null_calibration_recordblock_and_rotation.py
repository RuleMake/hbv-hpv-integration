"""Separate counting-unit and placement effects in spatial-null calibration.

The record-counting comparison retains the observed record unit while replacing
independent placements with a cluster-block null. One simulation uses the same
rigid-cluster generator and seeds as the baseline calibration, with pairing checked
against saved draws. A second rotates each chromosome's event track by an independent
uniform offset while leaving the real coordinate catalogue unchanged.

The g1 and g1b options select these simulations; smoke runs a small diagnostic case."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F                          # noqa: E402
import null_calibration_baseline as E               # noqa: E402

RES = Path(r"source_project\results")
PREREG = Path(r"source_project\docs\PREREG_G1_2026-10-03.md")
ALPHAS = (0.05, 0.01)             # 0.001 is not resolvable at n_perm = 500
ARMS_G1 = ("recordblock",)
ARMS_G1B = ("published", "dedup", "recordblock", "corrected")


def quiet(fn, *a, **k):
    """est_block prints a per-call overhang line; keep the log readable."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def pval(null: np.ndarray, obs: int, n_perm: int) -> float:
    return (int(np.sum(null >= obs)) + 1) / (n_perm + 1)


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return ctr - half, ctr + half


def summarize(d: pd.DataFrame, arms) -> pd.DataFrame:
    rows = []
    groups = [((w, v), g) for (w, v), g in d.groupby(["window_bp", "virus"])]
    groups.append((("pooled", "all"), d))
    for (w, v), g in groups:
        for arm in arms:
            p = g[f"p_{arm}"].to_numpy()
            r = dict(window_bp=w, virus=v, arm=arm, n_rep=len(g),
                     median_fold=round(float(g[f"fold_{arm}"].median()), 4),
                     p95_fold=round(float(g[f"fold_{arm}"].quantile(0.95)), 4),
                     median_p=round(float(np.median(p)), 3))
            for al in ALPHAS:
                k = int(np.sum(p <= al))
                lo, hi = wilson(k, len(p))
                r[f"FPR_{al}"] = round(k / len(p), 4)
                r[f"FPR_{al}_lo"] = round(lo, 4)
                r[f"FPR_{al}_hi"] = round(hi, 4)
            rows.append(r)
    return pd.DataFrame(rows)


def prereg_stamp() -> str:
    h = hashlib.sha256(PREREG.read_bytes()).hexdigest()
    m = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(PREREG.stat().st_mtime))
    print(f"[design] {PREREG.name}  sha256 {h}  mtime {m}  (local timestamp, not external)",
          flush=True)
    return h


# ---------------------------------------------------------------------------- stage g1
def stage_g1(a, out: Path) -> None:
    stored = pd.read_csv(RES / "null_calibration_detail.csv")
    stored = stored.set_index(["window_bp", "virus", "rep"])
    idx, _ = F.load_events()
    nul = F.Nuller()
    co = F.load_coordinates()
    strat, thin_pools = E.bin_lookup()

    rows = []
    n_pair_ok = 0
    for win in [int(w) for w in a.windows.split(",")]:
        for vc in a.viruses.split(","):
            sub = co[co["virus_class"] == vc].reset_index(drop=True)
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), win)
            gen = np.random.default_rng(a.seed + 1000)          # identical to the baseline calibration
            t0 = time.time()
            for rep in range(a.n_rep):
                syn = E.synth_set(gen, nul, sub, labels, strat, thin_pools)
                if len(syn) == 0:
                    continue
                slab = F.leader_clusters(syn["chrom"].to_numpy(str),
                                         syn["pos0"].to_numpy(np.int64), win)
                key = (win, vc, rep)
                if not a.smoke:
                    st = stored.loc[key]
                    n_cl = int(pd.Series(slab).nunique())
                    if int(st.n_coords) != len(syn) or int(st.n_clusters) != n_cl:
                        raise SystemExit(f"PAIRING FAILED at {key}: regenerated "
                                         f"({len(syn)}, {n_cl}) vs stored "
                                         f"({int(st.n_coords)}, {int(st.n_clusters)})")
                    n_pair_ok += 1

                rng = np.random.default_rng(a.seed + rep)
                o, n = quiet(F.est_block, idx, nul, rng, syn, win, a.n_perm, slab,
                             unit="pairs")[:2]
                row = dict(window_bp=win, virus=vc, rep=rep, n_coords=len(syn),
                           n_clusters=int(pd.Series(slab).nunique()),
                           fold_recordblock=o / float(n.mean()),
                           p_recordblock=pval(n, o, a.n_perm))

                if not a.smoke and rep < 5:       # recompute stored arms: same stream -> same p
                    rng = np.random.default_rng(a.seed + rep)
                    o1, n1 = F.est_independent(idx, nul, rng, syn, win, a.n_perm)
                    if abs(pval(n1, o1, a.n_perm) - float(st.p_published)) > 1e-12:
                        raise SystemExit(f"PAIRING FAILED (published p) at {key}")
                    if rep < 3:
                        rng = np.random.default_rng(a.seed + rep)
                        o3, n3 = quiet(F.est_block, idx, nul, rng, syn, win, a.n_perm, slab,
                                       unit="uniqSV")[:2]
                        if abs(pval(n3, o3, a.n_perm) - float(st.p_corrected)) > 1e-12:
                            raise SystemExit(f"PAIRING FAILED (corrected p) at {key}")
                    row["pairing_p_checked"] = True
                rows.append(row)
                if (rep + 1) % 10 == 0:
                    el = time.time() - t0
                    print(f"  [g1 {vc} {win // 1000} kb] {rep + 1}/{a.n_rep}  {el:.0f}s elapsed, "
                          f"{el / (rep + 1) * (a.n_rep - rep - 1):.0f}s left", flush=True)
            pd.DataFrame(rows).to_csv(out / "record_block_calibration_detail.csv", index=False)

    d = pd.DataFrame(rows)
    if not a.smoke:
        print(f"[pairing] n_coords / n_clusters identical to the baseline calibration in {n_pair_ok}/{len(d)} "
              f"replicates; published p (first 5 per cell) and corrected p (first 3 per cell) "
              f"reproduced exactly", flush=True)
        if n_pair_ok != len(d) or len(d) != 600:
            raise SystemExit(f"expected 600 paired replicates, got {n_pair_ok}/{len(d)}")
        # join the stored arms so that all four arms sit side by side
        d = d.join(stored[["fold_published", "p_published", "fold_dedup", "p_dedup",
                           "fold_corrected", "p_corrected"]], on=["window_bp", "virus", "rep"])
        arms = ("published", "dedup", "recordblock", "corrected")
    else:
        arms = ARMS_G1
    d.to_csv(out / "record_block_calibration_detail.csv", index=False)
    s = summarize(d, arms)
    s.to_csv(out / "record_block_calibration_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print("\n" + s[["window_bp", "virus", "arm", "n_rep", "median_fold", "FPR_0.05", "FPR_0.05_lo",
                    "FPR_0.05_hi", "FPR_0.01"]].to_string(index=False), flush=True)
    verdict_g1(s)


def verdict_g1(s: pd.DataFrame) -> None:
    r = s[s.arm == "recordblock"]
    pooled = float(r[r.window_bp.astype(str) == "pooled"]["FPR_0.05"].iloc[0])
    cells = r[r.window_bp.astype(str) != "pooled"]["FPR_0.05"]
    if pooled <= 0.08 and cells.max() <= 0.12:
        v = "CALIBRATED"
    elif pooled > 0.08:
        v = "MISCALIBRATED"
    else:
        v = "BORDERLINE"
    print(f"\n[Block-generator calibration] recordblock pooled FPR@0.05 = {pooled:.4f}, "
          f"worst cell {cells.max():.4f} -> {v}", flush=True)


# ---------------------------------------------------------------------------- stage g1b
def rotate_events(idx, nul, rng):
    """Circularly rotate each chromosome's TD-event track by an independent uniform offset.
    Events whose two endpoints straddle the wrap point are dropped.
    Returns (new idx, n dropped, n total)."""
    new, dropped, total = {}, 0, 0
    for c, t in idx.items():
        size = max(int(nul.chrom_size.get(c, 0)), int(t.right.max()) + 1)
        u = int(rng.integers(0, size))
        l2, r2 = t.left + u, t.right + u
        keep = (l2 // size) == (r2 // size)
        dropped += int((~keep).sum())
        total += len(keep)
        new[c] = F.ChromEvents((l2 % size)[keep], (r2 % size)[keep], t.eid[keep])
    return new, dropped, total


def stage_g1b(a, out: Path) -> None:
    idx0, _ = F.load_events()
    nul = F.Nuller()
    co = F.load_coordinates()
    rows = []
    for win in [int(w) for w in a.windows.split(",")]:
        for vc in a.viruses.split(","):
            sub = co[co["virus_class"] == vc].reset_index(drop=True)
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), win)
            thin = sub.loc[pd.Series(labels).drop_duplicates(keep="first").index]
            thin = thin.reset_index(drop=True)
            vcode = 0 if vc == "HBV" else 1
            t0 = time.time()
            for rep in range(a.n_rep):
                idx, dr, tot = rotate_events(
                    idx0, nul, np.random.default_rng([5000, rep, win, vcode]))
                row = dict(window_bp=win, virus=vc, rep=rep, n_coords=len(sub),
                           n_clusters=int(pd.Series(labels).nunique()),
                           events_dropped_frac=dr / tot)
                rng = np.random.default_rng(a.seed + rep)
                o, n = F.est_independent(idx, nul, rng, sub, win, a.n_perm)
                row.update(fold_published=o / float(n.mean()), p_published=pval(n, o, a.n_perm))
                rng = np.random.default_rng(a.seed + rep)
                o, n = F.est_independent(idx, nul, rng, thin, win, a.n_perm)
                row.update(fold_dedup=o / float(n.mean()), p_dedup=pval(n, o, a.n_perm))
                rng = np.random.default_rng(a.seed + rep)
                o, n = quiet(F.est_block, idx, nul, rng, sub, win, a.n_perm, labels,
                             unit="pairs")[:2]
                row.update(fold_recordblock=o / float(n.mean()),
                           p_recordblock=pval(n, o, a.n_perm))
                rng = np.random.default_rng(a.seed + rep)
                o, n = quiet(F.est_block, idx, nul, rng, sub, win, a.n_perm, labels,
                             unit="uniqSV")[:2]
                row.update(fold_corrected=o / float(n.mean()), p_corrected=pval(n, o, a.n_perm))
                rows.append(row)
                if (rep + 1) % 10 == 0:
                    el = time.time() - t0
                    print(f"  [g1b {vc} {win // 1000} kb] {rep + 1}/{a.n_rep}  {el:.0f}s elapsed, "
                          f"{el / (rep + 1) * (a.n_rep - rep - 1):.0f}s left", flush=True)
            pd.DataFrame(rows).to_csv(out / "chromosome_rotation_calibration_detail.csv", index=False)

    d = pd.DataFrame(rows)
    s = summarize(d, ARMS_G1B)
    s.to_csv(out / "chromosome_rotation_calibration_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print("\n" + s[["window_bp", "virus", "arm", "n_rep", "median_fold", "FPR_0.05", "FPR_0.05_lo",
                    "FPR_0.05_hi", "FPR_0.01"]].to_string(index=False), flush=True)
    print(f"\n[rotation] mean fraction of events dropped at the wrap point: "
          f"{d.events_dropped_frac.mean():.5f}", flush=True)
    p = s[s.window_bp.astype(str) == "pooled"].set_index("arm")["FPR_0.05"]
    for arm in ARMS_G1B:
        if arm == "published":
            v = "stays MISCALIBRATED" if p[arm] > 0.15 else "NOT miscalibrated (>0.15 not met)"
        else:
            v = "remains CALIBRATED" if p[arm] <= 0.08 else "GENERATOR-SPECIFIC / miscalibrated"
        print(f"[Rotation-generator calibration] {arm}: pooled FPR@0.05 = {p[arm]:.4f} "
              f"-> {v}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["g1", "g1b", "smoke"], required=True)
    ap.add_argument("--n-rep", type=int, default=150)
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--windows", default="100000,1000000")
    ap.add_argument("--viruses", default="HBV,HPV")
    a = ap.parse_args()
    a.smoke = a.stage == "smoke"
    prereg_stamp()
    if a.smoke:
        import tempfile
        out = Path(tempfile.mkdtemp(prefix="null_calibration_example_"))
        a.n_rep, a.n_perm = 3, 50
        a.windows, a.viruses = "100000", "HBV"
        print(f"[smoke] scratch output dir {out}", flush=True)
        stage_g1(a, out)
        stage_g1b(a, out)
        return
    {"g1": stage_g1, "g1b": stage_g1b}[a.stage](a, RES)


if __name__ == "__main__":
    main()
