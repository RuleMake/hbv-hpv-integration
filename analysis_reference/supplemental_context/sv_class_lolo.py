"""Computes leave-one-cluster-out sensitivity for class-specific enrichment.
"""

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
import discriminating_tests as ND  # noqa: E402

OUT = Path(r"source_project\results")
WINDOW = 100_000


def expected_folds() -> dict[tuple[str, str], float]:
    d = pd.read_csv(OUT / "sv_class_locus_level.csv")
    bc = d[d["estimator"] == "bc_uniqSV_blockNull"]
    return {(str(r.svclass), str(r.virus)): float(r.fold) for r in bc.itertuples()}


def est_block_uniqsv_percluster(idx, nul, rng, sub, win, n_perm, labels):
    """Byte-for-byte the uniqSV branch of F.est_block, retaining per-cluster nulls."""
    df = sub.copy()
    df["_cluster"] = labels
    obs_total = 0
    null_total = np.zeros(n_perm, dtype=np.int64)
    per_obs, per_null, meta = [], [], []
    for cid, g in df.groupby("_cluster", sort=True):
        c = str(g["chrom"].iloc[0])
        t = idx.get(c)
        if t is None:
            continue
        p = g["pos0"].to_numpy(np.int64)
        st = g["stratum"].to_numpy()
        off = g["within_bin_offset"].to_numpy()

        j = rng.integers(0, len(p), size=n_perm)
        delta = np.empty(n_perm, dtype=np.int64)
        for k in range(len(p)):
            mask = j == k
            nk = int(mask.sum())
            if nk:
                delta[mask] = nul.draw(rng, c, int(st[k]), int(off[k]), int(p[k]), nk) - int(p[k])
        lo_b, hi_b = int(p.min()) - win, int(p.max()) + win
        L = hi_b - lo_b
        o = int(t.count(L, np.array([lo_b]))[0])
        nl = t.count(L, lo_b + delta)

        obs_total += o
        null_total += nl
        per_obs.append(o)
        per_null.append(nl.astype(np.int32))
        meta.append({"cluster": int(cid), "chrom": c, "n_coords": len(p),
                     "start": lo_b, "end": hi_b, "obs": o,
                     "null_mean": float(nl.mean())})
    return obs_total, null_total, np.array(per_obs), per_null, pd.DataFrame(meta)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--top", type=int, default=12, help="clusters to report individually")
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), 100_000)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)
    co = F.load_coordinates()
    canon = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
    led = ND.build_ledger(set(canon["aliquot_metadata"].astype(str)))
    expected = expected_folds()

    per_rows, summary = [], []
    for cls in ND.CLASSES:
        idx, _ = ND.to_index(led[led.svclass == cls])
        for vc in ("HBV", "HPV", "EBV"):
            sub = ND.eligible(co, nul, mask, vc)
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), WINDOW)
            rng = np.random.default_rng(a.seed)
            o_tot, n_tot, per_obs, per_null, meta = est_block_uniqsv_percluster(
                idx, nul, rng, sub, WINDOW, a.n_perm, labels)

            full = F.summarize(o_tot, n_tot)
            exp_fold = expected[(cls, vc)]
            drift = abs(full["fold"] - exp_fold)
            flag = "OK" if drift <= 0.004 else "DRIFT"
            print(f"[{cls:3s} {vc}] full fold {full['fold']} (on-disk {exp_fold}) -> {flag}; "
                  f"obs={o_tot} clusters={len(per_obs)}", flush=True)
            if flag == "DRIFT":
                raise SystemExit(f"{cls}/{vc}: reproduction of sv_class_locus_level.csv failed")

            excess = per_obs - meta["null_mean"].to_numpy()
            meta = meta.assign(excess=excess,
                               excess_share=excess / max(excess.sum(), 1e-9))
            order = np.argsort(-excess)

            folds = np.empty(len(per_obs))
            for i in range(len(per_obs)):
                ni = n_tot - per_null[i].astype(np.int64)
                folds[i] = (o_tot - per_obs[i]) / ni.mean() if ni.mean() else np.nan
            meta = meta.assign(fold_without=folds)

            for i in order[:a.top]:
                oi = o_tot - per_obs[i]
                ni = n_tot - per_null[i].astype(np.int64)
                p_i = (int(np.sum(ni >= oi)) + 1) / (a.n_perm + 1)
                r = meta.iloc[i]
                per_rows.append({"svclass": cls, "virus": vc, "chrom": r["chrom"],
                                 "start": int(r["start"]), "end": int(r["end"]),
                                 "n_coords": int(r["n_coords"]), "obs": int(r["obs"]),
                                 "null_mean": round(float(r["null_mean"]), 2),
                                 "excess": round(float(r["excess"]), 2),
                                 "excess_share_pct": round(100 * float(r["excess_share"]), 2),
                                 "fold_without": round(float(r["fold_without"]), 4),
                                 "p_without": round(p_i, 6)})

            fw = meta["fold_without"].to_numpy()
            top10_share = 100 * np.sort(-excess)[:10].sum() / -max(excess.sum(), 1e-9)
            summary.append({
                "svclass": cls, "virus": vc, "n_clusters": len(per_obs), "n_perm": a.n_perm,
                "fold_full": full["fold"], "p_full": full["emp_p"],
                "fold_lolo_min": round(float(np.nanmin(fw)), 4),
                "fold_lolo_max": round(float(np.nanmax(fw)), 4),
                "worst_case_still_above_1": bool(np.nanmin(fw) > 1.0),
                "top10_excess_share_pct": round(float(top10_share), 1),
            })
            print(f"      LOLO fold range [{np.nanmin(fw):.4f}, {np.nanmax(fw):.4f}]; "
                  f"top-10 clusters carry {top10_share:.1f}% of the excess", flush=True)

    pd.DataFrame(per_rows).to_csv(OUT / "sv_class_lolo.csv", index=False)
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "sv_class_lolo_summary.csv", index=False)
    pd.set_option("display.width", 240)
    print("\n" + s.to_string(index=False))
    print(f"\nwrote {OUT / 'sv_class_lolo.csv'} and _summary.csv ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
