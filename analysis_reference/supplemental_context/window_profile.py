"""Computes class-enrichment profiles across spatial window sizes.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F          # noqa: E402
import joint_correction as J         # noqa: E402
import amplification_and_histology_control as A   # noqa: E402
import discriminating_tests as ND     # noqa: E402

OUT = Path(r"source_project\results")
COORDS = OUT / "htlv_coordinates.csv"
CLASSES = ["DUP", "DEL", "INV", "TRA"]
WINDOWS = (1_000, 5_000, 10_000, 25_000, 50_000, 100_000, 250_000, 500_000)
SEED = 7
Z_FLOOR = 0.05


def mde(fold: float, z: float) -> float:
    if not np.isfinite(fold) or not np.isfinite(z) or abs(z) < Z_FLOOR:
        return np.inf
    return 1.0 + 2.80 * abs((fold - 1.0) / z)


def load_all_coordinates(with_htlv: bool) -> pd.DataFrame:
    co = F.load_coordinates()[["virus_class", "chrom", "pos0"]].copy()
    if not with_htlv:
        return co
    h = pd.read_csv(COORDS)
    h = h.assign(virus_class="HTLV_all")[["virus_class", "chrom", "pos0"]]
    print(f"[htlv]    {len(h):,} HTLV-1 coordinates pooled as HTLV_all")
    return pd.concat([co, h], ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--viruses", default="HBV,HPV,EBV")
    ap.add_argument("--windows", type=int, nargs="*", default=list(WINDOWS))
    ap.add_argument("--out-suffix", default="")
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    viruses = tuple(a.viruses.split(","))
    allco = load_all_coordinates(any(v.startswith("HTLV") for v in viruses))

    canon = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
    led = ND.build_ledger(set(canon["aliquot_metadata"].astype(str)))
    idx_by_class = {c: ND.to_index(led[led.svclass == c]) for c in CLASSES}

    rows = []
    for window in a.windows:
        # The amplicon mask is built at the scoring window, exactly as every other script does:
        # a peak that matters at 500 kb does not matter at 1 kb, and using one mask for all
        # windows would make the small-window arms carry the large-window exclusion.
        mask = A.masked_bins(A.amp_peaks(), window)
        nul = J.JointNuller(feat, F.Nuller(), mask)
        for cls in CLASSES:
            idx, n_ev = idx_by_class[cls]
            for virus in viruses:
                sub = ND.eligible(allco, nul, mask, virus, win=window)
                labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                           sub["pos0"].to_numpy(np.int64), window)
                n_clu = int(len(np.unique(labels)))
                rng = np.random.default_rng(SEED)
                o, nl, _, so, sn = F.est_block(idx, nul, rng, sub, window, a.n_perm,
                                               labels, unit="uniqSV")
                for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                                 ("f_clusterSignTest", F.summarize(so, sn))):
                    res["mde_fold"] = round(mde(res["fold"], res["z"]), 4)
                    rows.append(dict(window_bp=window, svclass=cls, n_events=n_ev, virus=virus,
                                     estimator=est, n_coords=len(sub), n_clusters=n_clu, **res))
                b = rows[-2]
                print(f"  [{window//1000:>3} kb {cls:3s} {virus:10s}] "
                      f"n={b['n_coords']:6,}/{n_clu:6,}cl  fold={b['fold']:.4f} "
                      f"[{b['fold_ci_lo']:.4f},{b['fold_ci_hi']:.4f}] p={b['emp_p']:.4g} "
                      f"MDE={b['mde_fold']:.3f}", flush=True)

    d = pd.DataFrame(rows)
    fn = OUT / f"window_profile{a.out_suffix}.csv"
    d.to_csv(fn, index=False)

    piv = d[d.estimator == "bc_uniqSV_blockNull"].pivot_table(
        index=["virus", "svclass"], columns="window_bp", values="fold")
    pd.set_option("display.width", 240)
    print("\n=== fold by window (bc_uniqSV_blockNull) ===")
    print(piv.to_string(float_format=lambda x: f"{x:.4f}"))
    print(f"\nwrote {fn}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
