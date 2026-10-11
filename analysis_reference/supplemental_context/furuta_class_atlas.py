"""Computes class-specific spatial enrichment for the Furuta non-tumour coordinates.
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
import furuta_nontumour_test as FU  # noqa: E402

OUT = Path(r"source_project\results")
WINDOW = 100_000
SEED = 7


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    a = ap.parse_args()
    t0 = time.time()

    d = FU.load_furuta()
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)
    canon = pd.read_csv(F.BASE / "input_pooled_canonical_events.tsv", sep="\t")
    led = ND.build_ledger(set(canon["aliquot_metadata"].astype(str)))

    rows = []
    for cls in ND.CLASSES:
        idx, n_ev = ND.to_index(led[led.svclass == cls])
        for tag, sel in (("NON-TUMOUR (unselected)", d.tissue == "NT"),
                         ("TUMOUR (same donors)", d.tissue == "T")):
            sub = FU.prep(d[sel], nul, mask)
            if len(sub) < 30:
                continue
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), WINDOW)
            nul.n_calls = nul.n_fallback = 0
            rng = np.random.default_rng(SEED)
            o, nl, _, so, sn = F.est_block(idx, nul, rng, sub, WINDOW, a.n_perm,
                                           labels, unit="uniqSV")
            fb = round(100 * nul.n_fallback / max(nul.n_calls, 1), 2)
            for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                             ("f_clusterSignTest", F.summarize(so, sn))):
                rows.append(dict(svclass=cls, subset=tag, estimator=est, n_events=n_ev,
                                 n_coords=len(sub),
                                 n_clusters=int(len(np.unique(labels))),
                                 pct_pool_fallback=fb, **res))
            print(f"  [{cls:3s} {tag:24s}] bc={rows[-2]['fold']:.4f} "
                  f"[{rows[-2]['fold_ci_lo']:.4f},{rows[-2]['fold_ci_hi']:.4f}] "
                  f"p={rows[-2]['emp_p']:.5g} | sign={rows[-1]['fold']:.4f} "
                  f"p={rows[-1]['emp_p']:.5g}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "furuta_class_atlas.csv", index=False)

    # Check non-tumour duplication enrichment against the reference value of 1.0213.
    gate = df[(df.svclass == "DUP") & (df.subset.str.startswith("NON-TUMOUR"))
              & (df.estimator == "bc_uniqSV_blockNull")]["fold"]
    if len(gate):
        drift = abs(float(gate.iloc[0]) - 1.0213)
        print(f"\n[check] NT DUP {gate.iloc[0]:.4f} vs canonical-ledger 1.0213 "
              f"(drift {drift:.4f}) -> {'OK' if drift <= 0.01 else 'CHECK'}")
    pd.set_option("display.width", 240)
    print("\n" + df[["svclass", "subset", "estimator", "n_coords", "fold",
                     "fold_ci_lo", "fold_ci_hi", "emp_p"]].to_string(index=False))
    print(f"\nwrote {OUT / 'furuta_class_atlas.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
