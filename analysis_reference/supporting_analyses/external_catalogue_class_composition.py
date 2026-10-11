"""Estimate local class composition at CLCA and Sasa integration coordinates.

The coordinates are independent of the primary catalogue; the propensity map
remains the PCAWG integration-negative reference. Estimates, intervals and minimum
detectable displacements are reported together because the external catalogues
have fewer independent coordinate clusters than the primary analysis."""

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
from replication_catalogue_power_probe import prep_any    # noqa: E402

OUT = E.PROJ / "results" / "external_locus_class_propensity.csv"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--window", type=int, default=1_000_000)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t0 = time.time()

    cum, nb, p_global = E.propensity_map()
    half = a.window // E.BIN
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), J.WINDOW)
    joint = J.JointNuller(feat, F.Nuller(), mask)
    dens = E.DensityNuller(joint, cum, nb, half)

    import sasa_external_replication as S
    import clca_class_replication as C
    cats = [
        ("Sasa2025", "HPV", S.lift(S.load_sasa())[["chrom", "pos0"]]),
        ("CLCA", "HBV", C.load_integrations()[["chrom", "pos0"]]),
    ]

    own_i = [E.CLASSES.index(c) for c in E.AXIS[0]]
    oth_i = [E.CLASSES.index(c) for c in E.AXIS[1]]
    rows, axis_draws = [], {}

    for vi, (cat, virus, raw) in enumerate(cats):
        sub, labels, n_clu = prep_any(raw, joint, mask, f"{cat} {virus}")
        tier = dens.tier.reindex(pd.MultiIndex.from_arrays(
            [sub["chrom"].to_numpy(str),
             (sub["pos0"].to_numpy(np.int64) // E.BIN) * E.BIN])).to_numpy()
        sub["stratum"] = sub["stratum"].astype(int)

        for arm, nuller, dt in (("stratum", joint, None), ("density", dens, tier)):
            rng = np.random.default_rng([a.seed, vi])
            nuller.n_calls = nuller.n_fallback = 0
            obs, obs_d, null, null_d, _ = E.score(sub, labels, nuller, cum, nb, half,
                                                  p_global, rng, a.n_perm, density_tier=dt)
            base = dict(catalogue=cat, virus=virus, arm=arm, window_bp=a.window,
                        n_coords=len(sub), n_clusters=n_clu, n_perm=a.n_perm)
            for ci, cls in enumerate(E.CLASSES):
                rows.append({**base, **E.summarize(obs[ci], null[:, ci], f"prop_{cls}")})
            rows.append({**base, **E.summarize(obs_d, null_d, "breakpoint_density")})
            o_ax = obs[own_i].sum() - obs[oth_i].sum()
            n_ax = null[:, own_i].sum(axis=1) - null[:, oth_i].sum(axis=1)
            rows.append({**base, **E.summarize(o_ax, n_ax, "class_axis_TRAINV_minus_DELDUP",
                                               signed=True)})
            axis_draws[(arm, cat)] = (o_ax - n_ax.mean(), n_ax - n_ax.mean())
            print(f"  [{cat} {arm}] axis {o_ax - n_ax.mean():+.5f}  "
                  f"MDE {2.80*n_ax.std(ddof=1):.5f}", flush=True)

    # crossed check: CLCA (HBV side) minus Sasa (HPV side) -- both fully independent of Tian
    for arm in ("stratum", "density"):
        (eh, nh), (ep, npv) = axis_draws[(arm, "CLCA")], axis_draws[(arm, "Sasa2025")]
        r = float(np.corrcoef(nh, npv)[0, 1])
        rows.append({
            "catalogue": "CLCA-vs-Sasa2025", "virus": "HBV-vs-HPV", "arm": arm,
            "window_bp": a.window, "n_coords": np.nan, "n_clusters": np.nan,
            "n_perm": a.n_perm, "null_stream_corr": round(r, 4),
            **E.summarize(eh - ep, nh - npv, "interaction_axis_HBV_minus_HPV", signed=True)})
        print(f"  [{arm}] interaction {eh - ep:+.5f}  null-stream corr {r:+.4f}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    cols = ["catalogue", "virus", "arm", "statistic", "observed", "null_mean", "excess",
            "z", "emp_p_two_sided", "mde_excess", "clears_mde", "n_clusters"]
    key = df[df.statistic.str.contains("class_axis|interaction")]
    print("\n" + key[cols].to_string(index=False))
    print(f"\nwrote {OUT}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
