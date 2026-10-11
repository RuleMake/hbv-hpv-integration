"""Estimate the HTLV-1 class-composition contrast on the shared host reference.

The signed contrast is (TRA + INV) - (DEL + DUP), measured from donors without a
called viral integration. Joint-stratum matching and additional breakpoint-density
matching are evaluated separately. The reference viruses are recomputed using the
same propensity map, shrinkage, placement and summary routines."""

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
import locus_class_propensity as COMPOSITION  # noqa: E402

PROJ = Path(r"source_project")
OUT = PROJ / "results"
COORDS = OUT / "htlv_coordinates.csv"
REFERENCE = ("HBV", "HPV")
HTLV_ARMS = ("HTLV_all", "HTLV_ATL", "HTLV_nonmalignant", "HTLV_unknown",
             "HTLV_ATL_s35654946", "HTLV_nonmalig_s35654946")


def load_all_coordinates() -> pd.DataFrame:
    co = F.load_coordinates()[["virus_class", "chrom", "pos0"]].copy()
    h = pd.read_csv(COORDS)
    # PMID 35654946 is the only study contributing to both clinical groups, so it carries the
    # only malignant-vs-non-malignant contrast in this dataset that is not confounded with
    # study. Both of its groups are emitted as their own arms.
    # PRIMARY arm: every HTLV-1 coordinate, pooled. The premise being tested is "does a virus
    # that cuts its own site carry the host break grammar", and that question does not need the
    # clinical split at all -- every record here is integrase-chosen whatever the donor's disease
    # state. Pooling is therefore the primary analysis and the maximum-power one; the clinical
    # arms below are a sensitivity analysis on whether clonal selection matters, and they are
    # confounded with study (see htlv_build.py).
    allarm = h.copy()
    allarm["arm"] = "HTLV_all"

    s = h[h["pmid"] == 35654946]
    extra = []
    for grp, name in (("HTLV_ATL", "HTLV_ATL_s35654946"),
                      ("HTLV_nonmalignant", "HTLV_nonmalig_s35654946")):
        e = s[s["arm"] == grp].copy()
        e["arm"] = name
        extra.append(e)
    h = pd.concat([h, allarm] + extra, ignore_index=True)
    htlv = h.rename(columns={"arm": "virus_class"})[["virus_class", "chrom", "pos0"]]
    print(f"[htlv]    {len(htlv):,} HTLV-1 coordinates "
          f"{htlv['virus_class'].value_counts().to_dict()}")
    return pd.concat([co, htlv], ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--windows", type=int, nargs="*", default=[1_000_000])
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--viruses", default=",".join(REFERENCE + HTLV_ARMS))
    a = ap.parse_args()
    t0 = time.time()

    cum, nb, p_global = COMPOSITION.propensity_map()
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), J.WINDOW)
    allco = load_all_coordinates()
    joint = J.JointNuller(feat, F.Nuller(), mask)

    own_i = [COMPOSITION.CLASSES.index(c) for c in COMPOSITION.AXIS[0]]
    oth_i = [COMPOSITION.CLASSES.index(c) for c in COMPOSITION.AXIS[1]]
    viruses = tuple(a.viruses.split(","))
    axrows = []

    for window in a.windows:
        half = window // COMPOSITION.BIN
        dens = COMPOSITION.DensityNuller(joint, cum, nb, half)
        axis_draws = {}
        for vi, vc in enumerate(viruses):
            sub0, sub, tier, labels = COMPOSITION.prepare(allco, vc, joint, dens, mask)
            n_clu = int(len(np.unique(labels)))
            print(f"\n[{window//1000} kb | {vc}] {len(sub0):,} -> {len(sub):,} coordinates, "
                  f"{n_clu:,} clusters", flush=True)

            for arm, nuller, dt in (("stratum", joint, None), ("density", dens, tier)):
                # Use independent streams for each virus; a shared stream would
                # correlate the two viruses' nulls and make any crossed contrast
                # anti-conservative.
                rng = np.random.default_rng([a.seed, vi])
                nuller.n_calls = nuller.n_fallback = 0
                obs, obs_d, null, null_d, _ = COMPOSITION.score(sub, labels, nuller, cum, nb, half,
                                                        p_global, rng, a.n_perm,
                                                        density_tier=dt)
                base = dict(virus=vc, arm=arm, window_bp=window, n_coords=len(sub),
                            n_clusters=n_clu, n_perm=a.n_perm,
                            null_fallback_pct=round(100 * nuller.n_fallback
                                                    / max(nuller.n_calls, 1), 2))
                o_ax = obs[own_i].sum() - obs[oth_i].sum()
                n_ax = null[:, own_i].sum(axis=1) - null[:, oth_i].sum(axis=1)
                axrows.append({**base, **COMPOSITION.summarize(
                    o_ax, n_ax, "class_axis_TRAINV_minus_DELDUP", signed=True)})
                # Local breakpoint density is recorded for every arm because it is the
                # competing explanation for any residual class-atlas elevation: a coordinate
                # set that simply sits in break-dense chromatin will show small folds on the
                # abundant classes without carrying any grammar. The number belongs beside
                # the atlas, not in a footnote.
                axrows.append({**base, **COMPOSITION.summarize(obs_d, null_d, "breakpoint_density")})
                for ci, cls in enumerate(COMPOSITION.CLASSES):
                    axrows.append({**base, **COMPOSITION.summarize(obs[ci], null[:, ci],
                                                           f"prop_{cls}")})
                axis_draws[(arm, vc)] = (o_ax - n_ax.mean(), n_ax - n_ax.mean())
                print(f"  [{arm:7s}] " + "  ".join(
                    f"{c} {obs[i]:.4f}/{null[:, i].mean():.4f}"
                    for i, c in enumerate(COMPOSITION.CLASSES))
                    + f"   axis {o_ax - n_ax.mean():+.5f}"
                    + f"   fallback {base['null_fallback_pct']}%", flush=True)

        # every crossed contrast among the viruses present, not only HBV-vs-HPV
        names = [v for v in viruses]
        for arm in ("stratum", "density"):
            for i in range(len(names)):
                for k in range(i + 1, len(names)):
                    va, vb = names[i], names[k]
                    (ea, na), (eb, nbv) = axis_draws[(arm, va)], axis_draws[(arm, vb)]
                    r = float(np.corrcoef(na, nbv)[0, 1])
                    axrows.append({
                        "virus": f"{va}-vs-{vb}", "arm": arm, "window_bp": window,
                        "n_coords": np.nan, "n_clusters": np.nan, "n_perm": a.n_perm,
                        "null_fallback_pct": np.nan, "null_stream_corr": round(r, 4),
                        **COMPOSITION.summarize(ea - eb, na - nbv,
                                        f"interaction_axis_{va}_minus_{vb}", signed=True)})
                    print(f"  [{arm:7s}] interaction {va} - {vb}: {ea - eb:+.5f}   "
                          f"null-stream corr {r:+.4f}", flush=True)

    dfa = pd.DataFrame(axrows)
    dfa.to_csv(OUT / "integrase_positive_axis.csv", index=False)

    # Check the composition estimates against the reference table.
    pub = pd.read_csv(OUT / "locus_class_propensity_axes.csv")
    m = dfa.merge(pub[["virus", "arm", "window_bp", "statistic", "observed", "excess"]],
                  on=["virus", "arm", "window_bp", "statistic"], suffixes=("", "_pub"))
    m = m[m["virus"].isin(REFERENCE)]
    if len(m):
        dev = (m["observed"] - m["observed_pub"]).abs().max()
        print(f"\n[check]    {len(m)} reference rows; max |observed - published| = {dev:.6f}")
        if dev > 1e-4:
            raise SystemExit("reproduction gate FAILED on the axis table")
    else:
        print("\n[check]    No comparable reference rows; reference check unavailable")

    cols = ["window_bp", "virus", "arm", "statistic", "observed", "null_mean", "excess",
            "z", "emp_p_two_sided", "mde_excess", "clears_mde"]
    print("\n" + dfa[cols].to_string(index=False))
    print(f"\nwrote {OUT / 'integrase_positive_axis.csv'}   [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
