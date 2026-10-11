"""Estimate the resolution of external class-composition comparisons.

The script counts eligible coordinates and leader clusters in the Sasa and CLCA
catalogues after amplification masking and matched-stratum lookup. Minimum detectable
displacements are projected from the primary analysis's null standard deviation,
scaled by the number of independent clusters. Results are printed without writing files."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import corrected_estimator_baseline as F                    # noqa: E402
import joint_correction as J                   # noqa: E402
import amplification_and_histology_control as A    # noqa: E402
import locus_class_propensity as E              # noqa: E402

BIN = F.BIN_SIZE
WINDOW = 100_000


def prep_any(df: pd.DataFrame, joint, mask, label: str):
    """Amplicon mask + joint stratum + leader clustering, for any (chrom, pos0) frame."""
    n0 = len(df)
    df = df.dropna(subset=["chrom", "pos0"]).copy()
    df["chrom"] = df["chrom"].astype(str)
    df["pos0"] = df["pos0"].astype(np.int64)
    df = df.drop_duplicates(["chrom", "pos0"]).sort_values(["chrom", "pos0"]).reset_index(drop=True)
    n1 = len(df)
    in_mask = np.array([(c, int(p) // BIN * BIN) in mask
                        for c, p in zip(df["chrom"], df["pos0"])])
    df = df[~in_mask].reset_index(drop=True)
    key = pd.MultiIndex.from_arrays([df["chrom"].to_numpy(str),
                                     (df["pos0"].to_numpy(np.int64) // BIN) * BIN])
    df = df.assign(stratum=joint.lookup.reindex(key).to_numpy())
    df = df[df["stratum"].notna()].reset_index(drop=True)
    labels = F.leader_clusters(df["chrom"].to_numpy(str),
                               df["pos0"].to_numpy(np.int64), WINDOW)
    n_clu = int(len(np.unique(labels)))
    print(f"[{label}] {n0} raw -> {n1} unique -> {len(df)} after mask+stratum -> "
          f"{n_clu} clusters")
    return df, labels, n_clu


def main() -> None:
    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    joint = J.JointNuller(feat, F.Nuller(), mask)

    # Tian reference from the 1-Mb stratum-matched comparison
    ax = pd.read_csv(E.PROJ / "results" / "locus_class_propensity_axes.csv")
    ref = ax[(ax.window_bp == 1_000_000) & (ax.arm == "stratum")
             & ax.statistic.str.startswith("class_axis")]
    print("\n=== Tian reference (1 Mb, stratum-matched) ===")
    for _, r in ref.iterrows():
        sd = r.mde_excess / 2.80
        print(f"  {r.virus}: {int(r.n_clusters)} clusters, |excess| {abs(r.excess):.5f}, "
              f"null SD {sd:.5f}, MDE {r.mde_excess:.5f}, "
              f"excess/MDE = {abs(r.excess)/r.mde_excess:.2f}")
    tian = {r.virus: (int(r.n_clusters), abs(r.excess), r.mde_excess / 2.80)
            for _, r in ref.iterrows()}

    # --- Sasa 2025 (HPV)
    import sasa_external_replication as S
    sasa = S.lift(S.load_sasa())
    sasa_df, _, sasa_clu = prep_any(sasa[["chrom", "pos0"]], joint, mask, "Sasa HPV")

    # --- CLCA (HBV)
    import clca_class_replication as C
    clca = C.load_integrations()
    clca_df, _, clca_clu = prep_any(clca[["chrom", "pos0"]], joint, mask, "CLCA HBV")

    print("\n=== projected MDE, assuming the same per-cluster variance as Tian ===")
    print(f"{'catalogue':16s} {'clusters':>9s} {'scale':>7s} {'proj. MDE':>10s} "
          f"{'Tian effect':>12s} {'detectable?':>12s}")
    rows = [("Sasa 2025 HPV", sasa_clu, "HPV"), ("CLCA HBV", clca_clu, "HBV")]
    for name, n_clu, ref_v in rows:
        n_ref, eff, sd_ref = tian[ref_v]
        scale = np.sqrt(n_ref / n_clu)
        mde = 2.80 * sd_ref * scale
        ok = "YES" if mde < eff else "NO"
        print(f"{name:16s} {n_clu:9d} {scale:7.2f}x {mde:10.5f} {eff:12.5f} {ok:>12s}")

    print("\nA catalogue is only informative here if its projected MDE is BELOW the effect Tian\n"
          "measured. Otherwise a null means the instrument could not see it, not that it is absent.")


if __name__ == "__main__":
    main()
