"""Evaluate class-specific spatial enrichment in the CLCA HBV liver cohort.

CLCA integration coordinates and per-sample structural variants are scored with
the joint-matched cluster-block null in 100-kb windows, using 20,000 permutations
and seed 7. DUP, DEL and INV use intervals of 100 kb-7.5 Mb; TRA uses endpoints.

The analysis is repeated after removing events with endpoints within 0 or 10 kb
of an integration coordinate from the same sample. This checks sensitivity to
integration-associated calls. The CLCA background is HBV-positive HCC and therefore
differs from the integration-negative background used for the primary composition
analysis."""

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

SRC = Path(r"upstream_project")
OUT = Path(r"source_project\results")
BIN = F.BIN_SIZE
WINDOW = 100_000
SEED = 7
SPAN = (100_000, 7_500_000)
INTS = SRC / "data/interim/clca/integrations.csv"
SVS = SRC / "data/interim/clca/sv_api.csv"
CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}


def chr_norm(c) -> str:
    s = str(c)
    return s if s.startswith("chr") else f"chr{s}"


def load_integrations() -> pd.DataFrame:
    d = pd.read_csv(INTS)
    d["chrom"] = d["chrom"].astype(str).map(chr_norm)
    d["pos0"] = d["pos"].astype(np.int64)
    d = d[d["chrom"].isin(CHROMS)]
    print(f"[ints]    {len(d)} CLCA integration records, {d['sample'].nunique()} samples, "
          f"{d['chrom'].nunique()} chromosomes")
    return d


def build_events(j_excl: int, ints: pd.DataFrame) -> pd.DataFrame:
    sv = pd.read_csv(SVS, low_memory=False)
    sv["c1"] = sv["chra"].astype(str).map(chr_norm)
    sv["c2"] = sv["chrb"].astype(str).map(chr_norm)
    sv = sv[sv["c1"].isin(CHROMS) & sv["c2"].isin(CHROMS)]
    sv["p1"] = sv["posa"].to_numpy(np.int64)
    sv["p2"] = sv["posb"].to_numpy(np.int64)
    same = (sv["c1"] == sv["c2"]).to_numpy()
    sa = sv["strand_a"].astype(str).to_numpy()
    sb = sv["strand_b"].astype(str).to_numpy()
    span = np.abs(sv["p2"].to_numpy(np.int64) - sv["p1"].to_numpy(np.int64))
    cls = np.where(same,
                   np.where((sa == "-") & (sb == "+"), "DUP",
                            np.where((sa == "+") & (sb == "-"), "DEL",
                                     np.where(sa == sb, "INV", "OTHER"))),
                   "TRA")
    sv["svclass"] = cls
    ok = (cls == "TRA") | (same & (span > SPAN[0]) & (span <= SPAN[1]))
    sv = sv[ok].copy()

    # junction exclusion: an endpoint within J of an integration coordinate of the SAME sample
    # on the SAME chromosome (coordinates are per-chromosome, so the chrom key is mandatory)
    if j_excl > 0:
        per_samp_chrom = {(str(s), str(c)): np.sort(h["pos0"].to_numpy(np.int64))
                          for (s, c), h in ints.groupby(["sample", "chrom"])}

        def _near(arr: np.ndarray, x: int, jj: int) -> bool:
            d = np.searchsorted(arr, x)
            return (d < len(arr) and arr[d] - x <= jj) or \
                   (d > 0 and x - arr[d - 1] <= jj)

        keep = np.ones(len(sv), dtype=bool)
        c1n = sv["c1"].to_numpy()
        c2n = sv["c2"].to_numpy()
        sids = sv["caseid"].astype(str).to_numpy()
        p1n = sv["p1"].to_numpy(np.int64)
        p2n = sv["p2"].to_numpy(np.int64)
        for i in range(len(sv)):
            a1 = per_samp_chrom.get((sids[i], c1n[i]))
            if a1 is not None and _near(a1, int(p1n[i]), j_excl):
                keep[i] = False
                continue
            if c1n[i] == c2n[i]:
                a2 = per_samp_chrom.get((sids[i], c2n[i]))
                if a2 is not None and _near(a2, int(p2n[i]), j_excl):
                    keep[i] = False
        sv = sv[keep]

    rows = []
    for _, r in sv.iterrows():
        if r["svclass"] == "TRA":
            rows.append({"aliquot": r["caseid"], "svclass": "TRA", "chrom": r["c1"],
                         "left_pos0": r["p1"], "right_pos0": r["p1"]})
            rows.append({"aliquot": r["caseid"], "svclass": "TRA", "chrom": r["c2"],
                         "left_pos0": r["p2"], "right_pos0": r["p2"]})
        else:
            rows.append({"aliquot": r["caseid"], "svclass": r["svclass"], "chrom": r["c1"],
                         "left_pos0": min(r["p1"], r["p2"]),
                         "right_pos0": max(r["p1"], r["p2"])})
    ev = pd.DataFrame(rows)
    ev["unique_td_id"] = (ev["aliquot"].astype(str) + ":" + ev["svclass"] + ":" +
                          ev["chrom"] + ":" + ev["left_pos0"].astype(str) + "-" +
                          ev["right_pos0"].astype(str))
    ev = ev.drop_duplicates("unique_td_id")
    print(f"[sv]      J={j_excl // 1000:>2} kb -> {len(ev):,} events "
          f"{ev.svclass.value_counts().to_dict()}")
    return ev


def to_index(ev: pd.DataFrame):
    codes = pd.factorize(ev["unique_td_id"])[0].astype(np.int64)
    ev = ev.assign(_eid=codes)
    return {str(c): F.ChromEvents(g["left_pos0"].to_numpy(np.int64),
                                  g["right_pos0"].to_numpy(np.int64),
                                  g["_eid"].to_numpy(np.int64))
            for c, g in ev.groupby("chrom")}, len(ev)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    a = ap.parse_args()
    t0 = time.time()

    feat = J.build_joint_features()
    mask = A.masked_bins(A.amp_peaks(), WINDOW)
    base = F.Nuller()
    nul = J.JointNuller(feat, base, mask)

    ints = load_integrations()
    in_mask = np.array([(c, int(p) // BIN * BIN) in mask
                        for c, p in zip(ints["chrom"], ints["pos0"])])
    sub = ints[~in_mask].reset_index(drop=True)
    k = pd.MultiIndex.from_arrays([sub["chrom"].to_numpy(str),
                                   (sub["pos0"].to_numpy(np.int64) // BIN) * BIN])
    sub = sub.assign(stratum=nul.lookup.reindex(k).to_numpy())
    sub = sub[sub["stratum"].notna()].reset_index(drop=True)
    sub["stratum"] = sub["stratum"].astype(int)
    sub["within_bin_offset"] = sub["pos0"] - (sub["pos0"] // BIN) * BIN
    print(f"[coords]  {len(sub)} eligible integration records "
          f"({ints['sample'].nunique()} samples)")

    rows = []
    for j_excl in (0, 10_000):
        ev = build_events(j_excl, ints)
        for cls in ("DUP", "DEL", "INV", "TRA"):
            idx, n_ev = to_index(ev[ev.svclass == cls])
            labels = F.leader_clusters(sub["chrom"].to_numpy(str),
                                       sub["pos0"].to_numpy(np.int64), WINDOW)
            rng = np.random.default_rng(SEED)
            o, nl, _, so, sn = F.est_block(idx, nul, rng, sub, WINDOW, a.n_perm,
                                           labels, unit="uniqSV")
            for est, res in (("bc_uniqSV_blockNull", F.summarize(o, nl)),
                             ("f_clusterSignTest", F.summarize(so, sn))):
                rows.append(dict(svclass=cls, j_excl_kb=j_excl // 1000, estimator=est,
                                 n_events=n_ev, n_coords=len(sub),
                                 n_clusters=int(len(np.unique(labels))),
                                 multiplicity=round(len(sub) / len(np.unique(labels)), 2),
                                 **res))
            print(f"  [{cls:3s} J={j_excl//1000:>2}kb] n_ev={n_ev:6d} "
                  f"bc={rows[-2]['fold']:.4f} [{rows[-2]['fold_ci_lo']:.4f},"
                  f"{rows[-2]['fold_ci_hi']:.4f}] p={rows[-2]['emp_p']:.5g} | "
                  f"sign={rows[-1]['fold']:.4f} p={rows[-1]['emp_p']:.5g}", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "clca_class_replication.csv", index=False)
    pd.set_option("display.width", 240)
    print("\n" + d[["svclass", "j_excl_kb", "estimator", "n_events", "fold",
                    "fold_ci_lo", "fold_ci_hi", "emp_p"]].to_string(index=False))
    print(f"\nwrote {OUT / 'clca_class_replication.csv'} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
