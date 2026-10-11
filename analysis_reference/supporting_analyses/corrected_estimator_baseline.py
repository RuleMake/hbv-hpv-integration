"""Compare spatial counting units and coordinate-placement schemes on common inputs.

Reference observed counts are checked before estimation. Each interval event is
counted as left_hits + right_hits - both_hits, so an event with both endpoints
inside a query window is counted once.

The comparisons include coordinate-level counting with independent placements,
leader-thinned coordinates, record counts with rigid cluster placement, distinct
events within each cluster interval, cluster occupancy, and a cluster sign test.
Rigid placement preserves the spacing among coordinates in a cluster. For the
distinct-event estimator, observed and null counts use the same cluster interval.

Leader clustering uses a sorted sweep: every member lies within the specified
radius of its cluster leader, avoiding single-linkage chaining."""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SRC_ROOT = Path(r"upstream_project")
BASE = SRC_ROOT / "association_analysis" / "drafts" / "core_association" / "source_data"
FEATURES = SRC_ROOT / "data" / "interim" / "genome_features_50kb.parquet"
OUT = Path(r"source_project\results")

WINDOWS = (100_000, 1_000_000, 5_000_000)
BIN_SIZE = 50_000

# Published Fig. 3 observed coordinate-event pair counts, from
# drafts/core_association/source_data/tian_cumulative_summary.tsv.
# The reproduction arm must hit these exactly or the run aborts.
PUBLISHED_OBSERVED = {
    ("EBV", 100_000): 548, ("EBV", 1_000_000): 3843, ("EBV", 5_000_000): 14366,
    ("HBV", 100_000): 18653, ("HBV", 1_000_000): 125427, ("HBV", 5_000_000): 456116,
    ("HPV", 100_000): 25863, ("HPV", 1_000_000): 155423, ("HPV", 5_000_000): 524210,
}
PUBLISHED_FOLD = {
    ("EBV", 100_000): 1.0646, ("EBV", 1_000_000): 1.0654, ("EBV", 5_000_000): 1.0649,
    ("HBV", 100_000): 1.1589, ("HBV", 1_000_000): 1.0925, ("HBV", 5_000_000): 1.0338,
    ("HPV", 100_000): 1.4097, ("HPV", 1_000_000): 1.1926, ("HPV", 5_000_000): 1.0578,
}


# ===================================================================== event index
class ChromEvents:
    """Canonical TD events on one chromosome, indexed for fixed-length interval queries.

    `count(L, lo)` returns, for each start `lo`, the number of DISTINCT events with at least
    one endpoint in the closed interval [lo, lo+L] — the published rule
    `left_hits + right_hits - both_hits`, where `both_hits` corrects events whose two
    endpoints both land inside. Covering arrays are cached per interval length L.
    """

    def __init__(self, left: np.ndarray, right: np.ndarray, eid: np.ndarray):
        order = np.argsort(left, kind="stable")
        self.left = left[order]
        self.right = right[order]
        self.eid = eid[order]
        self.left_sorted = np.sort(self.left)
        self.right_sorted = np.sort(self.right)
        self.span = self.right - self.left
        self._cover: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        # endpoint stream for exact set-union queries (estimator b, diagnostics)
        pos = np.concatenate([self.left, self.right])
        ids = np.concatenate([self.eid, self.eid])
        o2 = np.argsort(pos, kind="stable")
        self.ep_pos = pos[o2]
        self.ep_eid = ids[o2]

    def _covering(self, L: int) -> tuple[np.ndarray, np.ndarray]:
        if L not in self._cover:
            ok = self.span <= L
            self._cover[L] = (np.sort(self.right[ok] - L), np.sort(self.left[ok]))
        return self._cover[L]

    def count(self, L: int, lo) -> np.ndarray:
        lo = np.asarray(lo, dtype=np.int64)
        hi = lo + L
        left_hits = (np.searchsorted(self.left_sorted, hi, "right")
                     - np.searchsorted(self.left_sorted, lo, "left"))
        right_hits = (np.searchsorted(self.right_sorted, hi, "right")
                      - np.searchsorted(self.right_sorted, lo, "left"))
        cover_start, cover_end = self._covering(L)
        both = (np.searchsorted(cover_start, lo, "right")
                - np.searchsorted(cover_end, lo, "left"))
        counts = left_hits + right_hits - both
        if (counts < 0).any() or (counts > len(self.left)).any():
            raise RuntimeError("event-pair count violated its bounds")
        return counts

    def ids_in(self, lo: int, hi: int) -> np.ndarray:
        i0 = np.searchsorted(self.ep_pos, lo, "left")
        i1 = np.searchsorted(self.ep_pos, hi, "right")
        return self.ep_eid[i0:i1]


def load_events() -> tuple[dict[str, ChromEvents], int]:
    ev = pd.read_csv(BASE / "input_pooled_canonical_events.tsv", sep="\t")
    ev = ev.drop_duplicates("unique_td_id")
    if (ev["left_pos0"] >= ev["right_pos0"]).any():
        raise SystemExit("canonical TD endpoints must be ordered")
    codes = pd.factorize(ev["unique_td_id"])[0].astype(np.int64)
    ev = ev.assign(_eid=codes)
    idx = {
        str(c): ChromEvents(g["left_pos0"].to_numpy(np.int64),
                            g["right_pos0"].to_numpy(np.int64),
                            g["_eid"].to_numpy(np.int64))
        for c, g in ev.groupby("chrom")
    }
    print(f"[events]  {len(ev):,} canonical TD events across {len(idx)} chromosomes")
    return idx, len(ev)


def load_coordinates() -> pd.DataFrame:
    co = pd.read_csv(BASE / "input_tian_coordinates.tsv", sep="\t")
    audit = pd.read_csv(BASE / "tian_matching_audit.tsv", sep="\t")
    eligible = set(audit.loc[audit["matching_status"] == "ELIGIBLE", "tian_coordinate_id"])
    co = co[co["tian_coordinate_id"].isin(eligible)].copy()
    counts = co["virus_class"].value_counts().to_dict()
    print(f"[coords]  {len(co):,} matched-eligible coordinates {counts}")
    if len(co) != 7846:
        raise SystemExit(f"expected 7,846 matched-eligible coordinates, got {len(co)}")
    return co


class EmptyMatchedPool(RuntimeError):
    """A matched pool contains no alternative position after excluding the source coordinate.
    
    The caller can catch this condition, omit the unscorable coordinate and record
    its exclusion before drawing spatial-null positions."""


class Nuller:
    """Exact-stratum matched pools, reproducing run_tian_matched_null's pool rule:
    same chromosome, same 18-way stratum, same within-bin offset, own position excluded."""

    def __init__(self) -> None:
        feat = pd.read_parquet(FEATURES)
        # the parquet stores bare chromosome names ("1"); prepare_feature_bins normalizes
        feat = feat.assign(chrom=feat["chrom"].astype(str).map(
            lambda c: c if c.startswith("chr") else f"chr{c}"))
        feat = feat[(feat["stratum"] >= 0) & (feat["chrom"] != "chrY")]
        self.bins: dict[tuple[str, int], np.ndarray] = {
            (str(c), int(s)): g["bin_start"].to_numpy(np.int64)
            for (c, s), g in feat.groupby(["chrom", "stratum"])
        }
        self.chrom_size = {
            str(c): int(g["bin_start"].max()) + BIN_SIZE
            for c, g in feat.groupby("chrom")
        }
        print(f"[null]    {len(feat):,} eligible 50-kb bins, {len(self.bins)} chrom x stratum pools")

    def pool(self, chrom: str, stratum: int, offset: int, pos0: int) -> np.ndarray:
        bins = self.bins.get((chrom, int(stratum)))
        if bins is None or bins.size == 0:
            raise RuntimeError(f"empty matched pool for {chrom} stratum {stratum}")
        cand = np.unique(bins + int(offset))
        out = cand[cand != int(pos0)]
        if out.size == 0:
            raise EmptyMatchedPool(
                f"no matched position for {chrom} stratum {stratum} at {pos0} "
                f"(offset {offset}); the pool holds {bins.size} bin(s)")
        return out

    def draw(self, rng, chrom, stratum, offset, pos0, n) -> np.ndarray:
        p = self.pool(chrom, stratum, offset, pos0)
        return p[rng.integers(0, p.size, size=n)]


def leader_clusters(chrom: np.ndarray, pos: np.ndarray, radius: int) -> np.ndarray:
    labels = np.empty(len(pos), dtype=np.int64)
    order = np.lexsort((pos, chrom))
    cur, lead_c, lead_p = -1, None, None
    for i in order:
        c, p = chrom[i], pos[i]
        if lead_c is None or c != lead_c or (p - lead_p) > radius:
            cur += 1
            lead_c, lead_p = c, p
        labels[i] = cur
    return labels


def summarize(obs, null: np.ndarray) -> dict:
    exp = float(null.mean())
    n = len(null)
    p = (int(np.sum(null >= obs)) + 1) / (n + 1)
    sd = float(null.std(ddof=1))
    lo, hi = np.percentile(null, [2.5, 97.5])
    return {
        "observed": int(obs), "expected": round(exp, 1),
        "fold": round(obs / exp, 4) if exp > 0 else np.nan,
        "fold_ci_lo": round(obs / hi, 4) if hi > 0 else np.nan,
        "fold_ci_hi": round(obs / lo, 4) if lo > 0 else np.nan,
        "z": round((obs - exp) / sd, 2) if sd > 0 else np.nan,
        "emp_p": p, "n_perm": n,
    }


# ===================================================================== estimators

def drop_unscorable(nul: "Nuller", sub: pd.DataFrame, *, label: str = "") -> tuple[pd.DataFrame, int]:
    """Remove and count coordinates without an eligible alternative null position.
    
    EmptyMatchedPool identifies a chromosome-stratum pool containing only the source
    coordinate. Such coordinates are excluded before scoring rather than compared
    with themselves."""
    keep = np.ones(len(sub), dtype=bool)
    for i, (c, p, st, off) in enumerate(zip(sub["chrom"], sub["pos0"],
                                            sub["stratum"], sub["within_bin_offset"])):
        try:
            nul.pool(str(c), int(st), int(off), int(p))
        except EmptyMatchedPool:
            keep[i] = False
    n_drop = int((~keep).sum())
    if n_drop:
        print(f"[unscorable] {label or 'coordinates'}: dropped {n_drop:,} of {len(sub):,} "
              f"with no matched position to draw against", flush=True)
    return sub.loc[keep].reset_index(drop=True), n_drop


def est_independent(idx, nul, rng, sub, win, n_perm):
    """Published design: per coordinate, unique TD events within +/- win; null repositions
    each coordinate independently within its exact matched stratum."""
    L = 2 * win
    obs = 0
    null = np.zeros(n_perm, dtype=np.int64)
    for c, p, st, off in zip(sub["chrom"], sub["pos0"], sub["stratum"], sub["within_bin_offset"]):
        t = idx.get(str(c))
        if t is None:
            continue
        obs += int(t.count(L, np.array([p - win]))[0])
        m = nul.draw(rng, str(c), st, off, p, n_perm)
        null += t.count(L, m - win)
    return obs, null


def est_block(idx, nul, rng, sub, win, n_perm, labels, *, unit: str):
    """Cluster-preserving null. One common offset per cluster per permutation, drawn as
    (matched draw for a random member) - (that member's position), clipped to the chromosome.

    unit='pairs'   observed identical to the published estimator; only the null changes.
    unit='uniqSV'  both sides count distinct events in the cluster's contiguous interval.
    unit='binary'  both sides count clusters containing >=1 event.

    Also returns a CALIBRATED cluster sign statistic: how many clusters exceed their own null
    mean. The reference proportion is NOT 0.5 — for skewed counts (many clusters expect well
    under one event) the null median sits below the null mean, so P(draw > mean) < 0.5 under
    H0 and a 0.5 reference manufactures spurious depletion. Instead the same statistic is
    recomputed on every null draw (against its leave-one-out mean), giving its exact null
    distribution.

    Returns (obs_total, null_total, per_cluster_frame, sign_obs, sign_null).
    """
    df = sub.copy()
    df["_cluster"] = labels
    obs_total = 0
    null_total = np.zeros(n_perm, dtype=np.int64)
    sign_obs = 0
    sign_null = np.zeros(n_perm, dtype=np.int64)
    n_overhang = 0
    n_draws = 0
    rows = []
    for cid, g in df.groupby("_cluster", sort=True):
        c = str(g["chrom"].iloc[0])
        t = idx.get(c)
        if t is None:
            continue
        p = g["pos0"].to_numpy(np.int64)
        st = g["stratum"].to_numpy()
        off = g["within_bin_offset"].to_numpy()
        size = nul.chrom_size[c]

        j = rng.integers(0, len(p), size=n_perm)
        delta = np.empty(n_perm, dtype=np.int64)
        for k in range(len(p)):
            mask = j == k
            nk = int(mask.sum())
            if nk:
                delta[mask] = nul.draw(rng, c, st[k], off[k], p[k], nk) - p[k]
        lo_b, hi_b = int(p.min()) - win, int(p.max()) + win
        # NO CLIP. An earlier version clipped delta so the window-PADDED block stayed wholly
        # on-chromosome, while the observed side (t.count(L, lo_b)) truncates freely at
        # telomeres — so observed and null scored windows of different effective length. The
        # published per-coordinate null (est_independent) applies no clip either, letting both
        # sides truncate identically. Blocks that overhang simply count fewer events, exactly
        # as an observed window at a telomere does.
        n_overhang += int(np.sum((lo_b + delta < 0) | (hi_b + delta > size)))
        n_draws += n_perm

        if unit == "pairs":
            o = 0
            nl = np.zeros(n_perm, dtype=np.int64)
            for pi in p:
                o += int(t.count(2 * win, np.array([pi - win]))[0])
                nl += t.count(2 * win, (pi - win) + delta)
        else:
            L = hi_b - lo_b
            o = int(t.count(L, np.array([lo_b]))[0])
            nl = t.count(L, lo_b + delta)
            if unit == "binary":
                o = int(o > 0)
                nl = (nl > 0).astype(np.int64)
        obs_total += o
        null_total += nl
        mean_all = float(nl.mean())
        sign_obs += int(o > mean_all)
        # leave-one-out reference so each null draw is scored the same way the observed is
        loo = (nl.sum() - nl) / max(n_perm - 1, 1)
        sign_null += (nl > loo).astype(np.int64)
        rows.append({"cluster": int(cid), "chrom": c, "n_coords": len(p),
                     "start": lo_b, "end": hi_b, "obs": int(o),
                     "null_mean": mean_all})
    if n_draws:
        print(f"        [{unit}] block overhang: "
              f"{100 * n_overhang / n_draws:.2f}% of null placements extend past a "
              f"chromosome end (observed windows truncate the same way)", flush=True)
    return obs_total, null_total, pd.DataFrame(rows), sign_obs, sign_null


def est_uniqsv_scatter(idx, nul, rng, sub, win, n_perm, labels):
    """Corrected observed (distinct events per cluster) against the PUBLISHED scattered null:
    cluster members repositioned independently, distinct events over the union of their windows.
    Mixes a corrected numerator with an uncorrected denominator, so it is a sensitivity arm
    rather than the estimator to believe."""
    df = sub.copy()
    df["_cluster"] = labels
    obs = 0
    null = np.zeros(n_perm, dtype=np.int64)
    for cid, g in df.groupby("_cluster", sort=True):
        c = str(g["chrom"].iloc[0])
        t = idx.get(c)
        if t is None:
            continue
        p = g["pos0"].to_numpy(np.int64)
        st = g["stratum"].to_numpy()
        off = g["within_bin_offset"].to_numpy()
        lo_b, hi_b = int(p.min()) - win, int(p.max()) + win
        obs += int(t.count(hi_b - lo_b, np.array([lo_b]))[0])
        m = len(p)
        if m == 1:
            x = nul.draw(rng, c, st[0], off[0], p[0], n_perm) - win
            null += t.count(2 * win, x)
        else:
            M = np.empty((m, n_perm), dtype=np.int64)
            for k in range(m):
                M[k] = nul.draw(rng, c, st[k], off[k], p[k], n_perm)
            for s in range(n_perm):
                cat = np.concatenate([t.ids_in(M[k, s] - win, M[k, s] + win) for k in range(m)])
                null[s] += 0 if cat.size == 0 else np.unique(cat).size
    return obs, null


# ===================================================================== driver
def run(n_perm: int, n_perm_scatter: int, seed: int, viruses: tuple[str, ...]) -> tuple:
    idx, n_events = load_events()
    co = load_coordinates()
    nul = Nuller()

    rows, cluster_frames = [], []
    for win in WINDOWS:
        for vc in viruses:
            sub = co[co["virus_class"] == vc].sort_values(
                ["chrom", "pos0"]).reset_index(drop=True)
            if sub.empty:
                continue
            labels = leader_clusters(sub["chrom"].to_numpy(str),
                                     sub["pos0"].to_numpy(np.int64), win)
            n_clu = int(len(np.unique(labels)))
            lead_idx = pd.Series(labels).drop_duplicates(keep="first").index
            thin = sub.loc[lead_idx].reset_index(drop=True)
            base = {"window_bp": win, "virus": vc, "n_coords": len(sub),
                    "n_clusters": n_clu, "multiplicity": round(len(sub) / n_clu, 3)}
            t0 = time.time()

            rng = np.random.default_rng(seed)
            o, nl = est_independent(idx, nul, rng, sub, win, n_perm)
            expected = PUBLISHED_OBSERVED[(vc, win)]
            if o != expected:
                raise SystemExit(
                    f"REPRODUCTION GATE FAILED: {vc} {win} bp observed {o}, "
                    f"published {expected}. Refusing to report corrections.")
            rows.append({**base, "estimator": "naive_published", "n_units": len(sub),
                         "published_fold": PUBLISHED_FOLD[(vc, win)], **summarize(o, nl)})

            rng = np.random.default_rng(seed)
            o, nl = est_independent(idx, nul, rng, thin, win, n_perm)
            rows.append({**base, "estimator": "a_dedup_coords", "n_units": len(thin),
                         "published_fold": np.nan, **summarize(o, nl)})

            rng = np.random.default_rng(seed)
            o, nl, _, _, _ = est_block(idx, nul, rng, sub, win, n_perm, labels, unit="pairs")
            rows.append({**base, "estimator": "c_clusterBlock_null", "n_units": n_clu,
                         "published_fold": np.nan, **summarize(o, nl)})

            rng = np.random.default_rng(seed)
            o, nl, per, s_obs, s_null = est_block(idx, nul, rng, sub, win, n_perm, labels,
                                                  unit="uniqSV")
            rows.append({**base, "estimator": "bc_uniqSV_blockNull", "n_units": n_clu,
                         "published_fold": np.nan, **summarize(o, nl)})
            per.insert(0, "virus", vc)
            per.insert(1, "window_bp", win)
            cluster_frames.append(per)

            # calibrated against its own null, NOT against 0.5 (see est_block docstring)
            rows.append({**base, "estimator": "f_clusterSignTest", "n_units": n_clu,
                         "published_fold": np.nan, **summarize(s_obs, s_null)})

            rng = np.random.default_rng(seed)
            o, nl, _, _, _ = est_block(idx, nul, rng, sub, win, n_perm, labels, unit="binary")
            sat = 100 * o / n_clu if n_clu else np.nan
            rows.append({**base, "estimator": "d_clusterBinary", "n_units": n_clu,
                         "published_fold": np.nan, "saturation_pct": round(sat, 1),
                         **summarize(o, nl)})

            if n_perm_scatter > 0:
                rng = np.random.default_rng(seed)
                o, nl = est_uniqsv_scatter(idx, nul, rng, sub, win, n_perm_scatter, labels)
                rows.append({**base, "estimator": "b_uniqSV_scatterNull", "n_units": n_clu,
                             "published_fold": np.nan, **summarize(o, nl)})

            print(f"  [{vc} {win // 1000} kb] {len(sub)} coords -> {n_clu} clusters "
                  f"(multiplicity {len(sub) / n_clu:.2f})  {time.time() - t0:.0f}s", flush=True)
            for r in rows[-7:]:
                print(f"      {r['estimator']:22s} obs={r['observed']:>8} "
                      f"exp={r['expected']:>10} fold={r['fold']} "
                      f"p={r['emp_p']:.5g}", flush=True)
    return pd.DataFrame(rows), pd.concat(cluster_frames, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=20000)
    ap.add_argument("--n-perm-scatter", type=int, default=0,
                    help="permutations for the b_uniqSV_scatterNull sensitivity arm (0 = skip)")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--viruses", default="HBV,HPV,EBV")
    # Default is the NO-CLIP name: the authoritative run is the one without the boundary
    # clip (see est_block). Defaulting to "corrected_estimator_baseline.csv" would silently
    # regenerate the superseded filename now sitting in results/_superseded/.
    ap.add_argument("--out", default="estimator_reference_unclipped.csv")
    a = ap.parse_args()
    res, clusters = run(a.n_perm, a.n_perm_scatter, a.seed, tuple(a.viruses.split(",")))
    OUT.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT / a.out, index=False)
    clusters.to_csv(OUT / a.out.replace(".csv", "_clusters.csv"), index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print("\n" + "=" * 140)
    print(res.to_string(index=False))
    print(f"\nwrote {OUT / a.out}")


if __name__ == "__main__":
    main()
