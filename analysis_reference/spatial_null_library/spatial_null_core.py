# Core spatial-null model: event counting, matched pools and cluster-preserving estimators.
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np
import pandas as pd
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WS = ROOT.parents[1]
DATA01 = WS / 'analysis_inputs/structural_variant_catalogue/donor_duplication_statistics/data'
sys.path.insert(0, str(WS / 'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
BIN_SIZE = 50_000
FEATURES = DATA01 / 'genome_features_50kb.parquet'
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

def est_block(idx, nul, rng, sub, win, n_perm, labels, *, unit: str):
    """Cluster-preserving null. One common offset per cluster per permutation, drawn as
    (matched draw for a random member) - (that member's position).

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

def load_coordinates():
    co = pd.read_csv(DATA01 / 'tian_coordinates.tsv', sep='\t')
    audit = pd.read_csv(DATA01 / 'tian_matching_audit.tsv', sep='\t')
    eligible = set(audit.loc[audit.matching_status == 'ELIGIBLE', 'tian_coordinate_id'])
    co = co[co.tian_coordinate_id.isin(eligible)].copy()
    assert len(co) == 7846
    return co
