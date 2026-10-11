"""Pure adapters/statistic summaries; no project paths or IO on import."""
import numpy as np
import pandas as pd
import hashlib


class TrackedRNG:
    """Record the exact draw stream while delegating unchanged to numpy's Generator."""
    def __init__(self, seed):
        self.generator = np.random.default_rng(seed)
        self.digest = hashlib.sha256()
        self.n_calls = 0

    def integers(self, *args, **kwargs):
        values = self.generator.integers(*args, **kwargs)
        array = np.asarray(values)
        self.digest.update(str(array.dtype).encode())
        self.digest.update(str(array.shape).encode())
        self.digest.update(array.tobytes())
        self.n_calls += 1
        return values

    def hexdigest(self):
        return self.digest.hexdigest()


def make_index(geometry, constructor):
    required = {"chrom", "left_pos0", "right_pos0", "event_id"}
    if not required.issubset(geometry.columns):
        raise ValueError("geometry lacks required fields")
    g = geometry.sort_values(["chrom", "left_pos0"]).copy()
    if g.duplicated(["chrom", "event_id"]).any():
        raise ValueError("multiple same-event geometry rows on one chromosome need an explicit union index")
    if (g.left_pos0 > g.right_pos0).any():
        raise ValueError("inverted interval geometry")
    g["_eid"] = pd.factorize(g.event_id)[0].astype(np.int64)
    index = {str(c): constructor(s.left_pos0.to_numpy(np.int64), s.right_pos0.to_numpy(np.int64),
                                 s._eid.to_numpy(np.int64)) for c, s in g.groupby("chrom")}
    return index, g.event_id.nunique()


def global_observed_union(index, clusters):
    """Diagnostic only: count an event once across the observed cluster catalogue."""
    ids = [index[str(r.chrom)].ids_in(int(r.start), int(r.end))
           for r in clusters.itertuples() if str(r.chrom) in index]
    return int(np.unique(np.concatenate(ids)).size) if ids else 0


def full_summary(observed, null):
    null = np.asarray(null)
    expected, sd = float(null.mean()), float(null.std(ddof=1))
    low, high = np.percentile(null, [2.5, 97.5])
    upper=float((np.count_nonzero(null >= observed)+1)/(len(null)+1))
    lower=float((np.count_nonzero(null <= observed)+1)/(len(null)+1))
    return dict(observed=int(observed), expected=expected, null_sd=sd,
        fold=float(observed/expected) if expected else np.nan,
        fold_null_interval_low=float(observed/high) if high else np.nan,
        fold_null_interval_high=float(observed/low) if low else np.nan,
        null_quantile_025=float(low), null_quantile_975=float(high),
        z=float((observed-expected)/sd) if sd else np.nan,
        emp_p=upper, emp_p_lower=lower, emp_p_two_sided_equal_tail=min(1.,2*min(upper,lower)), n_perm=len(null))


def compare_published(rebuilt, published):
    keys = ["svclass", "virus", "estimator"]
    merged = published.merge(rebuilt, on=keys, how="outer", suffixes=("_published", "_rebuilt"), indicator=True)
    rows = []
    for r in merged.to_dict("records"):
        for col in published.columns:
            if col in keys:
                continue
            a, b = r.get(col+"_published"), r.get(col+"_rebuilt")
            tolerance = 1e-15 if col == "emp_p" else 0
            ok = r["_merge"] == "both" and bool(np.isclose(a, b, atol=tolerance, rtol=0, equal_nan=True))
            rows.append({**{k:r[k] for k in keys}, "column":col, "published":a, "rebuilt":b,
                         "absolute_difference":float(abs(a-b)), "tolerance":tolerance, "passed":ok})
    return pd.DataFrame(rows)
