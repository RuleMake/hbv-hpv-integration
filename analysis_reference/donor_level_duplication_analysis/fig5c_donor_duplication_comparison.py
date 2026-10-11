"""Compare donor-level duplication proportions between integration-detection groups.

Inputs are verified against their manifest before analysis. With Python 3.12,
run --stage reproduce to check reference summaries, then --stage statistics.
NumPy, pandas, openpyxl and pyarrow are required."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "runtime_deps"))
import numpy as np
import pandas as pd

DATA, OUT = ROOT / "data", ROOT / "outputs"
VIRUSES = ("HBV", "HPV")
HIST = {"HBV": ("Liver-HCC",), "HPV": ("Cervix-SCC", "Cervix-AdenoCA", "Head-SCC")}
WINDOW = 100_000
BIN = 50_000
BASE_SEED = 20260908
N_BOOT, N_PERM = 50_000, 99_999


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(df, name):
    # Seventeen significant digits retain binary64 inputs across CSV round trips.
    df.to_csv(OUT / name, index=False, float_format="%.17g", na_rep="")


def verify_inputs():
    manifest = pd.read_csv(ROOT / "input_manifest.csv", encoding="utf-8-sig")
    checks = []
    for row in manifest.itertuples():
        p = ROOT / row.relative_path
        if not p.resolve().is_relative_to(ROOT):
            raise ValueError("Input path escapes this correction directory")
        ok = p.stat().st_size == row.bytes and sha(p) == row.sha256
        checks.append({"relative_path": row.relative_path, "passed": ok})
    write_csv(pd.DataFrame(checks), "input_hash_checks.csv")
    if not all(r["passed"] for r in checks):
        raise ValueError("Input checksum failed")


def norm_chrom(series):
    return series.astype(str).map(lambda s: s if s.startswith("chr") else "chr" + s)


def nearest_distance(positions, references):
    positions = np.asarray(positions, dtype=np.int64)
    references = np.sort(np.unique(np.asarray(references, dtype=np.int64)))
    if not len(references):
        return np.full(len(positions), np.inf)
    ix = np.searchsorted(references, positions)
    left = np.where(ix > 0, np.abs(positions - references[np.clip(ix - 1, 0, len(references)-1)]), np.inf)
    right = np.where(ix < len(references), np.abs(positions - references[np.clip(ix, 0, len(references)-1)]), np.inf)
    return np.minimum(left, right)


def score_events(events, coordinates, window=WINDOW):
    """One Boolean outcome per unique event, OR across its two endpoints."""
    out = events.copy()
    dleft = np.full(len(out), np.inf)
    dright = np.full(len(out), np.inf)
    for chrom, ix in out.groupby("chrom", sort=False).groups.items():
        ref = coordinates.loc[coordinates.chrom == chrom, "pos0"].to_numpy(np.int64)
        ii = np.asarray(ix, dtype=int)
        dleft[ii] = nearest_distance(out.loc[ii, "left_pos0"], ref)
        dright[ii] = nearest_distance(out.loc[ii, "right_pos0"], ref)
    out["nearest_endpoint_distance_bp"] = np.minimum(dleft, dright)
    out["left_hit"] = dleft <= window
    out["right_hit"] = dright <= window
    out["near"] = out.left_hit | out.right_hit
    return out


def load_coordinates():
    co = pd.read_csv(DATA / "tian_coordinates.tsv", sep="\t")
    aud = pd.read_csv(DATA / "tian_matching_audit.tsv", sep="\t")
    feat = pd.read_parquet(DATA / "genome_features_50kb.parquet")
    feat["chrom"] = norm_chrom(feat.chrom)
    feat = feat[(feat.stratum >= 0) & (feat.chrom != "chrY")]
    valid_bins = set(zip(feat.chrom, feat.bin_start.astype(int)))
    mask = set()
    peaks = []
    with gzip.open(DATA / "gistic_peaks.txt.gz", "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = header.index("Wide Peak Limits")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0].startswith("Amplification Peak"):
                m = re.match(r"chr([0-9XY]+):(\d+)-(\d+)", p[idx].strip())
                if m:
                    peaks.append(("chr" + m[1], int(m[2]), int(m[3])))
    for chrom, start, end in set(peaks):
        lo = max(0, (start-WINDOW)//BIN)*BIN
        hi = (end+WINDOW)//BIN*BIN
        mask.update((chrom, b) for b in range(lo, hi+BIN, BIN))
    eligible_ids = set(aud.loc[aud.matching_status == "ELIGIBLE", "tian_coordinate_id"])
    co["chrom"] = norm_chrom(co.chrom)
    keys = list(zip(co.chrom, (co.pos0 // BIN * BIN).astype(int)))
    co["legacy_eligible"] = co.tian_coordinate_id.isin(eligible_ids)
    co["masked_amplification_peak"] = [k in mask for k in keys]
    co["valid_feature_bin"] = [k in valid_bins for k in keys]
    co["joint_eligible"] = co.legacy_eligible & ~co.masked_amplification_peak & co.valid_feature_bin
    if co.legacy_eligible.sum() != 7846:
        raise ValueError("Original ELIGIBLE count differs from 7846")
    write_csv(co, "coordinate_eligibility.csv")
    result = co[co.joint_eligible].copy()
    counts = result.virus_class.value_counts().to_dict()
    if any(counts.get(v) != n for v, n in {"HBV": 3138, "HPV": 3712, "EBV": 107}.items()):
        raise ValueError(f"Joint coordinate count mismatch: {counts}")
    return result


def donor_metadata():
    h = pd.read_excel(DATA / "histology.xlsx")
    first = h.drop_duplicates("icgc_donor_id").set_index("icgc_donor_id")
    rows = []
    for donor, g in h.groupby("icgc_donor_id", sort=True):
        labels = sorted(g.histology_abbreviation.dropna().astype(str).unique())
        projects = sorted(g.project_code.dropna().astype(str).unique())
        rows.append({"donor": donor,
                     "histology_legacy_first_row": first.loc[donor, "histology_abbreviation"],
                     "histology_labels": "|".join(labels),
                     "histology_corrected": labels[0] if len(labels) == 1 else np.nan,
                     "histology_status": "consistent" if len(labels) == 1 else "conflicting" if labels else "missing",
                     "project": "|".join(projects)})
    meta = pd.DataFrame(rows).set_index("donor")
    write_csv(meta[meta.histology_status != "consistent"].reset_index(), "histology_issues.csv")
    cohort = pd.read_csv(DATA / "integration_cohort.csv")
    sets = {v: set(cohort.loc[cohort.virus_class == v, "donor"]) for v in VIRUSES}
    return meta, sets, set(cohort.donor)


def load_events():
    ev = pd.read_csv(DATA / "canonical_dup.tsv", sep="\t")
    key_fields = ["donor", "chrom", "left_pos0", "right_pos0", "left_strand", "right_strand"]
    for eid, g in ev[ev.unique_td_id.duplicated(False)].groupby("unique_td_id"):
        if len(g[key_fields].drop_duplicates()) != 1:
            raise ValueError(f"Conflicting event ID {eid}")
    ev = ev.drop_duplicates("unique_td_id").reset_index(drop=True)
    ev["chrom"] = norm_chrom(ev.chrom)
    if len(ev) != 31201 or ev.donor.nunique() != 1902:
        raise ValueError("Canonical event/donor counts differ")
    if not ((ev.left_strand == "-") & (ev.right_strand == "+")).all():
        raise ValueError("Non-DUP orientation in canonical ledger")
    return ev


def legacy_bootstrap(a, b, rng):
    # Exact historical order: two independent integer draws each replicate.
    ratios = np.empty(5000)
    for i in range(5000):
        ratios[i] = a[rng.integers(0, len(a), len(a))].mean() / b[rng.integers(0, len(b), len(b))].mean()
    return np.percentile(ratios, [2.5, 97.5])


def prepare():
    verify_inputs()
    co = load_coordinates()
    ev = load_events()
    meta, positives, all_positive = donor_metadata()
    roster = sorted(set(meta.index) | set(ev.donor) | all_positive)
    donor_tables, event_tables, checks, metadata_changes = [], [], [], []
    old = pd.read_csv(ROOT / "provenance/donor_proportion_reference.csv").set_index("virus")
    legacy_rng = np.random.default_rng(7)
    groups = []
    for v in VIRUSES:
        scored = score_events(ev, co[co.virus_class == v])
        scored["virus_coordinate_set"] = v
        event_tables.append(scored[["unique_td_id", "donor", "chrom", "left_pos0", "right_pos0", "virus_coordinate_set", "left_hit", "right_hit", "near", "nearest_endpoint_distance_bp"]])
        agg = scored.groupby("donor").near.agg(["sum", "count"]).rename(columns={"sum": "near_DUP", "count": "total_DUP"})
        d = pd.DataFrame(index=pd.Index(roster, name="donor")).join(meta).join(agg)
        d["virus_coordinate_set"] = v
        d["n_coordinates"] = int((co.virus_class == v).sum())
        d["window_bp"] = WINDOW
        d["integration_detected_for_virus"] = d.index.isin(positives[v])
        d["any_integration_detected_in_cohort"] = d.index.isin(all_positive)
        has_denom = d.total_DUP.notna() & (d.total_DUP > 0)
        d["fraction"] = d.near_DUP / d.total_DUP
        d["denominator_status"] = np.where(has_denom, "has_positive_canonical_DUP_denominator", "no_ledger_DUP_record_true_zero_vs_missing_unresolved")
        d["true_zero_fraction"] = has_denom & (d.near_DUP == 0)
        for name in ("legacy_first_row", "corrected"):
            in_tissue = d["histology_" + name].isin(HIST[v])
            d["included_" + name] = has_denom & in_tissue
            labels = np.where(d.integration_detected_for_virus, "corresponding_virus_integration_detected", "not_listed_as_corresponding_virus_integration_detected")
            d["group_" + name] = np.where(in_tissue, labels, "outside_target_histology_or_unresolved")
            a = d.loc[d["included_"+name] & d.integration_detected_for_virus, "fraction"]
            b = d.loc[d["included_"+name] & ~d.integration_detected_for_virus, "fraction"]
            for label, vals in (("detected", a), ("comparison", b)):
                groups.append({"virus": v, "coordinate_set": "joint_eligible_amplicon_masked", "n_coordinates": int((co.virus_class == v).sum()), "metadata_rule": name, "group": label, "n_donors": len(vals), "donor_mean_fraction": vals.mean(), "n_zero_fractions": int((vals==0).sum())})
            if name == "legacy_first_row":
                other = d.loc[has_denom & ~in_tissue & ~d.any_integration_detected_in_cohort & d.histology_legacy_first_row.notna(), "fraction"]
                lo, hi = legacy_bootstrap(a.to_numpy(), b.to_numpy(), legacy_rng)
                now = {"n_coords": int((co.virus_class == v).sum()), "n_exposed": len(a), "n_unexposed_same_tissue": len(b), "n_other_tissue": len(other), "rate_exposed": round(a.mean(), 4), "rate_unexposed_same_tissue": round(b.mean(), 4), "rate_other_tissue": round(other.mean(), 4), "rate_ratio_virus": round(a.mean()/b.mean(), 3), "rr_ci_lo": round(lo, 3), "rr_ci_hi": round(hi, 3)}
                for metric, value in now.items():
                    expected = old.loc[v, metric]
                    checks.append({"virus": v, "coordinate_set": "joint_eligible_amplicon_masked", "metric": metric, "legacy_value": expected, "reproduced_value": value, "difference": value-expected, "passed": bool(np.isclose(value, expected, rtol=0, atol=1e-12))})
        changed = d[d.included_legacy_first_row != d.included_corrected].copy()
        metadata_changes.append(changed.reset_index())
        donor_tables.append(d.reset_index())
    donors = pd.concat(donor_tables, ignore_index=True)
    write_csv(donors, "donor_scores.csv")
    write_csv(pd.concat(event_tables, ignore_index=True), "event_scores.csv")
    write_csv(pd.DataFrame(checks), "legacy_reproduction.csv")
    write_csv(pd.DataFrame(groups), "group_summaries.csv")
    write_csv(pd.concat(metadata_changes, ignore_index=True), "metadata_inclusion_changes.csv")
    write_csv(donors[~donors.included_corrected], "not_in_corrected_target.csv")
    included = donors[donors.included_corrected]
    support = included.groupby(["virus_coordinate_set", "histology_corrected", "project", "integration_detected_for_virus"], dropna=False).size().reset_index(name="n_donors")
    write_csv(support, "histology_project_support.csv")
    aliquots = ev.groupby("donor").aliquot_metadata.agg(lambda x: "|".join(sorted(set("|".join(x.astype(str)).split("|")))))
    write_csv(aliquots.rename("aliquot_metadata_values").reset_index(), "donor_aliquot_audit.csv")
    gate = {"passed": all(r["passed"] for r in checks), "protocol_sha256": sha(ROOT / "retrospective_correction_protocol.md"), "script_sha256": sha(__file__), "input_manifest_sha256": sha(ROOT / "input_manifest.csv"), "n_checks": len(checks)}
    (OUT / "reproduction_gate.json").write_text(json.dumps(gate, indent=2), encoding="utf-8")
    print(pd.DataFrame(groups).to_string(index=False))
    print("REFERENCE CHECK", gate["passed"])
    if not gate["passed"]:
        print(pd.DataFrame(checks).query("not passed").to_string(index=False))
        raise ValueError("Reference check failed; statistics were not executed")


def weighted_point(groups, weights):
    means_a = np.array([a.mean() for a, b in groups])
    means_b = np.array([b.mean() for a, b in groups])
    ma, mb = weights @ means_a, weights @ means_b
    variance = sum(w*w * (a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b)) for w, (a, b) in zip(weights, groups))
    diff = ma-mb
    t = diff / np.sqrt(variance) if variance > 0 else (0.0 if diff == 0 else np.copysign(np.inf, diff))
    return ma, mb, diff, ma/mb if mb != 0 else np.nan, t


def bootstrap_effects(groups, weights, rng, n_boot=N_BOOT):
    am, bm = np.zeros(n_boot), np.zeros(n_boot)
    for w, (a, b) in zip(weights, groups):
        for start in range(0, n_boot, 1024):
            stop = min(start+1024, n_boot)
            am[start:stop] += w*a[rng.integers(len(a), size=(stop-start, len(a)))].mean(axis=1)
            bm[start:stop] += w*b[rng.integers(len(b), size=(stop-start, len(b)))].mean(axis=1)
    diff = am-bm
    nonfinite = int((bm == 0).sum())
    ratios = np.divide(am, bm, out=np.full_like(am, np.nan), where=bm != 0)
    return np.percentile(diff, [2.5, 97.5]), (np.percentile(ratios, [2.5, 97.5]) if nonfinite == 0 else [np.nan, np.nan]), nonfinite


def permutation_test(groups, weights, rng, n_perm=N_PERM):
    obs = weighted_point(groups, weights)
    exceed_t, exceed_diff = 0, 0
    pools = [np.concatenate([a, b]) for a, b in groups]
    sizes = [len(a) for a, b in groups]
    for _ in range(n_perm):
        shuffled = []
        for pool, n in zip(pools, sizes):
            v = rng.permutation(pool)
            shuffled.append((v[:n], v[n:]))
        perm = weighted_point(shuffled, weights)
        exceed_t += abs(perm[4]) >= abs(obs[4])-1e-14
        exceed_diff += abs(perm[2]) >= abs(obs[2])-1e-14
    return (exceed_t+1)/(n_perm+1), (exceed_diff+1)/(n_perm+1)


def holm(values):
    p = np.asarray(values)
    order = np.argsort(p)
    adjusted = np.maximum.accumulate((len(p)-np.arange(len(p)))*p[order])
    out = np.empty(len(p))
    out[order] = np.minimum(adjusted, 1)
    return out


def run_analysis(label, virus, groups, weights, index, histologies, n_coords, interpretation):
    ma, mb, diff, ratio, t = weighted_point(groups, weights)
    rng_b = np.random.default_rng(np.random.SeedSequence([BASE_SEED, index, 1]))
    rng_p = np.random.default_rng(np.random.SeedSequence([BASE_SEED, index, 2]))
    ci_d, ci_r, bad = bootstrap_effects(groups, weights, rng_b)
    p_t, p_diff = permutation_test(groups, weights, rng_p)
    return {"analysis": label, "virus": virus, "coordinate_set": "joint_eligible_amplicon_masked", "n_coordinates": n_coords, "window_bp": WINDOW, "histologies": "|".join(histologies), "n_detected": sum(len(a) for a,b in groups), "n_comparison": sum(len(b) for a,b in groups), "mean_detected": ma, "mean_comparison": mb, "mean_difference": diff, "difference_percentage_points": 100*diff, "diff_ci_lo": ci_d[0], "diff_ci_hi": ci_d[1], "mean_ratio": ratio, "ratio_ci_lo": ci_r[0], "ratio_ci_hi": ci_r[1], "welch_statistic": t, "p_two_sided_studentized_permutation": p_t, "p_two_sided_raw_difference_permutation": p_diff, "p_holm_two_pooled_tests": np.nan, "n_bootstrap": N_BOOT, "n_permutations": N_PERM, "bootstrap_zero_comparison_means": bad, "seed_base": BASE_SEED, "arm_seed_index": index, "stratum_weights": "|".join(f"{w:.12g}" for w in weights), "old_binary_MDE_status": "withdrawn_incompatible_binary_sampling_model", "interpretation": interpretation}


def statistics():
    verify_inputs()
    gate = json.loads((OUT / "reproduction_gate.json").read_text(encoding="utf-8"))
    if not gate["passed"] or gate["protocol_sha256"] != sha(ROOT / "retrospective_correction_protocol.md"):
        raise ValueError("Reference check missing or failed, or protocol changed")
    d = pd.read_csv(OUT / "donor_scores.csv", float_precision="round_trip")
    d = d[d.included_corrected]
    rows = []
    for i, v in enumerate(VIRUSES):
        sub = d[d.virus_coordinate_set == v]
        a = sub.loc[sub.integration_detected_for_virus, "fraction"].to_numpy(float)
        b = sub.loc[~sub.integration_detected_for_virus, "fraction"].to_numpy(float)
        rows.append(run_analysis(f"{v}_pooled_donor_mean", v, [(a,b)], np.array([1.0]), i, HIST[v], int(sub.n_coordinates.iloc[0]), "same_histology_observational_association_no_causal_attribution" if v=="HBV" else "descriptive_pooled_histologies_not_histology_adjusted"))
        print("completed", rows[-1]["analysis"], flush=True)
    adjusted = holm([r["p_two_sided_studentized_permutation"] for r in rows])
    for r, p in zip(rows, adjusted):
        r["p_holm_two_pooled_tests"] = p
    hpv = d[d.virus_coordinate_set == "HPV"]
    support, strata = [], []
    for hist in HIST["HPV"]:
        h = hpv[hpv.histology_corrected == hist]
        a = h.loc[h.integration_detected_for_virus, "fraction"].to_numpy(float)
        b = h.loc[~h.integration_detected_for_virus, "fraction"].to_numpy(float)
        ok = len(a) >= 2 and len(b) >= 2
        support.append({"histology": hist, "n_detected": len(a), "n_comparison": len(b), "included_in_standardized_comparison": ok, "reason": "both_groups_at_least_two_for_variance_estimation" if ok else "no_both_group_support_or_variance_unestimable"})
        if ok:
            strata.append((hist, a, b))
    write_csv(pd.DataFrame(support), "HPV_standardization_support.csv")
    if strata:
        weights = np.array([len(a)+len(b) for h,a,b in strata], dtype=float)
        weights /= weights.sum()
        rows.append(run_analysis("HPV_common_support_histology_standardized", "HPV", [(a,b) for h,a,b in strata], weights, 2, [h for h,a,b in strata], 3712, "histology_standardized_overlap_population_sensitivity_no_causal_attribution"))
    result = pd.DataFrame(rows)
    write_csv(result, "corrected_statistics.csv")
    versions = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__}
    import pyarrow
    versions["pyarrow"] = pyarrow.__version__
    run = {"protocol_sha256": sha(ROOT / "retrospective_correction_protocol.md"), "script_sha256": sha(__file__), "input_manifest_sha256": sha(ROOT / "input_manifest.csv"), "donor_scores_sha256": sha(OUT / "donor_scores.csv"), "versions": versions, "scope": "PCAWG_DONOR_STATISTICS_joint_same_histology_only", "new_four_cell_analysis_executed": False, "CLCA_recomputed": False}
    (OUT / "statistical_run.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
    print(result[["analysis", "n_detected", "n_comparison", "mean_detected", "mean_comparison", "mean_difference", "diff_ci_lo", "diff_ci_hi", "mean_ratio", "ratio_ci_lo", "ratio_ci_hi", "p_two_sided_studentized_permutation"]].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", required=True, choices=("reproduce", "statistics"))
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    prepare() if args.stage == "reproduce" else statistics()
