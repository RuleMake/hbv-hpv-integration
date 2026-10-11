"""Build and check a paired-call catalogue from PCAWG BEDPE files.

The stage command copies source inputs; build reads manifest-verified copies under
data/. The verify-sources command checks the inputs without modifying them."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "runtime_deps"))
import pandas as pd
import pyarrow
from sv_event_definitions import (ANALYSIS_CHROMS, CLASS_MAP, atomic_aliquots, chrom, digest,
                        paired_signature, source_event_id, exclusion_reasons, legacy_records,
                        original_content_signature, donor_point_pair_key)

DATA, OUT = ROOT / "data", ROOT / "outputs"
SRC = Path(r"upstream_project")
WORKSPACE = ROOT.parents[2]
CANON_NAME = "input_pooled_canonical_events.tsv"
LEGACY_COLS = ["aliquot", "svclass", "chrom", "left_pos0", "right_pos0", "unique_td_id"]
CORE = {
    CANON_NAME: SRC / "association_analysis/drafts/core_association/source_data" / CANON_NAME,
    "cohort.csv": SRC / "data/processed/cohort.csv",
    "pcawg_sample_sheet.tsv": SRC / "data/raw/viral_integration/pcawg_sample_sheet.v1.4.2016-09-14.tsv",
    "legacy_sv_class_ledger.parquet": WORKSPACE / "reference_data/results/cache/sv_class_ledger.parquet",
    "legacy_build_ledger_source.py": WORKSPACE / "reference_data/analysis/discriminating_tests.py",
    "legacy_positive_filter_source.py": WORKSPACE / "reference_data/analysis/locus_class_propensity.py",
}


def safe(path: Path) -> Path:
    path = path.resolve()
    if not path.is_relative_to(ROOT):
        raise RuntimeError(f"refusing write outside task directory: {path}")
    return path


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = safe(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def write_csv(path, data):
    path = safe(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    data.to_csv(path, index=False, encoding="utf-8")


def read_text_table(path, sep=","):
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def stage():
    if (DATA / "input_manifest.csv").exists():
        raise RuntimeError("staged manifest already exists; verify it instead of overwriting baseline")
    safe(DATA).mkdir(parents=True, exist_ok=True)
    manifest = []

    def copy_one(source, dest, role):
        dest = safe(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        before = sha(source)
        shutil.copy2(source, dest)
        copied, after = sha(dest), sha(source)
        if before != copied or before != after:
            raise RuntimeError(f"source/copy hash changed: {source}")
        manifest.append(dict(role=role, source_path=str(source.resolve()),
                             copy_path=dest.relative_to(ROOT).as_posix(),
                             bytes=dest.stat().st_size, sha256=before,
                             source_sha256_after=after, copy_sha256=copied,
                             staged_utc=datetime.now(timezone.utc).isoformat()))

    for name, source in CORE.items():
        copy_one(source, DATA / name, "core")
    canon = read_text_table(DATA / CANON_NAME, "\t")
    strings = set(canon.aliquot_metadata)
    atoms = {a for value in strings for a in atomic_aliquots(value)}
    inventory = []
    # Preserve original directory traversal ordering for dictionary last-wins selection.
    for origin in ("icgc_open", "tcga_open"):
        directory = SRC / "data/raw/pcawg_sv" / origin
        for source in directory.rglob("*.bedpe.gz"):
            aliquot = source.name.split(".")[0]
            relative = Path("bedpe") / origin / source.relative_to(directory)
            selected = aliquot in atoms
            row = dict(discovery_order=len(inventory), aliquot=aliquot, origin=origin,
                       source_path=str(source.resolve()), source_bytes=source.stat().st_size,
                       canonical_atomic=selected, legacy_string_selected=aliquot in strings,
                       staged_relative=(DATA / relative).relative_to(ROOT).as_posix() if selected else "")
            inventory.append(row)
            if selected:
                copy_one(source, DATA / relative, "raw_bedpe")
    write_csv(DATA / "raw_source_inventory.csv", inventory)
    write_csv(DATA / "input_manifest.csv", manifest)
    write_json(DATA / "stage_summary.json", dict(
        source_files_discovered=len(inventory),
        source_compressed_bytes=sum(r["source_bytes"] for r in inventory),
        copied_bedpe_files=sum(r["role"] == "raw_bedpe" for r in manifest),
        copied_bedpe_bytes=sum(r["bytes"] for r in manifest if r["role"] == "raw_bedpe"),
        canonical_rows=len(canon), metadata_unique_strings=len(strings), atomic_aliquots=len(atoms),
        all_staged_hashes_match=True))
    print((DATA / "stage_summary.json").read_text(encoding="utf-8"), flush=True)


def verify_sources(include_sources=False):
    manifest = read_text_table(DATA / "input_manifest.csv")
    records = []
    for r in manifest.to_dict("records"):
        copy = safe(ROOT / r["copy_path"])
        copy_sha = sha(copy)
        source_sha = sha(Path(r["source_path"])) if include_sources else "not_rechecked"
        records.append(dict(copy_path=r["copy_path"], bytes=copy.stat().st_size,
                            expected_sha256=r["sha256"], actual_copy_sha256=copy_sha,
                            actual_source_sha256=source_sha,
                            copy_match=copy_sha == r["sha256"],
                            source_match=source_sha == r["sha256"] if include_sources else None))
    filename = "hash_verification_sources_and_copies.csv" if include_sources else "hash_verification_copies.csv"
    write_csv(OUT / filename, records)
    if any(not r["copy_match"] or (include_sources and not r["source_match"]) for r in records):
        raise RuntimeError("input hash verification failed")
    print(f"hash verification: {len(records)} copies match" + ("; all sources match" if include_sources else ""), flush=True)
    return len(records)


def donor_maps(canon, sample_sheet):
    explicit = {}
    for r in canon.to_dict("records"):
        for a in atomic_aliquots(r["aliquot_metadata"]):
            explicit.setdefault(a, set()).add(r["donor"])
    official = {}
    for r in sample_sheet.to_dict("records"):
        if r["aliquot_id"] and r["icgc_donor_id"]:
            official.setdefault(r["aliquot_id"], set()).add(r["icgc_donor_id"])
    rows, resolved, conflicts = [], {}, []
    for a in sorted(explicit):
        c, o = explicit[a], official.get(a, set())
        reasons = []
        if len(c) != 1:
            reasons.append("canonical_ambiguous")
        if len(o) > 1:
            reasons.append("sample_sheet_ambiguous")
        if o and c != o:
            reasons.append("canonical_sample_sheet_disagree")
        donor = next(iter(o or c)) if not reasons and len(o or c) == 1 else ""
        row = dict(aliquot=a, donor=donor, canonical_donors="|".join(sorted(c)),
                   sample_sheet_donors="|".join(sorted(o)),
                   mapping_source="sample_sheet_and_canonical" if o and c == o else "canonical_only",
                   conflict="|".join(reasons))
        rows.append(row)
        resolved[a] = donor
        if reasons:
            conflicts.append(row)
    write_csv(OUT / "aliquot_donor_mapping.csv", rows)
    write_csv(OUT / "mapping_conflicts.csv", pd.DataFrame(conflicts, columns=rows[0].keys()))
    if conflicts:
        raise RuntimeError(f"{len(conflicts)} sample-to-donor mapping conflicts")
    return resolved, pd.DataFrame(rows)


def parse_raw(chosen, dmap):
    rows = []
    schema_rows = []
    required = {f"{v}{s}" for v in ("chrom", "start", "end", "strand") for s in (1, 2)} | {"sv_id", "svclass"}
    for count, inv in enumerate(chosen.to_dict("records"), 1):
        path = safe(ROOT / inv["staged_relative"])
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            fields = reader.fieldnames or []
            if not required.issubset(fields):
                raise RuntimeError(f"incomplete paired BEDPE schema: {path}")
            n = 0
            for n, source in enumerate(reader, 1):
                row = dict(source)
                row["source_content_signature"] = original_content_signature(source)
                for key in ("start1", "end1", "start2", "end2"):
                    row[key] = int(row[key])
                row.update(aliquot=inv["aliquot"], donor=dmap.get(inv["aliquot"], ""),
                           source_origin=inv["origin"], source_copy=inv["staged_relative"],
                           source_data_row=n, source_line=n + 1,
                           legacy_string_selected=inv["legacy_string_selected"] == "True")
                row["source_row_uid"] = digest([row["source_copy"], n])
                row["event_id"] = source_event_id(row["aliquot"], row["sv_id"])
                row["paired_signature"] = digest(paired_signature(row))
                row["donor_event_exactsignature"] = digest([row["donor"], paired_signature(row)])
                row["analysis_class"] = CLASS_MAP.get(str(row["svclass"]).upper(), "OTHER")
                row["chrom1_normalized"], row["chrom2_normalized"] = chrom(row["chrom1"]), chrom(row["chrom2"])
                row["span_bp"] = abs(row["start2"] - row["start1"]) if chrom(row["chrom1"]) == chrom(row["chrom2"]) else None
                row["exclusion_reasons"] = "|".join(exclusion_reasons(row))
                row["analysis_eligible"] = row["exclusion_reasons"] == ""
                rows.append(row)
            schema_rows.append(dict(aliquot=inv["aliquot"], source_copy=inv["staged_relative"],
                                    n_rows=n, schema="|".join(fields)))
        if count % 400 == 0:
            print(f"parsed {count}/{len(chosen)} aliquot files, {len(rows):,} rows", flush=True)
    write_csv(OUT / "raw_file_schema_and_rows.csv", schema_rows)
    raw = pd.DataFrame(rows)
    # Compare all original fields before adding provenance; annotation conflicts are also fatal.
    conflicts = raw.groupby("event_id").source_content_signature.nunique()
    bad = set(conflicts[conflicts > 1].index)
    write_csv(OUT / "source_event_id_conflicts.csv", raw[raw.event_id.isin(bad)])
    if bad:
        raise RuntimeError(f"{len(bad)} source event IDs have contradictory original content")
    duplicate_rows = raw[raw.duplicated("event_id", keep=False)]
    write_csv(OUT / "duplicate_source_event_rows.csv", duplicate_rows)
    return raw, raw.drop_duplicates("event_id").copy()


def make_legacy(events, label):
    rows = [r for event in events.to_dict("records") for r in legacy_records(event)]
    provenance = pd.DataFrame(rows)
    provenance["donor"] = provenance.aliquot.map(events.drop_duplicates("aliquot").set_index("aliquot").donor)
    provenance.to_parquet(safe(OUT / f"{label}_row_provenance.parquet"), index=False)
    unique = provenance.drop_duplicates("unique_td_id")[LEGACY_COLS].copy()
    unique.to_parquet(safe(OUT / f"{label}.parquet"), index=False)
    g = provenance.groupby("unique_td_id", sort=False).agg(
        n_source_records=("source_event_id", "size"), n_source_events=("source_event_id", "nunique"),
        source_event_ids=("source_event_id", lambda x: "|".join(sorted(set(x)))))
    collisions = g[g.n_source_records > 1].reset_index().merge(unique, on="unique_td_id", how="left")
    write_csv(OUT / f"{label}_collapsed_records.csv", collisions)
    return unique, provenance, collisions


def compare_legacy(rebuilt):
    stored = pd.read_parquet(DATA / "legacy_sv_class_ledger.parquet")
    expected = stored[LEGACY_COLS].astype(str).sort_values("unique_td_id").reset_index(drop=True)
    actual = rebuilt[LEGACY_COLS].astype(str).sort_values("unique_td_id").reset_index(drop=True)
    expected_keys, actual_keys = set(expected.unique_td_id), set(actual.unique_td_id)
    missing, added = expected_keys - actual_keys, actual_keys - expected_keys
    write_csv(OUT / "legacy_cached_only_rows.csv", stored[stored.unique_td_id.isin(missing)])
    write_csv(OUT / "legacy_rebuilt_only_rows.csv", rebuilt[rebuilt.unique_td_id.isin(added)])
    joined = expected.merge(actual, on="unique_td_id", suffixes=("_cached", "_rebuilt"))
    bad = pd.Series(False, index=joined.index)
    for col in LEGACY_COLS[:-1]:
        bad |= joined[col + "_cached"] != joined[col + "_rebuilt"]
    write_csv(OUT / "legacy_value_differences.csv", joined[bad])
    result = dict(cache_rows=len(stored), rebuilt_rows=len(rebuilt),
                  cached_only_rows=len(missing), rebuilt_only_rows=len(added),
                  differing_values_rows=int(bad.sum()), exact_rowset_and_values=expected.equals(actual))
    write_json(OUT / "legacy_exact_comparison.json", result)
    if not result["exact_rowset_and_values"]:
        raise RuntimeError("legacy reconstruction does not exactly match submitted cache")
    return result


def split_breakends(events):
    frames = []
    for side in (1, 2):
        d = events[["event_id", "sv_id", "aliquot", "donor", "analysis_class", "analysis_eligible",
                    "exclusion_reasons", "donor_event_exactsignature", "source_copy", "source_data_row"]].copy()
        d["side"] = side
        d["breakend_id"] = d.event_id + f"::end{side}"
        for key in ("chrom", "start", "end", "strand"):
            d[key] = events[key + str(side)].to_numpy()
        d["chrom_normalized"] = events[f"chrom{side}_normalized"].to_numpy()
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def reconcile_canonical(canon, events):
    """Join known donor/chrom/start/strand pairs; canonical has no end1/end2."""
    dup = events[(events.analysis_class == "DUP") & (events.chrom1_normalized == events.chrom2_normalized) &
                 (events.span_bp > 100_000) & (events.span_bp <= 7_500_000)].copy()
    dup["point_pair_key"] = [donor_point_pair_key(r["donor"], "DUP",
        (r["chrom1"], r["start1"], r["strand1"]), (r["chrom2"], r["start2"], r["strand2"]))
        for r in dup.to_dict("records")]
    canon = canon.copy()
    canon["point_pair_key"] = [donor_point_pair_key(r["donor"], "DUP",
        (r["chrom"], r["left_pos0"], r["left_strand"]), (r["chrom"], r["right_pos0"], r["right_strand"]))
        for r in canon.to_dict("records")]
    g = dup.groupby("point_pair_key", sort=False).agg(
        donor=("donor", "first"), n_sample_calls=("event_id", "size"),
        n_aliquots=("aliquot", "nunique"), n_full_pair_signatures=("donor_event_exactsignature", "nunique"),
        strict_analysis_eligible=("analysis_eligible", "all"),
        aliquots=("aliquot", lambda x: "|".join(sorted(set(x)))),
        sample_event_ids=("event_id", lambda x: "|".join(sorted(set(x)))),
        donor_event_exactsignatures=("donor_event_exactsignature", lambda x: "|".join(sorted(set(x)))))
    x = canon.merge(g.reset_index(), on="point_pair_key", how="outer", suffixes=("_canonical", "_raw"), indicator=True)
    x["match_status"] = x.pop("_merge").astype(str)
    x.to_csv(safe(OUT / "canonical_DUP_point_pair_crosswalk.csv.gz"), index=False, compression="gzip")
    write_csv(OUT / "canonical_DUP_unmatched.csv", x[x.match_status != "both"])
    strict_excluded = x[(x.match_status == "both") & ~x.strict_analysis_eligible.fillna(False).astype(bool)]
    write_csv(OUT / "canonical_DUP_excluded_by_strict_chromosomes.csv", strict_excluded)
    source_keys, canonical_keys = set(dup.point_pair_key), set(canon.point_pair_key)
    strict_keys = set(dup.loc[dup.analysis_eligible, "point_pair_key"])
    result = dict(
        matching_definition="donor + unordered paired chrom/start0/strand; canonical has no endpoint end fields",
        all_raw_breakend_intervals_one_bp=bool(((events.end1-events.start1)==1).all() and ((events.end2-events.start2)==1).all()),
        canonical_events=len(canon), canonical_unique_point_pair_keys=canon.point_pair_key.nunique(),
        raw_large_DUP_sample_calls=len(dup), raw_unique_donor_point_pairs=len(source_keys),
        canonical_without_raw_pair=len(canonical_keys-source_keys), raw_pair_without_canonical=len(source_keys-canonical_keys),
        all_chromosome_pair_sets_identical=canonical_keys == source_keys,
        strict_chr1_22_X_donor_point_pairs=len(strict_keys),
        canonical_excluded_from_strict_set=len(canonical_keys-strict_keys),
        strict_set_not_in_canonical=len(strict_keys-canonical_keys),
        excluded_chromosomes=sorted(strict_excluded.chrom.unique().tolist()),
        raw_signature_ambiguities=int(g.n_full_pair_signatures.gt(1).sum()),
        metadata_ID_order_used_to_infer_pairing=False)
    write_json(OUT / "canonical_DUP_reconciliation.json", result)
    if not result["all_chromosome_pair_sets_identical"] or result["raw_signature_ambiguities"]:
        raise RuntimeError("canonical DUP pair reconciliation has unresolved differences")
    return result


def counts(frame, class_col="analysis_class", unit="sample_event"):
    rows = []
    for cls in ["ALL", "DEL", "DUP", "INV", "TRA"]:
        d = frame if cls == "ALL" else frame[frame[class_col] == cls]
        if unit == "legacy_mixed_record":
            nb = int((d[class_col] != "TRA").sum() * 2 + (d[class_col] == "TRA").sum())
            n = len(d)
        elif unit == "donor_exactsignature":
            n = d.donor_event_exactsignature.nunique()
            nb = 2 * n
        else:
            n, nb = len(d), 2 * len(d)
        rows.append(dict(svclass=cls, unit=unit, n_units=int(n), n_breakend_observations=int(nb),
                         n_aliquots=int(d.aliquot.nunique()), n_donors=int(d.loc[d.donor != "", "donor"].nunique())))
    return rows


def build():
    safe(OUT).mkdir(parents=True, exist_ok=True)
    verified_n = verify_sources()
    canon = read_text_table(DATA / CANON_NAME, "\t")
    cohort = read_text_table(DATA / "cohort.csv")
    sample_sheet = read_text_table(DATA / "pcawg_sample_sheet.tsv", "\t")
    inventory = read_text_table(DATA / "raw_source_inventory.csv")
    dmap, mapping = donor_maps(canon, sample_sheet)
    positive_donors = set(cohort.donor)
    positive_aliquots = {str(p).replace("\\", "/").split("/")[-1].split(".")[0] for p in cohort.sv_path if p}
    chosen = inventory[inventory.canonical_atomic == "True"].drop_duplicates("aliquot", keep="last").copy()
    duplicates = inventory[inventory.aliquot.isin(chosen.aliquot) & inventory.duplicated("aliquot", keep=False)]
    write_csv(OUT / "duplicate_raw_source_paths.csv", duplicates)
    # Every source path per aliquot is retained in inventory. Distinct uncompressed content is not silently chosen.
    for a, g in duplicates.groupby("aliquot"):
        content_hashes = {hashlib.sha256(gzip.decompress((ROOT / p).read_bytes())).hexdigest()
                          for p in g.staged_relative}
        if len(content_hashes) > 1:
            raise RuntimeError(f"multiple nonidentical raw files for aliquot {a}; resolve explicitly")
    strings = set(canon.aliquot_metadata)
    mapping["legacy_string_selected"] = mapping.aliquot.isin(strings)
    mapping["raw_file_available"] = mapping.aliquot.isin(chosen.aliquot)
    mapping["integration_positive_donor"] = mapping.donor.isin(positive_donors)
    mapping["legacy_positive_aliquot"] = mapping.aliquot.isin(positive_aliquots)
    write_csv(OUT / "aliquot_selection_audit.csv", mapping)
    raw, events = parse_raw(chosen, dmap)
    events["integration_positive_donor"] = events.donor.isin(positive_donors)
    events["legacy_positive_aliquot"] = events.aliquot.isin(positive_aliquots)
    events.to_parquet(safe(OUT / "normalized_events_all.parquet"), index=False)
    raw.to_parquet(safe(OUT / "source_rows_all.parquet"), index=False)
    selected_legacy = events[events.legacy_string_selected]
    legacy, legacy_prov, legacy_collisions = make_legacy(selected_legacy, "legacy_exact_rebuilt")
    comparison = compare_legacy(legacy)
    expanded, expanded_prov, expanded_collisions = make_legacy(events, "legacy_expanded_selection")
    eligible = events[events.analysis_eligible].copy()
    canonical_reconciliation = reconcile_canonical(canon, events)
    eligible.to_parquet(safe(OUT / "normalized_analysis_sample_events.parquet"), index=False)
    breakends = split_breakends(events)
    breakends.to_parquet(safe(OUT / "normalized_breakends_all.parquet"), index=False)
    breakends[breakends.analysis_eligible].to_parquet(safe(OUT / "normalized_analysis_breakends.parquet"), index=False)
    signature_groups = eligible.groupby("donor_event_exactsignature", sort=False).agg(
        donor=("donor", "first"), svclass=("analysis_class", "first"),
        n_sample_events=("event_id", "size"), n_aliquots=("aliquot", "nunique"),
        aliquots=("aliquot", lambda x: "|".join(sorted(set(x)))),
        event_ids=("event_id", lambda x: "|".join(sorted(set(x)))),
        source_sv_ids=("sv_id", lambda x: "|".join(sorted(set(x)))))
    signature_groups = signature_groups.reset_index()
    signature_groups.to_parquet(safe(OUT / "donor_event_exactsignature_groups.parquet"), index=False)
    write_csv(OUT / "donor_exactsignature_multiple_calls.csv", signature_groups[signature_groups.n_sample_events > 1])
    write_csv(OUT / "within_aliquot_same_pair_different_source_ids.csv",
              eligible[eligible.duplicated(["aliquot", "paired_signature"], keep=False)])
    reason_rows = []
    for reason in sorted({v for r in events.exclusion_reasons for v in r.split("|") if v}):
        d = events[events.exclusion_reasons.str.split("|").map(lambda x: reason in x)]
        reason_rows.append(dict(reason=reason, n_events=len(d), n_aliquots=d.aliquot.nunique(), n_donors=d.donor.nunique()))
    write_csv(OUT / "normalization_exclusion_reasons.csv", reason_rows)
    write_csv(OUT / "normalization_exclusion_combinations.csv", events.groupby("exclusion_reasons", dropna=False).agg(
        n_events=("event_id", "size"), n_aliquots=("aliquot", "nunique"), n_donors=("donor", "nunique")).reset_index())
    write_csv(OUT / "bedpe_interval_widths.csv", pd.concat([
        pd.DataFrame(dict(side=side, width=events[f"end{side}"] - events[f"start{side}"]))
        for side in (1, 2)]).groupby(["side", "width"]).size().rename("n_breakends").reset_index())
    summary_rows = []
    for name, frame, col, unit in [
        ("legacy_exact", legacy, "svclass", "legacy_mixed_record"),
        ("legacy_expanded", expanded, "svclass", "legacy_mixed_record"),
        ("normalized_sample_events", eligible, "analysis_class", "sample_event"),
        ("normalized_donor_exactsignature", eligible, "analysis_class", "donor_exactsignature"),
    ]:
        if "donor" not in frame:
            frame = frame.assign(donor=frame.aliquot.map(dmap))
        for rule, sub in [
            ("all", frame),
            ("legacy_aliquot_exclusion", frame[~frame.aliquot.isin(positive_aliquots)]),
            ("whole_positive_donor_exclusion", frame[~frame.donor.isin(positive_donors)]),
        ]:
            summary_rows.extend(dict(dataset=name, exclusion_rule=rule, **r) for r in counts(sub, col, unit))
    write_csv(OUT / "counts_by_unit_class_and_exclusion.csv", summary_rows)
    mapping["raw_event_count"] = mapping.aliquot.map(events.groupby("aliquot").size()).fillna(0).astype(int)
    mapping["eligible_event_count"] = mapping.aliquot.map(eligible.groupby("aliquot").size()).fillna(0).astype(int)
    mapping["legacy_record_count"] = mapping.aliquot.map(legacy.groupby("aliquot").size()).fillna(0).astype(int)
    mapping["expanded_legacy_record_count"] = mapping.aliquot.map(expanded.groupby("aliquot").size()).fillna(0).astype(int)
    mapping["selection_diagnosis"] = mapping.apply(lambda r:
        "missing_raw_file" if not r.raw_file_available else
        "omitted_by_unsplit_metadata" if not r.legacy_string_selected else
        "no_legacy_eligible_record" if r.legacy_record_count == 0 else "represented_in_legacy", axis=1)
    write_csv(OUT / "aliquot_selection_audit.csv", mapping)
    omitted = mapping[~mapping.legacy_string_selected]
    write_csv(OUT / "aliquots_omitted_by_unsplit_metadata.csv", omitted)
    added = set(expanded.unique_td_id) - set(legacy.unique_td_id)
    write_csv(OUT / "legacy_rows_added_by_metadata_expansion.csv", expanded[expanded.unique_td_id.isin(added)])
    half_tra = events[(events.analysis_class == "TRA") &
                      (events.chrom1_normalized.isin(ANALYSIS_CHROMS) ^ events.chrom2_normalized.isin(ANALYSIS_CHROMS))]
    write_csv(OUT / "paired_TRA_with_only_one_allowed_breakend.csv", half_tra)
    extra_positive = events[events.integration_positive_donor & ~events.legacy_positive_aliquot]
    write_csv(OUT / "positive_donor_extra_aliquot_events.csv", extra_positive)
    write_csv(OUT / "integration_positive_donor_inventory.csv", pd.DataFrame(dict(donor=sorted(positive_donors))).assign(
        canonical_atomic_aliquots=lambda d: d.donor.map(mapping.groupby("donor").aliquot.nunique()).fillna(0).astype(int),
        eligible_event_aliquots=lambda d: d.donor.map(eligible.groupby("donor").aliquot.nunique()).fillna(0).astype(int)))
    checks = []
    def check(name, ok, detail):
        checks.append(dict(check=name, passed=bool(ok), detail=str(detail)))
    check("legacy_exact_rowset_and_values", comparison["exact_rowset_and_values"], comparison)
    check("one_event_id_per_normalized_call", events.event_id.is_unique, len(events))
    check("two_breakends_per_source_event", len(breakends) == 2 * len(events) and breakends.groupby("event_id").size().eq(2).all(), len(breakends))
    check("unique_breakend_ids", breakends.breakend_id.is_unique, len(breakends))
    for side in (1, 2):
        restored = breakends[breakends.side == side].set_index("event_id").loc[events.event_id]
        ok = all(restored[k].astype(str).to_numpy().tolist() == events[f"{k}{side}"].astype(str).to_numpy().tolist()
                 for k in ("chrom", "start", "end", "strand"))
        check(f"endpoint_{side}_lossless_roundtrip", ok, "chrom/start/end/strand")
    check("eligible_both_ends_chr1_22_X", eligible.chrom1_normalized.isin(ANALYSIS_CHROMS).all() and eligible.chrom2_normalized.isin(ANALYSIS_CHROMS).all(), len(eligible))
    check("all_selected_aliquots_have_donor", mapping.donor.ne("").all(), len(mapping))
    check("expanded_legacy_contains_original", set(legacy.unique_td_id).issubset(set(expanded.unique_td_id)), len(added))
    check("source_row_count_preserved", len(raw) == len(events) + int(raw.duplicated("event_id").sum()), len(raw))
    write_csv(OUT / "validation_checks.csv", checks)
    if not all(r["passed"] for r in checks):
        raise RuntimeError("ledger validation check failed")
    snapshot = dict(
        verified_input_files=verified_n,
        canonical_rows=len(canon), canonical_unique_events=canon.unique_td_id.nunique(),
        canonical_metadata_strings=canon.aliquot_metadata.nunique(), canonical_atomic_aliquots=len(mapping),
        canonical_donors=canon.donor.nunique(), mapping_conflicts=0,
        raw_files_chosen=len(chosen), raw_source_rows=len(raw), normalized_all_events=len(events),
        exact_duplicate_source_rows=len(raw)-len(events), legacy_exact=comparison,
        canonical_DUP_reconciliation=canonical_reconciliation,
        expanded_legacy_rows=len(expanded), added_legacy_rows=len(added),
        omitted_atomic_aliquots=len(omitted), omitted_with_available_files=int(omitted.raw_file_available.sum()),
        omitted_with_filtered_events=int(omitted.eligible_event_count.gt(0).sum()),
        normalized_eligible_sample_events=len(eligible),
        donor_exactsignature_groups=len(signature_groups),
        donor_signatures_spanning_multiple_aliquots=int(signature_groups.n_aliquots.gt(1).sum()),
        same_donor_exactsignature_call_excess=int(len(eligible)-len(signature_groups)),
        legacy_collapsed_rows=int(len(legacy_prov)-len(legacy)),
        expanded_legacy_collapsed_rows=int(len(expanded_prov)-len(expanded)),
        legacy_collisions_involving_multiple_source_events=int(legacy_collisions.n_source_events.gt(1).sum()),
        TRA_source_events_with_one_end_outside=len(half_tra),
        positive_cohort_rows=len(cohort), positive_donors=len(positive_donors),
        positive_path_aliquots=len(positive_aliquots),
        extra_positive_donor_aliquots=int(extra_positive.aliquot.nunique()),
        extra_positive_donor_eligible_events=int(extra_positive.analysis_eligible.sum()),
        validation_passed=True,
        biological_main_statistics_recomputed=False,
        completed_utc=datetime.now(timezone.utc).isoformat(),
        environment=dict(python=sys.version, pandas=pd.__version__, pyarrow=pyarrow.__version__, platform=platform.platform()))
    write_json(OUT / "audit_summary.json", snapshot)
    print(json.dumps(snapshot, ensure_ascii=False, indent=2), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=["stage", "build", "verify-sources"])
    a = p.parse_args()
    if a.command == "stage":
        stage()
    elif a.command == "verify-sources":
        verify_sources(True)
    else:
        build()


if __name__ == "__main__":
    main()
