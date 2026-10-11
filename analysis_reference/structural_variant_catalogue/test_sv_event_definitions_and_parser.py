"""Counterexamples for event identity, endpoint preservation and denominators."""
import unittest
import csv
import gzip
import tempfile
from pathlib import Path
from unittest.mock import patch
import build_and_audit_sv_call_catalogue
from sv_event_definitions import (atomic_aliquots, paired_signature, source_event_id, original_content_signature,
                        donor_point_pair_key,
                        exclusion_reasons, legacy_records, negative_masks)


def row(**kw):
    d = dict(aliquot="a", donor="d", sv_id="SVMERGE1", event_id="a::SVMERGE1",
             chrom1="1", start1=10, end1=11, strand1="+",
             chrom2="2", start2=20, end2=21, strand2="-", svclass="TRA")
    d.update(kw)
    return d


class LedgerCounterexamples(unittest.TestCase):
    def test_atomic_metadata_is_not_a_sample(self):
        self.assertEqual(atomic_aliquots("a|b|a"), ("a", "b"))
        self.assertEqual(len(set(["a", "b", "a|b"])), 3)
        self.assertEqual(len({x for v in ["a", "b", "a|b"] for x in atomic_aliquots(v)}), 2)

    def test_shared_breakend_is_not_shared_event(self):
        a, b = row(), row(sv_id="SVMERGE2", event_id="a::SVMERGE2", chrom2="3")
        self.assertNotEqual(paired_signature(a), paired_signature(b))
        keys = {r["unique_td_id"] for e in (a, b) for r in legacy_records(e)}
        self.assertEqual(len(keys), 3)  # two true source events, four endpoint observations
        self.assertEqual(len(legacy_records(a)) + len(legacy_records(b)), 4)

    def test_local_ids_are_namespaced_by_aliquot(self):
        self.assertNotEqual(source_event_id("a", "SVMERGE1"), source_event_id("b", "SVMERGE1"))

    def test_endpoint_order_keeps_strands_attached(self):
        a = row()
        b = row(chrom1="2", start1=20, end1=21, strand1="-",
                chrom2="1", start2=10, end2=11, strand2="+")
        self.assertEqual(paired_signature(a), paired_signature(b))
        b["strand1"] = "+"
        self.assertNotEqual(paired_signature(a), paired_signature(b))

    def test_y_partner_is_not_dropped_from_paired_event(self):
        a = row(chrom2="Y")
        self.assertEqual(len(legacy_records(a)), 1)
        self.assertIn("one_or_both_breakends_outside_chr1_22_X", exclusion_reasons(a))

    def test_exact_span_boundaries(self):
        a = row(svclass="DUP", chrom2="1", start2=100010, end2=100011)
        self.assertIn("span_not_gt_100000", exclusion_reasons(a))
        a.update(start2=100011, end2=100012)
        self.assertEqual(exclusion_reasons(a), [])
        a.update(start2=7500010, end2=7500011)
        self.assertEqual(exclusion_reasons(a), [])
        a.update(start2=7500011, end2=7500012)
        self.assertIn("span_gt_7500000", exclusion_reasons(a))

    def test_donor_strategy_catches_extra_positive_aliquot(self):
        self.assertEqual(negative_masks("d", "b", {"d"}, {"a"}), (True, False))

    def test_same_donor_repeat_not_an_original_id_merge(self):
        a, b = row(), row(aliquot="b", event_id="b::SVMERGE1")
        self.assertEqual(paired_signature(a), paired_signature(b))
        self.assertNotEqual(a["event_id"], b["event_id"])

    def test_missing_source_id_rejected(self):
        with self.assertRaises(ValueError):
            source_event_id("a", "")

    def test_conflicting_source_id_has_different_content(self):
        a, b = row(), row(start2=21, end2=22)
        self.assertEqual(source_event_id(a["aliquot"], a["sv_id"]), source_event_id(b["aliquot"], b["sv_id"]))
        self.assertNotEqual(paired_signature(a), paired_signature(b))

    def test_identical_pair_does_not_hide_annotation_conflict(self):
        a, b = row(pe_support="10"), row(pe_support="12")
        self.assertEqual(paired_signature(a), paired_signature(b))
        self.assertNotEqual(original_content_signature(a), original_content_signature(b))

    def test_canonical_join_uses_donor_and_both_strands(self):
        key = donor_point_pair_key("d", "DUP", ("chr1", 10, "-"), ("chr1", 200010, "+"))
        self.assertEqual(key, donor_point_pair_key("d", "DUP", ("1", 200010, "+"), ("1", 10, "-")))
        self.assertNotEqual(key, donor_point_pair_key("other_donor", "DUP", ("1", 10, "-"), ("1", 200010, "+")))
        self.assertNotEqual(key, donor_point_pair_key("d", "DUP", ("1", 10, "+"), ("1", 200010, "-")))


class RealParserCounterexamples(unittest.TestCase):
    """Exercise the production parser with real gzipped BEDPE fixture files."""

    def setUp(self):
        self.task_root = Path(__file__).resolve().parent
        scratch_parent = self.task_root / "test_scratch"
        scratch_parent.mkdir(exist_ok=True)
        self.scratch = tempfile.TemporaryDirectory(prefix="parser_", dir=scratch_parent)
        self.scratch_path = Path(self.scratch.name).resolve()
        self.assertTrue(self.scratch_path.is_relative_to(self.task_root))
        self.out_patch = patch.object(run_audit, "OUT", self.scratch_path / "outputs")
        self.out_patch.start()

    def tearDown(self):
        self.out_patch.stop()
        # Verify the concrete recursive-cleanup target remains in this task's scratch folder.
        if not self.scratch_path.is_relative_to((self.task_root / "test_scratch").resolve()):
            raise RuntimeError("unsafe fixture cleanup target")
        self.scratch.cleanup()

    def parse(self, source_rows):
        fixture = self.scratch_path / "a.bedpe.gz"
        fields = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "sv_id",
                  "pe_support", "strand1", "strand2", "svclass", "svmethod"]
        with gzip.open(fixture, "wt", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for source in source_rows:
                writer.writerow({"pe_support": "10", "svmethod": "MERGE", **source})
        inventory = run_audit.pd.DataFrame([dict(aliquot="a", origin="fixture",
            staged_relative=fixture.relative_to(self.task_root).as_posix(), legacy_string_selected="True")])
        return run_audit.parse_raw(inventory, {"a": "d"})

    def test_parser_stops_for_same_id_same_pair_changed_annotation(self):
        for field, value in [("pe_support", "12"), ("svmethod", "DIFFERENT")]:
            with self.subTest(field=field):
                a, b = row(), row()
                b[field] = value
                with self.assertRaisesRegex(RuntimeError, "contradictory original content"):
                    self.parse([a, b])
                report = run_audit.pd.read_csv(self.scratch_path / "outputs/source_event_id_conflicts.csv")
                self.assertEqual(len(report), 2)
                self.assertEqual(report.paired_signature.nunique(), 1)
                self.assertEqual(report.source_content_signature.nunique(), 2)

    def test_parser_keeps_duplicate_provenance_without_false_content_conflict(self):
        raw, events = self.parse([row(), row()])
        self.assertEqual((len(raw), len(events)), (2, 1))
        self.assertEqual(raw.source_data_row.tolist(), [1, 2])
        self.assertEqual(raw.source_row_uid.nunique(), 2)
        self.assertEqual(raw.source_content_signature.nunique(), 1)

    def test_parser_rejects_missing_original_id(self):
        with self.assertRaises(ValueError):
            self.parse([row(sv_id="")])


if __name__ == "__main__":
    unittest.main(verbosity=2)
