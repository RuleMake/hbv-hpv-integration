"""Define event identities and interval geometry without external file access."""
from __future__ import annotations

import hashlib
import json

ANALYSIS_CHROMS = {str(i) for i in range(1, 23)} | {"X"}
CLASS_MAP = {"DUP": "DUP", "DEL": "DEL", "H2HINV": "INV", "T2TINV": "INV", "TRA": "TRA"}


def atomic_aliquots(value: str) -> tuple[str, ...]:
    return tuple(sorted({x.strip() for x in str(value).split("|") if x.strip()}))


def chrom(value: str) -> str:
    value = str(value)
    return value[3:] if value.startswith("chr") else value


def digest(parts) -> str:
    return hashlib.sha256(json.dumps(parts, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def endpoint_tuple(row: dict, side: int) -> tuple:
    return (chrom(row[f"chrom{side}"]), int(row[f"start{side}"]),
            int(row[f"end{side}"]), str(row[f"strand{side}"]))


def paired_signature(row: dict) -> tuple:
    return (str(row["svclass"]).upper(), *sorted((endpoint_tuple(row, 1), endpoint_tuple(row, 2))))


def original_content_signature(row: dict) -> str:
    return digest(sorted((str(k), str(v)) for k, v in row.items()))


def donor_point_pair_key(donor: str, svclass: str, first: tuple, second: tuple) -> str:
    """Exact donor + known chrom/start/strand pair; no invented interval ends."""
    ends = sorted(((chrom(first[0]), int(first[1]), str(first[2])),
                   (chrom(second[0]), int(second[1]), str(second[2]))))
    return digest([donor, str(svclass).upper(), ends])


def source_event_id(aliquot: str, sv_id: str) -> str:
    if not str(sv_id).strip():
        raise ValueError("missing source sv_id; cannot manufacture an event identity")
    return f"{aliquot}::{sv_id}"


def exclusion_reasons(row: dict) -> list[str]:
    reasons = []
    raw = str(row["svclass"]).upper()
    cls = CLASS_MAP.get(raw)
    if cls is None:
        reasons.append("unrepresented_svclass")
    c1, c2 = chrom(row["chrom1"]), chrom(row["chrom2"])
    if c1 not in ANALYSIS_CHROMS or c2 not in ANALYSIS_CHROMS:
        reasons.append("one_or_both_breakends_outside_chr1_22_X")
    for side in (1, 2):
        if int(row[f"start{side}"]) < 0 or int(row[f"end{side}"]) <= int(row[f"start{side}"]):
            reasons.append("invalid_bedpe_coordinate_interval")
            break
    if cls == "TRA":
        if c1 == c2:
            reasons.append("TRA_not_interchromosomal")
    elif cls is not None:
        if c1 != c2:
            reasons.append("interval_class_not_intrachromosomal")
        else:
            span = abs(int(row["start2"]) - int(row["start1"]))
            if span <= 100_000:
                reasons.append("span_not_gt_100000")
            elif span > 7_500_000:
                reasons.append("span_gt_7500000")
    if not row.get("donor"):
        reasons.append("donor_unmapped")
    return reasons


def legacy_records(row: dict) -> list[dict]:
    """Faithful representation of original build_ledger row logic (no old imports)."""
    cls = CLASS_MAP.get(str(row["svclass"]).upper())
    c1, c2 = chrom(row["chrom1"]), chrom(row["chrom2"])
    s1, s2 = int(row["start1"]), int(row["start2"])
    coords = []
    if cls == "TRA" and c1 == c2:
        return []
    if cls == "TRA":
        coords = [(c1, s1, s1, 1), (c2, s2, s2, 2)]
    elif cls is not None and c1 == c2:
        lo, hi = sorted((s1, s2))
        if 100_000 < hi - lo <= 7_500_000:
            coords = [(c1, lo, hi, 0)]
    out = []
    for c, lo, hi, side in coords:
        if c not in ANALYSIS_CHROMS:
            continue
        key = f"{row['aliquot']}:{cls}:chr{c}:{lo}-{hi}"
        out.append(dict(aliquot=row["aliquot"], svclass=cls, chrom="chr" + c,
                        left_pos0=lo, right_pos0=hi, unique_td_id=key,
                        source_event_id=row["event_id"], source_side=side))
    return out


def negative_masks(donor: str, aliquot: str, positive_donors: set, positive_aliquots: set):
    return aliquot not in positive_aliquots, donor not in positive_donors
