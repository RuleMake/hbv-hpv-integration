"""Prepare Figure 1 display tables without changing the saved analysis results."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
(HERE/'review').mkdir(exist_ok=True)
(HERE/'plotted_data').mkdir(exist_ok=True)
RAW = HERE / 'source_data'
OUT = HERE / 'plotted_data'
OUT.mkdir(exist_ok=True)

ledger = pd.read_csv(RAW / 'cohort_units.csv')
def cohort(rule):
    d = ledger[(ledger.dataset == 'normalized_sample_events') &
               (ledger.exclusion_rule == rule) & (ledger.svclass == 'ALL')]
    assert len(d) == 1
    return d.iloc[0]
full, neg = cohort('all'), cohort('whole_positive_donor_exclusion')
dup = pd.read_csv(HERE / 'supporting_data/Figure_1B_duplication_reference.csv').iloc[0]
pd.DataFrame([
    dict(panel='B', measure='Donors', unit='donor', total=int(full.n_donors),
         without_detected_integration=int(neg.n_donors), percentage=100*neg.n_donors/full.n_donors),
    dict(panel='B', measure='Aliquots', unit='aliquot', total=int(full.n_aliquots),
         without_detected_integration=np.nan, percentage=np.nan),
    dict(panel='B', measure='Structural variants', unit='paired call', total=int(full.n_units),
         without_detected_integration=int(neg.n_units), percentage=100*neg.n_units/full.n_units),
    dict(panel='B', measure='Large duplications', unit='record', total=int(dup.n_records),
         without_detected_integration=np.nan, percentage=np.nan),
]).to_csv(OUT / 'Figure_1B_legend_values.csv', index=False)

spacing = pd.read_csv(RAW / 'main/Figure_1C.csv')
rows = []
for (virus, kind), group in spacing.groupby(['virus', 'kind'], sort=False):
    v = group.spacing_bp.to_numpy()
    shown = np.maximum(v, 1)
    rows.append(dict(panel='C', virus=virus,
                     placement='Observed' if kind == 'observed' else 'Independent',
                     n_spacings=len(v), median_bp=float(np.median(v)),
                     p10_bp=float(np.percentile(v, 10)), p90_bp=float(np.percentile(v, 90)),
                     displayed_median_bp=float(np.median(shown)),
                     displayed_p10_bp=float(np.percentile(shown, 10)),
                     displayed_p90_bp=float(np.percentile(shown, 90)),
                     zeros_displayed_at_1bp=int((v == 0).sum()),
                     range_type='10th-90th percentiles of adjacent spacing'))
pd.DataFrame(rows).to_csv(OUT / 'Figure_1C_plotted_values.csv', index=False)

tracks = pd.read_csv(RAW / 'main/Figure_1D.csv', dtype={'chrom': str})
tracks = tracks[tracks.track.isin(['SV breakends', 'HBV', 'HPV'])].copy()
peaks = pd.read_csv(RAW / 'main/Figure_1D_peaks.csv', dtype={'chrom': str})
lengths = tracks.groupby('chrom').chrom_length.first().to_dict()
BIN, PAD = 50_000, 100_000
effective = []
for p in peaks.itertuples(index=False):
    # Reproduce the analysis mask: the inclusive last grid-bin start contributes
    # one whole 50-kb bin. These are filtering-bin edges, not adjusted peak calls.
    start = max(0, (int(p.start) - PAD) // BIN) * BIN
    end = ((int(p.end) + PAD) // BIN + 1) * BIN
    effective.append(dict(chrom=p.chrom, start=start, end=min(end, int(lengths[p.chrom]))))
merged = []
for chrom, group in pd.DataFrame(effective).groupby('chrom', sort=False):
    intervals = sorted(zip(group.start.astype(int), group.end.astype(int)))
    lo, hi = intervals[0]
    for s, e in intervals[1:]:
        if s <= hi:
            hi = max(hi, e)
        else:
            merged.append(dict(chrom=chrom, start=lo, end=hi))
            lo, hi = s, e
    merged.append(dict(chrom=chrom, start=lo, end=hi))
mask = pd.DataFrame(merged)
mask.to_csv(OUT / 'Figure_1D_query_exclusions.csv', index=False)

fractions = []
for row in tracks.itertuples(index=False):
    same = mask[mask.chrom == row.chrom]
    overlap = sum(max(0, min(int(row.bin_end), e) - max(int(row.bin_start), s))
                  for s, e in zip(same.start, same.end))
    fractions.append(overlap / (row.bin_end - row.bin_start))
tracks['query_excluded_fraction'] = fractions
tracks['query_exclusion_applies'] = tracks.track != 'SV breakends'
tracks['count_basis'] = np.where(tracks.track == 'SV breakends',
    'All eligible paired-call breakends', 'Retained integration query coordinates')
tracks.to_csv(OUT / 'Figure_1D_plotted_values.csv', index=False)

complete = tracks[(tracks.track == 'SV breakends') & (tracks.query_excluded_fraction == 1)]
assert len(complete) == 20 and complete['count'].sum() == 9592
assert tracks[(tracks.track != 'SV breakends') & (tracks.query_excluded_fraction == 1)]['count'].sum() == 0
report = dict(
    new_biological_analyses=False,
    mask_grid_bp=BIN, mask_padding_bp=PAD, merged_exclusion_intervals=len(mask),
    fully_excluded_display_bins=len(complete), sv_breakends_in_these_bins=int(complete['count'].sum()),
    count_totals=tracks.groupby('track')['count'].sum().to_dict(),
    source_hashes={str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sorted(RAW.rglob('*.csv'))},
    transformation_notes={
      'B':'Counts are selected by explicit cohort and unit. Percentages are numerator/denominator times 100. Blank fields are not plotted.',
      'C':'Percentiles and ECDF are derived from the saved spacing records; only zero spacing positions are displayed at 1 bp on the log axis.',
      'D':'Counts unchanged. Query exclusions reproduce the saved 100-kb padding and 50-kb filtering grid; mask intervals are merged and clipped to chromosome length for display.',
    })
(HERE / 'review/source_preparation.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ['source_hashes','transformation_notes']}, indent=2))
