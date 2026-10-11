"""Load structural-variant catalogue representations for the shared spatial estimators."""
from __future__ import annotations
import hashlib
import json
import sys
from functools import lru_cache
from pathlib import Path, PureWindowsPath

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WS = ROOT.parents[1]
CATALOGUE = WS / 'analysis_inputs/structural_variant_catalogue'
LEDGER = CATALOGUE / 'event_ledger/outputs'
DATA01 = CATALOGUE / 'donor_duplication_statistics/data'
sys.path.insert(0, str(CATALOGUE / 'event_ledger/runtime_deps'))
import numpy as np
import pandas as pd
import spatial_null_core as F
import joint_covariate_matching as J
import amplification_peak_mask as A

VARIANTS = ('legacy_exact', 'legacy_expanded', 'paired_sample', 'donor_signature')
CLASSES = ('DEL', 'DUP', 'INV', 'TRA')
VIRUSES = ('HBV', 'HPV', 'EBV')
FILES = {
    'legacy_exact': LEDGER / 'legacy_exact_rebuilt.parquet',
    'legacy_expanded': LEDGER / 'legacy_expanded_selection.parquet',
    'paired_sample': LEDGER / 'normalized_analysis_sample_events.parquet',
    'donor_signature': LEDGER / 'donor_event_exactsignature_groups.parquet',
}

def load_modules():
    return F, J, A

@lru_cache(None)
def _metadata():
    mapping = pd.read_csv(LEDGER / 'aliquot_donor_mapping.csv').set_index('aliquot')['donor']
    cohort = pd.read_csv(DATA01 / 'integration_cohort.csv')
    positive_donors = set(cohort.donor.astype(str))
    positive_aliquots = {PureWindowsPath(str(p)).name.split('.')[0] for p in cohort.sv_path}
    return mapping, positive_donors, positive_aliquots

@lru_cache(None)
def _events(name):
    if name not in VARIANTS:
        raise ValueError(name)
    mapping, positive_donors, positive_aliquots = _metadata()
    if name.startswith('legacy_'):
        d = pd.read_parquet(FILES[name]).rename(columns={'unique_td_id':'event_id'})
        d['donor'] = d.aliquot.map(mapping)
        assert d.donor.notna().all()
    else:
        d = pd.read_parquet(FILES['paired_sample'])
        d['source_svclass'] = d.svclass
        d['svclass'] = d.analysis_class
        if name == 'donor_signature':
            groups = pd.read_parquet(FILES[name])
            d = d.drop_duplicates('donor_event_exactsignature').copy()
            assert set(d.donor_event_exactsignature) == set(groups.donor_event_exactsignature)
            d['event_id'] = d.donor_event_exactsignature
            # Metadata only: geometry comes from an exact paired signature representative.
            d['aliquot'] = d.donor_event_exactsignature.map(groups.set_index('donor_event_exactsignature').aliquots)
    d['integration_positive_donor'] = d.donor.isin(positive_donors)
    d['legacy_positive_aliquot'] = d.aliquot.isin(positive_aliquots)
    assert d.event_id.is_unique
    return d

def load_variant_events(name, exclude_positive=False):
    d = _events(name)
    if exclude_positive:
        column = 'legacy_positive_aliquot' if name.startswith('legacy_') else 'integration_positive_donor'
        d = d.loc[~d[column]]
    return d.copy()

def get_variant_geometry(name, svclass):
    """Return (per-chromosome endpoint geometry, globally unique ID count).

    F.ChromEvents counts endpoint hits, NOT overlaps of the interval's interior.
    TRA has two chromosome-specific point rows with the same globally unique event ID.
    """
    d = load_variant_events(name)
    d = d[d.svclass == svclass].copy()
    common = ['event_id', 'svclass', 'donor', 'aliquot']
    if name.startswith('legacy_'):
        geometry = d[['chrom', 'left_pos0', 'right_pos0'] + common].copy()
    elif svclass == 'TRA':
        rows = []
        for side in (1, 2):
            g = d[common].copy()
            g['chrom'] = 'chr' + d[f'chrom{side}_normalized'].astype(str)
            g['left_pos0'] = d[f'start{side}'].to_numpy(np.int64)
            g['right_pos0'] = g.left_pos0
            rows.append(g)
        geometry = pd.concat(rows, ignore_index=True)
    else:
        geometry = d[common].copy()
        assert (d.chrom1_normalized == d.chrom2_normalized).all()
        geometry['chrom'] = 'chr' + d.chrom1_normalized.astype(str)
        geometry['left_pos0'] = np.minimum(d.start1.to_numpy(np.int64), d.start2.to_numpy(np.int64))
        geometry['right_pos0'] = np.maximum(d.start1.to_numpy(np.int64), d.start2.to_numpy(np.int64))
    assert not geometry.duplicated(['chrom', 'event_id']).any()
    return geometry.reset_index(drop=True), int(d.event_id.nunique())

def load_breakpoints(name, exclude_positive=True):
    d = load_variant_events(name, exclude_positive=exclude_positive)
    common = ['event_id','donor','aliquot','svclass']
    rows = []
    for side in (1, 2):
        if name.startswith('legacy_'):
            sub = d if side == 1 else d[d.svclass != 'TRA']
            g = sub[common].copy()
            g['chrom'] = sub.chrom
            g['pos0'] = sub['left_pos0' if side == 1 else 'right_pos0'].to_numpy(np.int64)
        else:
            g = d[common].copy()
            g['chrom'] = 'chr' + d[f'chrom{side}_normalized'].astype(str)
            g['pos0'] = d[f'start{side}'].to_numpy(np.int64)
        g['endpoint_index'] = side
        rows.append(g)
    out = pd.concat(rows, ignore_index=True)
    assert not out.duplicated(['event_id','endpoint_index']).any()
    if not name.startswith('legacy_'):
        assert len(out) == 2 * len(d)
    return out[['chrom','pos0','event_id','donor','aliquot','svclass','endpoint_index']]

def eligible_coordinates(co, nul, mask, virus):
    sub = co[co.virus_class == virus].sort_values(['chrom','pos0']).reset_index(drop=True)
    in_mask = np.array([(str(c), int(p)//F.BIN_SIZE*F.BIN_SIZE) in mask
                        for c,p in zip(sub.chrom, sub.pos0)])
    sub = sub[~in_mask].reset_index(drop=True)
    k = pd.MultiIndex.from_arrays([sub.chrom.astype(str).to_numpy(),
                                   sub.pos0.to_numpy(np.int64)//F.BIN_SIZE*F.BIN_SIZE])
    sub = sub.assign(stratum=nul.lookup.reindex(k).to_numpy())
    sub = sub[sub.stratum.notna()].reset_index(drop=True)
    sub['stratum'] = sub.stratum.astype(int)
    sub['within_bin_offset'] = sub.pos0 - sub.pos0//F.BIN_SIZE*F.BIN_SIZE
    return sub

def load_context():
    feat = pd.read_parquet(ROOT / 'data/joint_features.parquet')
    mask = A.masked_bins(A.amp_peaks(), 100_000)
    nul = J.JointNuller(feat, F.Nuller(), mask)
    co = F.load_coordinates()  # Original pre-joint table; filter identically per virus.
    co_joint = pd.concat([eligible_coordinates(co, nul, mask, v) for v in VIRUSES], ignore_index=True)
    assert co_joint.virus_class.value_counts().to_dict() == {'HPV':3712, 'HBV':3138, 'EBV':107}
    return dict(F=F, J=J, A=A, feat=feat, mask=mask, nul=nul, co=co, co_joint=co_joint)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for c in iter(lambda:f.read(1024*1024), b''):
            h.update(c)
    return h.hexdigest()

def main():
    out = ROOT / 'data'
    rows = []
    for name in VARIANTS:
        for exclude in (False, True):
            d = load_variant_events(name, exclude_positive=exclude)
            bp = load_breakpoints(name, exclude_positive=exclude)
            for cls in CLASSES:
                sub = d[d.svclass == cls]
                rows.append(dict(variant=name, exclude_positive=exclude, svclass=cls,
                                 n_units=len(sub), n_endpoints=int((bp.svclass == cls).sum()),
                                 n_donors=sub.donor.nunique()))
    pd.DataFrame(rows).to_csv(out / 'adapter_unit_counts.csv', index=False)
    context = load_context()
    context['co_joint'].to_parquet(out / 'joint_eligible_coordinates.parquet', index=False)
    context['feat'].to_csv(out / 'joint_feature_audit.csv', index=False, float_format='%.17g')
    for v in VIRUSES:
        c = context['co_joint'].query('virus_class == @v')
        assert len(np.unique(F.leader_clusters(c.chrom.to_numpy(str),c.pos0.to_numpy(np.int64),100_000))) == {'HBV':2421,'HPV':2390,'EBV':93}[v]
    paths = list(FILES.values()) + [LEDGER/'aliquot_donor_mapping.csv']
    paths += [DATA01/n for n in ['integration_cohort.csv','genome_features_50kb.parquet',
              'tian_coordinates.tsv','tian_matching_audit.tsv','gistic_peaks.txt.gz']]
    paths += [ROOT/'data/joint_features.parquet', ROOT/'data/wgEncodeCrgMapabilityAlign100mer.bigWig']
    pd.DataFrame([dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in paths]).to_csv(out/'input_manifest.csv',index=False)
    print(pd.DataFrame(rows).to_string(index=False))
    print('Shared adapters and context validated.')

if __name__ == '__main__':
    main()
