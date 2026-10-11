"""Load fixed coordinates and variant geometry for replication-timing sensitivity.

The imported helpers define data structures and estimators without running their
analysis entry points. No file I/O occurs on import; run with Python -B."""
from __future__ import annotations
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
WS = HERE.parents[2]
COMPOSITION_ANALYSIS = WS / 'analysis_inputs/class_composition'
OUT = HERE / 'outputs'
sys.path.insert(0, str(WS / 'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
sys.path.insert(0, str(COMPOSITION_ANALYSIS / 'shared'))
import numpy as np
import pandas as pd
import shared_analysis_context as common

def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

E = _module('composition_reference', COMPOSITION_ANALYSIS / 'axis/axis_recompute.py')
AT = _module('atlas_reference', COMPOSITION_ANALYSIS / 'atlas/atlas_core.py')
F, J, A = common.load_modules()
VARIANTS = ('paired_sample', 'donor_signature')
VIRUSES = ('HBV', 'HPV')
CLASSES = ('DEL', 'DUP', 'INV', 'TRA')
BIN = 50_000
N_PERM = 20_000

def guarded(path):
    path = Path(path).resolve()
    if not path.is_relative_to(HERE):
        raise ValueError(f'Output outside the configured baseline directory: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    return path

def write_csv(df, path):
    df.to_csv(guarded(path), index=False, float_format='%.17g', na_rep='')

def write_json(value, path):
    guarded(path).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def load_reference_context():
    """Load the consensus timing reference with its original fallback matching rule."""
    return common.load_context()

def load_coordinates(virus):
    if virus not in VIRUSES:
        raise ValueError(virus)
    return pd.read_csv(OUT / f'coords_{virus}.csv', float_precision='round_trip')

def load_endpoints(variant):
    if variant not in VARIANTS:
        raise ValueError(variant)
    return pd.read_parquet(OUT / f'composition_breakpoints_{variant}.parquet')

def load_geometry(variant, svclass=None):
    if variant not in VARIANTS or (svclass is not None and svclass not in CLASSES):
        raise ValueError((variant, svclass))
    d = pd.read_parquet(OUT / f'atlas_geometry_{variant}.parquet')
    return d if svclass is None else d[d.svclass == svclass].reset_index(drop=True)

def make_index(geometry):
    return AT.make_index(geometry, F.ChromEvents)

def endpoint_maps(endpoints):
    return E.propensity_map(endpoints)

def score_composition(coordinates, nuller, maps, virus, n_perm=N_PERM):
    """Original leader-anchor seed [7, virus index], equal-cluster mean, 41 bins."""
    if virus not in VIRUSES:
        raise ValueError(virus)
    cum, nb, pg = maps
    labels = coordinates.original_cluster.to_numpy(np.int64)
    return E.score(coordinates, labels, nuller, cum, nb, pg,
                   np.random.default_rng([7, VIRUSES.index(virus)]), n_perm,
                   verify_lookup=True)

def score_atlas(coordinates, nuller, index, n_perm=N_PERM):
    """Original random-member anchor, seed7; endpoint hits per cluster, sum clusters."""
    rng = AT.TrackedRNG(7)
    result = F.est_block(index, nuller, rng, coordinates, 100_000, n_perm,
                         coordinates.original_cluster.to_numpy(np.int64), unit='uniqSV')
    return result, rng.hexdigest()
