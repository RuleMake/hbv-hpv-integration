"""Load fixed inputs and evaluate the HPV duplication spatial statistic.

This module provides the input and estimator interface for gene-length matching.
It performs no file writes on import."""
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
WS=ROOT.parents[2]
COMPOSITION_ANALYSIS=WS/'analysis_inputs/class_composition'
sys.path.insert(0,str(WS/'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
sys.path.insert(0,str(COMPOSITION_ANALYSIS/'shared'))
sys.path.insert(0,str(COMPOSITION_ANALYSIS/'atlas'))
import numpy as np
import pandas as pd
import shared_analysis_context as common
from atlas_core import make_index as _make_index,full_summary,TrackedRNG

VARIANTS=('canonical_legacy','class_legacy','donor_signature','paired_sample')
SUBSETS=('ALL','HOTSPOT','NON-HOTSPOT')
WINDOW=100_000
N_PERM=20_000
SEED=7

def output(relative):
 p=(ROOT/relative).resolve()
 if not p.is_relative_to(ROOT):raise ValueError('Output outside the configured baseline directory')
 p.parent.mkdir(parents=True,exist_ok=True)
 return p

def load_context():
 """Load the shared genomic context, preserving matching pools and fallback behavior."""
 return common.load_context()

def load_coordinates(subset='NON-HOTSPOT'):
 if subset not in SUBSETS:raise ValueError(subset)
 co=pd.read_csv(ROOT/'outputs/coordinates_HPV.csv',float_precision='round_trip',keep_default_na=False)
 if subset!='ALL':co=co[co.hotspot_recurrence==(subset=='HOTSPOT')].copy()
 co=co.reset_index(drop=True)
 co['subset_cluster']=common.F.leader_clusters(co.chrom.to_numpy(str),co.pos0.to_numpy(np.int64),WINDOW)
 return co

def load_geometry(variant='donor_signature'):
 if variant not in VARIANTS:raise ValueError(variant)
 return pd.read_parquet(ROOT/f'outputs/geometry_{variant}.parquet')

def make_index(geometry):return _make_index(geometry,common.F.ChromEvents)

def run_statistic(index,nuller,coordinates,seed=SEED,n_perm=N_PERM):
 """Evaluate the random-member rigid-block estimator.
 
 Returns observed counts, total_null[draw], per-cluster observed counts and null
 means, sign statistics, and random-stream metadata. Each subset starts with
 random seed 7. The nuller supplies chromosome sizes, candidate draws and fallback
 counts. Placement retains complete clusters without boundary clipping."""
 rng=TrackedRNG(seed)
 labels=common.F.leader_clusters(coordinates.chrom.to_numpy(str),coordinates.pos0.to_numpy(np.int64),WINDOW)
 nuller.n_calls=nuller.n_fallback=0
 o,n,cl,s,sn=common.F.est_block(index,nuller,rng,coordinates,WINDOW,n_perm,labels,unit='uniqSV')
 return o,n,cl,s,sn,{'seed':seed,'rng_sha256':rng.hexdigest(),'rng_integer_calls':rng.n_calls,
  'pool_calls':nuller.n_calls,'pool_fallbacks':nuller.n_fallback,
  'pct_pool_fallback':100*nuller.n_fallback/max(nuller.n_calls,1)}
