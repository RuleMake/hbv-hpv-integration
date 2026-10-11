"""Reuse spatial counting and placement routines for replication-timing sensitivity."""
from pathlib import Path
import importlib.util
import numpy as np

SOURCE=Path(__file__).resolve().parents[2]/'gene_length_matching/analysis/matching_core.py'
spec=importlib.util.spec_from_file_location('gene_length_matching_core',SOURCE)
_previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(_previous)
joint_displacements=_previous.joint_displacements
summarize_counts=_previous.summarize_spatial
calibrated_sign=_previous.calibrated_sign

def axis_of(composition):
    p=np.asarray(composition,dtype=float)
    if p.ndim not in (1,2) or p.shape[-1]!=4 or not np.isfinite(p).all() or np.any(p<0) or not np.allclose(p.sum(axis=-1),1.,atol=1e-12,rtol=0):
        raise ValueError('Expected finite DEL,DUP,INV,TRA composition rows summing to one')
    return p[...,2:].sum(axis=-1)-p[...,:2].sum(axis=-1)

def signed_summary(observed,null):
    o=np.asarray(observed,float);n=np.asarray(null,float)
    if o.ndim!=0 or not np.isfinite(o) or n.ndim!=1 or len(n)<2 or not np.isfinite(n).all():
        raise ValueError('Invalid signed observed/null values')
    o=float(o);mu=float(n.mean());q=np.percentile(n,[2.5,97.5]);N=len(n)
    pu=(np.count_nonzero(n>=o)+1)/(N+1);pl=(np.count_nonzero(n<=o)+1)/(N+1)
    return dict(observed=o,null_mean=mu,excess=o-mu,null_sd=float(n.std(ddof=1)),
                null_reference_q025=float(q[0]),null_reference_q975=float(q[1]),
                excess_reference_low=float(o-q[1]),excess_reference_high=float(o-q[0]),
                p_upper=float(pu),p_lower=float(pl),p_two_sided=float(min(1.,2*min(pu,pl))),
                n_permutations=N,interval_type='spatial-null reference; not donor sampling CI')

def interaction(obs_hbv,obs_hpv,null_hbv,null_hpv):
    h=np.asarray(null_hbv,float);p=np.asarray(null_hpv,float)
    if h.ndim!=1 or h.shape!=p.shape or len(h)<2 or not np.isfinite(h).all() or not np.isfinite(p).all():
        raise ValueError('Interaction requires aligned finite virus null streams')
    if np.ndim(obs_hbv)!=0 or np.ndim(obs_hpv)!=0 or not np.isfinite(obs_hbv) or not np.isfinite(obs_hpv):
        raise ValueError('Interaction observations must be finite scalars')
    observed=(float(obs_hbv)-h.mean())-(float(obs_hpv)-p.mean())
    centered=(h-h.mean())-(p-p.mean())
    return observed,centered
