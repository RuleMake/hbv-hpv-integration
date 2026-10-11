"""Define gene-length candidate matching and paired spatial-null operations without file I/O."""
import numpy as np


def checked_counts(values):
    a=np.asarray(values)
    if a.dtype.kind not in 'iuf' or not np.all(np.isfinite(a)) or np.any(a<0) or np.any(a!=np.floor(a)) or np.any(a>=2**53):
        raise ValueError('counts must be finite nonnegative exact integers')
    return a.astype(np.int64)


def overlap_category(counts):
    counts=np.asarray(counts)
    if not np.issubdtype(counts.dtype,np.integer) or np.any(counts<0):
        raise ValueError('gene overlap counts must be nonnegative integers')
    return np.minimum(counts,2)


def length_candidate_masks(observed_counts,observed_logspan,candidate_counts,candidate_logspan,caliper=.20):
    oc=np.asarray(observed_counts); ol=np.asarray(observed_logspan,dtype=float)
    cc=np.asarray(candidate_counts); cl=np.asarray(candidate_logspan,dtype=float)
    if oc.ndim!=1 or ol.shape!=oc.shape or cc.ndim!=2 or cc.shape[1]!=oc.size or cl.shape!=cc.shape:
        raise ValueError('annotation dimensions disagree')
    if not np.isfinite(caliper) or caliper<0:
        raise ValueError('invalid caliper')
    og=oc>0; cg=cc>0
    if np.any(~np.isfinite(ol[og])) or np.any(~np.isnan(ol[~og])):
        raise ValueError('observed genic length must be finite; intergenic length must be NaN')
    if np.any(~np.isfinite(cl[cg])) or np.any(~np.isnan(cl[~cg])):
        raise ValueError('candidate genic length must be finite; intergenic length must be NaN')
    status=np.all(overlap_category(cc)==overlap_category(oc)[None,:],axis=1)
    close=np.all(np.where(og[None,:],np.abs(cl-ol[None,:])<=caliper,True),axis=1)
    return status,status & close


def joint_displacements(pools,n_draws,seed):
    """pools[anchor][arm], each sorted unique nonzero integer displacement list.

    Original anchor distribution is preserved; no anchor is removed. Common uniform
    quantiles pair candidate lists, which is a computational coupling, not paired patients.
    """
    if not isinstance(n_draws,(int,np.integer)) or n_draws<2 or not pools:
        raise ValueError('need >=2 draws and >=1 anchor')
    n_arms=len(pools[0])
    if n_arms<1 or any(len(x)!=n_arms for x in pools):
        raise ValueError('arms inconsistent across anchors')
    pools=[[np.asarray(p) for p in anchor] for anchor in pools]
    for anchor in pools:
        for p in anchor:
            p=np.asarray(p)
            if p.ndim!=1 or not len(p) or not np.issubdtype(p.dtype,np.integer) or np.any(p==0) or np.any(np.diff(p)<=0):
                raise ValueError('candidate displacements must be nonempty, sorted, unique, nonzero integers')
    rng=np.random.default_rng(seed)
    anchor=rng.integers(0,len(pools),size=n_draws)
    u=rng.random(n_draws)
    result=np.empty((n_arms,n_draws),dtype=np.int64)
    for k,p in enumerate(pools):
        mask=anchor==k
        for a,values in enumerate(p):
            result[a,mask]=values[np.floor(u[mask]*len(values)).astype(np.int64)]
    return result,anchor,u


def summarize_spatial(observed,null):
    observed=checked_counts(observed)
    if observed.ndim!=0:raise ValueError('observed count must be scalar')
    observed=int(observed)
    null=checked_counts(null).astype(float)
    if null.ndim!=1 or len(null)<2:
        raise ValueError('invalid null count stream')
    mu=float(null.mean()); lo,hi=np.percentile(null,[2.5,97.5])
    upper=(int(np.sum(null>=observed))+1)/(len(null)+1)
    lower=(int(np.sum(null<=observed))+1)/(len(null)+1)
    return dict(observed=int(observed),expected=mu,fold=observed/mu if mu>0 else None,
                null_quantile_025=float(lo),null_quantile_975=float(hi),
                fold_null_low=observed/hi if hi>0 else None,
                fold_null_high=observed/lo if lo>0 else None,
                p_upper=upper,p_lower=lower,p_two_sided=min(1.,2*min(upper,lower)),
                n_permutations=len(null),interval_type='observed/null spatial reference, not sampling CI')


def calibrated_sign(observed_by_cluster,null_by_cluster):
    obs=checked_counts(observed_by_cluster)
    null=checked_counts(null_by_cluster)
    if null.ndim!=2 or obs.shape!=(null.shape[0],) or null.shape[0]<1 or null.shape[1]<2:
        raise ValueError('cluster/null dimensions invalid')
    if np.max(null)>np.iinfo(np.int64).max//null.shape[1]:raise ValueError('count multiplication would overflow')
    observed=int(np.sum(obs>null.mean(axis=1)))
    # Equivalent to x > leave-one-out mean without a precision-dependent division.
    null_sign=np.sum(null*null.shape[1]>null.sum(axis=1)[:,None],axis=0,dtype=np.int64)
    return observed,null_sign
