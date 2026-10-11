"""Estimate the four-class HBV-HPV contrast from paired 100-kb spatial randomizations."""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path
sys.dont_write_bytecode=True
from local_sv_class_composition import ROOT, OUT, SHARED, VARIANTS, CLASSES, SEED, N_PERM, sha, write_csv
import numpy as np
import pandas as pd

WINDOW=100_000
VIRUSES=('HBV','HPV')
CROSSED=(('crossed_all',('DEL','DUP'),('TRA','INV')),('crossed_concordant',('DUP',),('TRA','INV')))

def log_contrast(oh,nh,ob,nb):
    eh,eb=nh.mean(),nb.mean()
    if min(oh,ob,eh,eb)<=0:
        raise ValueError('Observed or expected zero makes original log contrast undefined')
    return float(np.log(oh/eh)-np.log(ob/eb)), np.log(np.maximum(nh,1)/eh)-np.log(np.maximum(nb,1)/eb)

def summarize(dobs,dnull):
    sd=float(dnull.std(ddof=1))
    mean=float(dnull.mean())
    ge,le=int((dnull>=dobs).sum()),int((dnull<=dobs).sum())
    hi,lo=(ge+1)/(len(dnull)+1),(le+1)/(len(dnull)+1)
    qlo,qhi=np.percentile(dnull,[2.5,97.5])
    return dict(log_ratio=dobs,ratio_of_folds=float(np.exp(dobs)),null_mean=mean,
        null_reference_q025=float(qlo),null_reference_q975=float(qhi),null_sd=sd,
        interval_type='randomization_reference_not_sampling_CI',
        z_legacy_zero_centered=dobs/sd if sd>0 else np.nan,
        z_null_mean_centered_diagnostic=(dobs-mean)/sd if sd>0 else np.nan,
        n_null_ge_observed=ge,n_null_le_observed=le,emp_p_greater=hi,emp_p_less=lo,
        emp_p_two_sided=min(1.,2*min(hi,lo)),mde_log_ratio_approx_2p80sd=2.80*sd,
        clears_legacy_mde=bool(abs(dobs)>=2.80*sd),
        mde_interpretation='approximate_sensitivity_not_effect_upper_bound')

def compare_legacy(table):
    ref=pd.read_csv(ROOT/'provenance/legacy_class_interaction_interaction.csv',float_precision='round_trip')
    mapping={'z':'z_legacy_zero_centered','null_ci_lo':'null_reference_q025','null_ci_hi':'null_reference_q975',
        'mde_log_ratio':'mde_log_ratio_approx_2p80sd','clears_mde':'clears_legacy_mde'}
    rounding={'fold_HPV':4,'fold_HBV':4,'log_ratio':5,'ratio_of_folds':4,
              'null_ci_lo':5,'null_ci_hi':5,'z':2,'null_stream_corr':4}
    checks=[]
    for row in ref.to_dict('records'):
        got=table[(table.estimator==row['estimator'])&(table.statistic==row['statistic'])]
        assert len(got)==1
        got=got.iloc[0]
        for key,expected in row.items():
            if pd.isna(expected):
                continue
            actual=got[mapping.get(key,key)]
            if key=='mde_log_ratio':
                actual=round(round(2.486*float(got.null_sd),5)*(2.80/2.486),5)
            elif key in rounding:
                actual=round(float(actual),rounding[key])
            passed=bool(np.isclose(actual,expected,rtol=0,atol=1e-13)) if isinstance(expected,(int,float)) else actual==expected
            checks.append(dict(estimator=row['estimator'],statistic=row['statistic'],metric=key,
                legacy=expected,recomputed_at_legacy_precision=actual,passed=passed))
    return pd.DataFrame(checks)

def run_variant(common,ctx,variant):
    start=time.time()
    out=OUT/variant
    out.mkdir(parents=True,exist_ok=True)
    F,nul=ctx['F'],ctx['nul']
    draws,arrays,counts,provenance={},{},[],[]
    for cls in CLASSES:
        geometry,n_events=common.get_variant_geometry(variant,cls)
        idx={str(c):F.ChromEvents(g.left_pos0.to_numpy(np.int64),g.right_pos0.to_numpy(np.int64),g.event_id.to_numpy(str)) for c,g in geometry.groupby('chrom')}
        geometry_hash=__import__('hashlib').sha256(geometry.to_csv(index=False,lineterminator='\n').encode()).hexdigest()
        for vi,virus in enumerate(VIRUSES):
            sub=ctx['co_joint'][ctx['co_joint'].virus_class==virus].sort_values(['chrom','pos0']).reset_index(drop=True)
            labels=F.leader_clusters(sub.chrom.to_numpy(str),sub.pos0.to_numpy(np.int64),WINDOW)
            nul.n_calls=nul.n_fallback=0
            o,nl,per_cluster,so,sn=F.est_block(idx,nul,np.random.default_rng([SEED,vi]),sub,WINDOW,N_PERM,labels,unit='uniqSV')
            fallback=100*nul.n_fallback/max(nul.n_calls,1)
            for est,observed,null in [('bc',float(o),np.asarray(nl,float)),('sign',float(so),np.asarray(sn,float))]:
                draws[(est,cls,virus)]=(observed,null)
                arrays[f'{est}__{cls}__{virus}__null_counts']=null
                arrays[f'{est}__{cls}__{virus}__observed_count']=np.asarray(observed)
                counts.append(dict(variant=variant,estimator=est,svclass=cls,virus=virus,
                    n_events=n_events,n_coords=len(sub),n_clusters=len(np.unique(labels)),observed=observed,
                    null_mean=float(null.mean()),fold=observed/null.mean(),
                    n_null_zero=int((null==0).sum()),n_null_below_one=int((null<1).sum()),
                    null_min=float(null.min()),null_max=float(null.max()),null_fallback_pct=fallback))
            write_csv(per_cluster,out/f'class_interaction_{cls}_{virus}_clusters.csv')
            provenance.append(dict(variant=variant,svclass=cls,virus=virus,seed=f'[{SEED},{vi}]',
                geometry_sha256=geometry_hash,n_geometry_rows=len(geometry),n_global_events=n_events,
                coordinate_sha256=__import__('hashlib').sha256(sub.to_csv(index=False,lineterminator='\n').encode()).hexdigest(),
                n_pool_calls=nul.n_calls,n_pool_fallback=nul.n_fallback))
            print(f'CLASS_INTERACTION {variant} {cls} {virus}: bc={o}/{np.mean(nl):.8f} sign={so}/{np.mean(sn):.8f} ({time.time()-start:.0f}s)',flush=True)
    rows=[]
    for est in ('bc','sign'):
        for cls in CLASSES:
            oh,nh=draws[(est,cls,'HPV')]
            ob,nb=draws[(est,cls,'HBV')]
            observed,null=log_contrast(oh,nh,ob,nb)
            statistic=f'HPV_vs_HBV_{cls}'
            arrays[f'{est}__{statistic}__log_null']=null
            rows.append(dict(variant=variant,estimator=est,svclass=cls,statistic=statistic,n_perm=N_PERM,
                window_bp=WINDOW,fold_HPV=oh/nh.mean(),fold_HBV=ob/nb.mean(),
                null_stream_corr=float(np.corrcoef(nh,nb)[0,1]),
                n_floor_applications=int((nh<1).sum()+(nb<1).sum()),**summarize(observed,null)))
        for name,hpv,hbv in CROSSED:
            observed,null=0.,np.zeros(N_PERM)
            floor=0
            for cls,sign in [(c,1) for c in hpv]+[(c,-1) for c in hbv]:
                oh,nh=draws[(est,cls,'HPV')]
                ob,nb=draws[(est,cls,'HBV')]
                do,dn=log_contrast(oh,nh,ob,nb)
                observed+=sign*do
                null+=sign*dn
                floor+=int((nh<1).sum()+(nb<1).sum())
            arrays[f'{est}__{name}__log_null']=null
            rows.append(dict(variant=variant,estimator=est,svclass=f"({'+'.join(hpv)})-({'+'.join(hbv)})",
                statistic=name,n_perm=N_PERM,window_bp=WINDOW,fold_HPV=np.nan,fold_HBV=np.nan,
                null_stream_corr=np.nan,n_floor_applications=floor,**summarize(observed,null)))
    table=pd.DataFrame(rows)
    write_csv(table,out/'class_interaction_interaction.csv')
    write_csv(pd.DataFrame(counts),out/'class_interaction_virus_counts.csv')
    write_csv(pd.DataFrame(provenance),out/'class_interaction_call_provenance.csv')
    np.savez_compressed(out/'class_interaction_null_draws.npz',**arrays)
    if variant=='legacy_exact':
        audit=compare_legacy(table)
        write_csv(audit,out/'class_interaction_legacy_reproduction.csv')
        passed=bool(audit.passed.all())
        (out/'class_interaction_legacy_gate.json').write_text(json.dumps(dict(passed=passed,n_checks=len(audit),n_passed=int(audit.passed.sum())),indent=2))
        if not passed:
            print(audit[~audit.passed].to_string(index=False),flush=True)
            raise AssertionError('Virus-contrast reference check failed')
    meta=dict(variant=variant,seed=SEED,n_perm=N_PERM,window_bp=WINDOW,
        script_sha256=sha(__file__),protocol_sha256=sha(ROOT/'PROTOCOL.md'),shared_adapter_sha256=sha(SHARED/'shared_analysis_context.py'),
        rng_rule='reset default_rng([7,vi]) for each class; HBV vi0, HPV vi1',
        geometry_exclusion='all canonical-largeDUP roster donors retained; no negative-reference exclusion in CLASS_INTERACTION')
    (out/'class_interaction_run.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    return table

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--variant',choices=VARIANTS+('all',),default='legacy_exact')
    a=p.parse_args()
    sys.path.insert(0,str(SHARED))
    import shared_analysis_context as common
    ctx=common.load_context()
    variants=VARIANTS if a.variant=='all' else (a.variant,)
    if variants[0]!='legacy_exact':
        gate=OUT/'legacy_exact/class_interaction_legacy_gate.json'
        if not gate.exists() or not json.loads(gate.read_text())['passed']:
            raise ValueError('Virus-contrast reference check must pass first')
    for variant in variants:
        run_variant(common,ctx,variant)

if __name__=='__main__':
    main()
