"""Evaluate gene-length sensitivity using the prespecified matching design."""
from pathlib import Path
import sys,json,hashlib,time
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'baseline'))
import load_inputs_and_original_null as B
import numpy as np
import pandas as pd
from matching_core import joint_displacements,summarize_spatial,calibrated_sign

OUT=ROOT/'analysis/outputs'
SCOPES=('protein_coding','all_gene')
ARMS=('joint_geometry','gene_status','gene_length')
VARIANTS=('donor_signature','paired_sample')
N=20000
SEED=20260910
WINDOW=100000

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writejson(name,x):(OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')

def check_freeze():
    freeze=json.loads((ROOT/'analysis/freeze_record.json').read_text(encoding='utf-8'))
    bad=[p for p,h in freeze['files'].items() if sha(ROOT/p)!=h]
    if bad:raise RuntimeError('Frozen files changed: '+str(bad))
    if not freeze['baseline_gate_passed'] or not freeze['design_review_passed']:raise RuntimeError('Required reference or design checks are missing')
    return freeze

def main():
    freeze=check_freeze()
    if (OUT/'run_metadata.json').exists():raise RuntimeError('Finished outcome run exists; do not overwrite')
    OUT.mkdir(exist_ok=True)
    geometries={v:B.load_geometry(v) for v in VARIANTS}
    indexes={v:B.make_index(g)[0] for v,g in geometries.items()}
    estimates=[];per=[];differences=[];support_effect=[];placement_diagnostics=[];old_common=[];start=time.time()
    for si,scope in enumerate(SCOPES):
        co=pd.read_csv(ROOT/f'design/common_coordinates_{scope}.csv',float_precision='round_trip')
        support=json.loads((ROOT/f'design/support_{scope}.json').read_text(encoding='utf-8'))
        npz=np.load(ROOT/f'design/pools_{scope}.npz')
        # Both representations and all three matching arms use identical placements.
        totals={v:np.zeros((3,N),dtype=np.int64) for v in VARIANTS}
        signs={v:np.zeros((3,N),dtype=np.int64) for v in VARIANTS}
        observed={v:0 for v in VARIANTS};observed_signs={v:np.zeros(3,dtype=np.int64) for v in VARIANTS}
        digest=hashlib.sha256();observed_union={v:[] for v in VARIANTS}
        for number,(cid,g) in enumerate(co.groupby('subset_cluster',sort=True)):
            p=g.pos0.to_numpy(np.int64);chrom=str(g.chrom.iloc[0]);lo=int(p.min())-WINDOW;L=int(p.max()-p.min())+2*WINDOW
            pools=[[npz[f'c{int(cid)}_k{k}_{arm}'] for arm in ARMS] for k in range(len(p))]
            # Stable original cluster identity, not retained-list order, identifies the stream.
            # Scope shares the same underlying uniforms as well, but may have different support.
            delta,anchor,u=joint_displacements(pools,N,np.random.SeedSequence([SEED,int(cid)]))
            digest.update(np.array([cid],np.int64).tobytes());digest.update(delta.tobytes())
            for a,arm in enumerate(ARMS):
                placement_diagnostics.append(dict(scope=scope,cluster=int(cid),arm=arm,
                    n_anchors=len(p),observed_window_start=lo,observed_window_end=lo+L,
                    sampled_distinct_displacements=len(np.unique(delta[a])),
                    minimum_displacement=int(delta[a].min()),maximum_displacement=int(delta[a].max())))
            for variant in VARIANTS:
                idx=indexes[variant].get(chrom)
                if idx is None:raise ValueError('No matching chromosome event index')
                o=int(idx.count(L,np.array([lo]))[0]);observed[variant]+=o
                observed_union[variant].extend(idx.ids_in(lo,lo+L).tolist())
                nl=idx.count(L,lo+delta.reshape(-1)).reshape(3,N)
                totals[variant]+=nl
                mus=nl.mean(axis=1)
                scored=[calibrated_sign(np.array([o]),nl[a:a+1]) for a in range(3)]
                sg=np.array([s[0] for s in scored],dtype=np.int64)
                sn=np.array([s[1] for s in scored],dtype=np.int64)
                observed_signs[variant]+=sg;signs[variant]+=sn
                for a,arm in enumerate(ARMS):
                    per.append(dict(scope=scope,variant=variant,cluster=int(cid),chrom=chrom,n_coordinates=len(p),arm=arm,
                        observed=o,expected=float(mus[a]),null_sd=float(nl[a].std(ddof=1)),
                        observed_cluster_sign=int(sg[a]),expected_cluster_sign=float(sn[a].mean()),
                        window_start=lo,window_end=lo+L))
            if number%400==0:print(f'{scope} outcomes {number}/{len(support["common_clusters"])} clusters elapsed {time.time()-start:.1f}s',flush=True)
        npz.close()
        for variant in VARIANTS:
            np.savez_compressed(OUT/f'null_{scope}_{variant}.npz',block=totals[variant],sign=signs[variant],
                                observed_block=np.array(observed[variant],dtype=np.int64),observed_sign=observed_signs[variant],arms=np.array(ARMS))
            for a,arm in enumerate(ARMS):
                for estimator,o,nl in [('cluster_event_incidence',observed[variant],totals[variant][a]),
                                        ('calibrated_cluster_sign',observed_signs[variant][a],signs[variant][a])]:
                    row=dict(scope=scope,variant=variant,arm=arm,estimator=estimator,n_coordinates=len(co),
                             n_clusters=co.subset_cluster.nunique(),**summarize_spatial(o,nl))
                    row['p_two_sided_bonferroni2']=min(1.,2*row['p_two_sided']) if (scope=='protein_coding' and variant=='donor_signature' and arm=='gene_length') else None
                    estimates.append(row)
            for before,after in [(0,1),(1,2),(0,2)]:
                for estimator,os,nls in [('cluster_event_incidence',np.repeat(observed[variant],3),totals[variant]),
                                          ('calibrated_cluster_sign',observed_signs[variant],signs[variant])]:
                    m=nls.mean(axis=1);fold=os/m
                    count_diff=nls[after]-nls[before]
                    # Delta-method Monte Carlo error of fixed observed / estimated null means.
                    influence=-os[after]/m[after]**2*(nls[after]-m[after])+os[before]/m[before]**2*(nls[before]-m[before])
                    differences.append(dict(scope=scope,variant=variant,estimator=estimator,
                        before=ARMS[before],after=ARMS[after],fold_before=float(fold[before]),fold_after=float(fold[after]),
                        fold_change=float(fold[after]-fold[before]),log_fold_change=float(np.log(fold[after]/fold[before])),
                        observed_minus_expected_change=float((os[after]-m[after])-(os[before]-m[before])),
                        expected_change=float(m[after]-m[before]),expected_change_MC_SE=float(count_diff.std(ddof=1)/np.sqrt(N)) if estimator=='cluster_event_incidence' else None,
                        fold_change_MC_SE=float(influence.std(ddof=1)/np.sqrt(N)) if estimator=='cluster_event_incidence' else None,
                        uncertainty_label=('Monte Carlo precision only; not effect CI; no attenuation hypothesis test' if estimator=='cluster_event_incidence'
                            else 'Sign MC precision not estimated because realized sign classifications and null signs depend on simulated reference means; point changes only, no effect CI')))
            support_effect.append(dict(scope=scope,variant=variant,n_coordinates=len(co),n_clusters=co.subset_cluster.nunique(),
                 global_unique_events_observed=len(set(observed_union[variant])),
                 cross_cluster_event_repeat_incidences=int(observed[variant]-len(set(observed_union[variant]))),
                 observed_cluster_event_incidences=observed[variant],placement_sha256=digest.hexdigest()))
            # Separate support selection from removal of fallback/full-block mask correction.
            # This preserves the original pool behavior on the selected cluster set; its RNG
            # resets to historical seed7 and is not paired with new rule A/B/C.
            ctx=B.load_context()
            o,nl,cl,so,sn,rngmeta=B.run_statistic(indexes[variant],ctx['nul'],co,seed=7,n_perm=N)
            np.savez_compressed(OUT/f'old_null_common_{scope}_{variant}.npz',block=nl,sign=sn,
                                observed_block=np.array(o),observed_sign=np.array(so))
            cl.to_csv(OUT/f'old_null_common_clusters_{scope}_{variant}.csv',index=False,float_format='%.17g')
            for estimator,value,stream in [('cluster_event_incidence',o,nl),('calibrated_cluster_sign',so,sn)]:
                old_common.append(dict(scope=scope,variant=variant,arm='original_null_on_common_support',estimator=estimator,
                    n_coordinates=len(co),n_clusters=co.subset_cluster.nunique(),**summarize_spatial(value,stream),**rngmeta))
    pd.DataFrame(estimates).to_csv(OUT/'estimates.csv',index=False,float_format='%.17g')
    pd.DataFrame(per).to_csv(OUT/'per_cluster_results.csv',index=False,float_format='%.17g')
    pd.DataFrame(differences).to_csv(OUT/'paired_changes.csv',index=False,float_format='%.17g')
    pd.DataFrame(support_effect).to_csv(OUT/'event_union_diagnostics.csv',index=False)
    pd.DataFrame(placement_diagnostics).to_csv(OUT/'placement_diagnostics.csv',index=False)
    pd.DataFrame(old_common).to_csv(OUT/'original_null_common_support.csv',index=False,float_format='%.17g')
    writejson('run_metadata.json',dict(completed=True,seed=SEED,seed_rule='SeedSequence([20260910, original_nonhotspot_cluster_id])',
        n_permutations=N,scopes=SCOPES,arms=ARMS,variants=VARIANTS,window_bp=WINDOW,
        script_sha256=sha(__file__),freeze_sha256=sha(ROOT/'analysis/freeze_record.json'),
        elapsed_seconds=time.time()-start,source_writes=False,original_fourcell_reopened=False,
        new_effect_difference_hypothesis_tests=0,interval_type='spatial null reference; no patient sampling confidence intervals'))
    print('Finished all frozen outcome arms',flush=True)

if __name__=='__main__':main()
