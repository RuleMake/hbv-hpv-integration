"""Verify replication-timing sensitivity estimates from source inputs and saved outputs.

Production statistics and index-building functions are not imported. Run with
Python -B; statistical review summaries are written to the configured output path."""
from pathlib import Path
import sys, json, hashlib, time, ast
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
WS = STAGE.parents[1]
sys.path.insert(0, str(WS / 'analysis_inputs/structural_variant_catalogue/event_ledger/runtime_deps'))
import numpy as np
import pandas as pd

TRACKS = ('consensus', 'helas3', 'hepg2')
VIRUSES = ('HBV', 'HPV')
VARIANTS = ('paired_sample', 'donor_signature')
CLASSES = ('DEL', 'DUP', 'INV', 'TRA')
N, SEED, BIN = 20000, 20260911, 50000
CHECKS, SUMMARY_ROWS, DIFFERENCE_ROWS, STREAMS = [], [], [], []
COUNTERS = dict(composition_candidate_rows=0, atlas_candidate_class_counts=0,
                replayed_cluster_track_draws=0, brute_event_set_queries=0)
MAX_ERROR = dict(composition_candidates=0., composition_aggregate=0., atlas_counts=0., summaries=0.)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()


def readcsv(path): return pd.read_csv(path, float_precision='round_trip')


def check(label, actual, expected, atol=0., family='exact'):
    a, e = np.asarray(actual), np.asarray(expected)
    shape = a.shape == e.shape
    if shape and a.dtype.kind in 'biufc' and e.dtype.kind in 'biufc':
        equal = bool(np.allclose(a, e, rtol=0, atol=atol, equal_nan=True))
        finite = np.isfinite(a) & np.isfinite(e)
        error = float(np.max(np.abs(a[finite].astype(float)-e[finite].astype(float)))) if finite.any() else 0.
    else:
        equal = shape and bool(np.array_equal(a, e)); error = 0. if equal else None
    CHECKS.append(dict(check=label, family=family, passed=equal, max_abs_error=error, tolerance=atol))
    if family in MAX_ERROR and error is not None: MAX_ERROR[family] = max(MAX_ERROR[family], error)
    if not equal:
        print('FAILED', label, 'shape', a.shape, e.shape, 'error', error, flush=True)
    return equal


def axis(p):
    p = np.asarray(p)
    return p[..., 2] + p[..., 3] - (p[..., 0] + p[..., 1])


def summarize(obs, null, signed):
    x = np.asarray(null, float)
    mean = np.mean(x); qlo, qhi = np.quantile(x, [0.025, .975], method='linear')
    up = (np.sum(x >= obs) + 1) / (len(x)+1)
    low = (np.sum(x <= obs) + 1) / (len(x)+1)
    row = dict(observed=float(obs), p_upper=up, p_lower=low,
               p_two_sided=min(1., 2*min(up, low)), n_permutations=len(x))
    if signed:
        row.update(null_mean=mean, excess=obs-mean, null_sd=np.std(x, ddof=1),
                   null_reference_q025=qlo, null_reference_q975=qhi,
                   excess_reference_low=obs-qhi, excess_reference_high=obs-qlo,
                   interval_type='spatial-null reference; not donor sampling CI')
    else:
        row.update(expected=mean, fold=obs/mean if mean else np.nan,
                   null_quantile_025=qlo, null_quantile_975=qhi,
                   fold_null_low=obs/qhi if qhi else np.nan,
                   fold_null_high=obs/qlo if qlo else np.nan,
                   interval_type='observed/null spatial reference, not sampling CI')
    return row


def poolname(virus, cid, k, track): return f'{virus}_c{cid}_k{k}_{track}'


def replay(est, virus, cid, K, pools):
    rng = np.random.default_rng([SEED, ('composition', 'atlas').index(est), VIRUSES.index(virus), int(cid)])
    an = rng.integers(0, 1 if est == 'composition' else K, size=N)
    u = rng.random(N)
    ds = np.empty((3, N), dtype=np.int64)
    for k in range(1 if est == 'composition' else K):
        rows = np.flatnonzero(an == k)
        for ti, track in enumerate(TRACKS):
            p = pools[poolname(virus, cid, k, track)]
            ds[ti, rows] = p[(u[rows]*len(p)).astype(np.int64)]
    h = hashlib.sha256()
    for a in (ds, an, u):
        h.update(str(a.dtype).encode()); h.update(str(a.shape).encode()); h.update(a.tobytes())
    COUNTERS['replayed_cluster_track_draws'] += 3*N
    return ds, h.hexdigest()


def load_pools(virus):
    with np.load(STAGE/f'design/pools_{virus}.npz') as z:
        return {k:z[k] for k in z.files}


def compare_table(actual, expected, keys, label, tolerance=1e-12):
    a = actual.set_index(keys).sort_index(); e = pd.DataFrame(expected).set_index(keys).sort_index()
    check(label+' unique keys', not a.index.has_duplicates and not e.index.has_duplicates, True)
    check(label+' keys', a.index.to_numpy(), e.index.to_numpy())
    for col in e:
        if col not in a:
            check(label+' missing '+col, False, True); continue
        check(label+' '+col, a[col].to_numpy(), e[col].to_numpy(),
              tolerance if e[col].dtype.kind in 'f' else 0., 'summaries')


def endpoint_oracle(variant):
    ep = pd.read_parquet(STAGE/f'baseline/outputs/composition_breakpoints_{variant}.parquet')
    check(variant+' endpoint valid positions', (ep.pos0 >= 0).all(), True)
    check(variant+' endpoint class domain', sorted(ep.svclass.unique()), sorted(CLASSES))
    global_counts = np.array([(ep.svclass == c).sum() for c in CLASSES])
    pg = global_counts / len(ep)
    # Each physical query bin is scored by search over endpoint-bin multisets.
    # No cumulative histogram or production window implementation is used.
    maps = {}
    feat = pd.read_parquet(STAGE/'design/features.parquet')
    for chrom, fg in feat.groupby('chrom', sort=True):
        maxbin = int(fg.bin_start.max())//BIN + 2
        centers = np.arange(maxbin)
        counts = np.zeros((maxbin, 4), dtype=np.int64)
        group = ep[ep.chrom == chrom]
        for ci, cl in enumerate(CLASSES):
            bb = np.sort(group.loc[group.svclass == cl, 'pos0'].to_numpy(np.int64)//BIN)
            counts[:, ci] = np.searchsorted(bb, centers+20, side='right') - np.searchsorted(bb, centers-20, side='left')
        maps[str(chrom)] = (counts+30*pg)/(counts.sum(axis=1)[:,None]+30.)
    return ep, pg, maps


def endpoint_scores(chrom, points, delta, maps):
    bins = (np.asarray(delta)[:,None]+points[None,:])//BIN
    if np.any(bins < 0) or np.any(bins >= len(maps[chrom])):
        raise ValueError('Unqualified query outside independent physical-bin array')
    return maps[chrom][bins].mean(axis=1)


def run_composition(variant, pools, meta):
    ep, pg, maps = endpoint_oracle(variant)
    out = STAGE/'outputs'
    expected_detail = []
    observed = {}; null_axis = {}; observed_composition = {}; null_composition = {}
    with np.load(out/f'composition_candidate_compositions_{variant}.npz') as ca, np.load(out/f'composition_null_{variant}.npz') as ag:
        for virus in VIRUSES:
            co = readcsv(STAGE/f'design/common_coordinates_composition_{virus}.csv')
            C = co.original_cluster.nunique()
            ob = np.zeros(4); nl = np.zeros((3,N,4))
            m = meta[(meta.variant == variant)&(meta.estimand == 'composition')&(meta.virus == virus)].iloc[0]
            check(f'{variant} {virus} full endpoint denominator', m.n_reference_breakends, len(ep))
            check(f'{variant} {virus} fixed prior', ast.literal_eval(m.global_prior), pg)
            check(f'{variant} {virus} composition common support', [m.n_coordinates,m.n_clusters], [len(co),C])
            for j, (cid,g) in enumerate(co.groupby('original_cluster',sort=True)):
                cid = int(cid); chrom = str(g.chrom.iloc[0]); p = g.pos0.to_numpy(np.int64)
                oc = endpoint_scores(chrom,p,np.array([0]),maps)[0]
                ob += oc/C
                ds,digest = replay('composition',virus,cid,len(p),pools[virus])
                STREAMS.append(dict(estimand='composition',variant=variant,virus=virus,cluster=cid,draw_sha256=digest))
                for ti,t in enumerate(TRACKS):
                    d = pools[virus][poolname(virus,cid,0,t)]
                    pc = endpoint_scores(chrom,p,d,maps)
                    label = f'composition {variant} {virus} {cid} {t}'
                    check(label+' all candidate compositions',ca[f'{virus}_c{cid}_{t}'],pc,1e-14,'composition_candidates')
                    COUNTERS['composition_candidate_rows'] += len(pc)
                    nl[ti] += pc[np.searchsorted(d,ds[ti])]/C
                    expected_detail.append(dict(variant=variant,virus=virus,cluster=cid,n_members=len(p),track=t,
                        observed_axis=axis(oc), candidate_expected_axis=np.mean(axis(pc)),
                        **{f'observed_{cl}':oc[ci] for ci,cl in enumerate(CLASSES)}))
                if j%1000 == 0: print('oracle composition',variant,virus,j,flush=True)
            observed[virus] = float(axis(ob)); null_axis[virus] = axis(nl)
            observed_composition[virus],null_composition[virus] = ob,nl
            check(f'{variant} {virus} all observed composition',ag[virus+'_observed_composition'],ob,1e-13,'composition_aggregate')
            check(f'{variant} {virus} all 20k composition',ag[virus+'_null_composition'],nl,1e-13,'composition_aggregate')
            check(f'{variant} {virus} all 20k axis',ag[virus+'_null_axis'],null_axis[virus],1e-13,'composition_aggregate')
            for ti,t in enumerate(TRACKS):
                SUMMARY_ROWS.append(dict(table='composition',variant=variant,track=t,contrast=virus,
                    n_clusters=C,n_coordinates=len(co),p_bonferroni_3=np.nan,
                    **summarize(observed[virus],null_axis[virus][ti],True)))
        for ti,t in enumerate(TRACKS):
            h,p = null_axis['HBV'][ti],null_axis['HPV'][ti]
            oi = (observed['HBV']-h.mean())-(observed['HPV']-p.mean())
            ni = (h-h.mean())-(p-p.mean())
            check(f'{variant} {t} interaction observed',ag[t+'_interaction_observed'],oi,1e-13,'composition_aggregate')
            check(f'{variant} {t} all 20k centered interaction',ag[t+'_interaction_null'],ni,1e-13,'composition_aggregate')
            check(f'{variant} {t} interaction centered',float(np.mean(ni)),0.,1e-14)
            row = summarize(oi,ni,True)
            SUMMARY_ROWS.append(dict(table='composition',variant=variant,track=t,contrast='HBV-minus-HPV',
                n_clusters=np.nan,n_coordinates=np.nan,p_bonferroni_3=min(1.,3*row['p_two_sided']) if variant=='paired_sample' else np.nan,**row))
    compare_table(readcsv(out/f'composition_clusters_{variant}.csv'),expected_detail,
                  ['variant','virus','cluster','track'],variant+' composition cluster table')
    for contrast in (*VIRUSES,'HBV-minus-HPV'):
        for ti,t in enumerate(TRACKS[1:],1):
            if contrast in VIRUSES:
                n0,nt = null_axis[contrast][0],null_axis[contrast][ti]
            else:
                n0 = null_axis['HBV'][0]-null_axis['HPV'][0]
                nt = null_axis['HBV'][ti]-null_axis['HPV'][ti]
            # Observed axis/contrast is identical on all three tracks: track change
            # is the negative change of expected null axis, with paired covariance.
            covariance = np.cov(np.vstack([nt,n0]),ddof=1)
            mc = np.sqrt((covariance[0,0]+covariance[1,1]-2*covariance[0,1])/N)
            DIFFERENCE_ROWS.append(dict(estimand='composition',variant=variant,contrast=contrast,track=t,
                reference='consensus',metric='signed_excess',difference=n0.mean()-nt.mean(),
                monte_carlo_se=mc,uncertainty_type='MC precision of null means; not effect-change CI',
                svclass='',estimator=''))


def prepare_geometry(variant):
    data = pd.read_parquet(STAGE/f'baseline/outputs/atlas_geometry_{variant}.parquet')
    check(variant+' geometry unique chromosome event',data.duplicated(['chrom','event_id']).sum(),0)
    check(variant+' geometry ordered endpoints',bool((data.left_pos0 <= data.right_pos0).all()),True)
    idx = {}
    for (chrom,cl),g in data.groupby(['chrom','svclass'],sort=True):
        idx[(str(chrom),str(cl))] = (g.left_pos0.to_numpy(np.int64),g.right_pos0.to_numpy(np.int64),g.event_id.to_numpy())
    return idx


def event_counts(chrom,cl,points,delta,geometry):
    """Independent closed endpoint-hit union by range updates over query starts.

    Each endpoint votes for all window starts in [endpoint-L,endpoint].
    Both-endpoint-hit intersection votes are subtracted once per event.
    This is not span overlap; one SV ID contributes at most one per cluster.
    """
    q = int(np.min(points))-100000+np.asarray(delta,np.int64)
    L = int(np.max(points))-int(np.min(points))+200000
    if np.any(np.diff(q)<0): raise ValueError('Queries must be sorted')
    left,right,_ = geometry.get((chrom,cl),(np.empty(0,np.int64),np.empty(0,np.int64),[]))
    change = np.zeros(len(q)+1,np.int64)
    for endpoint in (left,right):
        lower = np.searchsorted(q,endpoint-L,side='left')
        upper = np.searchsorted(q,endpoint,side='right')
        np.add.at(change,lower,1); np.add.at(change,upper,-1)
    both = right-left <= L
    lower = np.searchsorted(q,right[both]-L,side='left')
    upper = np.searchsorted(q,left[both],side='right')
    np.add.at(change,lower,-1); np.add.at(change,upper,1)
    return np.cumsum(change)[:-1]


def brute_counts(chrom,cl,points,delta,geometry):
    left,right,ids = geometry.get((chrom,cl),(np.empty(0,np.int64),np.empty(0,np.int64),np.array([])))
    result = []
    for d in delta:
        low,high = int(np.min(points))-100000+int(d),int(np.max(points))+100000+int(d)
        mask = ((left>=low)&(left<=high))|((right>=low)&(right<=high))
        result.append(len(set(ids[mask].tolist())))
    COUNTERS['brute_event_set_queries'] += len(delta)
    return np.array(result,np.int64)


def run_atlas(variant,pools,meta):
    geometry = prepare_geometry(variant); out = STAGE/'outputs'; expected_detail = []
    for virus in VIRUSES:
        co = readcsv(STAGE/f'design/common_coordinates_atlas_{virus}.csv')
        C = co.original_cluster.nunique()
        m = meta[(meta.variant==variant)&(meta.estimand=='atlas')&(meta.virus==virus)].iloc[0]
        check(f'{variant} {virus} atlas common support',[m.n_coordinates,m.n_clusters],[len(co),C])
        ob = np.zeros(4,np.int64); total = np.zeros((3,4,N),np.int64)
        sob = np.zeros((3,4),np.int64); snull = np.zeros((3,4,N),np.int64)
        clusterids,clusterobs = [],[]; brute_seen=set()
        with np.load(out/f'atlas_candidate_counts_{variant}_{virus}.npz') as ca, np.load(out/f'atlas_null_{variant}_{virus}.npz') as ag:
            for j,(cid,g) in enumerate(co.groupby('original_cluster',sort=True)):
                cid=int(cid); chrom=str(g.chrom.iloc[0]); points=g.pos0.to_numpy(np.int64)
                oc=np.array([event_counts(chrom,c,points,np.array([0]),geometry)[0] for c in CLASSES])
                ob+=oc;clusterids.append(cid);clusterobs.append(oc)
                ds,digest=replay('atlas',virus,cid,len(points),pools[virus])
                STREAMS.append(dict(estimand='atlas',variant=variant,virus=virus,cluster=cid,draw_sha256=digest))
                do_brute = chrom not in brute_seen or len(points)>1 and (chrom,'multi') not in brute_seen
                if do_brute:
                    brute_seen.add(chrom)
                    if len(points)>1:brute_seen.add((chrom,'multi'))
                    for ci,c in enumerate(CLASSES):
                        check(f'atlas {variant} {virus} c{cid} {c} observed ID-set',oc[ci],brute_counts(chrom,c,points,[0],geometry)[0])
                for ti,t in enumerate(TRACKS):
                    d=np.unique(np.concatenate([pools[virus][poolname(virus,cid,k,t)] for k in range(len(points))]))
                    counts=np.stack([event_counts(chrom,c,points,d,geometry) for c in CLASSES],axis=1)
                    check(f'atlas {variant} {virus} c{cid} {t} union displacement',ca[f'c{cid}_{t}_delta'],d)
                    check(f'atlas {variant} {virus} c{cid} {t} all candidate counts',ca[f'c{cid}_{t}_counts'],counts,0.,'atlas_counts')
                    COUNTERS['atlas_candidate_class_counts']+=counts.size
                    if do_brute:
                        selected=np.unique([0,len(d)//2,len(d)-1])
                        for ci,c in enumerate(CLASSES):
                            check(f'atlas {variant} {virus} c{cid} {t} {c} candidate ID-set',counts[selected,ci],
                                  brute_counts(chrom,c,points,d[selected],geometry))
                    local=counts[np.searchsorted(d,ds[ti])].T
                    total[ti]+=local
                    sums=local.sum(axis=1)
                    so=oc*N>sums
                    # x > (sum-x)/(N-1), evaluated in exact integer arithmetic.
                    sn=local*(N-1)>(sums[:,None]-local)
                    sob[ti]+=so;snull[ti]+=sn
                    for ci,c in enumerate(CLASSES):
                        expected_detail.append(dict(variant=variant,virus=virus,cluster=cid,n_members=len(points),track=t,
                            svclass=c,observed=int(oc[ci]),null_mean=local[ci].mean(),null_sd=local[ci].std(ddof=1),
                            observed_sign=int(so[ci]),expected_sign=sn[ci].mean()))
                if j%800==0:print('oracle atlas',variant,virus,j,flush=True)
            expected_arrays=dict(observed=ob,null_counts=total,sign_observed=sob,sign_null=snull,
                                 original_cluster=np.array(clusterids),observed_by_cluster=np.array(clusterobs))
            for name,arr in expected_arrays.items():check(f'{variant} {virus} all atlas '+name,ag[name],arr,0.,'atlas_counts')
        for ti,t in enumerate(TRACKS):
            for ci,c in enumerate(CLASSES):
                for estimator,o,n in [('event_incidence',ob[ci],total[ti,ci]),('calibrated_sign',sob[ti,ci],snull[ti,ci])]:
                    r=summarize(o,n,False)
                    SUMMARY_ROWS.append(dict(table='atlas',variant=variant,virus=virus,track=t,svclass=c,estimator=estimator,
                        n_clusters=C,n_coordinates=len(co),p_bonferroni_48=min(1.,48*r['p_two_sided']) if variant=='paired_sample' else np.nan,**r))
        for ci,c in enumerate(CLASSES):
            for ti,t in enumerate(TRACKS[1:],1):
                for estimator in ('event_incidence','calibrated_sign'):
                    if estimator=='event_incidence':
                        n0,nt=total[0,ci].astype(float),total[ti,ci].astype(float);o=float(ob[ci])
                        m0,mt=n0.mean(),nt.mean();a,b=-o/mt**2,o/m0**2
                        cov=np.cov(np.vstack([nt,n0]),ddof=1)
                        mc=np.sqrt((a*a*cov[0,0]+b*b*cov[1,1]+2*a*b*cov[0,1])/N)
                        diff=o/mt-o/m0
                    else:
                        diff=sob[ti,ci]/snull[ti,ci].mean()-sob[0,ci]/snull[0,ci].mean();mc=np.nan
                    DIFFERENCE_ROWS.append(dict(estimand='atlas',variant=variant,contrast=virus,svclass=c,estimator=estimator,track=t,
                        reference='consensus',metric='fold',difference=diff,monte_carlo_se=mc,
                        uncertainty_type='delta-method MC precision' if estimator=='event_incidence' else 'not estimated: MC-calibrated sign is discontinuous'))
    compare_table(readcsv(out/f'atlas_clusters_{variant}.csv'),expected_detail,
                  ['variant','virus','cluster','track','svclass'],variant+' atlas cluster table')


def main():
    start=time.time()
    complete=json.loads((STAGE/'outputs/analysis_complete.json').read_text())
    check('production completed 20k',complete['n_permutations'],N)
    check('production completed fixed seed',complete['seed'],SEED)
    check('production completed track order',complete['tracks'],list(TRACKS))
    check('production completed class order',complete['class_order'],list(CLASSES))
    check('production completed freeze hash',complete['freeze_sha256'],sha(STAGE/'analysis/freeze_record.json'))
    files=[p for p in (STAGE/'outputs').iterdir() if p.is_file()]
    bound=[dict(path=str(p.relative_to(STAGE)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
    pd.DataFrame(bound).to_csv(HERE/'real_bound_scientific_outputs.csv',index=False)
    freeze=json.loads((STAGE/'analysis/freeze_record.json').read_text(encoding='utf-8'))
    for rel,h in freeze['sha256'].items():check('frozen '+rel,sha(STAGE/rel),h)
    pools={v:load_pools(v) for v in VIRUSES}
    meta=readcsv(STAGE/'outputs/reference_metadata.csv')
    for variant in VARIANTS:
        run_composition(variant,pools,meta)
        run_atlas(variant,pools,meta)
    e=[{k:v for k,v in r.items() if k!='table'} for r in SUMMARY_ROWS if r['table']=='composition']
    a=[{k:v for k,v in r.items() if k!='table'} for r in SUMMARY_ROWS if r['table']=='atlas']
    compare_table(readcsv(STAGE/'outputs/composition_summary.csv'),e,['variant','track','contrast'],'composition complete summary')
    compare_table(readcsv(STAGE/'outputs/atlas_summary.csv'),a,['variant','virus','track','svclass','estimator'],'atlas complete summary')
    actual_diff=readcsv(STAGE/'outputs/track_differences.csv').fillna({'svclass':'','estimator':''})
    compare_table(actual_diff,DIFFERENCE_ROWS,['estimand','variant','contrast','track','svclass','estimator'],'all track differences')
    compare_table(readcsv(STAGE/'outputs/draw_streams.csv'),STREAMS,['estimand','variant','virus','cluster'],'all actual draw digests')
    streams=pd.DataFrame(STREAMS)
    check('all cluster track draws shared across representations',streams.groupby(['estimand','virus','cluster']).draw_sha256.nunique().max(),1)
    check('primary composition correction exactly three',pd.DataFrame(e).p_bonferroni_3.notna().sum(),3)
    check('primary atlas correction exactly 48',pd.DataFrame(a).p_bonferroni_48.notna().sum(),48)
    for row in bound:check('output unchanged '+row['path'],sha(STAGE/row['path']),row['sha256'])
    pd.DataFrame(e).to_csv(HERE/'real_oracle_composition_summary.csv',index=False,float_format='%.17g')
    pd.DataFrame(a).to_csv(HERE/'real_oracle_atlas_summary.csv',index=False,float_format='%.17g')
    pd.DataFrame(DIFFERENCE_ROWS).to_csv(HERE/'real_oracle_track_differences.csv',index=False,float_format='%.17g')
    checks=pd.DataFrame(CHECKS);checks.to_csv(HERE/'real_checks.csv',index=False,float_format='%.17g')
    summary=dict(status='PASS' if checks.passed.all() else 'FAIL',checks=len(checks),passed=int(checks.passed.sum()),
        failed=int((~checks.passed).sum()),seconds=time.time()-start,n_permutations=N,seed=SEED,
        counters=COUNTERS,max_absolute_error=MAX_ERROR,scientific_files_bound=len(bound),
        bound_output_manifest_sha256=sha(HERE/'real_bound_scientific_outputs.csv'),
        analysis_freeze_sha256=sha(STAGE/'analysis/freeze_record.json'),
        independent_algorithm='endpoint-bin sorted-search; event-ID endpoint range-vote union and ID-set spotchecks; complete random-stream replay; integer leave-one-out sign; independent covariance MCSE',
        no_production_statistics_imported=True,originals_or_frozen_stages_written=False)
    (HERE/'real_review_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
    if not checks.passed.all():raise SystemExit(1)


if __name__=='__main__':main()
