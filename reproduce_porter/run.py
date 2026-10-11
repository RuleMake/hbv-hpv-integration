"""Reproduce record-, event- and patient-weighted summaries of Porter et al. (2025)."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'outputs'
METHODS = ['breakpoint_record', 'author_event', 'patient_equal']
CONTRASTS = [('event_minus_record', 1, 0), ('patient_minus_record', 2, 0),
             ('patient_minus_event', 2, 1)]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def aggregates(frame, categories, include_fusion):
    required = ['patient', 'event_id', 'category', 'record_count']
    if any(c not in frame for c in required):
        raise ValueError('Missing normalized event fields')
    if frame[required].isna().any().any() or any(frame[c].astype(str).str.strip().eq('').any() for c in required[:3]):
        raise ValueError('Missing event identity/category')
    if frame.duplicated(['patient', 'event_id']).any():
        raise ValueError('Event keys must be qualified by patient and unique')
    m = pd.to_numeric(frame.record_count, errors='raise').to_numpy(float)
    if not np.isfinite(m).all() or (m<1).any() or not np.equal(m, np.floor(m)).all():
        raise ValueError('Every common-support event needs a positive integer record count')
    if set(frame.category)!=set(categories):
        raise ValueError('Observed categories differ from frozen configuration')
    y = np.column_stack([(frame.category==c).to_numpy(float) for c in categories])
    if not np.equal(y.sum(axis=1), 1).all():
        raise ValueError('Author categories must be an exhaustive single label per event')
    measures = ['structure:'+c for c in categories]
    if include_fusion:
        if 'fusion_transcript_detected' not in frame:
            raise ValueError('Required explicit published fusion label missing')
        labels = frame.fusion_transcript_detected.astype(str)
        if not labels.isin(['yes', 'no']).all():
            raise ValueError('Fusion comparison requires explicit yes/no on entire common support')
        y = np.column_stack([y, labels.eq('yes').to_numpy(float)])
        measures.append('published_fusion_transcript_detected')
    patients = sorted(frame.patient.unique())
    index = pd.Index(patients).get_indexer(frame.patient)
    P, K = len(patients), y.shape[1]
    record_n = np.zeros(P)
    event_n = np.zeros(P)
    record_y = np.zeros((P, K))
    event_y = np.zeros((P, K))
    np.add.at(record_n, index, m)
    np.add.at(event_n, index, 1)
    np.add.at(record_y, index, m[:, None]*y)
    np.add.at(event_y, index, y)
    patient_y = event_y/event_n[:, None]
    return dict(patients=patients, measures=measures, record_n=record_n,
                event_n=event_n, record_y=record_y, event_y=event_y, patient_y=patient_y)


def estimate(a, patient_counts):
    """Counts retain repeated bootstrap draws; output is draw × method × measure."""
    counts = np.asarray(patient_counts, dtype=float)
    if counts.ndim not in (1, 2):
        raise ValueError('Patient multiplicities must be a vector or draw-by-patient matrix')
    if counts.ndim==1:
        counts = counts[None, :]
    if (counts.shape[1]!=len(a['patients']) or not np.isfinite(counts).all()
            or (counts<0).any() or not np.equal(counts, np.floor(counts)).all()
            or (counts.sum(axis=1)<=0).any()):
        raise ValueError('Invalid patient multiplicities')
    record = counts@a['record_y']/(counts@a['record_n'])[:, None]
    event = counts@a['event_y']/(counts@a['event_n'])[:, None]
    patient = counts@a['patient_y']/counts.sum(axis=1)[:, None]
    return np.stack([record, event, patient], axis=1)


def write(frame, filename):
    frame.to_csv(OUT/filename, index=False, float_format='%.17g')


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'outputs',
                        help='Output directory; defaults to outputs beside this script.')
    parser.add_argument('--verify', action='store_true',
                        help='Compare reproduced arrays and tables with the supplied reference results.')
    args = parser.parse_args()
    OUT = args.output_dir.resolve()
    protected = [ROOT/'inputs', ROOT/'expected']
    if OUT == ROOT or any(OUT == p or p in OUT.parents for p in protected):
        parser.error('Choose an output directory other than the module root, inputs, or expected.')
    config_path = ROOT/'run_config.json'
    cfg = json.loads(config_path.read_text(encoding='utf-8'))
    source = ROOT/cfg['normalized_events_path']
    if sha(source)!=cfg['normalized_events_sha256']:
        raise RuntimeError('Input checksum mismatch')
    frame = pd.read_csv(source, keep_default_na=False)
    a = aggregates(frame, cfg['category_order'], cfg['include_published_fusion_secondary'])
    P = len(a['patients'])
    if [len(frame), P, int(a['record_n'].sum())]!=cfg['expected_event_patient_record_counts']:
        raise ValueError('Common-support denominator changed')
    OUT.mkdir(parents=True, exist_ok=True)
    point = estimate(a, np.ones(P))[0]
    rng = np.random.default_rng(cfg['seed'])
    draws = rng.multinomial(P, np.full(P, 1/P), size=cfg['bootstrap_replicates'])
    boot = estimate(a, draws)
    lo, hi = np.quantile(boot, [.025, .975], axis=0, method='linear')
    estimates = []
    for mi, method in enumerate(METHODS):
        for ki, measure in enumerate(a['measures']):
            estimates.append(dict(measure=measure, weighting=method, estimate=point[mi, ki],
                                  ci_low=lo[mi, ki], ci_high=hi[mi, ki], patients=P,
                                  events=len(frame), breakpoint_records=int(a['record_n'].sum()),
                                  interval='pointwise patient-bootstrap percentile'))
    write(pd.DataFrame(estimates), 'weight_estimates.csv')

    differences = []
    tv_rows = []
    n_struct = len(cfg['category_order'])
    loo = estimate(a, np.ones((P, P))-np.eye(P)) if P>1 else None
    for name, left, right in CONTRASTS:
        delta = boot[:, left]-boot[:, right]
        dlo, dhi = np.quantile(delta, [.025, .975], axis=0, method='linear')
        ldelta = loo[:, left]-loo[:, right] if loo is not None else None
        for ki, measure in enumerate(a['measures']):
            differences.append(dict(measure=measure, contrast=name,
                                    difference=point[left, ki]-point[right, ki],
                                    ci_low=dlo[ki], ci_high=dhi[ki],
                                    loo_min=np.min(ldelta[:, ki]) if loo is not None else np.nan,
                                    loo_max=np.max(ldelta[:, ki]) if loo is not None else np.nan,
                                    interval='paired pointwise patient-bootstrap percentile'))
        tv = .5*np.abs(point[left, :n_struct]-point[right, :n_struct]).sum()
        tv_rows.append(dict(contrast=name, total_variation=tv,
                            interpretation='descriptive redistribution across author structure labels; no null test'))
    write(pd.DataFrame(differences), 'weight_differences.csv')
    write(pd.DataFrame(tv_rows), 'composition_redistribution.csv')

    weighted = frame.copy()
    ep = frame.patient.map(frame.groupby('patient').size()).to_numpy(float)
    weighted['weight_record'] = frame.record_count.to_numpy(float)/a['record_n'].sum()
    weighted['weight_event'] = 1/len(frame)
    weighted['weight_patient'] = 1/(P*ep)
    write(weighted, 'event_weights.csv')
    rows = []
    for pi, patient in enumerate(a['patients']):
        row = dict(patient=patient, events=int(a['event_n'][pi]), records=int(a['record_n'][pi]),
                   record_weight=a['record_n'][pi]/a['record_n'].sum(),
                   event_weight=a['event_n'][pi]/a['event_n'].sum(), patient_weight=1/P)
        for ki, measure in enumerate(a['measures']):
            row[measure+':events'] = int(a['event_y'][pi, ki])
            row[measure+':records'] = int(a['record_y'][pi, ki])
        rows.append(row)
    write(pd.DataFrame(rows), 'patient_contributions.csv')
    loo_rows = []
    if loo is not None:
        for pi, patient in enumerate(a['patients']):
            for name, left, right in CONTRASTS:
                for ki, measure in enumerate(a['measures']):
                    loo_rows.append(dict(omitted_patient=patient, measure=measure, contrast=name,
                                         difference=loo[pi, left, ki]-loo[pi, right, ki]))
    write(pd.DataFrame(loo_rows), 'leave_one_patient_out.csv')
    np.save(OUT/'patient_bootstrap_counts.npy', draws.astype(np.uint16), allow_pickle=False)
    np.savez_compressed(OUT/'joint_bootstrap_estimates.npz', estimates=boot,
                        methods=np.array(METHODS), measures=np.array(a['measures']),
                        patients=np.array(a['patients']))
    run = dict(input_sha256=sha(source), config_sha256=sha(config_path), script_sha256=sha(Path(__file__)),
               seed=cfg['seed'], bootstrap_replicates=cfg['bootstrap_replicates'],
               numpy=np.__version__, pandas=pd.__version__, patients=P, events=len(frame),
               records=int(a['record_n'].sum()), categories=cfg['category_order'],
               primary_null_tests_run=0,
               interpretation='Retrospective descriptive weighting of published author-defined events and labels')
    (OUT/'run_metadata.json').write_text(json.dumps(run, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(run))
    if args.verify:
        from verify import compare_outputs
        verification = compare_outputs(OUT, ROOT/'expected')
        (OUT/'verification.json').write_text(json.dumps(verification, indent=2) + '\n', encoding='utf-8')
        print(json.dumps(verification))
        if not verification['passed']:
            raise SystemExit('Reproduction check failed; see verification.json')


if __name__=='__main__':
    main()
