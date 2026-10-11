"""Check released Porter outputs against the supplied numerical reference results."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare_outputs(output_dir, expected_dir):
    output_dir, expected_dir = Path(output_dir), Path(expected_dir)
    manifest = json.loads((expected_dir/'manifest.json').read_text(encoding='utf-8'))
    results = []
    for item in manifest:
        reference = expected_dir/item['file']
        observed = output_dir/item['file']
        row = {'file': item['file'], 'passed': False}
        try:
            if sha(reference) != item['sha256']:
                raise ValueError('Reference-file checksum differs from its manifest.')
            row['byte_identical'] = sha(observed) == item['sha256']
            if reference.suffix == '.csv':
                original = pd.read_csv(reference, keep_default_na=False)
                reproduced = pd.read_csv(observed, keep_default_na=False)
                try:
                    pd.testing.assert_frame_equal(original, reproduced, check_exact=True)
                    row['values_identical'] = True
                except AssertionError:
                    pd.testing.assert_frame_equal(original, reproduced, check_exact=False,
                                                  rtol=1e-12, atol=1e-12)
                    row['values_identical'] = False
                row['rows'], row['columns'] = original.shape
            elif reference.suffix == '.npy':
                original = np.load(reference, allow_pickle=False)
                reproduced = np.load(observed, allow_pickle=False)
                np.testing.assert_array_equal(original, reproduced)
                if original.dtype != reproduced.dtype:
                    raise ValueError('Bootstrap-count data type differs.')
                row['values_identical'] = True
                row['shape'] = list(original.shape)
            elif reference.suffix == '.npz':
                with np.load(reference, allow_pickle=False) as original, np.load(observed, allow_pickle=False) as reproduced:
                    if original.files != reproduced.files:
                        raise ValueError('Array names or order differ.')
                    row['values_identical'] = True
                    for key in original.files:
                        if original[key].dtype != reproduced[key].dtype:
                            raise ValueError('Array data type differs: '+key)
                        if np.issubdtype(original[key].dtype, np.number):
                            np.testing.assert_allclose(original[key], reproduced[key],
                                                       rtol=1e-12, atol=1e-12)
                        else:
                            np.testing.assert_array_equal(original[key], reproduced[key])
                        row['values_identical'] &= bool(np.array_equal(original[key], reproduced[key]))
                    row['array_shapes'] = {key: list(original[key].shape) for key in original.files}
            else:
                raise ValueError('Unexpected reference format.')
            row['passed'] = True
        except (OSError, ValueError, AssertionError, KeyError) as exc:
            row['error'] = str(exc)
        results.append(row)
    passed = all(row['passed'] for row in results)
    return {'passed': passed,
            'numeric_tolerance': {'relative': 1e-12, 'absolute': 1e-12},
            'table_and_array_order_preserved': passed,
            'excluded_from_reference_comparison': 'run_metadata.json records the current script, configuration and package versions.',
            'files': results}
