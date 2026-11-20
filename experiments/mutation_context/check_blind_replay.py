"""Verify immutable solver artifacts and a byte-identical prediction replay."""
import argparse
import json
from pathlib import Path

import benchmark as b


def check(frozen, replay, entrypoint, python_version, requirements='requirements.txt'):
    frozen, replay = Path(frozen), Path(replay)
    manifest = json.loads((frozen/'submission_manifest.json').read_text())
    for name, expected in manifest['source_and_report_files'].items():
        if b.sha256(frozen/'agent_submission'/name) != expected:
            raise ValueError(f'Frozen file changed: {name}')
    for name, expected in manifest['inputs'].items():
        if b.sha256(replay/name) != expected:
            raise ValueError(f'Replay input changed: {name}')
    actual = b.sha256(replay/'predictions.tsv')
    if actual != manifest['prediction_sha256']:
        raise ValueError('Replay predictions differ from frozen submission')
    result = {'replay_directory': str(replay.resolve()), 'prediction_sha256': actual,
              'byte_identical_predictions': True, 'frozen_artifact_hashes_verified': True,
              'public_input_hashes_verified': True, 'entrypoint': entrypoint,
              'python': python_version, 'requirements': f'agent_submission/{requirements}',
              'scope': 'Reproduces the submitted estimator; not a new agent discovery trial'}
    b.write_json(frozen/'replay_check.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('frozen', type=Path)
    parser.add_argument('replay', type=Path)
    parser.add_argument('--entrypoint', required=True)
    parser.add_argument('--python-version', required=True)
    parser.add_argument('--requirements', default='requirements.txt')
    args = parser.parse_args()
    print(check(args.frozen, args.replay, args.entrypoint, args.python_version, args.requirements))
