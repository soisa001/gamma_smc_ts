"""Author-side verification of a replay against the frozen blind submission."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('replay', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    frozen = json.loads((root/'submission_manifest.json').read_text())
    for name, expected in frozen['source_and_report_files'].items():
        assert digest(root/'agent_submission'/name) == expected, name
    for name, expected in frozen['inputs'].items():
        assert digest(args.replay/name) == expected, name
    actual = digest(args.replay/'predictions.tsv')
    assert actual == frozen['prediction_sha256'], 'Prediction replay differs'
    result = {'replay_directory': str(args.replay.resolve()),
              'prediction_sha256': actual, 'byte_identical_predictions': True,
              'frozen_source_and_report_hashes_verified': True,
              'public_input_hashes_verified': True,
              'python': '3.13.15', 'requirements': 'agent_submission/requirements.txt',
              'entrypoint': 'solve.py, followed by verify.py',
              'scope': 'Reproduces the frozen estimator; not a repeat of agent discovery'}
    (root/'replay_check.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
