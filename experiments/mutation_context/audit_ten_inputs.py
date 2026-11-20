"""Verify that the two blind releases differ only in nested calibration size."""
import argparse
from pathlib import Path

import benchmark as b


def audit(small, large, output):
    small, large = Path(small), Path(large)
    for root in (small, large):
        if not b.manifest_valid(root, 'samples.manifest.json'):
            raise ValueError(f'Invalid data manifest: {root}')
    shared = ('public/TASK.md', 'public/metadata.json', 'public/new_parent.fasta',
              'public/queries.tsv', 'private/new_offspring.fasta',
              'private/kernel.npy', 'private/kernel_encoding.json')
    hashes = {}
    for name in shared:
        hashes[name] = b.sha256(small/name)
        if hashes[name] != b.sha256(large/name):
            raise ValueError(f'Differing shared input: {name}')
    prefixes = {}
    for name in ('calibration_parent.fasta', 'calibration_offspring.fasta'):
        data = (small/'public'/name).read_bytes()
        with (large/'public'/name).open('rb') as handle:
            if handle.read(len(data)) != data:
                raise ValueError(f'Calibration prefix differs: {name}')
        prefixes[name] = {'prefix_bytes': len(data), 'prefix_sha256': b.sha256(small/'public'/name)}
    result = {'small_task': str(small.resolve()), 'large_task': str(large.resolve()),
              'shared_files_sha256': hashes, 'exact_calibration_prefixes': prefixes,
              'source_sha256': b.sha256(__file__)}
    b.write_json(output, result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('small', type=Path)
    parser.add_argument('large', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(audit(args.small, args.large, args.output))
