"""Package a frozen public inheritance task and a separate private Harbor verifier."""
import argparse
import csv
from pathlib import Path
import shutil

import benchmark as b
import numpy as np
from release import PUBLIC_FILES
from ten_task import circular_contexts

BASE_IMAGE = 'python:3.12-slim-bookworm@sha256:392307d22300de8b5986851a12d9176dfc0fc073e65bf6523ebd7dcbeb23564e'


def package(task, destination):
    task, destination = Path(task), Path(destination)
    if not b.manifest_valid(task, 'samples.manifest.json'):
        raise ValueError('Corrupt or missing input manifest')
    root = destination/'tasks/offspring-alleles'
    env, tests = root/'environment', root/'tests'
    public = env/'input'
    public.mkdir(parents=True, exist_ok=True)
    tests.mkdir(parents=True, exist_ok=True)
    for name in PUBLIC_FILES:
        shutil.copyfile(task/'public'/name, public/name)
    parents = b.read_fasta(task/'public/new_parent.fasta')
    children = b.read_fasta(task/'private/new_offspring.fasta')
    kernel = np.load(task/'private/kernel.npy', allow_pickle=False)
    if kernel.shape != (4**10, 4):
        raise ValueError('Expected the frozen ten-base kernel')
    ids = {name: circular_contexts(seq) for name, seq in parents.items()}
    with (public/'queries.tsv').open() as f:
        queries = list(csv.DictReader(f, delimiter='\t'))
    codes = np.array([ids[q['contig']][int(q['position'])-1] for q in queries])
    parent = np.array([parents[q['contig']][int(q['position'])-1] for q in queries], dtype=np.int8)
    child = np.array([children[q['contig']][int(q['position'])-1] for q in queries], dtype=np.int8)
    np.savez_compressed(tests/'oracle.npz', probabilities=kernel[codes], parent=parent, child=child)
    shutil.copyfile(public/'queries.tsv', tests/'queries.tsv')
    shutil.copyfile(Path(__file__).with_name('platform_grade.py'), tests/'grade.py')
    (tests/'test.sh').write_text('#!/bin/sh\nset -eu\nexec python -I /tests/grade.py\n')
    (env/'Dockerfile').write_text(f'''FROM {BASE_IMAGE}
ENV OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 PYTHONHASHSEED=0 MPLBACKEND=Agg
RUN pip install --no-cache-dir uv==0.11.16 && uv pip install --system --no-cache numpy==2.4.6 scipy==1.17.1 pandas==3.1.0 scikit-learn==1.9.1 matplotlib==3.10.8 edlib==1.3.9.post1
RUN useradd -m engineer && mkdir -p /app/output && chown engineer /app/output
COPY input/ /app/input/
WORKDIR /app/output
''')
    (tests/'Dockerfile').write_text(f'''FROM {BASE_IMAGE}
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0
RUN pip install --no-cache-dir uv==0.11.16 && uv pip install --system --no-cache numpy==2.4.6
COPY . /tests/
RUN chmod 700 /tests && chmod +x /tests/test.sh
WORKDIR /tests
''')
    (root/'task.toml').write_text('''schema_version = "1.4"
name = "offspring-alleles"
reward_type = "bounded_continuous"
submit_paths = ["/app/output/predictions.tsv", "/app/output/summary.json", "/app/output/REPORT.md", "/app/output/decision_log.md"]
artifacts = ["/app/output/predictions.tsv", "/app/output/summary.json", "/app/output/REPORT.md", "/app/output/decision_log.md"]

[metadata]
category = "statistical-genetics"
tags = ["genetics", "probabilistic-inference"]
score_direction = "higher_is_better"
platforms = ["linux/amd64"]

[agent]
timeout_sec = 3600
user = "engineer"
network_mode = "no-network"

[environment]
cpus = 4
memory_mb = 26624
storage_mb = 10240
build_timeout_sec = 1200
network_mode = "no-network"

[verifier]
timeout_sec = 600
environment_mode = "separate"
user = "root"
network_mode = "no-network"

[verifier.environment]
cpus = 1
memory_mb = 8192
storage_mb = 10240
build_timeout_sec = 1200
network_mode = "no-network"
''')
    (root/'instruction.md').write_text('''Predict offspring alleles using the six supplied files in `/app/input`. Read `/app/input/TASK.md` for the scientific task and output format. Independently develop and validate your estimates from these inputs.

Write `/app/output/predictions.tsv` in exact query order, `/app/output/summary.json` with the expected number of changed alleles, and `/app/output/REPORT.md` explaining your hypotheses, alternatives, validation, final method, and uncertainty. Keep a chronological `/app/output/decision_log.md` as you work. Keep any code you develop in `/app/output` as well. Do not modify the supplied inputs.

The workspace is offline. Python, uv, NumPy, SciPy, pandas, scikit-learn, Matplotlib and edlib are installed. You have four CPUs and 26 GiB memory. Seed random sampling with 0. Do not display plots inline. No development offspring outcomes or score feedback are available for the new parent.

The primary evaluation is a log-loss-based score for the four allele probabilities, including unchanged alleles. Reported total mutation burden is supplementary. The maximum prediction file size is 300 MB. A probability row must sum to one within 1e-9.
''')
    (destination/'.gitattributes').write_text('* -text\n*.py text eol=lf\n*.sh text eol=lf\n*.toml text eol=lf\n*.md text eol=lf\n**/Dockerfile text eol=lf\n')
    (destination/'.gitignore').write_text('__pycache__/\n.venv/\nlocal-checks/\n')
    (destination/'README.md').write_text('''# Inheritance probability evaluation

Harbor task `tasks/offspring-alleles` supplies 10 million calibration bases and one million target queries. The solver image contains only the six public inputs and general-purpose scientific packages. The private oracle and grader run in a separate verifier image; they must never be copied into the solver environment.

The requested evaluation is three independent Claude Opus 5.5 rollouts on the same fixed data. The primary raw metric is excess expected log loss per expected mutation. Platform reward is `exp(-excess_log_loss)`; zero excess loss yields reward one. Invalid or missing predictions receive zero reward. Aggregate mutation burden is reported separately. Window discovery is assessed from the submitted report and agent trajectory, not inferred from a reward threshold.

The simulation, packaging source, and earlier experiments are maintained in `soisa001/gamma_smc_ts`, under `experiments/mutation_context`. This repository holds the immutable platform task. Seed is 0, and all scaffolds are circular with complete shuffled query coverage.
''')
    cfg = {'package_version': 1, 'sample_manifest_sha256': b.sha256(task/'samples.manifest.json'),
           'source_sha256': b.sha256(__file__), 'grader_sha256': b.sha256(tests/'grade.py'),
           'public_input_sha256': {name: b.sha256(public/name) for name in PUBLIC_FILES},
           'base_image': BASE_IMAGE, 'oracle_private_only': True}
    files = [destination/name for name in ('.gitattributes', '.gitignore', 'README.md')]
    files += [root/'instruction.md', root/'task.toml', env/'Dockerfile']
    files += [public/name for name in PUBLIC_FILES]
    files += [tests/name for name in ('Dockerfile', 'grade.py', 'test.sh', 'queries.tsv', 'oracle.npz')]
    b.seal(destination, 'package.manifest.json', cfg, files)
    return {'destination': str(destination), 'task': 'tasks/offspring-alleles', 'config': cfg}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    print(package(args.task, args.destination))
