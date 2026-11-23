# Ten-million-base platform experiment

On November 22, 2026 (Chicago time), the launch retry started **three instances** in Hyperfocal run `run_01a205f5-54b2-7132-911f-9fbcc93a7223`. The start receipt confirmed the requested model and pinned task commit. The environment run list was empty immediately before launch. At the initial status check the worker was running and the problem was pending.

## Frozen protocol

- Environment: `mew-popgen-inheritance-10mb`; task: `offspring-alleles`.
- Repository: https://github.com/Hyperfocal-AI-Environments/mew-popgen-inheritance-10mb
- Task commit: `e2b25ea2f1f35c5407fa687e45ae6e10215bfa3c`.
- Requested agent: `claude-code`, model `claude-opus-5-5`, thinking level `high`; three independent rollouts on the same fixed dataset, concurrency three. These are agent replicates, not simulation replicates.
- Compute: Modal, four reserved CPUs and 26,624 MB RAM per rollout; agent timeout 3,600 seconds, verifier timeout 600 seconds, maximum platform runtime two hours. Shared Claude subscription billing, with no API-key fallback.
- Calibration: 10,000,000 bases, 49,683 observed mutations. Target: 1,000,000 queries. Seed 0; 50,000-base circular scaffolds; complete shuffled target query coverage.
- Retains the preceding arbitrary 10-mer truth, offsets -5 through +4, mean mutation probability 0.005 and log-rate standard deviation 1.25. This is an accelerated mock mutation process, not a human germline rate. The target and kernel are identical to the preceding 4M/64M comparisons.
- Public release: `../../releases/task_007.zip`, SHA256 `bb7f8fb2ef8775d12560870abb5cb8316f7bc418173bc7e92c6621f6608b1194`.
- Sample manifest SHA256: `e81baec8d32510e33d8d7d50c5c5c31462038d142b33956b3e6ecf5af9f1e0d0`.

The solver image contains only six public inputs and scientific software. The prompt does not mention context dependence or window length. Private truth and grading code are mounted only in the separate verifier environment. The retained template task is not selected.

Primary metric: expected excess log loss divided by the expected number of mutations; lower is better. Platform reward is exp(-this metric). Report observed log loss and burden error separately. Assess discovery from each report and trajectory: exact asymmetric window, inclusion of the true window within a larger model, or no discovery. Do not infer discovery solely from reward. Preserve all three outcomes, including invalid submissions and timeouts, without selecting the best run.

## Validation and reference results

Five relevant grader tests passed. A full one-million-row submission from the earlier constant-rate solver was regraded against this task's private target: KL per expected mutation 0.7851081811910477 and reward 0.45607036088196673, matching the previous grader. Public input hashes and package manifests were checked. The platform read-only launch plan resolves the pinned commit and requested compute. Docker images have not been built locally: the available WSL Docker daemon requires permissions unavailable to the current user.

The author-side reference search is **not a blind agent rollout**. Its configuration and complete metrics are saved alongside this document. At 10M bases there are 9.54 observations and 0.0474 mutations per possible 10-mer on average; 98.4757% of the three alternative-base channels per context have zero observed mutations. Only 0.0060% of contexts are entirely unseen, so missing mutation events, rather than missing contexts, dominate sparsity.

| Reference model | Oracle KL per expected mutation |
| --- | ---: |
| Parent base only | 0.785323 |
| Centered 5-mer | 0.785323 |
| Centered 9-mer | 0.771465 |
| True 10-mer (-5..+4) | 0.751480 |
| Shifted 10-mer (-4..+5) | 0.771526 |
| 11-mer (-5..+5), shrunk toward true 10-mer | 0.751235 |

The true 10-mer reduces excess log loss by about 4.31% versus the parent-base model. Validation nevertheless selects the 11-mer, with a small 5.51-nat validation advantage over the true window. This fit uses hierarchical shrinkage; the 11-mer inherits the fitted 10-mer. These results show some detectable signal but do not establish that a blind agent will identify the minimal window. An arbitrary independent 10-mer table does not guarantee that a 5-mer will be the best compromise under limited data.

## Launch record

The start command below documents the launched configuration. To check this experiment use `hyperfocal run status run_01a205f5-54b2-7132-911f-9fbcc93a7223 --json`.

```bash
if [ ! -d mew-popgen-inheritance-10mb/.git ]; then
  git clone https://github.com/Hyperfocal-AI-Environments/mew-popgen-inheritance-10mb.git
fi
cd mew-popgen-inheritance-10mb
git pull --ff-only origin main
hyperfocal run start mew-popgen-inheritance-10mb \
  --task offspring-alleles --agent-type claude-code \
  --model claude-opus-5-5 --thinking-level high \
  --rollouts 3 --concurrency 3 \
  --sandbox-provider modal --cpus 4 --memory-mb 26624 --cpu-mode fixed \
  --max-hours 2 \
  --expected-commit e2b25ea2f1f35c5407fa687e45ae6e10215bfa3c \
  --yes --json
```

Grades, verifier scores and submitted reports were pending at launch.
