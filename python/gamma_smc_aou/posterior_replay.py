"""Replay native hard-call summaries after one unchanged HMM pass."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import traceback

import numpy as np

from .decoder import run_within_decoder
from .fresh_power import atomic_json, canonical_hash, digest, read_profile

RULES = {"median": ("median", 0.5), "prob80": ("prob", 0.8), "prob90": ("prob", 0.9)}


def checked_source(item, name):
    path = Path(item["directory"]) / name
    spec = item["artifacts"][name]
    if path.stat().st_size != spec["bytes"] or digest(path) != spec["sha256"]:
        raise ValueError(f"Corrupt simulation input: {path}")
    return path


def posterior_identity(item, decoder_sha):
    return canonical_hash(
        dict(
            tree=item["artifacts"]["decoded_input.trees"]["sha256"],
            pairs=digest(Path(item["root"]) / "pairs.tsv"),
            positions=digest(Path(item["root"]) / "positions.txt"),
            decoder=decoder_sha,
            parameters={
                k: v
                for k, v in item["cfg"].items()
                if k.startswith("decoder_")
                or k
                in (
                    "mutation_rate",
                    "generation_time_years",
                    "tmrca_cutoffs_years",
                    "haplotype_pairs",
                )
            },
            input_transform="one_based",
            no_recent_probability=True,
            schema="recent-call-posteriors/v1",
        )
    )


def decode_identity(item, rule, decoder_sha):
    return canonical_hash(
        dict(
            posterior=posterior_identity(item, decoder_sha),
            rule=RULES[rule],
            replay_helper=item["replay_helper_sha256"],
            schema="recent-call-comparison/v2",
        )
    )


def cached_decode(directory, fingerprint, cfg):
    receipt = directory / "complete.json"
    if not receipt.exists():
        return False
    saved = json.loads(receipt.read_text())
    if saved["fingerprint"] != fingerprint:
        raise ValueError(f"Different decoding inputs: {directory}")
    for name, spec in saved["outputs"].items():
        path = directory / name
        if path.stat().st_size != spec["bytes"] or digest(path) != spec["sha256"]:
            raise ValueError(f"Corrupt saved decoding: {path}")
    read_profile(directory / "frac_recent.tsv", cfg)
    return True


def decode_bundle(payload):
    item, out, decoder, decoder_sha = payload
    out = Path(out)
    stage = out / "posteriors" / item["key"]
    stage.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    try:
        tree = checked_source(item, "decoded_input.trees")
        cfg = item["cfg"]
        identities = {rule: decode_identity(item, rule, decoder_sha) for rule in RULES}
        if all(
            cached_decode(out / "decoded" / rule / item["key"], identity, cfg)
            for rule, identity in identities.items()
        ):
            return dict(key=item["key"], status="reused", summaries=3)
        if (stage / "failed.json").exists():
            history = stage / "failure_history"
            history.mkdir(exist_ok=True)
            os.replace(
                stage / "failed.json", history / f"failure_{time.time_ns()}.json"
            )
        fingerprint = posterior_identity(item, decoder_sha)
        raw = stage / "posteriors.zst"
        raw_receipt = stage / "posterior_complete.json"
        if raw_receipt.exists():
            saved = json.loads(raw_receipt.read_text())
            if saved["fingerprint"] != fingerprint:
                raise ValueError("Saved posterior inputs differ")
            for name, spec in saved["outputs"].items():
                path = stage / name
                if (
                    path.stat().st_size != spec["bytes"]
                    or digest(path) != spec["sha256"]
                ):
                    raise ValueError(f"Corrupt saved posterior: {path}")
        else:
            result = run_within_decoder(
                decoder,
                tree,
                stage / "native_mean.tsv",
                raw_output=raw,
                scaled_mutation_rate=cfg["decoder_scaled_mutation_rate"],
                recombination_to_mutation_ratio=cfg[
                    "decoder_recombination_to_mutation_ratio"
                ],
                mutation_rate=cfg["mutation_rate"],
                threshold_years=cfg["tmrca_cutoffs_years"],
                generation_time=cfg["generation_time_years"],
                output_at_stride=-1,
                output_at_hets=False,
                only_within=False,
                output_positions_file=Path(item["root"]) / "positions.txt",
                pairs_file=Path(item["root"]) / "pairs.tsv",
                recent_call="mean",
                threads=1,
                cache_size=cfg["decoder_cache_size"],
                pair_block=cfg["decoder_pair_block"],
                exp10=cfg["decoder_exp10"],
                backward_alignment=cfg["decoder_backward_alignment"],
                vcf_position_transform="one_based",
                extra_args=["--no_recent_probability"],
            )
            for channel in ("stdout", "stderr"):
                (stage / f"decoder.{channel}.log").write_text(result.pop(channel))
            saved = dict(
                fingerprint=fingerprint,
                decode=result,
                outputs={
                    name: dict(
                        bytes=(stage / name).stat().st_size, sha256=digest(stage / name)
                    )
                    for name in (
                        "posteriors.zst",
                        "posteriors.zst.meta",
                        "native_mean.tsv",
                    )
                },
            )
            atomic_json(raw_receipt, saved)
        metadata = json.loads(raw.with_name(raw.name + ".meta").read_text())
        positions = np.loadtxt(Path(item["root"]) / "positions.txt", dtype=int)
        pairs = np.loadtxt(Path(item["root"]) / "pairs.tsv", dtype=int)
        if (
            metadata["num_pairs"] != cfg["haplotype_pairs"]
            or metadata["chunk_size"] != 8
            or metadata["sequence_length"] != len(positions)
            or not np.array_equal(metadata["output_positions"], positions)
            or not np.array_equal(metadata["pairs"], pairs)
        ):
            raise ValueError("Raw posterior pair/position manifest mismatch")
        helper = Path(decoder).with_name("summarize_recent_rules")
        command = [
            str(helper),
            str(raw),
            str(len(pairs)),
            str(Path(item["root"]) / "positions.txt"),
            str(cfg["decoder_scaled_mutation_rate"]),
            str(cfg["mutation_rate"]),
            str(cfg["generation_time_years"]),
            ",".join(map(str, cfg["tmrca_cutoffs_years"])),
            str(stage / "summaries"),
        ]
        tick = time.perf_counter()
        replay = subprocess.run(command, check=True, text=True, capture_output=True)
        baseline = (
            read_profile(checked_source(item, "frac_recent.tsv"), cfg)
            .iloc[:, 1:]
            .to_numpy()
        )
        native = read_profile(stage / "native_mean.tsv", cfg).iloc[:, 1:].to_numpy()
        recovered = (
            read_profile(stage / "summaries/mean.tsv", cfg).iloc[:, 1:].to_numpy()
        )
        for actual in (native, recovered):
            if not np.array_equal(
                np.rint(actual * len(pairs)), np.rint(baseline * len(pairs))
            ):
                raise ValueError(
                    "Posterior replay does not reproduce original mean counts"
                )
        extraction = dict(
            command=command,
            stdout=replay.stdout,
            stderr=replay.stderr,
            seconds=time.perf_counter() - tick,
            helper_sha256=digest(helper),
            mean_parity="Exact integer counts agree with native and original baseline",
            posterior_receipt_sha256=digest(raw_receipt),
        )
        for rule, identity in identities.items():
            directory = out / "decoded" / rule / item["key"]
            directory.mkdir(parents=True, exist_ok=True)
            read_profile(stage / f"summaries/{rule}.tsv", cfg)
            temporary = directory / "frac_recent.tsv.tmp"
            shutil.copyfile(stage / f"summaries/{rule}.tsv", temporary)
            os.replace(temporary, directory / "frac_recent.tsv")
            path = directory / "frac_recent.tsv"
            atomic_json(
                directory / "complete.json",
                dict(
                    fingerprint=identity,
                    input_receipt_sha256=item["receipt_sha256"],
                    rule=rule,
                    recent_call=RULES[rule][0],
                    probability=RULES[rule][1],
                    extraction=extraction,
                    outputs={
                        path.name: dict(bytes=path.stat().st_size, sha256=digest(path))
                    },
                ),
            )
        atomic_json(stage / "replay_complete.json", extraction)
        return dict(
            key=item["key"],
            status="completed",
            summaries=3,
            seconds=time.perf_counter() - started,
        )
    except Exception:
        error = traceback.format_exc()
        atomic_json(stage / "failed.json", dict(error=error))
        return dict(key=item["key"], status="failed", error=error)
