import json

import msprime
import pandas as pd
import tszip

from gamma_smc_aou.simulation import SimulationConfig, simulate_replicates, within_individual_truth
from gamma_smc_aou.tree_sequence import diploid_individuals, tree_sequence_to_vcf


def diploid_ts(length=50_000):
    ancestry = msprime.sim_ancestry(
        samples=[msprime.SampleSet(4, ploidy=2)],
        population_size=1_000,
        sequence_length=length,
        recombination_rate=1e-8,
        model=msprime.StandardCoalescent(),
        random_seed=42,
    )
    return msprime.sim_mutations(
        ancestry, rate=2e-7, model=msprime.BinaryMutationModel(), random_seed=43
    )


def test_tree_sequence_diploids_and_vcf_only_segregating_sites(tmp_path):
    ts = diploid_ts()
    assert len(diploid_individuals(ts)) == 4
    source = tmp_path / "input.trees"
    vcf = tmp_path / "input.vcf"
    ts.dump(source)
    tree_sequence_to_vcf(source, vcf)
    records = [line for line in vcf.read_text().splitlines() if not line.startswith("#")]
    assert len(records) == ts.num_sites
    assert all(len(record.split("\t")) == 9 + 4 for record in records)


def test_tsz_conversion_uses_tszip_loader(tmp_path):
    ts = diploid_ts()
    source = tmp_path / "input.tsz"
    vcf = tmp_path / "input.vcf"
    tszip.compress(ts, source)
    tree_sequence_to_vcf(source, vcf)
    records = [line for line in vcf.read_text().splitlines() if not line.startswith("#")]
    assert len(records) == ts.num_sites


def test_truth_summary_has_one_row_per_segregating_site():
    ts = diploid_ts()
    summary = within_individual_truth(ts, threshold_generations=150)
    assert len(summary) == ts.num_sites
    assert summary["n_pairs"].eq(4).all()
    assert summary["mean_p_tmrca_lt_threshold"].between(0, 1).all()


def test_simulation_records_fixed_standard_coalescent(tmp_path):
    config = SimulationConfig(
        n_replicates=2, n_diploids=4, sequence_length=20_000,
        effective_size=1_000, mutation_rate=2e-7, seed=7,
    )
    manifest = simulate_replicates(config, tmp_path, workers=2)
    assert manifest["model"].eq("StandardCoalescent").all()
    assert len(list((tmp_path / "truth_summaries").glob("*.tsv"))) == 2
    saved = json.loads((tmp_path / "config.json").read_text())
    assert saved["ancestry_model"] == "StandardCoalescent"
    assert saved["workers"] == 2
    assert (tmp_path / "plots" / "simulation_qc.png").exists()
    assert all(len(pd.read_csv(path, sep="\t")) > 0 for path in (tmp_path / "truth_summaries").glob("*.tsv"))
