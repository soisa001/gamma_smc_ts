import pandas as pd

from gamma_smc_aou.plotting import plot_truth_tmrca_relationship


def test_truth_vs_simulated_tmrca_plot_and_metrics(tmp_path):
    summaries = []
    for replicate, (tmrca, recent) in enumerate([(18_000, 0.009), (20_000, 0.007), (22_000, 0.006)]):
        path = tmp_path / f"replicate_{replicate:05d}.tsv"
        pd.DataFrame({
            "position_0based": [49_900, 90_000],
            "position_1based": [49_901, 90_001],
            "n_pairs": [2000, 2000],
            "mean_p_tmrca_lt_threshold": [recent, recent / 2],
            "mean_tmrca_generations": [tmrca, tmrca * 1.2],
        }).to_csv(path, sep="\t", index=False)
        summaries.append(path)
    output = tmp_path / "plots"
    metrics = plot_truth_tmrca_relationship(
        summaries, sequence_length=100_000, effective_size=10_000,
        threshold_generations=150, output_dir=output,
    )
    assert metrics["n_replicates"] == 3
    assert metrics["random_locus_expected_mean_tmrca_generations"] == 20_000
    assert abs(metrics["simulated_mean_tmrca_generations"] - 20_000) < 1e-12
    assert (output / "truth_vs_simulated_tmrca.png").exists()
    assert (output / "truth_tmrca_relationship.tsv").exists()
