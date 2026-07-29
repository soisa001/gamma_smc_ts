import re
from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_missing_genotypes_are_not_xored_as_alleles():
    # Genotypes are bit-packed, so the guard is a missing-plane bit test rather
    # than the old `alleles[i] < 0 || alleles[j] < 0`. What must hold is that a
    # missing call at either haplotype yields HOM_STRETCH and never reaches the
    # allele comparison.
    source = (ROOT / "src" / "gamma_smc.h").read_text()
    assert "SiteMatrix::test_bit(missing, i) || SiteMatrix::test_bit(missing, j)" in source
    assert "cur_ptr[lane] = HOM_STRETCH;" in source


def test_missing_genotypes_are_recorded_when_reading():
    # The missing plane is only meaningful if the reader actually fills it.
    source = (ROOT / "src" / "io.h").read_text()
    assert "bcf_gt_is_missing(ptr[j])" in source
    assert "SiteMatrix::set_bit(scratch_missing.data(), haplotype_base + j)" in source


def test_tree_converter_does_not_use_shell_command():
    source = (ROOT / "src" / "io.h").read_text()
    assert "std::system" not in source
    assert "tskit.load" in source
    assert "execlp" in source


def test_flow_field_hot_path_has_no_shared_scratch():
    # at_flat_vectorized used to write class members, which made a single
    # FlowFieldCache unusable from more than one thread. Marking it const is
    # what lets the compiler prove the fix holds.
    source = (ROOT / "src" / "flow_field.h").read_text()
    assert "bool forward) const {" in source
    assert "const float* const* tables" in source


def test_masks_are_merged_before_intersection():
    # Overlapping BED intervals would double-count n_called, drive n_missing
    # negative, and index the cached tables before their base.
    source = (ROOT / "src" / "io.h").read_text()
    assert "std::sort(global_mask.begin(), global_mask.end());" in source
    assert "global_mask.swap(merged);" in source


def test_segment_intersection_does_not_read_past_the_last_segment():
    source = (ROOT / "src" / "data_processor.h").read_text()
    assert "if (k >= _n_segments) {" in source


def test_openmp_is_enabled_in_the_build():
    makefile = (ROOT / "Makefile").read_text()
    assert "-fopenmp" in makefile
    assert "$(OPENMP_FLAGS)" in makefile


def test_numerical_corrections_are_on_by_default():
    # The C++ defaults, the Python wrapper defaults and the documentation all
    # have to agree; they drifted apart once already when the defaults flipped.
    source = (ROOT / "src" / "gamma_smc.cpp").read_text()
    assert '("exp10", "accurate (default)' in source
    assert '("backward_alignment", "fixed (default)' in source

    header = (ROOT / "src" / "gamma_smc.h").read_text()
    assert "bool _accurate_exp10 = true;" in header
    assert "bool _fix_backward_alignment = true;" in header

    cli = (ROOT / "python" / "gamma_smc_aou" / "cli.py").read_text()
    assert 'choices=["accurate", "fast"], default="accurate"' in cli
    assert 'choices=["fixed", "legacy"], default="fixed"' in cli

    decoder = (ROOT / "python" / "gamma_smc_aou" / "decoder.py").read_text()
    assert 'exp10: str = "accurate"' in decoder
    assert 'backward_alignment: str = "fixed"' in decoder


def test_docs_do_not_call_the_corrections_opt_in():
    for name in ("README.md", "AOU_WORKFLOW.md"):
        text = (ROOT / name).read_text()
        assert "opt-in numerical correction" not in text, name
        # The flag was renamed; the old spelling would silently do nothing.
        assert "--accurate_exp10" not in text, name
        assert "--accurate-exp10" not in text, name


def test_manifest_header_fields_are_documented():
    # The manifest is the resume path; a field added to the writer without a
    # doc entry is the kind of thing that only surfaces months later.
    writer = (ROOT / "src" / "pair_manifest.h").read_text()
    # The C++ source contains a literal backslash-t, not a tab character.
    emitted = set(re.findall(r'"# ([a-z_]+)\\t', writer))
    assert len(emitted) >= 10, f"regex stopped matching the writer: {emitted}"
    documented = (ROOT / "AOU_WORKFLOW.md").read_text()
    for field in emitted:
        assert f"# {field}\t" in documented, field


def test_generation_time_defaults_to_25_everywhere():
    # 4500 years is 180 generations at 25 years/generation, which is the
    # statistic this branch exists to compute. A surface still defaulting to 30
    # would silently make it 150 generations instead, and nulls decoded through
    # one surface would not be comparable with observed data decoded through
    # another.
    expected = [
        ("src/gamma_smc.cpp", '("generation_time", "Generation time in years", '
                              'cxxopts::value<double>()->default_value("25"))'),
        ("python/gamma_smc_aou/decoder.py", "generation_time: float = 25,"),
        ("python/gamma_smc_aou/simulation.py", "generation_time: float = 25"),
        ("python/gamma_smc_aou/container_decoder.py", "generation_time: float = 25,"),
    ]
    for relative, needle in expected:
        assert needle in (ROOT / relative).read_text(), relative

    cli = (ROOT / "python" / "gamma_smc_aou" / "cli.py").read_text()
    defaults = re.findall(r'"--generation-time", type=float, default=(\d+)', cli)
    assert defaults, "the generation-time arguments moved; this check went blind"
    assert set(defaults) == {"25"}, defaults


def test_reference_rates_are_the_defaults_everywhere():
    # Fixed rates, not per-file estimates: coalescent time is measured in units
    # of 2Ne = theta/(2*mu), so a per-file theta rescales the time axis and makes
    # P(T<t) incomparable between datasets and between observed and null.
    # theta/rho match the Gamma-SMC paper's 1000 Genomes analysis.
    source = (ROOT / "src" / "gamma_smc.cpp").read_text()
    assert "const float default_scaled_mutation_rate = 0.00075f;" in source
    assert "const float default_scaled_recombination_rate = 0.0006f;" in source
    assert "const double default_unscaled_mutation_rate = 1.29e-8;" in source
    # rho/theta must stay self-consistent with the ratio the wrappers pass.
    assert abs(0.0006 / 0.00075 - 0.8) < 1e-12

    # Attaching cxxopts default_value to these would break presence detection,
    # because count() stays 0 for an option that only received its default.
    for option in ("scaled_mutation_rate", "scaled_recombination_rate",
                   "unscaled_mutation_rate"):
        line = next(l for l in source.splitlines() if f'("{option}"' in l
                    or f',{option}"' in l or f'"m,{option}"' in l
                    or f'"r,{option}"' in l)
        assert "default_value" not in line, option

    decoder = (ROOT / "python" / "gamma_smc_aou" / "decoder.py").read_text()
    assert "scaled_mutation_rate: float = 0.00075," in decoder
    assert "recombination_to_mutation_ratio: float = 0.8," in decoder
    assert "mutation_rate: float = 1.29e-8," in decoder

    cli = (ROOT / "python" / "gamma_smc_aou" / "cli.py").read_text()
    assert '"--theta", type=float, required=True' not in cli
    assert '"--mutation-rate", type=float, required=True' not in cli
