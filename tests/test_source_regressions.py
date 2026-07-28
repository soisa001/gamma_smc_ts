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


def test_mask_option_typo_is_fixed():
    source = (ROOT / "src" / "gamma_smc.cpp").read_text()
    assert 'vm.count("masks_per_filename")' not in source
    assert 'vm.count("masks_per_sample")' in source


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
