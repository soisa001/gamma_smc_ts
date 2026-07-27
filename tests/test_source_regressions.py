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
