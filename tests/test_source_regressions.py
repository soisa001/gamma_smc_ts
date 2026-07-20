from pathlib import Path


ROOT = Path(__file__).parents[1]


def test_missing_genotypes_are_not_xored_as_alleles():
    source = (ROOT / "src" / "gamma_smc.h").read_text()
    assert "alleles[i] < 0 || alleles[j] < 0" in source


def test_mask_option_typo_is_fixed():
    source = (ROOT / "src" / "gamma_smc.cpp").read_text()
    assert 'vm.count("masks_per_filename")' not in source
    assert 'vm.count("masks_per_sample")' in source


def test_tree_converter_does_not_use_shell_command():
    source = (ROOT / "src" / "io.h").read_text()
    assert "std::system" not in source
    assert "tskit.load" in source
    assert "execlp" in source
