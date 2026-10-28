"""Repeated allele keys must preserve the panel and never invent source phase."""
import csv
import gzip
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import phase_singletons as ph
from test_phase_singletons import fixture
from test_singleton_chunks import execute


def duplicate(path, pos, transform=lambda line: line):
    lines = path.read_text().splitlines(keepends=True)
    index = next(i for i, line in enumerate(lines) if line.startswith(f"chr1\t{pos}\t"))
    lines.insert(index + 1, transform(lines[index]))
    path.write_text("".join(lines))


def audit_rows(path):
    with gzip.open(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


@pytest.mark.parametrize('count', [1, 4])
def test_repeated_sv_keys_ignored_in_source_and_preserved_in_panel(tmp_path, count):
    scaffold, hiphase, manifest = fixture(tmp_path)
    duplicate(hiphase, 230, lambda line: line.replace('END=290', 'END=310'))
    duplicate(scaffold, 230, lambda line: line.replace('END=290', 'END=310').replace('original_sv', 'other_sv'))
    duplicate(scaffold, 230)
    result, added = execute(manifest, tmp_path / 'out', count=count)
    assert {key[1] for key in added} == {150, 160, 210}
    assert result['preservation']['preserved_records'] == 7
    with ph.text_input(result['output']) as handle:
        svs = [raw for row, raw in ph.VCF(handle, 'chr1').records() if row[1] == '230']
    with ph.text_input(scaffold) as handle:
        assert svs == [raw for row, raw in ph.VCF(handle, 'chr1').records() if row[1] == '230']
    duplicates = audit_rows(tmp_path / 'out/chr1/duplicate_sites.tsv.gz')
    assert {r['action'] for r in duplicates} == {'PANEL_ROWS_PRESERVED', 'IGNORED_NON_SNV'}


@pytest.mark.parametrize('count', [1, 4])
def test_concordant_copies_add_once_and_do_not_inflate_anchor_votes(tmp_path, count):
    scaffold, hiphase, manifest = fixture(tmp_path)
    duplicate(hiphase, 150, lambda line: line.replace('v150', 'repeat150'))
    duplicate(hiphase, 100)
    duplicate(scaffold, 100)
    result, added = execute(manifest, tmp_path / 'out', count=count)
    assert result['final_counts'] == {'physical': 3, 'excluded': 3}
    target = next(row for key, row in added.items() if key[1] == 150)
    assert ph.info_value(target, 'SP_ANCHORS') == '0,2'
    assert result['preservation']['preserved_records'] == 6


@pytest.mark.parametrize('mutation', [
    lambda line: line.replace('1|0:10:10', '0|1:10:10'),
    lambda line: line.replace('1|0:10:10', '1|0:99:10'),
    lambda line: line.replace('1|0:10:10', '0|0:.:10').replace('0|0:.:13', '0|1:40:13'),
    lambda line: line.replace('1|0:10:10', '0|0:.:10'),
])
def test_conflicting_singleton_copies_excluded_without_picking_carrier(tmp_path, mutation):
    _, hiphase, manifest = fixture(tmp_path)
    duplicate(hiphase, 150, mutation)
    result, added = execute(manifest, tmp_path / 'out', count=4)
    assert {key[1] for key in added} == {160, 210}
    assert result['final_counts'] == {'physical': 2, 'excluded': 4}
    excluded = next(r for r in audit_rows(tmp_path / 'out/chr1/excluded_singletons.tsv.gz') if r['pos'] == '150')
    assert excluded['reason'] == 'CONFLICTING_DUPLICATE_SNV'
    assert excluded['sample'] == excluded['source_ps'] == '.'


@pytest.mark.parametrize('input_name,mutation', [
    ('hiphase', lambda line: line.replace('1|0:10:10', '0|1:10:10')),
    ('hiphase', lambda line: line.replace('1|0:10:10', '1|0:99:10')),
    ('scaffold', lambda line: line.replace('0|1:1:10', '1|0:1:10')),
])
def test_conflicting_nearest_duplicate_anchor_is_not_skipped(tmp_path, input_name, mutation):
    scaffold, hiphase, manifest = fixture(tmp_path)
    duplicate({'scaffold': scaffold, 'hiphase': hiphase}[input_name], 100, mutation)
    _, added = execute(manifest, tmp_path / 'out', count=4)
    assert {key[1] for key in added} == {160, 210}
    rejected = next(r for r in audit_rows(tmp_path / 'out/chr1/excluded_singletons.tsv.gz') if r['pos'] == '150')
    assert rejected['reason'] == 'CONFLICTING_NEAREST_ANCHORS'
    assert rejected['left_pos'] == '100'


def test_audit_rejects_missing_original_copy_and_repeated_addition(tmp_path):
    scaffold, _, manifest = fixture(tmp_path)
    duplicate(scaffold, 230)
    result, _ = execute(manifest, tmp_path / 'out', count=1)
    output = Path(result['output'])
    with gzip.open(output, 'rt') as handle:
        original = handle.readlines()
    damaged = list(original)
    del damaged[next(i for i, line in enumerate(damaged) if line.startswith('chr1\t230\t'))]
    with gzip.open(output, 'wt') as handle:
        handle.writelines(damaged)
    with pytest.raises(ValueError, match='Scaffold records missing/reordered'):
        ph.audit_preservation(scaffold, output, 'chr1')
    damaged = list(original)
    index = next(i for i, line in enumerate(damaged) if line.startswith('chr1\t150\t'))
    damaged.insert(index, damaged[index])
    with gzip.open(output, 'wt') as handle:
        handle.writelines(damaged)
    with pytest.raises(ValueError, match='Duplicate singleton addition'):
        ph.audit_preservation(scaffold, output, 'chr1')


def test_filtered_duplicate_does_not_override_passing_source_evidence(tmp_path):
    _, hiphase, manifest = fixture(tmp_path)
    duplicate(hiphase, 150, lambda line: line.replace('\tPASS\t', '\tLowQual\t').replace('1|0:10:10', '0|1:99:10'))
    result, added = execute(manifest, tmp_path / 'out', count=4)
    assert {key[1] for key in added} == {150, 160, 210}
    assert result['final_counts'] == {'physical': 3, 'excluded': 3}
