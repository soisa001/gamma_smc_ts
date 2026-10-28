"""Source annotations with conflicting schemas must not alter panel records."""
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import phase_singletons as ph
from test_phase_singletons import fixture
from test_singleton_chunks import execute, shrink_contig


@pytest.mark.parametrize("count", [1, 4])
@pytest.mark.parametrize("source_schema", ['Number=1,Type=Integer', 'Number=A,Type=Float'])
def test_conflicting_ac_keeps_both_schemas_and_original_panel(tmp_path, count, source_schema):
    scaffold, hiphase, manifest = fixture(tmp_path)
    shrink_contig(scaffold, hiphase)
    hiphase.write_text(hiphase.read_text().replace('ID=AC,Number=A,Type=Integer', f'ID=AC,{source_schema}'))
    result, records = execute(manifest, tmp_path / "out", count=count)
    assert result['final_counts'] == {'physical': 3, 'excluded': 3}
    assert result['preservation']['exact_record_preservation']
    bypos = {key[1]: row for key, row in records.items()}
    assert ph.info_value(bypos[150], 'SP_HIPHASE_AC') == '1'
    assert ph.info_value(bypos[210], 'SP_HIPHASE_AC') == '7'
    assert ph.info_value(bypos[150], 'AC') is None
    with ph.text_input(result['output']) as handle:
        vcf = ph.VCF(handle, 'chr1')
        assert any(line.startswith('##INFO=<ID=AC,Number=A,Type=Integer,') for line in vcf.meta)
        assert any(line.startswith(f'##INFO=<ID=SP_HIPHASE_AC,{source_schema},') for line in vcf.meta)
    directory = tmp_path / 'out/chr1'
    plan = directory / 'header_merge.json' if count == 1 else directory / 'chunk_000/header_merge.json'
    assert json.loads(plan.read_text())['INFO']['AC'] == 'SP_HIPHASE_AC'


def test_renaming_avoids_existing_source_and_panel_names(tmp_path):
    scaffold, hiphase, manifest = fixture(tmp_path)
    extra = '##INFO=<ID=SP_HIPHASE_AC,Number=1,Type=Integer,Description="Existing tag">\n'
    for path in (scaffold, hiphase):
        path.write_text(path.read_text().replace('#CHROM\t', extra + '#CHROM\t'))
    hiphase.write_text(hiphase.read_text().replace('ID=AC,Number=A,Type=Integer', 'ID=AC,Number=1,Type=Integer').replace('\tAC=1\t', '\tAC=1;SP_HIPHASE_AC=99\t'))
    _, records = execute(manifest, tmp_path / 'out', count=1)
    added = next(row for key, row in records.items() if key[1] == 150)
    assert ph.info_value(added, 'SP_HIPHASE_AC_2') == '1'
    assert ph.info_value(added, 'SP_HIPHASE_AC') == '99'


def test_non_gt_format_conflict_renamed_without_touching_genotypes(tmp_path):
    scaffold, hiphase, manifest = fixture(tmp_path)
    hiphase.write_text(hiphase.read_text().replace('ID=DP,Number=1,Type=Integer', 'ID=DP,Number=1,Type=Float'))
    result, records = execute(manifest, tmp_path / 'out', count=1)
    added = next(row for key, row in records.items() if key[1] == 150)
    assert added[8] == 'GT:SP_HIPHASE_DP'
    assert added[9] == '0|1:10'
    assert result['preservation']['exact_record_preservation']
