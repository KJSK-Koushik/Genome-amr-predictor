from pathlib import Path

from genome_amr_predictor.data.file_import import read_ast_file

FIXTURES = Path(__file__).parent / "fixtures"


def test_read_ast_browser_style_csv():
    records = read_ast_file(FIXTURES / "ast_export.csv")

    assert len(records) == 3
    assert records[0].biosample_accession == "SAMN10000001"
    assert records[0].antibiotic == "ampicillin"
    assert records[0].binary_label == 1
    assert records[1].binary_label == 0
    assert records[2].phenotype_normalized == "intermediate"
