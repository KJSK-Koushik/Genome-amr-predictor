from pathlib import Path

from genome_amr_predictor.data.file_import import read_ast_file

FIXTURES = Path(__file__).parent / "fixtures"


def test_read_ast_browser_style_csv():
    records = read_ast_file(FIXTURES / "ast_export.csv")

    assert len(records) == 4
    assert records[0].biosample_accession == "SAMN10000001"
    assert records[0].antibiotic == "ampicillin"
    assert records[0].binary_label == 1
    assert records[1].binary_label == 0
    assert records[2].phenotype_normalized == "intermediate"


def test_hash_prefixed_accession_species_filter_and_malformed_quote_repair(tmp_path):
    path = tmp_path / "download.csv"
    path.write_text(
        "#BioSample,Scientific name,Isolation source,Isolate,Antibiotic,"
        "Resistance phenotype,Measurement sign,MIC (mg/L),Testing standard\n"
        'SAMN1,"Escherichia coli",""Food\\, retail"",EC-1,ampicillin,'
        "resistant,==,32,CLSI\n"
        "SAMN2,Salmonella enterica,food,SAL-1,ampicillin,susceptible,==,2,CLSI\n",
        encoding="utf-8",
    )

    records = read_ast_file(path, species="Escherichia coli")

    assert len(records) == 1
    assert records[0].biosample_accession == "SAMN1"
    assert records[0].isolation_source == "Food\\, retail"
    assert records[0].antibiotic == "ampicillin"
