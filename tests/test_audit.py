from pathlib import Path

from genome_amr_predictor.config import AuditConfig, SelectionThresholds
from genome_amr_predictor.data.audit import audit_records, write_audit_outputs
from genome_amr_predictor.data.ncbi_biosample import parse_biosample_xml

FIXTURES = Path(__file__).parent / "fixtures"


def _config() -> AuditConfig:
    return AuditConfig(
        species="Escherichia coli",
        entrez_query="test query",
        candidate_antibiotics=("ampicillin", "ciprofloxacin", "meropenem"),
        selection_thresholds=SelectionThresholds(
            minimum_binary_records=2,
            minimum_resistant_records=1,
            minimum_susceptible_records=1,
            minimum_minority_fraction=0.25,
        ),
        batch_size=100,
    )


def test_audit_recommends_candidates_with_both_classes():
    records = parse_biosample_xml((FIXTURES / "biosample_antibiogram.xml").read_bytes())
    bundle = audit_records(records, _config())

    decisions = {row["antibiotic"]: row["recommendation"] for row in bundle.candidate_rows}
    assert decisions == {
        "ampicillin": "retain",
        "ciprofloxacin": "retain",
        "meropenem": "defer",
    }
    assert bundle.summary["unique_biosamples"] == 2
    assert bundle.summary["binary_records"] == 4
    assert bundle.summary["excluded_from_binary"] == 2


def test_write_audit_outputs(tmp_path):
    records = parse_biosample_xml((FIXTURES / "biosample_antibiogram.xml").read_bytes())
    bundle = audit_records(records, _config())
    manifest = {
        "query": "test query",
        "ncbi_reported_biosamples": 2,
        "processed_biosamples": 2,
    }

    write_audit_outputs(tmp_path, bundle, manifest, records=records)

    expected = {
        "audit_summary.json",
        "data_manifest.json",
        "antibiotic_summary.csv",
        "phenotype_summary.csv",
        "missingness_summary.csv",
        "candidate_recommendations.csv",
        "audit_report.md",
        "records.csv",
    }
    assert expected == {path.name for path in tmp_path.iterdir()}
