import json
from pathlib import Path

from genome_amr_predictor.cli import main

FIXTURES = Path(__file__).parent / "fixtures"


def test_audit_file_cli(tmp_path):
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "species": "Escherichia coli",
                "entrez_query": "test query",
                "candidate_antibiotics": ["ampicillin", "ciprofloxacin"],
                "selection_thresholds": {
                    "minimum_binary_records": 2,
                    "minimum_resistant_records": 1,
                    "minimum_susceptible_records": 1,
                    "minimum_minority_fraction": 0.25,
                },
                "batch_size": 100,
            }
        ),
        encoding="utf-8",
    )
    output_dir = tmp_path / "audit"

    exit_code = main(
        [
            "audit-file",
            "--input",
            str(FIXTURES / "ast_export.csv"),
            "--config",
            str(config_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    assert exit_code == 0
    assert (output_dir / "audit_report.md").is_file()
    summary = json.loads((output_dir / "audit_summary.json").read_text(encoding="utf-8"))
    assert summary["ast_records"] == 3
    manifest = json.loads((output_dir / "data_manifest.json").read_text(encoding="utf-8"))
    assert manifest["species_filter"] == "Escherichia coli"
    assert len(manifest["source_file_sha256"]) == 64
