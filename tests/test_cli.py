import csv
import json
from pathlib import Path

from genome_amr_predictor.cli import main
from genome_amr_predictor.data.fasta import inspect_fasta

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


def test_download_genomes_writes_provenance_manifest(tmp_path, monkeypatch):
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("GCF_1\n", encoding="utf-8")
    executable = tmp_path / "datasets.exe"
    executable.write_text("placeholder", encoding="utf-8")
    output = tmp_path / "genomes"

    monkeypatch.setattr(
        "genome_amr_predictor.cli.datasets_version", lambda _path: "datasets 18.35.0"
    )
    monkeypatch.setattr(
        "genome_amr_predictor.cli.download_genome_shards",
        lambda *args, **kwargs: [{"assembly_accession": "GCF_1", "download_status": "verified"}],
    )
    output.mkdir()

    exit_code = main(
        [
            "download-genomes",
            "--accessions",
            str(accessions),
            "--output-dir",
            str(output),
            "--datasets-executable",
            str(executable),
        ]
    )

    assert exit_code == 0
    manifest = json.loads((output / "download_manifest.json").read_text(encoding="utf-8"))
    assert manifest["ncbi_datasets_version"] == "datasets 18.35.0"
    assert manifest["genomes_verified"] == 1
    assert len(manifest["accession_file_sha256"]) == 64


def test_finalize_cli_reverifies_files_and_records_output_fingerprints(tmp_path):
    fasta = tmp_path / "sample.fna"
    fasta.write_text(">software-test-fixture\nACGTACGT\n", encoding="utf-8")
    files = {
        "cohort": [
            {
                "biosample_accession": "SAMN1",
                "assembly_accession": "GCF_1",
                "antibiotic": "ampicillin",
                "binary_label": "1",
                "split": "train",
            }
        ],
        "assemblies": [{"biosample_accession": "SAMN1", "assembly_accession": "GCF_1"}],
        "genomes": [
            {
                "assembly_accession": "GCF_1",
                "download_status": "verified",
                **inspect_fasta(fasta).as_dict(),
                "path": "sample.fna",
            }
        ],
    }
    for name, rows in files.items():
        with (tmp_path / f"{name}.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    config = json.loads(Path("configs/phase2_cohort.json").read_text(encoding="utf-8"))
    config["assembly_qc"].update(minimum_length=4, maximum_length=100, minimum_contig_n50=1)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    output = tmp_path / "final"
    exit_code = main(
        [
            "finalize-cohort",
            "--cohort",
            str(tmp_path / "cohort.csv"),
            "--assemblies",
            str(tmp_path / "assemblies.csv"),
            "--genome-manifest",
            str(tmp_path / "genomes.csv"),
            "--genome-dir",
            str(tmp_path),
            "--config",
            str(config_path),
            "--output-dir",
            str(output),
        ]
    )
    assert exit_code == 0
    manifest = json.loads((output / "final_manifest.json").read_text(encoding="utf-8"))
    assert manifest["summary"]["final_unique_genomes"] == 1
    assert manifest["summary"]["sequence_files_reverified"] is True
    assert len(manifest["outputs"]["final_cohort.csv"]["sha256"]) == 64
