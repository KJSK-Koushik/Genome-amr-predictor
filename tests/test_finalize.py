import csv
from pathlib import Path

import pytest

from genome_amr_predictor.cohort_config import Phase2Config
from genome_amr_predictor.data.fasta import inspect_fasta
from genome_amr_predictor.data.finalize import finalize_sequence_cohort


def _write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_finalize_removes_exact_sequences_and_reports_label_conflict(tmp_path):
    cohort = tmp_path / "cohort.csv"
    assemblies = tmp_path / "assemblies.csv"
    genomes = tmp_path / "genomes.csv"
    output = tmp_path / "final"
    _write_csv(
        cohort,
        [
            {
                "biosample_accession": "SAMN1",
                "antibiotic": "ampicillin",
                "binary_label": "R",
                "split": "train",
            },
            {
                "biosample_accession": "SAMN2",
                "antibiotic": "ampicillin",
                "binary_label": "S",
                "split": "test",
            },
        ],
    )
    _write_csv(
        assemblies,
        [
            {"biosample_accession": "SAMN1", "assembly_accession": "GCF_1"},
            {"biosample_accession": "SAMN2", "assembly_accession": "GCF_2"},
        ],
    )
    _write_csv(
        genomes,
        [
            {
                "assembly_accession": "GCF_1",
                "download_status": "verified",
                "sequence_sha256": "same",
                "path": "shard-0001/GCF_1.fna.gz",
            },
            {
                "assembly_accession": "GCF_2",
                "download_status": "verified",
                "sequence_sha256": "same",
                "path": "shard-0001/GCF_2.fna.gz",
            },
        ],
    )

    summary = finalize_sequence_cohort(cohort, assemblies, genomes, output)

    assert summary["final_unique_genomes"] == 1
    assert summary["duplicate_sequence_label_conflicts"] == 1
    assert summary["exclusion_reasons"] == {"exact_duplicate_sequence": 1}
    assert summary["final_cohort_rows"] == 0
    assert "SAMN2" in (output / "genome_exclusions.csv").read_text(encoding="utf-8")


def test_finalize_excludes_missing_and_invalid_genomes(tmp_path):
    cohort = tmp_path / "cohort.csv"
    assemblies = tmp_path / "assemblies.csv"
    genomes = tmp_path / "genomes.csv"
    _write_csv(
        cohort,
        [
            {
                "biosample_accession": "SAMN1",
                "antibiotic": "ampicillin",
                "binary_label": "R",
                "split": "train",
            },
            {
                "biosample_accession": "SAMN2",
                "antibiotic": "ampicillin",
                "binary_label": "S",
                "split": "test",
            },
        ],
    )
    _write_csv(
        assemblies,
        [
            {"biosample_accession": "SAMN1", "assembly_accession": "GCF_1"},
            {"biosample_accession": "SAMN2", "assembly_accession": "GCF_2"},
        ],
    )
    _write_csv(
        genomes,
        [
            {
                "assembly_accession": "GCF_2",
                "download_status": "invalid_fasta",
                "sequence_sha256": "",
                "path": "",
            }
        ],
    )

    summary = finalize_sequence_cohort(cohort, assemblies, genomes, tmp_path / "final")

    assert summary["final_unique_genomes"] == 0
    assert summary["exclusion_reasons"] == {
        "genome_not_downloaded": 1,
        "invalid_fasta": 1,
    }


def test_verified_genome_without_sequence_hash_is_rejected(tmp_path):
    cohort = tmp_path / "cohort.csv"
    assemblies = tmp_path / "assemblies.csv"
    genomes = tmp_path / "genomes.csv"
    _write_csv(
        cohort,
        [
            {
                "biosample_accession": "SAMN1",
                "antibiotic": "ampicillin",
                "binary_label": "1",
                "split": "train",
            }
        ],
    )
    _write_csv(
        assemblies,
        [{"biosample_accession": "SAMN1", "assembly_accession": "GCF_1"}],
    )
    _write_csv(
        genomes,
        [
            {
                "assembly_accession": "GCF_1",
                "download_status": "verified",
                "sequence_sha256": "",
                "path": "x.fna.gz",
            }
        ],
    )

    summary = finalize_sequence_cohort(cohort, assemblies, genomes, tmp_path / "final")

    assert summary["exclusion_reasons"] == {"missing_sequence_sha256": 1}


@pytest.mark.parametrize("damage", [False, True])
def test_finalization_rechecks_downloaded_sequence_and_qc(tmp_path, damage):
    fasta = tmp_path / "sample.fna"
    fasta.write_text(">sample\nACGT\n", encoding="utf-8")
    genome = {
        "assembly_accession": "GCF_1",
        "download_status": "verified",
        **inspect_fasta(fasta).as_dict(),
        "path": "sample.fna",
    }
    if damage:
        fasta.write_text(">sample\nAAAA\n", encoding="utf-8")
    cohort = tmp_path / "cohort.csv"
    assemblies = tmp_path / "assemblies.csv"
    manifest = tmp_path / "genomes.csv"
    _write_csv(
        cohort,
        [
            {
                "biosample_accession": "SAMN1",
                "antibiotic": "ampicillin",
                "binary_label": "1",
                "split": "train",
            }
        ],
    )
    _write_csv(assemblies, [{"biosample_accession": "SAMN1", "assembly_accession": "GCF_1"}])
    _write_csv(manifest, [genome])
    config = Phase2Config.load(Path("configs/phase2_cohort.json"))
    summary = finalize_sequence_cohort(
        cohort,
        assemblies,
        manifest,
        tmp_path / "final",
        genome_dir=tmp_path,
        qc=config.assembly_qc,
    )
    assert summary["final_unique_genomes"] == 0
    assert summary["sequence_files_reverified"] is True
    expected = "fasta_hash_mismatch" if damage else "fasta_genome_too_short"
    assert expected in summary["exclusion_reasons"]
