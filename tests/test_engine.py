from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome_amr_predictor import engine
from genome_amr_predictor.cli import main
from genome_amr_predictor.detection import AMRFinderRuntime, parse_amrfinder, sha256, write_evidence


@pytest.fixture
def setup(tmp_path, monkeypatch):
    config_dir = tmp_path / "configs"
    config_dir.mkdir()
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    (model_dir / "manifest.json").write_text("registry")
    executable = tmp_path / "amrfinder"
    executable.write_text("fixture")
    database = tmp_path / "database"
    database.mkdir()
    (database / "version.txt").write_text("2026-08-07.1\n")
    fasta = tmp_path / "input.fna"
    fasta.write_text(">contig1\nACGTACGTACGTACGT\n")
    config = {
        "schema_version": "genome-amr-engine-config-v1",
        "species": "Escherichia coli",
        "model_dir": "models",
        "registry_sha256": sha256(model_dir / "manifest.json"),
        "detector": {
            "executable": "amrfinder",
            "database": "database",
            "executable_sha256": sha256(executable),
            "software_version": "4.2.7",
            "database_version": "2026-08-07.1",
            "threads": 1,
            "timeout_seconds": 60,
        },
    }
    config_path = config_dir / "engine.json"
    config_path.write_text(json.dumps(config))

    def prediction(path, model_dir):
        return {
            "input_sha256": sha256(path),
            "warnings": ["External model"],
            "predictions": [
                {
                    "antibiotic": "ampicillin",
                    "prediction": "R",
                    "resistance_score": 0.9,
                    "status": "predicted",
                    "features_present": 1,
                    "features_total": 1000,
                }
            ],
        }

    monkeypatch.setattr(engine, "predict_pretrained", prediction)
    cache = tmp_path / "cache"
    cache.mkdir()
    fixture = Path("tests/fixtures/amrfinder_v4.2.7_dna.expected.tsv")
    (cache / "amrfinder.tsv").write_bytes(fixture.read_bytes())
    hits = parse_amrfinder(cache / "amrfinder.tsv", "original_sample")
    manifest = {
        "status": "completed",
        "schema_version": "amrfinder-evidence-v1",
        "sample_id": "original_sample",
        "input": {"file_sha256": sha256(fasta)},
        "software_version": "4.2.7",
        "database_version_report": "Database version: 2026-08-07.1\n",
        "organism": "Escherichia",
        "mode": "nucleotide",
        "plus": False,
        "adapter_sha256": sha256(Path(engine.__file__).with_name("detection.py")),
        "executable_sha256": sha256(executable),
        "database_files_sha256": {"version.txt": sha256(database / "version.txt")},
        "raw_output_sha256": sha256(cache / "amrfinder.tsv"),
    }
    write_evidence(cache / "detections.json", hits, manifest)
    manifest["total_hits"] = len(hits)
    manifest["evidence_sha256"] = sha256(cache / "detections.json")
    (cache / "run_manifest.json").write_text(json.dumps(manifest))
    return fasta, config_path, cache, AMRFinderRuntime(executable, database)


def test_verified_reuse_and_exports(setup, tmp_path):
    fasta, config, cache, _ = setup
    result = engine.analyze_genome(
        fasta, tmp_path / "run", "new_sample", config, detection_cache=cache
    )
    assert result["status"] == "completed"
    assert result["detection"]["reused"] is True
    assert result["detection"]["evidence"]["hits"][0]["sample_id"] == "original_sample"
    for name, digest in result["artifacts_sha256"].items():
        assert sha256(tmp_path / "run" / name) == digest
    assert "## Phenotype predictions" in (tmp_path / "run/report.md").read_text()
    assert (tmp_path / "run/input.fna").read_bytes() == fasta.read_bytes()


@pytest.mark.parametrize("change", ["input", "evidence", "raw", "database", "failed"])
def test_bad_detection_never_becomes_absence(setup, tmp_path, change):
    fasta, config, cache, runtime = setup
    if change == "input":
        fasta.write_text(">contig1\nAAAAAAAAAAAAAAA\n")
    elif change == "database":
        (runtime.database / "version.txt").write_text("wrong")
    elif change == "failed":
        path = cache / "run_manifest.json"
        value = json.loads(path.read_text())
        value["status"] = "failed"
        path.write_text(json.dumps(value))
    else:
        (cache / ("detections.json" if change == "evidence" else "amrfinder.tsv")).write_text(
            "altered"
        )
    result = engine.analyze_genome(fasta, tmp_path / "run", "sample", config, detection_cache=cache)
    assert result["status"] == "partial"
    assert result["detection"]["status"] == "failed"
    assert result["detection"]["evidence"] is None
    assert result["phenotype_prediction"]["status"] == "completed"
    assert len((tmp_path / "run/detections.csv").read_text().splitlines()) == 1


def test_model_failure_retains_detection(setup, tmp_path, monkeypatch):
    fasta, config, cache, _ = setup

    def fail(*args):
        raise ValueError("broken model")

    monkeypatch.setattr(engine, "predict_pretrained", fail)
    result = engine.analyze_genome(fasta, tmp_path / "run", "sample", config, detection_cache=cache)
    assert result["status"] == "partial"
    assert result["detection"]["status"] == "completed"
    assert result["phenotype_prediction"]["result"] is None


def test_invalid_input_skips_all_stages(setup, tmp_path, monkeypatch):
    fasta, config, _, _ = setup
    fasta.write_text(">duplicate\nACGT\n>duplicate\nACGT\n")

    def unexpected(*args, **kwargs):
        pytest.fail("Scientific stage should not execute for invalid input")

    monkeypatch.setattr(engine, "predict_pretrained", unexpected)
    monkeypatch.setattr(engine, "detect_genome", unexpected)
    result = engine.analyze_genome(fasta, tmp_path / "run", "sample", config)
    assert result["status"] == "failed"
    assert result["input"] is None


def test_existing_run_not_overwritten(setup, tmp_path):
    fasta, config, cache, _ = setup
    engine.analyze_genome(fasta, tmp_path / "run", "sample", config, detection_cache=cache)
    before = sha256(tmp_path / "run/analysis.json")
    with pytest.raises(ValueError, match="empty"):
        engine.analyze_genome(fasta, tmp_path / "run", "sample", config)
    assert sha256(tmp_path / "run/analysis.json") == before


def test_cli_success(setup, tmp_path):
    fasta, config, cache, _ = setup
    assert (
        main(
            [
                "analyze",
                "--fasta",
                str(fasta),
                "--sample-id",
                "sample",
                "--config",
                str(config),
                "--output-dir",
                str(tmp_path / "run"),
                "--detection-cache",
                str(cache),
            ]
        )
        == 0
    )
