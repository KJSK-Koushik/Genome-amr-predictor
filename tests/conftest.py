"""Shared fixtures for integrated analysis tests."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome_amr_predictor import engine
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
