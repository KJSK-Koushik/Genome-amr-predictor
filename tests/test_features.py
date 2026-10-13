from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome_amr_predictor.detection import parse_amrfinder
from genome_amr_predictor.features import build_features, detector_features


def test_features_exclude_non_amr_and_do_not_invent_mechanisms():
    hits = parse_amrfinder(
        Path(__file__).parent / "fixtures/amrfinder_v4.2.7_dna.expected.tsv", "fixture"
    )
    features = detector_features(hits)
    assert "mutation:pmrB_C84R" in features
    assert "gene:blaTEM-156" in features
    assert not any("stx" in key or "mechanism" in key for key in features)
    assert not any(key.startswith(("sample", "project", "label", "laboratory")) for key in features)
    assert detector_features([]) == {}


def test_duplicate_hits_are_not_duplicate_loci():
    hits = parse_amrfinder(
        Path(__file__).parent / "fixtures/amrfinder_v4.2.7_dna.expected.tsv", "fixture"
    )
    hit = next(hit for hit in hits if hit["element_symbol"] == "blaTEM-156")
    features = detector_features([hit, hit])
    assert features["gene:blaTEM-156"] == 1
    assert features["loci:blaTEM-156"] == 1
    other = {**hit, "start": hit["start"] + 1000, "stop": hit["stop"] + 1000}
    assert detector_features([hit, other])["loci:blaTEM-156"] == 2


def test_missing_detection_is_not_an_absent_gene(tmp_path):
    assemblies = tmp_path / "assemblies.csv"
    assemblies.write_text("biosample_accession\nSAMN_fixture\n")
    with pytest.raises(FileNotFoundError):
        build_features(assemblies, tmp_path / "detections", tmp_path / "features")


def test_failed_detection_blocks_features(tmp_path):
    assemblies = tmp_path / "assemblies.csv"
    assemblies.write_text("biosample_accession\nSAMN_fixture\n")
    run = tmp_path / "detections/SAMN_fixture"
    run.mkdir(parents=True)
    (run / "run_manifest.json").write_text(json.dumps({"status": "failed"}))
    (run / "detections.json").write_text(json.dumps({"status": "failed"}))
    with pytest.raises(ValueError, match="Unsuccessful"):
        build_features(assemblies, run.parent, tmp_path / "features")
