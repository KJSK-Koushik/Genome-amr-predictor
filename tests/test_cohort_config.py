import json

import pytest

from genome_amr_predictor.cohort_config import Phase2Config


def _config() -> dict:
    return {
        "species": "Escherichia coli",
        "antibiotics": ["ampicillin"],
        "assembly_request_batch_size": 50,
        "assembly_qc": {
            "minimum_length": 10,
            "maximum_length": 100,
            "maximum_contigs": 10,
            "minimum_contig_n50": 5,
            "minimum_checkm_completeness": 95.0,
            "maximum_checkm_contamination": 5.0,
            "require_checkm": False,
            "reject_ani_mismatch": True,
        },
        "splits": {
            "seed": "test",
            "train_fraction": 0.7,
            "validation_fraction": 0.15,
            "test_fraction": 0.15,
        },
    }


def test_phase2_config_loads_and_hashes(tmp_path):
    path = tmp_path / "phase2.json"
    path.write_text(json.dumps(_config()), encoding="utf-8")

    config = Phase2Config.load(path)

    assert config.antibiotics == ("ampicillin",)
    assert len(config.sha256()) == 64


def test_phase2_config_rejects_invalid_split(tmp_path):
    value = _config()
    value["splits"]["train_fraction"] = 0.8
    path = tmp_path / "phase2.json"
    path.write_text(json.dumps(value), encoding="utf-8")

    with pytest.raises(ValueError, match="sum to 1"):
        Phase2Config.load(path)


def test_phase2_config_rejects_impossible_checkm_threshold(tmp_path):
    value = _config()
    value["assembly_qc"]["minimum_checkm_completeness"] = 101
    path = tmp_path / "phase2.json"
    path.write_text(json.dumps(value), encoding="utf-8")

    with pytest.raises(ValueError, match="between 0 and 100"):
        Phase2Config.load(path)
