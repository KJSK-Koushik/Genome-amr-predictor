from __future__ import annotations

import json

import pytest

from genome_amr_predictor.detection import sha256
from genome_amr_predictor.pretrained import (
    canonical,
    extract_presence,
    model_score,
    predict_pretrained,
)


def test_strands_boundaries_ambiguity_and_line_wrapping(tmp_path):
    path = tmp_path / "input.fna"
    word = "AAAACCCCGGGGT"
    reverse = "ACCCCGGGGTTTT"
    path.write_text(f">a\n{word[:7]}\n{word[7:]}\n>b\n{reverse}\n>c\nAAAAANCCCCGGGGT\n")
    found, qc = extract_presence(path, {word, reverse, "GGGGTACCCCGGG"})
    assert found == {canonical(word)}
    assert qc["contigs"] == 3
    assert qc["ambiguous_bases"] == 1


@pytest.mark.parametrize("text", ["ACGT", ">a\nXYZ", ">a\nACGT\n>a\nACGT", ">a\n"])
def test_bad_fasta_rejected(tmp_path, text):
    path = tmp_path / "bad.fna"
    path.write_text(text)
    with pytest.raises(ValueError):
        extract_presence(path, {"A" * 13})


def test_tree_scores_and_logistic():
    tree = {
        "left": [1, -1, -1],
        "right": [2, -1, -1],
        "feature": [0, -2, -2],
        "threshold": [0.5, -2, -2],
        "p_r": [0.5, 0.1, 0.9],
    }
    model = {"kind": "random_forest", "trees": [tree]}
    assert model_score(model, [0]) == pytest.approx(0.1)
    assert model_score(model, [1]) == pytest.approx(0.9)
    assert (
        model_score({"kind": "logistic_regression", "weights": [1000], "intercept": -500}, [0])
        < 1e-100
    )


def test_qc_abstention_unsupported_and_checksum(tmp_path):
    path = tmp_path / "tiny.fna"
    path.write_text(">a\n" + "A" * 20)
    artifact = tmp_path / "ampicillin.json"
    artifact.write_text(
        json.dumps(
            {
                "antibiotic": "ampicillin",
                "kind": "logistic_regression",
                "weights": [1],
                "intercept": 0,
                "vocabulary": ["A" * 13],
            }
        )
    )
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "schema": "hzi-pretrained-v1",
                "artifacts": {artifact.name: sha256(artifact)},
                "limitations": ["external"],
            }
        )
    )
    result = predict_pretrained(path, tmp_path)
    assert result["predictions"][0]["status"] == "withheld_qc"
    assert result["predictions"][0]["prediction"] is None
    assert all(row["status"] == "unsupported" for row in result["predictions"][1:])
    artifact.write_text("changed")
    with pytest.raises(ValueError, match="checksum"):
        predict_pretrained(path, tmp_path)
