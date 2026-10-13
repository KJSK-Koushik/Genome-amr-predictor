"""Verify own saved models against frozen endpoints and saved test probabilities; never refit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import joblib
import numpy as np

from genome_amr_predictor.detection import sha256
from genome_amr_predictor.modeling import metrics, predict_bundle


def verify(args):
    manifest = json.loads((args.report_dir / "training_manifest.json").read_text())
    if manifest["status"] != "completed":
        raise ValueError("Model run incomplete")
    for key, path in (
        ("features_sha256", args.features),
        ("cohort_sha256", args.cohort),
        ("similarity_sha256", args.similarity_groups),
    ):
        if manifest[key] != sha256(path):
            raise ValueError("Frozen model input changed")
    for name, digest in manifest["report_outputs"].items():
        if sha256(args.report_dir / name) != digest:
            raise ValueError(f"Changed report output: {name}")
    decision_path = args.report_dir / "model_decisions.json"
    if sha256(decision_path) != manifest["frozen_decisions_sha256"]:
        raise ValueError("Model choices changed after test evaluation")
    decisions = json.loads(decision_path.read_text())["decisions"]
    features = {
        row["sample_id"]: row
        for line in args.features.read_text().splitlines()
        if line
        for row in [json.loads(line)]
    }
    with args.cohort.open(newline="", encoding="utf-8") as stream:
        labels = {
            (row["biosample_accession"], row["antibiotic"]): row for row in csv.DictReader(stream)
        }
    saved_metrics = json.loads((args.report_dir / "test_metrics.json").read_text())
    for antibiotic, digest in manifest["models"].items():
        artifact = args.model_dir / f"{antibiotic}.joblib"
        if sha256(artifact) != digest or digest != decisions[antibiotic]["model_sha256"]:
            raise ValueError("Model artifact changed")
        bundle = joblib.load(artifact)
        if bundle["threshold"] != decisions[antibiotic]["threshold"]:
            raise ValueError("Model threshold changed")
        with (args.report_dir / f"{antibiotic}-test-predictions.csv").open(newline="") as stream:
            predictions = list(csv.DictReader(stream))
        expected = {
            sample
            for (sample, drug), row in labels.items()
            if drug == antibiotic and row["split"] == "test"
        }
        samples = [row["biosample_accession"] for row in predictions]
        if set(samples) != expected or len(samples) != len(expected):
            raise ValueError("Saved test endpoint coverage differs")
        y = [int(labels[(sample, antibiotic)]["binary_label"]) for sample in samples]
        if y != [int(row["binary_label"]) for row in predictions]:
            raise ValueError("Saved test labels changed")
        p = predict_bundle(bundle, [features[sample][bundle["feature_set"]] for sample in samples])
        np.testing.assert_allclose(
            p, [float(row["resistance_probability"]) for row in predictions], rtol=1e-12, atol=1e-12
        )
        if metrics(y, p, bundle["threshold"]) != saved_metrics[antibiotic]["test"]:
            raise ValueError("Saved test metrics differ from reproducible inference")
    return {
        "status": "passed",
        "models": len(manifest["models"]),
        "checks": [
            "model/input/report hashes",
            "frozen choices",
            "test coverage and labels",
            "persisted model inference",
            "recomputed saved metrics",
        ],
        "refitting": False,
        "test_retuning": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("features", "cohort", "similarity-groups", "model-dir", "report-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = verify(args)
    (args.report_dir / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
