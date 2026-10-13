from __future__ import annotations

import csv
import json
from types import SimpleNamespace

import joblib
import numpy as np
import pytest

pytest.importorskip("sklearn")

from scripts.verify_phase4 import verify

from genome_amr_predictor.detection import sha256
from genome_amr_predictor.features import FEATURE_VERSION
from genome_amr_predictor.modeling import (
    fit_antibiotic,
    metrics,
    predict_bundle,
    select_threshold,
    train_models,
)


@pytest.fixture
def config():
    return {
        "seed": 1729,
        "feature_sets": ["detector", "detector_quality"],
        "candidates": {
            "logistic_regression_C": [1.0],
            "random_forest_max_depth": [3],
            "random_forest_estimators": 8,
        },
        "calibration_folds": 3,
        "mash": {"distance_thresholds": [0.005, 0.01]},
    }


def make_rows(split, groups):
    # Constructed software fixtures are never used in the real evaluation cohort.
    return [
        {
            "sample_id": f"{split}_{group}_{i}",
            "group_id": f"{split}_{group}",
            "split": split,
            "y": i % 2,
            "detector": {"gene:fixture_signal": float(i % 2)},
            "detector_quality": {
                "gene:fixture_signal": float(i % 2),
                "quality:length_mb": 5.0 + 0.1 * group,
            },
        }
        for group in range(groups)
        for i in range(6)
    ]


def test_train_only_vocabulary_and_calibration_groups(config):
    train = make_rows("train", 6)
    validation = make_rows("validation", 2)
    for row in validation:
        row["detector"]["validation_only_feature"] = 1
        row["detector_quality"]["validation_only_feature"] = 1
    bundle, decision, comparisons = fit_antibiotic(train, validation, config)
    names = bundle["pipeline"].named_steps["vectorizer"].get_feature_names_out()
    assert "validation_only_feature" not in names
    assert all(fold["project_overlap"] == 0 for fold in decision["calibration_folds"])
    assert len(comparisons) == 4
    p = predict_bundle(bundle, [row[bundle["feature_set"]] for row in validation])
    assert np.isfinite(p).all() and np.all((p >= 0) & (p <= 1))


def test_overlap_and_single_class_are_rejected(config):
    train, validation = make_rows("train", 6), make_rows("validation", 2)
    validation[0]["group_id"] = train[0]["group_id"]
    with pytest.raises(ValueError, match="overlap"):
        fit_antibiotic(train, validation, config)
    validation = make_rows("validation", 2)
    for row in validation:
        row["y"] = 0
    with pytest.raises(ValueError, match="both phenotype"):
        fit_antibiotic(train, validation, config)


def test_threshold_and_unestimable_similarity_subset():
    assert select_threshold(np.array([0, 0, 1, 1]), np.array([0.1, 0.2, 0.8, 0.9])) in {0.5, 0.8}
    assert metrics([], np.array([]))["status"] == "unestimable_empty_subset"
    assert metrics([1, 1], [0.7, 0.9])["roc_auc"] is None
    with pytest.raises(ValueError, match="finite"):
        metrics([0, 1], [float("nan"), 0.5])


def write_inputs(root, flip_test=False):
    rows = make_rows("train", 6) + make_rows("validation", 2) + make_rows("test", 2)
    feature_path = root / "features.jsonl"
    feature_path.write_text("".join(json.dumps(row) + "\n" for row in rows))
    (root / "feature_manifest.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "schema_version": FEATURE_VERSION,
                "features_sha256": sha256(feature_path),
            }
        )
    )
    cohort_path = root / "cohort.csv"
    with cohort_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "biosample_accession",
                "antibiotic",
                "binary_label",
                "split",
                "group_id",
                "testing_standard",
                "laboratory_platform",
            ],
        )
        writer.writeheader()
        for row in rows:
            y = 1 - row["y"] if flip_test and row["split"] == "test" else row["y"]
            writer.writerow(
                {
                    "biosample_accession": row["sample_id"],
                    "antibiotic": "fixture-drug",
                    "binary_label": y,
                    "split": row["split"],
                    "group_id": row["group_id"],
                    "testing_standard": "fixture",
                    "laboratory_platform": "fixture",
                }
            )
    similarity_path = root / "similarity.csv"
    with similarity_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=["biosample_accession", "mash_0.005", "mash_0.01"]
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "biosample_accession": row["sample_id"],
                    "mash_0.005": row["sample_id"],
                    "mash_0.01": "all-connected",
                }
            )
    (root / "similarity_manifest.json").write_text(
        json.dumps({"status": "completed", "groups_sha256": sha256(similarity_path)})
    )
    return feature_path, cohort_path, similarity_path


def test_test_labels_cannot_change_model_choices_or_weights(tmp_path, config):
    decisions = []
    bundles = []
    for flip in (False, True):
        root = tmp_path / str(flip)
        root.mkdir()
        paths = write_inputs(root, flip)
        config_path = root / "config.json"
        config_path.write_text(json.dumps(config))
        model_dir, reports = root / "models", root / "reports"
        manifest = train_models(*paths, config_path, model_dir, reports)
        assert manifest["status"] == "completed"
        verification = verify(
            SimpleNamespace(
                features=paths[0],
                cohort=paths[1],
                similarity_groups=paths[2],
                model_dir=model_dir,
                report_dir=reports,
            )
        )
        assert verification["status"] == "passed" and verification["refitting"] is False
        decision = json.loads((reports / "model_decisions.json").read_text())["decisions"][
            "fixture-drug"
        ]
        decision.pop("model_sha256")
        decisions.append(decision)
        bundle = joblib.load(model_dir / "fixture-drug.joblib")
        bundles.append(bundle)
        metrics_file = json.loads((reports / "test_metrics.json").read_text())["fixture-drug"]
        assert metrics_file["genomic_similarity_sensitivity"]["0.01"]["n"] == 0
        assert (reports / "fixture-drug-model-card.md").exists()
    assert decisions[0] == decisions[1]
    probe = [{"gene:fixture_signal": 0}, {"gene:fixture_signal": 1}]
    np.testing.assert_allclose(predict_bundle(bundles[0], probe), predict_bundle(bundles[1], probe))
