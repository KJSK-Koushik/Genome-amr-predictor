"""Frozen-split research models; fitting and decisions never use test labels."""

from __future__ import annotations

import csv
import importlib.metadata
import json
from pathlib import Path

import joblib
import numpy as np
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    log_loss,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from genome_amr_predictor.detection import sha256
from genome_amr_predictor.features import FEATURE_VERSION


def metrics(y, probabilities, threshold=0.5):
    y, probabilities = np.asarray(y), np.asarray(probabilities, dtype=float)
    if len(y) == 0:
        return {"n": 0, "status": "unestimable_empty_subset"}
    if not np.isfinite(probabilities).all() or np.any((probabilities < 0) | (probabilities > 1)):
        raise ValueError("Probabilities must be finite and within 0..1")
    predicted = probabilities >= threshold
    tn, fp, fn, tp = confusion_matrix(y, predicted, labels=[0, 1]).ravel()
    two_classes = len(set(y)) == 2
    return {
        "n": len(y),
        "resistant": int(sum(y)),
        "susceptible": int(len(y) - sum(y)),
        "threshold": float(threshold),
        "confusion_matrix": [[int(tn), int(fp)], [int(fn), int(tp)]],
        "precision": float(precision_score(y, predicted, zero_division=0)),
        "recall_sensitivity": float(recall_score(y, predicted, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if tn + fp else None,
        "f1": float(f1_score(y, predicted, zero_division=0)),
        "mcc": float(matthews_corrcoef(y, predicted)) if two_classes else None,
        "balanced_accuracy": float(balanced_accuracy_score(y, predicted)) if two_classes else None,
        "roc_auc": float(roc_auc_score(y, probabilities)) if two_classes else None,
        "average_precision": float(average_precision_score(y, probabilities))
        if two_classes
        else None,
        "brier": float(brier_score_loss(y, probabilities)),
        "log_loss": float(log_loss(y, probabilities, labels=[0, 1])),
    }


def reliability(y, probabilities):
    y, probabilities = np.asarray(y), np.asarray(probabilities)
    bins = np.minimum((probabilities * 10).astype(int), 9)
    return [
        {
            "bin": i,
            "n": int(sum(bins == i)),
            "mean_probability": float(probabilities[bins == i].mean()),
            "observed_resistance": float(y[bins == i].mean()),
        }
        for i in range(10)
        if np.any(bins == i)
    ]


def select_threshold(y, probabilities):
    if len(set(y)) != 2:
        raise ValueError("Threshold selection requires both classes in validation")
    candidates = np.unique(np.r_[0.0, 0.5, 1.0, probabilities])
    return float(
        max(
            candidates,
            key=lambda threshold: (
                matthews_corrcoef(y, probabilities >= threshold),
                -abs(threshold - 0.5),
                threshold,
            ),
        )
    )


def _logit(probabilities):
    values = np.clip(probabilities, 1e-6, 1 - 1e-6)
    return np.log(values / (1 - values)).reshape(-1, 1)


def predict_bundle(bundle, dictionaries):
    raw = bundle["pipeline"].predict_proba(dictionaries)[:, 1]
    if bundle["calibrator"] is None:
        return raw
    return bundle["calibrator"].predict_proba(_logit(raw))[:, 1]


def candidates(config):
    for value in config["candidates"]["logistic_regression_C"]:
        yield (
            f"logistic_C_{value}",
            Pipeline(
                [
                    ("vectorizer", DictVectorizer(sparse=True)),
                    ("scale", StandardScaler(with_mean=False)),
                    (
                        "classifier",
                        LogisticRegression(
                            C=value,
                            max_iter=3000,
                            class_weight="balanced",
                            random_state=config["seed"],
                        ),
                    ),
                ]
            ),
        )
    for depth in config["candidates"]["random_forest_max_depth"]:
        yield (
            f"random_forest_depth_{depth}",
            Pipeline(
                [
                    ("vectorizer", DictVectorizer(sparse=True)),
                    (
                        "classifier",
                        RandomForestClassifier(
                            n_estimators=config["candidates"]["random_forest_estimators"],
                            max_depth=depth,
                            min_samples_leaf=2,
                            class_weight="balanced_subsample",
                        n_jobs=1,
                            random_state=config["seed"],
                        ),
                    ),
                ]
            ),
        )


def fit_antibiotic(train, validation, config):
    """Only training/validation records enter this function; test records are not accepted."""
    train_y = np.array([row["y"] for row in train])
    validation_y = np.array([row["y"] for row in validation])
    groups = np.array([row["group_id"] for row in train])
    if len(set(train_y)) != 2 or len(set(validation_y)) != 2:
        raise ValueError("Training and validation must each contain both phenotype classes")
    if set(groups) & {row["group_id"] for row in validation}:
        raise ValueError("Training and validation project groups overlap")
    comparisons = []
    fitted = []
    for feature_set in config["feature_sets"]:
        train_x = [row[feature_set] for row in train]
        validation_x = [row[feature_set] for row in validation]
        for name, pipeline in candidates(config):
            pipeline.fit(train_x, train_y)
            probabilities = pipeline.predict_proba(validation_x)[:, 1]
            result = metrics(validation_y, probabilities)
            comparison = {"candidate": name, "feature_set": feature_set, "validation": result}
            comparisons.append(comparison)
            fitted.append((pipeline, comparison))
    selected, decision = max(
        fitted,
        key=lambda item: (
            item[1]["validation"]["average_precision"],
            -item[1]["validation"]["brier"],
        ),
    )
    feature_set = decision["feature_set"]
    train_x = [row[feature_set] for row in train]
    validation_x = [row[feature_set] for row in validation]
    n_folds = min(config["calibration_folds"], len(set(groups)))
    if n_folds < 2:
        raise ValueError("Not enough independent training projects for calibration")
    splitter = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=config["seed"])
    oof = np.full(len(train), np.nan)
    fold_summary = []
    for fit_indices, predict_indices in splitter.split(np.zeros(len(train)), train_y, groups):
        if len(set(train_y[fit_indices])) != 2:
            raise ValueError("A calibration training fold lacks both classes")
        if set(groups[fit_indices]) & set(groups[predict_indices]):
            raise ValueError("Calibration groups overlap")
        fold_model = clone(selected).fit([train_x[i] for i in fit_indices], train_y[fit_indices])
        oof[predict_indices] = fold_model.predict_proba([train_x[i] for i in predict_indices])[:, 1]
        fold_summary.append(
            {
                "fit_rows": len(fit_indices),
                "predict_rows": len(predict_indices),
                "project_overlap": 0,
            }
        )
    if not np.isfinite(oof).all():
        raise ValueError("Missing out-of-fold calibration scores")
    calibrator = LogisticRegression(C=1e6, max_iter=1000, random_state=config["seed"])
    calibrator.fit(_logit(oof), train_y)
    raw_validation = selected.predict_proba(validation_x)[:, 1]
    calibrated_validation = calibrator.predict_proba(_logit(raw_validation))[:, 1]
    use_calibration = calibrator.coef_[0, 0] >= 0 and brier_score_loss(
        validation_y, calibrated_validation
    ) < brier_score_loss(validation_y, raw_validation)
    final_probabilities = calibrated_validation if use_calibration else raw_validation
    threshold = select_threshold(validation_y, final_probabilities)
    bundle = {
        "schema_version": FEATURE_VERSION,
        "pipeline": selected,
        "calibrator": calibrator if use_calibration else None,
        "threshold": threshold,
        "feature_set": feature_set,
        "training_prevalence": float(train_y.mean()),
        "status": "research_candidate_not_approved_for_user_confidence",
    }
    decision = {
        **decision,
        "threshold": threshold,
        "calibration": "training_group_oof_sigmoid" if use_calibration else "uncalibrated",
        "calibration_folds": fold_summary,
        "validation_raw": metrics(validation_y, raw_validation, threshold),
        "validation_sigmoid": metrics(validation_y, calibrated_validation, threshold),
        "validation_final": metrics(validation_y, final_probabilities, threshold),
        "validation_reliability": reliability(validation_y, final_probabilities),
        "training_rows": len(train),
        "validation_rows": len(validation),
    }
    return bundle, decision, comparisons


def grouped_intervals(y, probabilities, groups, threshold, seed, repeats=200):
    unique_groups = sorted(set(groups))
    if len(unique_groups) < 2:
        return {"status": "unestimable_fewer_than_two_projects"}
    rng = np.random.default_rng(seed)
    samples = {name: [] for name in ("roc_auc", "average_precision", "brier", "mcc")}
    lookup = {group: np.flatnonzero(np.asarray(groups) == group) for group in unique_groups}
    for _ in range(repeats):
        chosen = rng.choice(unique_groups, len(unique_groups), replace=True)
        indices = np.concatenate([lookup[group] for group in chosen])
        result = metrics(np.asarray(y)[indices], probabilities[indices], threshold)
        for name in samples:
            if result[name] is not None:
                samples[name].append(result[name])
    return {
        name: {
            "lower_95": float(np.quantile(values, 0.025)),
            "upper_95": float(np.quantile(values, 0.975)),
            "valid_resamples": len(values),
        }
        for name, values in samples.items()
        if values
    }


def train_models(
    features_path: Path,
    cohort_path: Path,
    similarity_path: Path,
    config_path: Path,
    model_dir: Path,
    report_dir: Path,
):
    config = json.loads(config_path.read_text())
    feature_manifest = json.loads((features_path.parent / "feature_manifest.json").read_text())
    similarity_manifest = json.loads(
        (similarity_path.parent / "similarity_manifest.json").read_text()
    )
    if (
        feature_manifest["status"] != "completed"
        or feature_manifest["schema_version"] != FEATURE_VERSION
        or feature_manifest["features_sha256"] != sha256(features_path)
        or similarity_manifest["status"] != "completed"
        or similarity_manifest["groups_sha256"] != sha256(similarity_path)
    ):
        raise ValueError("Incomplete or changed feature/similarity artifacts")
    features = [json.loads(line) for line in features_path.read_text().splitlines() if line]
    by_sample = {row["sample_id"]: row for row in features}
    if len(by_sample) != len(features):
        raise ValueError("Duplicate feature BioSample")
    with cohort_path.open(newline="", encoding="utf-8") as stream:
        labels = list(csv.DictReader(stream))
    with similarity_path.open(newline="", encoding="utf-8") as stream:
        similarity = {row["biosample_accession"]: row for row in csv.DictReader(stream)}
    if set(similarity) != set(by_sample):
        raise ValueError("Similarity groups and feature genomes differ")
    seen = set()
    endpoints = {}
    for label in labels:
        pair = (label["biosample_accession"], label["antibiotic"])
        if pair in seen or label["binary_label"] not in {"0", "1"}:
            raise ValueError("Duplicate/nonbinary endpoint")
        seen.add(pair)
        feature = by_sample[pair[0]]
        if label["split"] != feature["split"] or label["group_id"] != feature["group_id"]:
            raise ValueError("Frozen endpoint and feature split/group differ")
        endpoints.setdefault(label["antibiotic"], []).append(
            {**feature, "y": int(label["binary_label"]), "label_metadata": label}
        )
    if {row["biosample_accession"] for row in labels} != set(by_sample):
        raise ValueError("Feature and labeled genome coverage differ")
    project_splits = {}
    for feature in features:
        project_splits.setdefault(feature["group_id"], set()).add(feature["split"])
    if any(len(values) != 1 for values in project_splits.values()):
        raise ValueError("Project overlap between frozen splits")
    for path in (model_dir, report_dir):
        path.mkdir(parents=True, exist_ok=True)
        if any(path.iterdir()):
            raise ValueError("Model/report run directories must be empty")
    manifest = {
        "status": "fitting",
        "config": config,
        "config_sha256": sha256(config_path),
        "features_sha256": sha256(features_path),
        "cohort_sha256": sha256(cohort_path),
        "similarity_sha256": sha256(similarity_path),
        "source_sha256": sha256(Path(__file__)),
        "feature_provenance": feature_manifest,
        "similarity_provenance": similarity_manifest,
        "versions": {
            name: importlib.metadata.version(name)
            for name in ("scikit-learn", "numpy", "scipy", "joblib")
        },
    }
    bundles, decisions, comparisons = {}, {}, {}
    for antibiotic, rows in sorted(endpoints.items()):
        train = [row for row in rows if row["split"] == "train"]
        validation = [row for row in rows if row["split"] == "validation"]
        bundle, decision, comparison = fit_antibiotic(train, validation, config)
        bundle["antibiotic"] = antibiotic
        bundle["provenance"] = {
            key: manifest[key]
            for key in (
                "config_sha256",
                "features_sha256",
                "cohort_sha256",
                "similarity_sha256",
                "source_sha256",
                "versions",
                "feature_provenance",
                "similarity_provenance",
            )
        }
        bundles[antibiotic], decisions[antibiotic], comparisons[antibiotic] = (
            bundle,
            decision,
            comparison,
        )
        artifact = model_dir / f"{antibiotic}.joblib"
        joblib.dump(bundle, artifact)
        decision["model_sha256"] = sha256(artifact)
        print(
            f"Fitted {antibiotic}: {decision['candidate']} / {decision['feature_set']}", flush=True
        )
    # Freeze every decision before evaluating test endpoints.
    decision_path = report_dir / "model_decisions.json"
    decision_path.write_text(
        json.dumps(
            {"manifest": manifest, "decisions": decisions, "comparisons": comparisons}, indent=2
        )
        + "\n"
    )
    manifest["frozen_decisions_sha256"] = sha256(decision_path)
    evaluations = {}
    for antibiotic, bundle in bundles.items():
        rows = [row for row in endpoints[antibiotic] if row["split"] == "test"]
        y = np.array([row["y"] for row in rows])
        p = predict_bundle(bundle, [row[bundle["feature_set"]] for row in rows])
        baseline = np.full(len(y), bundle["training_prevalence"])
        result = {
            "test": metrics(y, p, bundle["threshold"]),
            "prevalence_baseline": metrics(y, baseline),
            "reliability": reliability(y, p),
            "project_bootstrap_95": grouped_intervals(
                y, p, [row["group_id"] for row in rows], bundle["threshold"], config["seed"]
            ),
            "genomic_similarity_sensitivity": {},
            "laboratory_strata": {},
        }
        for threshold in config["mash"]["distance_thresholds"]:
            field = f"mash_{threshold}"
            development = {
                similarity[sample][field]
                for sample, feature in by_sample.items()
                if feature["split"] in {"train", "validation"}
            }
            indices = [
                i
                for i, row in enumerate(rows)
                if similarity[row["sample_id"]][field] not in development
            ]
            result["genomic_similarity_sensitivity"][str(threshold)] = metrics(
                y[indices], p[indices], bundle["threshold"]
            )
        for field in ("testing_standard", "laboratory_platform"):
            strata = sorted({row["label_metadata"][field] or "not_reported" for row in rows})
            result["laboratory_strata"][field] = {
                value: metrics(y[indices], p[indices], bundle["threshold"])
                for value in strata
                for indices in [
                    [
                        i
                        for i, row in enumerate(rows)
                        if (row["label_metadata"][field] or "not_reported") == value
                    ]
                ]
            }
        result["outperforms_prevalence_baseline"] = (
            result["test"]["average_precision"] > result["prevalence_baseline"]["average_precision"]
            and result["test"]["brier"] < result["prevalence_baseline"]["brier"]
        )
        evaluations[antibiotic] = result
        with (report_dir / f"{antibiotic}-test-predictions.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "biosample_accession",
                    "binary_label",
                    "resistance_probability",
                    "threshold",
                    "predicted_label",
                ]
            )
            writer.writerows(
                (
                    row["sample_id"],
                    int(label),
                    float(probability),
                    bundle["threshold"],
                    int(probability >= bundle["threshold"]),
                )
                for row, label, probability in zip(rows, y, p, strict=True)
            )
        classifier = bundle["pipeline"].named_steps["classifier"]
        names = bundle["pipeline"].named_steps["vectorizer"].get_feature_names_out()
        weights = (
            classifier.coef_[0] if hasattr(classifier, "coef_") else classifier.feature_importances_
        )
        top = sorted(zip(names, weights, strict=True), key=lambda pair: abs(pair[1]), reverse=True)[
            :20
        ]
        card = [
            f"# {antibiotic} research model card",
            "",
            "Status: research candidate; not approved for user-facing confidence or clinical use.",
            "",
            f"Candidate: {decisions[antibiotic]['candidate']}",
            f"Feature set: {bundle['feature_set']}",
            f"Threshold: {bundle['threshold']:.6f}",
            f"Calibration: {decisions[antibiotic]['calibration']}",
            "",
            "Training and vocabulary fitting use training genomes only. Sigmoid calibration uses",
            "project-separated training out-of-fold scores. Selection uses validation only.",
            "All model choices were frozen before test evaluation. Test results are not retuned.",
            "",
            "## Test results",
            "",
            "| Measure | Value |",
            "|---|---:|",
            *[
                f"| {name} | {value} |"
                for name, value in result["test"].items()
                if name != "confusion_matrix"
            ],
            "",
            f"Confusion matrix (true/predicted S,R): {result['test']['confusion_matrix']}",
            "",
            f"Outperforms prevalence baseline on AP and Brier: "
            f"{result['outperforms_prevalence_baseline']}",
            "",
            "Full calibration, confidence intervals, similarity sensitivity, and laboratory strata",
            "are saved in [test_metrics.json](test_metrics.json). Validation candidate and quality",
            "ablation comparisons are in [model_decisions.json](model_decisions.json).",
            "Provenance/hashes: [training_manifest.json](training_manifest.json).",
            "",
            "## Influential features",
            "",
            "| Feature | Weight/importance |",
            "|---|---:|",
        ]
        card.extend(f"| {name} | {float(weight):.6g} |" for name, weight in top)
        card.extend(
            [
                "",
                "## Limits",
                "",
                "Laboratory methods/standard years are missing; standards and platforms are mixed.",
                "Drug class evidence does not establish a phenotype. Locus counts do not prove",
                "functional copy number. Separate curated mechanism features are unavailable.",
                "Similarity analyses exclude test components connected to development genomes;",
                "small or empty remaining subsets cannot substantiate broad generalization.",
                "Quality ablation is recorded in validation comparisons. Coefficients and",
                "forest importance show association, not causation. Project bootstrap intervals",
                "can be wide with few projects. Baseline failure must prevent confidence claims.",
            ]
        )
        (report_dir / f"{antibiotic}-model-card.md").write_text("\n".join(card) + "\n")
    if sha256(decision_path) != manifest["frozen_decisions_sha256"]:
        raise ValueError("Model decisions changed during test evaluation")
    (report_dir / "test_metrics.json").write_text(json.dumps(evaluations, indent=2) + "\n")
    manifest["status"] = "completed"
    manifest["models"] = {name: sha256(model_dir / f"{name}.joblib") for name in bundles}
    manifest["report_outputs"] = {
        path.name: sha256(path) for path in report_dir.iterdir() if path.is_file()
    }
    (report_dir / "training_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
