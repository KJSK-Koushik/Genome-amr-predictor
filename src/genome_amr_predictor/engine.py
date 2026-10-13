"""Single-genome research workflow with independent detection and prediction stages."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

from genome_amr_predictor import __version__
from genome_amr_predictor.data.fasta import inspect_fasta
from genome_amr_predictor.detection import (
    SCHEMA_VERSION,
    AMRFinderRuntime,
    detect_genome,
    parse_amrfinder,
    sha256,
)
from genome_amr_predictor.pretrained import predict_pretrained

PREDICTION_FIELDS = [
    "antibiotic",
    "prediction",
    "resistance_score",
    "status",
    "features_present",
    "features_total",
]
HIT_FIELDS = [
    "sample_id",
    "element_symbol",
    "evidence_kind",
    "element_type",
    "drug_class",
    "drug_subclass",
    "contig_id",
    "start",
    "stop",
    "strand",
    "method",
    "identity_percent",
    "coverage_percent",
    "warnings",
]


def _json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: ";".join(value) if isinstance(value, list) else value
                    for key, value in row.items()
                }
            )


def _runtime(config: dict, root: Path) -> AMRFinderRuntime:
    detector = config["detector"]
    runtime = AMRFinderRuntime(
        root / detector["executable"], root / detector["database"], detector.get("wsl_distribution")
    )
    if sha256(runtime.executable) != detector["executable_sha256"]:
        raise ValueError("Detector executable checksum differs from the frozen configuration")
    if (runtime.database / "version.txt").read_text().strip() != detector["database_version"]:
        raise ValueError("Detector database version differs from the frozen configuration")
    return runtime


def verified_detection(
    directory: Path, fasta: Path, runtime: AMRFinderRuntime, expected: dict
) -> dict:
    """Reject incomplete, altered, mismatched or outdated detection evidence."""
    manifest = json.loads((directory / "run_manifest.json").read_text(encoding="utf-8"))
    evidence_path = directory / "detections.json"
    raw_path = directory / "amrfinder.tsv"
    if manifest["status"] != "completed" or manifest["schema_version"] != SCHEMA_VERSION:
        raise ValueError("Detection result is incomplete or has an unsupported schema")
    if manifest["input"]["file_sha256"] != sha256(fasta):
        raise ValueError("Detection result belongs to a different input genome")
    if manifest["software_version"] != expected["software_version"]:
        raise ValueError("Detection software version mismatch")
    if (
        f"Database version: {expected['database_version']}\n"
        not in manifest["database_version_report"]
    ):
        raise ValueError("Detection database version mismatch")
    if (
        manifest["organism"] != "Escherichia"
        or manifest["mode"] != "nucleotide"
        or manifest["plus"]
    ):
        raise ValueError("Detection mode differs from the supported research contract")
    if manifest["adapter_sha256"] != sha256(Path(__file__).with_name("detection.py")):
        raise ValueError("Detection adapter version mismatch")
    if manifest["executable_sha256"] != sha256(runtime.executable):
        raise ValueError("Detection executable mismatch")
    db_hashes = {
        p.relative_to(runtime.database).as_posix(): sha256(p)
        for p in sorted(runtime.database.rglob("*"))
        if p.is_file()
    }
    if not db_hashes or db_hashes != manifest["database_files_sha256"]:
        raise ValueError("Detection reference database checksum mismatch")
    if sha256(evidence_path) != manifest["evidence_sha256"]:
        raise ValueError("Detection evidence checksum mismatch")
    if sha256(raw_path) != manifest["raw_output_sha256"]:
        raise ValueError("Detection raw output checksum mismatch")
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if evidence["status"] != "completed" or evidence["schema_version"] != SCHEMA_VERSION:
        raise ValueError("Detection evidence is not complete")
    if evidence["hits"] != parse_amrfinder(raw_path, manifest["sample_id"]):
        raise ValueError("Normalized detection evidence does not match the raw output")
    evidence["verified_manifest"] = manifest
    return evidence


def analyze_genome(
    fasta: Path,
    output_dir: Path,
    sample_id: str,
    config_path: Path,
    *,
    detection_cache: Path | None = None,
) -> dict:
    """Persist a report even if one scientific stage fails; never infer absent genes on failure."""
    if not sample_id.strip():
        raise ValueError("sample_id must be nonempty")
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise ValueError("Analysis output directory must be empty; choose a new run directory")
    result = {
        "schema_version": "genome-amr-analysis-v1",
        "package_version": __version__,
        "sample_id": sample_id,
        "started_at_utc": datetime.now(UTC).isoformat(),
        "status": "running",
        "stage": "input_validation",
        "errors": [],
        "input": None,
        "detection": {"status": "not_run", "evidence": None},
        "phenotype_prediction": {"status": "not_run", "result": None},
        "warnings": [
            "Research and educational use only; laboratory confirmation is required.",
            "Known determinant detection and predicted phenotypes are separate outputs.",
            "No detected determinant does not establish susceptibility.",
        ],
    }
    progress = output_dir / "run_manifest.json"
    _json(progress, result)
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
        if config["schema_version"] != "genome-amr-engine-config-v1":
            raise ValueError("Unsupported engine configuration")
        if config["species"] != "Escherichia coli":
            raise ValueError("Only Escherichia coli is supported")
        root = config_path.resolve().parent.parent
        result["config_sha256"] = sha256(config_path)
        result["engine_sha256"] = sha256(Path(__file__))
        result["prediction_code_sha256"] = sha256(Path(__file__).with_name("pretrained.py"))
        metrics = inspect_fasta(fasta)
        opener = gzip.open if fasta.suffix.lower() == ".gz" else open
        with opener(fasta, "rt", encoding="ascii") as stream:
            identifiers = [
                line[1:].split()[0] if line[1:].split() else ""
                for line in stream
                if line.startswith(">")
            ]
        if "" in identifiers or len(identifiers) != len(set(identifiers)):
            raise ValueError("FASTA contig identifiers must be nonempty and unique")
        snapshot = output_dir / ("input.fna.gz" if fasta.suffix.lower() == ".gz" else "input.fna")
        shutil.copyfile(fasta, snapshot)
        if sha256(snapshot) != metrics.file_sha256:
            raise ValueError("Input changed while taking the analysis snapshot")
        result["input"] = metrics.as_dict()
        result["input"]["snapshot"] = str(snapshot.resolve())
        result["input"]["species"] = config["species"]
        result["warnings"].append(
            "E. coli species is user-supplied; automatic identification is unavailable."
        )
        result["stage"] = "phenotype_prediction"
        _json(progress, result)
        try:
            model_dir = root / config["model_dir"]
            if sha256(model_dir / "manifest.json") != config["registry_sha256"]:
                raise ValueError("Pretrained model registry differs from the frozen configuration")
            prediction = predict_pretrained(snapshot, model_dir)
            if prediction["input_sha256"] != metrics.file_sha256:
                raise ValueError("Prediction input changed during analysis")
            result["phenotype_prediction"] = {"status": "completed", "result": prediction}
            result["warnings"].extend(prediction["warnings"])
        except Exception as exc:
            result["phenotype_prediction"] = {"status": "failed", "result": None, "error": str(exc)}
            result["errors"].append({"stage": "phenotype_prediction", "message": str(exc)})
        result["stage"] = "detection"
        _json(progress, result)
        try:
            runtime = _runtime(config, root)
            if detection_cache is None:
                directory = output_dir / "detection"
                detect_genome(
                    snapshot,
                    directory,
                    sample_id,
                    runtime,
                    threads=config["detector"]["threads"],
                    timeout=config["detector"]["timeout_seconds"],
                )
            else:
                directory = detection_cache
            evidence = verified_detection(directory, snapshot, runtime, config["detector"])
            result["detection"] = {
                "status": "completed",
                "evidence": evidence,
                "reused": detection_cache is not None,
                "source_directory": str(directory.resolve()),
            }
        except Exception as exc:
            result["detection"] = {"status": "failed", "evidence": None, "error": str(exc)}
            result["errors"].append({"stage": "detection", "message": str(exc)})
            result["warnings"].append("Detection failed; gene presence/absence is unknown.")
        if sha256(snapshot) != metrics.file_sha256:
            raise ValueError("Analysis input snapshot changed during execution")
        count = sum(
            result[key]["status"] == "completed" for key in ("detection", "phenotype_prediction")
        )
        result["status"] = "completed" if count == 2 else "partial" if count == 1 else "failed"
    except Exception as exc:
        result["status"] = "failed"
        result["errors"].append({"stage": result["stage"], "message": str(exc)})
        # Any input integrity failure invalidates both outputs.
        result["detection"] = {"status": "not_available", "evidence": None}
        result["phenotype_prediction"] = {"status": "not_available", "result": None}
    result["stage"] = "finished"
    result["finished_at_utc"] = datetime.now(UTC).isoformat()
    _json(output_dir / "analysis.json", result)
    prediction = result["phenotype_prediction"]["result"]
    _csv(
        output_dir / "predictions.csv",
        prediction["predictions"] if prediction else [],
        PREDICTION_FIELDS,
    )
    evidence = result["detection"]["evidence"]
    _csv(output_dir / "detections.csv", evidence["hits"] if evidence else [], HIT_FIELDS)
    note = [
        f"# Genome AMR analysis: {sample_id}",
        "",
        f"Run status: **{result['status']}**",
        "",
        f"Detection status: {result['detection']['status']}",
        f"Prediction status: {result['phenotype_prediction']['status']}",
        "",
        "Empty CSV files from failed stages do not indicate absent genes or susceptibility.",
        "",
        "## Limitations",
        "",
        *[f"- {warning}" for warning in result["warnings"]],
    ]
    if prediction:
        note += [
            "",
            "## Phenotype predictions",
            "",
            "| Antibiotic | Prediction | Resistance score | Status |",
            "| --- | --- | --- | --- |",
        ]
        for row in prediction["predictions"]:
            score = f"{row['resistance_score']:.4f}" if row["resistance_score"] is not None else "—"
            note.append(
                f"| {row['antibiotic']} | {row['prediction'] or '—'} | {score} | {row['status']} |"
            )
        note += ["", "Scores are external model outputs, not locally calibrated confidence."]
    if evidence:
        note += [
            "",
            "## Known determinant evidence",
            "",
            f"Total detected elements: {evidence['summary']['total_hits']}; "
            f"AMR elements: {evidence['summary']['amr_hits']}.",
            "",
            "| Element | Type | Evidence | Class | Method |",
            "| --- | --- | --- | --- | --- |",
        ]
        for hit in evidence["hits"]:
            note.append(
                f"| {hit['element_symbol']} | {hit['element_type']} | "
                f"{hit['evidence_kind']} | {hit['drug_class']} | {hit['method']} |"
            )
    note += ["", "## Errors", "", *[f"- {e['stage']}: {e['message']}" for e in result["errors"]]]
    (output_dir / "report.md").write_text("\n".join(note) + "\n", encoding="utf-8")
    result["artifacts_sha256"] = {
        name: sha256(output_dir / name)
        for name in ("analysis.json", "predictions.csv", "detections.csv", "report.md")
    }
    _json(progress, result)
    return result
