"""Independently verify saved Phase 5 live, reuse, and mismatched-cache analysis runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from genome_amr_predictor.detection import parse_amrfinder, sha256

ARTIFACTS = ("analysis.json", "predictions.csv", "detections.csv", "report.md")


def _run(directory: Path) -> dict:
    manifest = json.loads((directory / "run_manifest.json").read_text(encoding="utf-8"))
    analysis = json.loads((directory / "analysis.json").read_text(encoding="utf-8"))
    if manifest["stage"] != "finished" or analysis["status"] != manifest["status"]:
        raise ValueError(f"{directory.name}: run did not finish consistently")
    for name in ARTIFACTS:
        if sha256(directory / name) != manifest["artifacts_sha256"][name]:
            raise ValueError(f"{directory.name}: {name} changed after the run")
    source = Path(manifest["input"]["path"])
    if source.exists() and sha256(source) != manifest["input"]["file_sha256"]:
        raise ValueError(f"{directory.name}: source genome changed")
    snapshot = Path(manifest["input"]["snapshot"])
    if snapshot.exists() and sha256(snapshot) != manifest["input"]["file_sha256"]:
        raise ValueError(f"{directory.name}: input snapshot changed")
    return manifest


def _detection_rows(directory: Path) -> int:
    return len((directory / "detections.csv").read_text(encoding="utf-8").splitlines()) - 1


def verify(report_dir: Path) -> dict:
    live = _run(report_dir / "live")
    reuse = _run(report_dir / "reuse")
    mismatch = _run(report_dir / "mismatched-cache")

    if live["status"] != "completed" or live["errors"] or live["detection"]["reused"]:
        raise ValueError("Live run must complete with fresh detection and no errors")
    source = Path(live["detection"]["source_directory"])
    evidence = live["detection"]["evidence"]
    if parse_amrfinder(source / "amrfinder.tsv", live["sample_id"]) != evidence["hits"]:
        raise ValueError("Live structured evidence differs from raw detector output")
    if _detection_rows(report_dir / "live") != evidence["summary"]["total_hits"]:
        raise ValueError("Live detection export row count differs from evidence")
    predictions = live["phenotype_prediction"]["result"]["predictions"]
    if live["phenotype_prediction"]["result"]["input_sha256"] != live["input"]["file_sha256"]:
        raise ValueError("Predictions were computed on a different input")

    if reuse["status"] != "completed" or not reuse["detection"]["reused"]:
        raise ValueError("Reuse run must complete using cached detection")
    if Path(reuse["detection"]["source_directory"]).resolve() != source.resolve():
        raise ValueError("Reuse run did not read the live detection directory")
    if reuse["input"]["file_sha256"] != live["input"]["file_sha256"]:
        raise ValueError("Reuse run analysed a different genome")
    for name in ("predictions.csv", "detections.csv"):
        if sha256(report_dir / "reuse" / name) != sha256(report_dir / "live" / name):
            raise ValueError(f"Reuse {name} differs from the live run")

    if mismatch["status"] != "partial" or mismatch["detection"]["status"] != "failed":
        raise ValueError("Mismatched cache must fail detection and leave the run partial")
    if mismatch["detection"]["evidence"] is not None or _detection_rows(
        report_dir / "mismatched-cache"
    ):
        raise ValueError("Rejected detection must not export any hits")
    if mismatch["input"]["file_sha256"] == live["input"]["file_sha256"]:
        raise ValueError("Mismatch check must use a different genome")
    if mismatch["phenotype_prediction"]["status"] != "completed":
        raise ValueError("Prediction must survive a detection failure")
    if "presence/absence is unknown" not in " ".join(mismatch["warnings"]):
        raise ValueError("Detection failure warning missing")

    return {
        "status": "passed",
        "runs": ["live", "reuse", "mismatched-cache"],
        "live_sample": live["sample_id"],
        "live_total_hits": evidence["summary"]["total_hits"],
        "live_amr_hits": evidence["summary"]["amr_hits"],
        "predicted_antibiotics": sum(row["prediction"] is not None for row in predictions),
        "unavailable_antibiotics": sum(row["prediction"] is None for row in predictions),
        "checks": [
            "artifact hashes",
            "input and snapshot hashes",
            "raw/structured detection equality",
            "prediction input binding",
            "reuse equals live exports",
            "mismatched cache rejected without hits",
            "prediction retained on detection failure",
        ],
        "limitation": "Integration validation on one genome; no accuracy claim.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.report_dir)
    (args.report_dir / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
