"""Independently verify saved live detection smoke artifacts against raw TSV and frozen inputs."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from genome_amr_predictor.detection import parse_amrfinder, sha256


def verify(report_dir: Path) -> dict:
    summary = json.loads((report_dir / "smoke_summary.json").read_text())
    if summary["status"] != "passed" or len(summary["results"]) != 3:
        raise ValueError("Incomplete smoke test")
    if {row["split"] for row in summary["results"]} != {"train", "validation", "test"}:
        raise ValueError("Smoke samples do not cover the frozen splits")
    for row in summary["results"]:
        run = report_dir / row["biosample_accession"]
        manifest = json.loads((run / "run_manifest.json").read_text())
        evidence = json.loads((run / "detections.json").read_text())
        if row["status"] != manifest["status"] or manifest["status"] != "completed":
            raise ValueError("Unsuccessful detection run")
        if row["evidence_sha256"] != sha256(run / "detections.json"):
            raise ValueError("Evidence changed after smoke test")
        if row["run_manifest_sha256"] != sha256(run / "run_manifest.json"):
            raise ValueError("Run manifest changed after smoke test")
        if manifest["raw_output_sha256"] != sha256(run / "amrfinder.tsv"):
            raise ValueError("Raw TSV fingerprint mismatch")
        if manifest["evidence_sha256"] != sha256(run / "detections.json"):
            raise ValueError("Run evidence fingerprint mismatch")
        parsed = parse_amrfinder(run / "amrfinder.tsv", row["biosample_accession"])
        if parsed != evidence["hits"]:
            raise ValueError("Structured evidence differs from raw detector output")
        counts = {
            "total_hits": len(parsed),
            "amr_hits": sum(hit["element_type"] == "AMR" for hit in parsed),
            "mutation_hits": sum(hit["evidence_kind"] == "mutation" for hit in parsed),
            "methods": dict(sorted(Counter(hit["method"] for hit in parsed).items())),
        }
        if counts != row["summary"] or counts != evidence["summary"]:
            raise ValueError("Detector summary differs from retained evidence")
        if not manifest["software_version"] or not manifest["database_version_report"]:
            raise ValueError("Detector versions missing")
        if sha256(Path(manifest["executable"])) != manifest["executable_sha256"]:
            raise ValueError("Installed executable changed")
        db = Path(manifest["database"])
        for name, digest in manifest["database_files_sha256"].items():
            path = (db / name).resolve()
            if not path.is_relative_to(db.resolve()) or sha256(path) != digest:
                raise ValueError("Installed database changed")
        input_path = Path(manifest["input"]["path"])
        if sha256(input_path) != manifest["input"]["file_sha256"]:
            raise ValueError("Smoke input changed")
    return {
        "status": "passed",
        "live_genomes": 3,
        "total_hits": sum(row["summary"]["total_hits"] for row in summary["results"]),
        "amr_hits": sum(row["summary"]["amr_hits"] for row in summary["results"]),
        "checks": [
            "output hashes",
            "raw/structured equality",
            "summary counts",
            "software/database provenance",
            "installed database hashes",
            "input hashes",
        ],
        "limitation": "Integration validation only; no phenotype or population accuracy claim.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.report_dir)
    (args.report_dir / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
