"""Versioned detector features without phenotype, project, location, or lineage shortcuts."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from genome_amr_predictor.detection import SCHEMA_VERSION, parse_amrfinder, sha256

FEATURE_VERSION = "amr-detector-features-v1"
QUALITY_FIELDS = {
    "fasta_total_length": ("quality:length_mb", 1_000_000),
    "fasta_contig_count": ("quality:contigs_500", 500),
    "fasta_contig_n50": ("quality:n50_mb", 1_000_000),
    "fasta_gc_percent": ("quality:gc_fraction", 100),
    "fasta_ambiguous_percent": ("quality:ambiguity_fraction", 100),
}


def detector_features(hits: list[dict]) -> dict[str, float]:
    features: dict[str, float] = {}
    loci = defaultdict(set)
    for hit in hits:
        if hit["element_type"] != "AMR":
            continue
        kind = "mutation" if hit["evidence_kind"] == "mutation" else "gene"
        symbol = hit["element_symbol"]
        features[f"{kind}:{symbol}"] = 1.0
        # A detected assembly locus is not a functional gene copy-number estimate.
        loci[f"loci:{symbol}"].add((hit["contig_id"], hit["start"], hit["stop"]))
        for node in hit["hierarchy_nodes"]:
            features[f"family:{node}"] = 1.0
        for field, prefix in (
            ("drug_class", "class"),
            ("drug_subclass", "subclass"),
            ("method", "method"),
            ("scope", "scope"),
        ):
            if hit[field] not in {"", "NA"}:
                features[f"{prefix}:{hit[field]}"] = 1.0
        for warning in hit["warnings"]:
            features[f"warning:{warning}"] = 1.0
    features.update({name: float(len(locations)) for name, locations in loci.items()})
    return dict(sorted(features.items()))


def build_features(assemblies: Path, detection_dir: Path, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise ValueError("Feature output directory must be empty")
    with assemblies.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    seen = set()
    outputs = []
    fingerprints = []
    detector_identity = None
    for row in rows:
        sample = row["biosample_accession"]
        if sample in seen:
            raise ValueError("Duplicate assembly BioSample")
        seen.add(sample)
        run = detection_dir / sample
        manifest = json.loads((run / "run_manifest.json").read_text())
        evidence = json.loads((run / "detections.json").read_text())
        if manifest["status"] != "completed" or evidence["status"] != "completed":
            raise ValueError(f"Unsuccessful detection: {sample}")
        if evidence["schema_version"] != SCHEMA_VERSION:
            raise ValueError("Unsupported detector evidence schema")
        if (
            manifest["sample_id"] != sample
            or manifest["input"]["file_sha256"] != row["fasta_file_sha256"]
            or manifest["input"]["sequence_sha256"] != row["sequence_sha256"]
            or sha256(run / "detections.json") != manifest["evidence_sha256"]
            or sha256(run / "amrfinder.tsv") != manifest["raw_output_sha256"]
        ):
            raise ValueError(f"Detection/frozen cohort provenance mismatch: {sample}")
        identity = {
            "software_version": manifest["software_version"],
            "executable_sha256": manifest["executable_sha256"],
            "database_files_sha256": manifest["database_files_sha256"],
            "adapter_sha256": manifest["adapter_sha256"],
            "organism": manifest["organism"],
            "mode": manifest["mode"],
            "plus": manifest["plus"],
        }
        if detector_identity is None:
            detector_identity = identity
        elif identity != detector_identity:
            raise ValueError("Cohort detector versions/settings differ")
        hits = parse_amrfinder(run / "amrfinder.tsv", sample)
        if hits != evidence["hits"]:
            raise ValueError("Raw TSV and structured evidence differ")
        detector = detector_features(hits)
        quality = {
            name: float(row[field]) / divisor for field, (name, divisor) in QUALITY_FIELDS.items()
        }
        outputs.append(
            {
                "sample_id": sample,
                "split": row["split"],
                "group_id": row["group_id"],
                "detector": detector,
                "detector_quality": {**detector, **quality},
            }
        )
        fingerprints.append(
            {
                "sample_id": sample,
                "evidence_sha256": manifest["evidence_sha256"],
                "manifest_sha256": sha256(run / "run_manifest.json"),
            }
        )
    if not outputs:
        raise ValueError("No genomes available for feature building")
    feature_path = output_dir / "features.jsonl"
    with feature_path.open("w", encoding="utf-8") as stream:
        for row in outputs:
            stream.write(json.dumps(row, sort_keys=True) + "\n")
    summary = {
        "status": "completed",
        "schema_version": FEATURE_VERSION,
        "genomes": len(outputs),
        "assemblies_sha256": sha256(assemblies),
        "features_sha256": sha256(feature_path),
        "detector_identity": detector_identity,
        "detector_results_sha256": hashlib.sha256(
            json.dumps(fingerprints, sort_keys=True).encode()
        ).hexdigest(),
        "builder_sha256": sha256(Path(__file__)),
        "mechanisms": "No mechanism feature inferred when curated annotations are absent",
        "vocabulary_policy": "Fit vocabulary only on each antibiotic training set",
    }
    (output_dir / "feature_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
