"""Run one real, hash-checked frozen genome per split; this is not a model evaluation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from genome_amr_predictor.detection import AMRFinderRuntime, detect_genome, sha256


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--assemblies",
        type=Path,
        default=Path("reports/phase2/2026-08-18/final/final_assemblies.csv"),
    )
    parser.add_argument("--genome-dir", type=Path, default=Path("data/raw/phase2-genomes"))
    parser.add_argument(
        "--executable",
        type=Path,
        default=Path("//wsl.localhost/Ubuntu/opt/genome-amr/amrfinderplus/4.2.7/amrfinder"),
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=Path("//wsl.localhost/Ubuntu/opt/genome-amr/amrfinderplus/4.2.7/data/2026-08-07.1"),
    )
    parser.add_argument("--wsl-distribution", default="Ubuntu")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    runtime = AMRFinderRuntime(args.executable, args.database, args.wsl_distribution)
    with args.assemblies.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    results = []
    root = args.genome_dir.resolve()
    for split in ("train", "validation", "test"):
        row = next(row for row in rows if row["split"] == split)
        fasta = (root / row["fasta_path"]).resolve()
        if not fasta.is_relative_to(root) or sha256(fasta) != row["fasta_file_sha256"]:
            raise ValueError(f"Frozen FASTA fingerprint mismatch: {row['assembly_accession']}")
        output = args.output_dir / row["biosample_accession"]
        print(f"Detecting {row['biosample_accession']} ({split})", flush=True)
        manifest = detect_genome(fasta, output, row["biosample_accession"], runtime)
        if manifest["input"]["sequence_sha256"] != row["sequence_sha256"]:
            raise ValueError("Frozen sequence fingerprint mismatch")
        evidence = json.loads((output / "detections.json").read_text())
        results.append(
            {
                "biosample_accession": row["biosample_accession"],
                "assembly_accession": row["assembly_accession"],
                "split": split,
                "status": manifest["status"],
                "summary": evidence["summary"],
                "evidence_sha256": sha256(output / "detections.json"),
                "run_manifest_sha256": sha256(output / "run_manifest.json"),
                "software_version": manifest["software_version"],
                "database_version_report": manifest["database_version_report"],
            }
        )
        print(f"Completed: {manifest['total_hits']} hits", flush=True)
    summary = {
        "status": "passed",
        "assembly_table_sha256": sha256(args.assemblies),
        "purpose": "Integration smoke test; no phenotype performance claims",
        "results": results,
    }
    (args.output_dir / "smoke_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
