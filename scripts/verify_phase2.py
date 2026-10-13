"""Check the final cohort against its frozen source and provenance manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def verify(report_dir: Path) -> dict[str, object]:
    final = report_dir / "final"
    manifest = json.loads((final / "final_manifest.json").read_text(encoding="utf-8"))
    for name in ("cohort", "assemblies"):
        metadata = manifest["inputs"][name]
        with (report_dir / metadata["name"]).open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != metadata["sha256"]:
            raise ValueError(f"Frozen source fingerprint changed: {name}")
    for name, metadata in manifest["outputs"].items():
        with (final / name).open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if digest != metadata["sha256"]:
            raise ValueError(f"Output fingerprint changed: {name}")

    original = read_csv(report_dir / "cohort.csv")
    rows = read_csv(final / "final_cohort.csv")
    assemblies = read_csv(final / "final_assemblies.csv")
    if not rows or not assemblies:
        raise ValueError("The final cohort is empty")
    originals = {(r["biosample_accession"], r["antibiotic"]): r for r in original}
    assembly_by_sample = {r["biosample_accession"]: r for r in assemblies}
    if len(assembly_by_sample) != len(assemblies):
        raise ValueError("Final assemblies contain duplicate BioSamples")
    if len({r["assembly_accession"] for r in assemblies}) != len(assemblies):
        raise ValueError("Final assemblies contain duplicate accessions")
    if len({r["sequence_sha256"] for r in assemblies}) != len(assemblies):
        raise ValueError("Exact duplicate genomes remain in the final cohort")

    qc = manifest["config"]["assembly_qc"]
    project_splits: defaultdict[str, set[str]] = defaultdict(set)
    for assembly in assemblies:
        project_splits[assembly["group_id"]].add(assembly["split"])
        length = int(assembly["fasta_total_length"])
        if not qc["minimum_length"] <= length <= qc["maximum_length"]:
            raise ValueError("An assembly fails the configured length thresholds")
        if int(assembly["fasta_contig_count"]) > qc["maximum_contigs"]:
            raise ValueError("An assembly fails the configured contig threshold")
        if int(assembly["fasta_contig_n50"]) < qc["minimum_contig_n50"]:
            raise ValueError("An assembly fails the configured N50 threshold")
    if any(len(splits) != 1 for splits in project_splits.values()):
        raise ValueError("A BioProject appears in multiple splits")

    conflicts = {
        (r["sequence_sha256"], r["antibiotic"])
        for r in read_csv(final / "duplicate_sequence_label_conflicts.csv")
    }
    seen = set()
    for row in rows:
        key = (row["biosample_accession"], row["antibiotic"])
        if key in seen or key not in originals:
            raise ValueError("Final cohort contains a repeated or invented endpoint")
        seen.add(key)
        source = originals[key]
        if any(row[field] != value for field, value in source.items()):
            raise ValueError("A frozen source label, split or metadata field changed")
        assembly = assembly_by_sample[row["biosample_accession"]]
        for field in ("assembly_accession", "sequence_sha256", "fasta_path", "split"):
            if row[field] != assembly[field]:
                raise ValueError(f"Cohort/assembly join disagrees on {field}")
        if (row["sequence_sha256"], row["antibiotic"]) in conflicts:
            raise ValueError("A conflicting duplicate phenotype remains")
        if row["binary_label"] not in {"0", "1"}:
            raise ValueError("A nonbinary phenotype remains")

    counts = Counter((r["antibiotic"], r["split"], r["binary_label"]) for r in rows)
    summary = manifest["summary"]
    if not summary["sequence_files_reverified"] or not summary["sequence_qc_applied"]:
        raise ValueError("Final cohort requires FASTA reverification and sequence QC")
    if (
        len(rows) != summary["final_cohort_rows"]
        or len(assemblies) != summary["final_unique_genomes"]
    ):
        raise ValueError("Final row counts disagree with the manifest")
    return {
        "status": "passed",
        "unique_genomes": len(assemblies),
        "phenotype_rows": len(rows),
        "project_groups": len(project_splits),
        "exact_sequence_duplicates": 0,
        "project_split_overlap": 0,
        "changed_frozen_endpoints": 0,
        "retained_conflicting_endpoints": 0,
        "both_classes_per_antibiotic_split": all(
            counts[(drug, split, label)] > 0
            for drug in manifest["config"]["antibiotics"]
            for split in ("train", "validation", "test")
            for label in ("0", "1")
        ),
        "limitation": "Genomic-distance/lineage leakage is not evaluated by these checks.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.report_dir)
    print(json.dumps(result, indent=2))
    (args.report_dir / "final" / "verification.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
