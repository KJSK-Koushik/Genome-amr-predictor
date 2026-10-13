from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from genome_amr_predictor.cohort_config import AssemblyQCConfig
from genome_amr_predictor.data.fasta import inspect_fasta


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0])
    for row in rows[1:]:
        fields.extend(field for field in row if field not in fields)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _verify_sequence(task: tuple[dict[str, str], Path]) -> tuple[dict[str, Any], str]:
    genome, genome_dir = task
    fasta = (genome_dir / genome.get("path", "")).resolve()
    if not fasta.is_relative_to(genome_dir.resolve()):
        return {}, "unsafe_fasta_path"
    try:
        actual = inspect_fasta(fasta).as_dict()
        if any(str(actual[key]) != genome.get(key) for key in ("file_sha256", "sequence_sha256")):
            return actual, "fasta_hash_mismatch"
        return actual, ""
    except (EOFError, OSError, ValueError) as exc:
        return {}, f"fasta_verification_failed:{type(exc).__name__}"


def finalize_sequence_cohort(
    cohort_path: Path,
    assemblies_path: Path,
    genome_manifest_path: Path,
    output_dir: Path,
    *,
    genome_dir: Path | None = None,
    qc: AssemblyQCConfig | None = None,
    max_workers: int = 1,
    progress: Callable[[int, int], None] | None = None,
) -> dict[str, Any]:
    if not 1 <= max_workers <= 8:
        raise ValueError("Finalization max-workers must be between 1 and 8")
    cohort_rows = _read_csv(cohort_path)
    assembly_rows = _read_csv(assemblies_path)
    genome_rows = _read_csv(genome_manifest_path)
    genome_by_assembly = {row["assembly_accession"]: row for row in genome_rows}
    if len(genome_by_assembly) != len(genome_rows):
        raise ValueError("Genome manifest contains duplicate assembly accessions")
    verification_results = {}
    if genome_dir is not None:
        candidates = [row for row in genome_rows if row.get("download_status") == "verified"]
        tasks = [(row, genome_dir) for row in candidates]
        if max_workers == 1:
            results = map(_verify_sequence, tasks)
            for index, (row, result) in enumerate(zip(candidates, results, strict=True), start=1):
                verification_results[row["assembly_accession"]] = result
                if progress and (index % 250 == 0 or index == len(candidates)):
                    progress(index, len(candidates))
        else:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                results = executor.map(_verify_sequence, tasks)
                for index, (row, result) in enumerate(
                    zip(candidates, results, strict=True), start=1
                ):
                    verification_results[row["assembly_accession"]] = result
                    if progress and (index % 250 == 0 or index == len(candidates)):
                        progress(index, len(candidates))
    labels_by_sample: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cohort_rows:
        labels_by_sample[row["biosample_accession"]].append(row)

    exclusions: list[dict[str, Any]] = []
    eligible_assemblies: dict[str, dict[str, str]] = {}
    for assembly in assembly_rows:
        accession = assembly["assembly_accession"]
        genome = genome_by_assembly.get(accession)
        if genome is None:
            exclusions.append(
                {
                    "biosample_accession": assembly["biosample_accession"],
                    "assembly_accession": accession,
                    "reason": "genome_not_downloaded",
                    "canonical_biosample": "",
                    "canonical_assembly": "",
                }
            )
        elif genome.get("download_status") != "verified":
            exclusions.append(
                {
                    "biosample_accession": assembly["biosample_accession"],
                    "assembly_accession": accession,
                    "reason": str(genome.get("download_status") or "genome_not_verified"),
                    "canonical_biosample": "",
                    "canonical_assembly": "",
                }
            )
        elif not genome.get("sequence_sha256"):
            exclusions.append(
                {
                    "biosample_accession": assembly["biosample_accession"],
                    "assembly_accession": accession,
                    "reason": "missing_sequence_sha256",
                    "canonical_biosample": "",
                    "canonical_assembly": "",
                }
            )
        else:
            reasons = []
            if genome_dir is not None:
                actual, error = verification_results[accession]
                if error:
                    reasons.append(error)
                genome.update({key: str(value) for key, value in actual.items() if key != "path"})
            if qc is not None and not reasons:
                try:
                    length = int(genome["total_length"])
                    contigs = int(genome["contig_count"])
                    n50 = int(genome["contig_n50"])
                    if length < qc.minimum_length:
                        reasons.append("fasta_genome_too_short")
                    if length > qc.maximum_length:
                        reasons.append("fasta_genome_too_long")
                    if contigs > qc.maximum_contigs:
                        reasons.append("fasta_too_many_contigs")
                    if n50 < qc.minimum_contig_n50:
                        reasons.append("fasta_contig_n50_too_low")
                except (KeyError, ValueError):
                    reasons.append("missing_or_invalid_fasta_metrics")
            if reasons:
                for reason in reasons:
                    exclusions.append(
                        {
                            "biosample_accession": assembly["biosample_accession"],
                            "assembly_accession": accession,
                            "reason": reason,
                            "canonical_biosample": "",
                            "canonical_assembly": "",
                        }
                    )
                continue
            eligible_assemblies[accession] = assembly

    by_digest: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    for accession, assembly in eligible_assemblies.items():
        digest = genome_by_assembly[accession]["sequence_sha256"]
        by_digest[digest].append(assembly)

    retained_biosamples: set[str] = set()
    duplicate_conflicts: list[dict[str, Any]] = []
    conflicting_endpoints: set[tuple[str, str]] = set()
    for digest, assemblies in sorted(by_digest.items()):
        representative = min(
            assemblies,
            key=lambda row: (
                -len(labels_by_sample[row["biosample_accession"]]),
                row["biosample_accession"],
            ),
        )
        canonical_biosample = representative["biosample_accession"]
        retained_biosamples.add(canonical_biosample)
        if len(assemblies) == 1:
            continue
        labels_by_antibiotic: defaultdict[str, set[str]] = defaultdict(set)
        samples_by_antibiotic: defaultdict[str, set[str]] = defaultdict(set)
        for assembly in assemblies:
            biosample = assembly["biosample_accession"]
            for row in labels_by_sample[biosample]:
                labels_by_antibiotic[row["antibiotic"]].add(row["binary_label"])
                samples_by_antibiotic[row["antibiotic"]].add(biosample)
            if biosample != canonical_biosample:
                exclusions.append(
                    {
                        "biosample_accession": biosample,
                        "assembly_accession": assembly["assembly_accession"],
                        "reason": "exact_duplicate_sequence",
                        "canonical_biosample": canonical_biosample,
                        "canonical_assembly": representative["assembly_accession"],
                    }
                )
        for antibiotic, labels in labels_by_antibiotic.items():
            if len(labels) > 1:
                conflicting_endpoints.add((digest, antibiotic))
                duplicate_conflicts.append(
                    {
                        "sequence_sha256": digest,
                        "antibiotic": antibiotic,
                        "labels": ";".join(sorted(labels)),
                        "biosamples": ";".join(sorted(samples_by_antibiotic[antibiotic])),
                    }
                )

    final_cohort: list[dict[str, Any]] = []
    final_assemblies: list[dict[str, Any]] = []
    for assembly in assembly_rows:
        biosample = assembly["biosample_accession"]
        accession = assembly["assembly_accession"]
        if biosample not in retained_biosamples or accession not in eligible_assemblies:
            continue
        genome = genome_by_assembly[accession]
        merged_assembly: dict[str, Any] = {
            **assembly,
            "fasta_path": genome.get("path", ""),
            "fasta_file_sha256": genome.get("file_sha256", ""),
            "sequence_sha256": genome.get("sequence_sha256", ""),
            "fasta_total_length": genome.get("total_length", ""),
            "fasta_contig_count": genome.get("contig_count", ""),
            "fasta_contig_n50": genome.get("contig_n50", ""),
            "fasta_gc_percent": genome.get("gc_percent", ""),
            "fasta_ambiguous_percent": genome.get("ambiguous_percent", ""),
        }
        final_assemblies.append(merged_assembly)
        for row in labels_by_sample[biosample]:
            if (genome["sequence_sha256"], row["antibiotic"]) in conflicting_endpoints:
                continue
            final_cohort.append(
                {
                    **row,
                    "sequence_sha256": genome.get("sequence_sha256", ""),
                    "fasta_path": genome.get("path", ""),
                }
            )

    exclusion_counts = Counter(str(row["reason"]) for row in exclusions)
    split_counts = Counter(
        (row["antibiotic"], row["split"], row["binary_label"]) for row in final_cohort
    )
    split_rows = [
        {"antibiotic": antibiotic, "split": split, "binary_label": label, "rows": count}
        for (antibiotic, split, label), count in sorted(split_counts.items())
    ]
    summary = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "metadata_qc_assemblies": len(assembly_rows),
        "verified_genomes": sum(row.get("download_status") == "verified" for row in genome_rows),
        "final_unique_genomes": len(final_assemblies),
        "final_cohort_rows": len(final_cohort),
        "exact_duplicate_sequence_groups": sum(len(values) > 1 for values in by_digest.values()),
        "duplicate_sequence_label_conflicts": len(duplicate_conflicts),
        "sequence_files_reverified": genome_dir is not None,
        "sequence_qc_applied": qc is not None,
        "conflicting_duplicate_endpoints_excluded": len(conflicting_endpoints),
        "ast_metadata_counts": {
            field: dict(
                sorted(Counter(row.get(field) or "not_reported" for row in final_cohort).items())
            )
            for field in ("laboratory_method", "laboratory_platform", "testing_standard")
        },
        "exclusion_reasons": dict(sorted(exclusion_counts.items())),
        "split_biosamples": {
            split: len(
                {row["biosample_accession"] for row in final_cohort if row["split"] == split}
            )
            for split in ("train", "validation", "test")
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "final_cohort.csv", final_cohort)
    _write_csv(output_dir / "final_assemblies.csv", final_assemblies)
    _write_csv(output_dir / "split_summary.csv", split_rows)
    _write_csv(output_dir / "genome_exclusions.csv", exclusions)
    _write_csv(output_dir / "duplicate_sequence_label_conflicts.csv", duplicate_conflicts)
    (output_dir / "final_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "final_report.md").write_text(_render_report(summary), encoding="utf-8")
    return summary


def _render_report(summary: dict[str, Any]) -> str:
    ast_metadata = "\n".join(
        f"| {field} | {value} | {count:,} |"
        for field, counts in summary["ast_metadata_counts"].items()
        for value, count in counts.items()
    )
    exclusions = (
        "\n".join(
            f"- {reason}: {count:,}" for reason, count in summary["exclusion_reasons"].items()
        )
        or "- None"
    )
    splits = "\n".join(
        f"- {split}: {count:,} BioSamples" for split, count in summary["split_biosamples"].items()
    )
    return f"""# Phase 2 Sequence-Verified Cohort Report

Generated: {summary["generated_at_utc"]}

## Outcome

- Metadata-QC assemblies: {summary["metadata_qc_assemblies"]:,}
- Verified downloaded genomes: {summary["verified_genomes"]:,}
- Final unique genomes: {summary["final_unique_genomes"]:,}
- Final isolate-antibiotic rows: {summary["final_cohort_rows"]:,}
- Exact duplicate sequence groups: {summary["exact_duplicate_sequence_groups"]:,}
- Duplicate-sequence phenotype conflicts flagged: {summary["duplicate_sequence_label_conflicts"]:,}
- Downloaded sequence files reread: {summary["sequence_files_reverified"]}
- Sequence quality thresholds applied: {summary["sequence_qc_applied"]}

## Genome exclusions

{exclusions}

## Frozen splits after exact-sequence deduplication

{splits}

## AST metadata coverage

| Field | Value | Isolate-antibiotic rows |
|---|---|---:|
{ast_metadata}

Missing method or testing-standard details remain explicit. These counts do not establish that
laboratory breakpoints and methods are comparable across all submitters; model development must
review these strata and the available measurement metadata.

The final cohort uses verified genome manifest entries. Exact duplicate sequences retain
one deterministic BioSample representative; conflicting antibiotic labels among duplicate sequences
are excluded and reported separately for scientific review. Sequence quality thresholds are applied
when configured. Closely related strains may remain across splits; genomic-distance validation is
still required before reporting model performance.
"""
