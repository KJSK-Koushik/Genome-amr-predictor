from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from genome_amr_predictor.cohort_config import AssemblyQCConfig, Phase2Config, SplitConfig
from genome_amr_predictor.data.ncbi_datasets import AssemblyMetadata
from genome_amr_predictor.data.schema import ASTRecord


@dataclass(frozen=True, slots=True)
class NormalizedLabels:
    records: tuple[ASTRecord, ...]
    duplicate_rows_removed: int
    conflicting_pairs_removed: int
    rows_without_biosample_removed: int


@dataclass(frozen=True, slots=True)
class CohortBundle:
    summary: dict[str, Any]
    cohort_rows: list[dict[str, Any]]
    assembly_rows: list[dict[str, Any]]
    excluded_rows: list[dict[str, Any]]
    split_rows: list[dict[str, Any]]


def normalize_training_labels(
    records: list[ASTRecord], antibiotics: tuple[str, ...]
) -> NormalizedLabels:
    allowed = set(antibiotics)
    grouped: defaultdict[tuple[str, str], list[ASTRecord]] = defaultdict(list)
    rows_without_biosample = 0
    for record in records:
        if record.antibiotic not in allowed or record.binary_label is None:
            continue
        if not record.biosample_accession:
            rows_without_biosample += 1
            continue
        grouped[(record.biosample_accession, record.antibiotic)].append(record)

    normalized: list[ASTRecord] = []
    duplicate_rows = 0
    conflicts = 0
    for key in sorted(grouped):
        values = grouped[key]
        labels = {record.binary_label for record in values}
        duplicate_rows += max(0, len(values) - 1)
        if len(labels) != 1:
            conflicts += 1
            continue
        normalized.append(values[0])
    return NormalizedLabels(
        records=tuple(normalized),
        duplicate_rows_removed=duplicate_rows,
        conflicting_pairs_removed=conflicts,
        rows_without_biosample_removed=rows_without_biosample,
    )


def assembly_qc_reasons(
    assembly: AssemblyMetadata, species: str, config: AssemblyQCConfig
) -> tuple[list[str], list[str]]:
    reasons: list[str] = []
    warnings: list[str] = []
    organism = " ".join(assembly.organism.casefold().split())
    expected = " ".join(species.casefold().split())
    if not (organism == expected or organism.startswith(f"{expected} ")):
        reasons.append("organism_mismatch")

    if assembly.total_sequence_length is None:
        reasons.append("missing_total_sequence_length")
    elif assembly.total_sequence_length < config.minimum_length:
        reasons.append("genome_too_short")
    elif assembly.total_sequence_length > config.maximum_length:
        reasons.append("genome_too_long")

    if assembly.number_of_contigs is None:
        reasons.append("missing_contig_count")
    elif assembly.number_of_contigs > config.maximum_contigs:
        reasons.append("too_many_contigs")

    if assembly.contig_n50 is None:
        reasons.append("missing_contig_n50")
    elif assembly.contig_n50 < config.minimum_contig_n50:
        reasons.append("contig_n50_too_low")

    if assembly.checkm_completeness is None or assembly.checkm_contamination is None:
        if config.require_checkm:
            reasons.append("missing_checkm")
        else:
            warnings.append("missing_checkm")
    else:
        if assembly.checkm_completeness < config.minimum_checkm_completeness:
            reasons.append("checkm_completeness_too_low")
        if assembly.checkm_contamination > config.maximum_checkm_contamination:
            reasons.append("checkm_contamination_too_high")

    if assembly.ani_match_status:
        if config.reject_ani_mismatch and assembly.ani_match_status != "species_match":
            reasons.append("ani_species_mismatch")
    else:
        warnings.append("missing_ani")
    return reasons, warnings


def _stable_group_order(group_id: str, seed: str) -> int:
    digest = hashlib.sha256(f"{seed}|{group_id}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def _assign_group_splits(
    labels_by_sample: dict[str, list[ASTRecord]],
    accepted: dict[str, tuple[AssemblyMetadata, list[str]]],
    config: SplitConfig,
) -> dict[str, str]:
    split_fractions = {
        "train": config.train_fraction,
        "validation": config.validation_fraction,
        "test": config.test_fraction,
    }
    samples_by_group: defaultdict[str, list[str]] = defaultdict(list)
    for biosample, (assembly, _) in accepted.items():
        samples_by_group[assembly.bioproject_accession or biosample].append(biosample)

    features_by_group: dict[str, Counter[str]] = {}
    total_features: Counter[str] = Counter()
    for group_id, biosamples in samples_by_group.items():
        features: Counter[str] = Counter({"samples": len(biosamples)})
        for biosample in biosamples:
            for record in labels_by_sample[biosample]:
                features[f"{record.antibiotic}|{record.binary_label}"] += 1
        features_by_group[group_id] = features
        total_features.update(features)

    targets = {
        split: {feature: total * fraction for feature, total in total_features.items()}
        for split, fraction in split_fractions.items()
    }
    current = {split: Counter() for split in split_fractions}
    assignments: dict[str, str] = {}
    ordered_groups = sorted(
        samples_by_group,
        key=lambda group_id: (
            -features_by_group[group_id]["samples"],
            _stable_group_order(group_id, config.seed),
        ),
    )

    def assignment_cost(candidate_split: str, group_id: str) -> float:
        cost = 0.0
        for split in split_fractions:
            for feature, target in targets[split].items():
                value = current[split][feature]
                if split == candidate_split:
                    value += features_by_group[group_id][feature]
                cost += ((value - target) / max(target, 1.0)) ** 2
        return cost

    for group_id in ordered_groups:
        selected_split = min(
            split_fractions,
            key=lambda split: (
                assignment_cost(split, group_id),
                _stable_group_order(f"{group_id}|{split}", config.seed),
            ),
        )
        assignments[group_id] = selected_split
        current[selected_split].update(features_by_group[group_id])
    return assignments


def _assembly_row(
    assembly: AssemblyMetadata, *, split: str, group_id: str, warnings: list[str]
) -> dict[str, Any]:
    row = assembly.as_dict()
    row.update(
        {
            "group_id": group_id,
            "split": split,
            "qc_warnings": ";".join(warnings),
        }
    )
    return row


def build_cohort(
    normalized: NormalizedLabels,
    assemblies: dict[str, AssemblyMetadata | None],
    config: Phase2Config,
) -> CohortBundle:
    labels_by_sample: defaultdict[str, list[ASTRecord]] = defaultdict(list)
    for record in normalized.records:
        labels_by_sample[record.biosample_accession].append(record)

    accepted: dict[str, tuple[AssemblyMetadata, list[str]]] = {}
    excluded: list[dict[str, Any]] = []
    for biosample in sorted(labels_by_sample):
        assembly = assemblies.get(biosample)
        if assembly is None:
            excluded.append(
                {
                    "biosample_accession": biosample,
                    "assembly_accession": "",
                    "reason": "no_assembly",
                }
            )
            continue
        reasons, warnings = assembly_qc_reasons(assembly, config.species, config.assembly_qc)
        if reasons:
            excluded.append(
                {
                    "biosample_accession": biosample,
                    "assembly_accession": assembly.assembly_accession,
                    "reason": ";".join(reasons),
                }
            )
            continue
        accepted[biosample] = (assembly, warnings)

    by_assembly: defaultdict[str, list[str]] = defaultdict(list)
    for biosample, (assembly, _) in accepted.items():
        by_assembly[assembly.assembly_accession].append(biosample)
    for assembly_accession, biosamples in sorted(by_assembly.items()):
        for duplicate in sorted(biosamples)[1:]:
            accepted.pop(duplicate)
            excluded.append(
                {
                    "biosample_accession": duplicate,
                    "assembly_accession": assembly_accession,
                    "reason": "duplicate_assembly_accession",
                }
            )

    group_splits = _assign_group_splits(labels_by_sample, accepted, config.splits)
    assembly_rows: list[dict[str, Any]] = []
    cohort_rows: list[dict[str, Any]] = []
    for biosample in sorted(accepted):
        assembly, warnings = accepted[biosample]
        group_id = assembly.bioproject_accession or biosample
        split = group_splits[group_id]
        assembly_rows.append(
            _assembly_row(assembly, split=split, group_id=group_id, warnings=warnings)
        )
        for record in sorted(labels_by_sample[biosample], key=lambda value: value.antibiotic):
            cohort_rows.append(
                {
                    "biosample_accession": biosample,
                    "assembly_accession": assembly.assembly_accession,
                    "bioproject_accession": assembly.bioproject_accession,
                    "group_id": group_id,
                    "split": split,
                    "antibiotic": record.antibiotic,
                    "binary_label": record.binary_label,
                    "phenotype": record.phenotype_normalized,
                    "measurement_sign": record.measurement_sign,
                    "measurement": record.measurement,
                    "measurement_units": record.measurement_units,
                    "laboratory_method": record.laboratory_method,
                    "laboratory_platform": record.laboratory_platform,
                    "testing_standard": record.testing_standard,
                    "collection_date": record.collection_date,
                    "geo_loc_name": record.geo_loc_name,
                    "isolation_source": record.isolation_source,
                    "host": record.host,
                }
            )

    split_rows = _split_summary(cohort_rows)
    group_counts = Counter(str(row["group_id"]) for row in assembly_rows)
    split_counts = {
        split: sum(row["split"] == split for row in assembly_rows)
        for split in ("train", "validation", "test")
    }
    exclusions = Counter(
        reason for row in excluded for reason in str(row["reason"]).split(";") if reason
    )
    summary = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "species": config.species,
        "normalized_binary_labels": len(normalized.records),
        "duplicate_label_rows_removed": normalized.duplicate_rows_removed,
        "conflicting_label_pairs_removed": normalized.conflicting_pairs_removed,
        "rows_without_biosample_removed": normalized.rows_without_biosample_removed,
        "biosamples_with_binary_labels": len(labels_by_sample),
        "biosamples_with_selected_assemblies": sum(
            assembly is not None for assembly in assemblies.values()
        ),
        "biosamples_passing_metadata_qc": len(assembly_rows),
        "frozen_cohort_rows": len(cohort_rows),
        "exclusion_reasons": dict(sorted(exclusions.items())),
        "split_biosamples": split_counts,
        "split_actual_fractions": {
            split: round(count / max(len(assembly_rows), 1), 6)
            for split, count in split_counts.items()
        },
        "split_target_fractions": {
            "train": config.splits.train_fraction,
            "validation": config.splits.validation_fraction,
            "test": config.splits.test_fraction,
        },
        "group_count": len(group_counts),
        "largest_groups": [
            {"group_id": group_id, "biosamples": count}
            for group_id, count in group_counts.most_common(5)
        ],
    }
    return CohortBundle(
        summary=summary,
        cohort_rows=cohort_rows,
        assembly_rows=assembly_rows,
        excluded_rows=excluded,
        split_rows=split_rows,
    )


def _split_summary(cohort_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    counts: defaultdict[tuple[str, str], Counter[int]] = defaultdict(Counter)
    samples: defaultdict[tuple[str, str], set[str]] = defaultdict(set)
    for row in cohort_rows:
        key = (str(row["antibiotic"]), str(row["split"]))
        counts[key][int(row["binary_label"])] += 1
        samples[key].add(str(row["biosample_accession"]))
    return [
        {
            "antibiotic": antibiotic,
            "split": split,
            "unique_biosamples": len(samples[(antibiotic, split)]),
            "resistant": counts[(antibiotic, split)][1],
            "susceptible": counts[(antibiotic, split)][0],
            "total": sum(counts[(antibiotic, split)].values()),
        }
        for antibiotic, split in sorted(counts)
    ]


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_cohort_outputs(output_dir: Path, bundle: CohortBundle, manifest: dict[str, Any]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "cohort_summary.json").write_text(
        json.dumps(bundle.summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "data_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_csv(output_dir / "cohort.csv", bundle.cohort_rows)
    _write_csv(output_dir / "assemblies.csv", bundle.assembly_rows)
    _write_csv(output_dir / "excluded_biosamples.csv", bundle.excluded_rows)
    _write_csv(output_dir / "split_summary.csv", bundle.split_rows)
    accessions = sorted({str(row["assembly_accession"]) for row in bundle.assembly_rows})
    (output_dir / "assembly_accessions.txt").write_text(
        "\n".join(accessions) + "\n", encoding="utf-8"
    )
    (output_dir / "cohort_report.md").write_text(
        render_cohort_report(bundle, manifest), encoding="utf-8"
    )


def render_cohort_report(bundle: CohortBundle, manifest: dict[str, Any]) -> str:
    summary = bundle.summary
    exclusions = (
        "\n".join(
            f"- {reason}: {count:,}" for reason, count in summary["exclusion_reasons"].items()
        )
        or "- None"
    )
    split_lines = "\n".join(
        f"- {split}: {count:,} BioSamples "
        f"({100 * summary['split_actual_fractions'][split]:.1f}%; "
        f"target {100 * summary['split_target_fractions'][split]:.1f}%)"
        for split, count in summary["split_biosamples"].items()
    )
    largest_groups = "\n".join(
        f"- {row['group_id']}: {row['biosamples']:,} BioSamples"
        for row in summary["largest_groups"]
    )
    return f"""# Phase 2 Assembly-Linked Cohort Report

Generated: {summary["generated_at_utc"]}

## Outcome

- Species: *{summary["species"]}*
- Normalized binary AST labels: {summary["normalized_binary_labels"]:,}
- BioSamples with usable labels: {summary["biosamples_with_binary_labels"]:,}
- BioSamples with a selected assembly: {summary["biosamples_with_selected_assemblies"]:,}
- BioSamples passing metadata QC: {summary["biosamples_passing_metadata_qc"]:,}
- Frozen isolate-antibiotic rows: {summary["frozen_cohort_rows"]:,}

## Label cleaning

- Duplicate label rows removed: {summary["duplicate_label_rows_removed"]:,}
- Conflicting isolate-antibiotic pairs removed: {summary["conflicting_label_pairs_removed"]:,}
- Binary rows without BioSample accessions removed: {summary["rows_without_biosample_removed"]:,}

## Metadata-QC exclusions

{exclusions}

## Frozen group-aware split

{split_lines}

The cohort contains {summary["group_count"]:,} indivisible BioProject/BioSample groups. The five
largest are:

{largest_groups}

All rows from the same BioProject are assigned to one split. If BioProject is unavailable, the
BioSample accession is used as the group. A deterministic allocator uses the versioned seed while
balancing sample counts and antibiotic/class strata as closely as whole-project grouping permits.

## Important boundary

This report freezes the assembly-linked, metadata-QC cohort. The genome download stage verifies
FASTA integrity and exact sequence hashes. No resistance detector or phenotype model is trained
in Phase 2.

## Reproducibility

See `data_manifest.json` for source hashes, NCBI endpoint, retrieval time, configuration hash,
software version, and cohort counts. The raw AST and genome files remain outside version control.
"""
