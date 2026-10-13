from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from genome_amr_predictor.config import AuditConfig
from genome_amr_predictor.data.schema import ASTRecord

AUDITED_MISSINGNESS_FIELDS = (
    "biosample_accession",
    "organism",
    "bioproject_accession",
    "sra_accession",
    "isolate",
    "collection_date",
    "geo_loc_name",
    "isolation_source",
    "host",
    "measurement",
    "laboratory_method",
    "laboratory_platform",
    "testing_standard",
)


@dataclass(frozen=True, slots=True)
class AuditBundle:
    summary: dict[str, Any]
    antibiotic_rows: list[dict[str, Any]]
    phenotype_rows: list[dict[str, Any]]
    missingness_rows: list[dict[str, Any]]
    candidate_rows: list[dict[str, Any]]


def _percent(numerator: int, denominator: int) -> float:
    return round((numerator / denominator * 100) if denominator else 0.0, 2)


def audit_records(records: list[ASTRecord], config: AuditConfig) -> AuditBundle:
    phenotype_counts = Counter(record.phenotype_normalized for record in records)
    antibiotic_counts: dict[str, Counter[str]] = defaultdict(Counter)
    biosamples_by_antibiotic: dict[str, set[str]] = defaultdict(set)
    for record in records:
        antibiotic_counts[record.antibiotic][record.phenotype_normalized] += 1
        if record.biosample_accession:
            biosamples_by_antibiotic[record.antibiotic].add(record.biosample_accession)

    antibiotic_rows: list[dict[str, Any]] = []
    for antibiotic, counts in sorted(antibiotic_counts.items()):
        resistant = counts["resistant"]
        susceptible = counts["susceptible"]
        intermediate = counts["intermediate"]
        excluded = counts["excluded"]
        binary_total = resistant + susceptible
        total = binary_total + intermediate + excluded
        minority = min(resistant, susceptible)
        antibiotic_rows.append(
            {
                "antibiotic": antibiotic,
                "unique_biosamples": len(biosamples_by_antibiotic[antibiotic]),
                "total_records": total,
                "resistant": resistant,
                "susceptible": susceptible,
                "intermediate": intermediate,
                "excluded": excluded,
                "binary_total": binary_total,
                "minority_fraction": round(minority / binary_total, 4) if binary_total else 0.0,
            }
        )

    phenotype_rows = [
        {
            "phenotype": phenotype,
            "records": count,
            "percent": _percent(count, len(records)),
        }
        for phenotype, count in sorted(phenotype_counts.items())
    ]

    missingness_rows: list[dict[str, Any]] = []
    for field in AUDITED_MISSINGNESS_FIELDS:
        missing = sum(not str(getattr(record, field) or "").strip() for record in records)
        missingness_rows.append(
            {
                "field": field,
                "missing_records": missing,
                "total_records": len(records),
                "missing_percent": _percent(missing, len(records)),
            }
        )

    rows_by_antibiotic = {row["antibiotic"]: row for row in antibiotic_rows}
    thresholds = config.selection_thresholds
    candidate_rows: list[dict[str, Any]] = []
    for antibiotic in config.candidate_antibiotics:
        values = rows_by_antibiotic.get(
            antibiotic,
            {
                "unique_biosamples": 0,
                "total_records": 0,
                "resistant": 0,
                "susceptible": 0,
                "intermediate": 0,
                "excluded": 0,
                "binary_total": 0,
                "minority_fraction": 0.0,
            },
        )
        reasons: list[str] = []
        if values["binary_total"] < thresholds.minimum_binary_records:
            reasons.append("insufficient binary records")
        if values["resistant"] < thresholds.minimum_resistant_records:
            reasons.append("insufficient resistant records")
        if values["susceptible"] < thresholds.minimum_susceptible_records:
            reasons.append("insufficient susceptible records")
        if values["minority_fraction"] < thresholds.minimum_minority_fraction:
            reasons.append("minority fraction below threshold")
        candidate_rows.append(
            {
                "antibiotic": antibiotic,
                **values,
                "recommendation": "retain" if not reasons else "defer",
                "reason": "; ".join(reasons) if reasons else "meets Phase 1 thresholds",
            }
        )

    unique_biosamples = {
        record.biosample_accession for record in records if record.biosample_accession
    }
    duplicate_keys = Counter((record.biosample_accession, record.antibiotic) for record in records)
    duplicate_pairs = sum(count > 1 for count in duplicate_keys.values())
    conflicts = 0
    labels_by_key: dict[tuple[str, str], set[int]] = defaultdict(set)
    for record in records:
        if record.binary_label is not None:
            labels_by_key[(record.biosample_accession, record.antibiotic)].add(record.binary_label)
    conflicts = sum(len(labels) > 1 for labels in labels_by_key.values())

    summary = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "species": config.species,
        "ast_records": len(records),
        "unique_biosamples": len(unique_biosamples),
        "unique_antibiotics": len(antibiotic_counts),
        "binary_records": phenotype_counts["resistant"] + phenotype_counts["susceptible"],
        "excluded_from_binary": phenotype_counts["intermediate"] + phenotype_counts["excluded"],
        "duplicate_isolate_antibiotic_pairs": duplicate_pairs,
        "conflicting_binary_label_pairs": conflicts,
        "candidate_antibiotics_retained": [
            row["antibiotic"] for row in candidate_rows if row["recommendation"] == "retain"
        ],
        "candidate_antibiotics_deferred": [
            row["antibiotic"] for row in candidate_rows if row["recommendation"] == "defer"
        ],
    }
    return AuditBundle(
        summary=summary,
        antibiotic_rows=antibiotic_rows,
        phenotype_rows=phenotype_rows,
        missingness_rows=missingness_rows,
        candidate_rows=candidate_rows,
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    output = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        output.append("| " + " | ".join(str(value) for value in row) + " |")
    return "\n".join(output)


def render_markdown_report(bundle: AuditBundle, manifest: dict[str, Any]) -> str:
    summary = bundle.summary
    candidate_table = _markdown_table(
        ["Antibiotic", "R", "S", "Binary total", "Minority fraction", "Decision"],
        [
            [
                row["antibiotic"],
                row["resistant"],
                row["susceptible"],
                row["binary_total"],
                f"{row['minority_fraction']:.3f}",
                row["recommendation"],
            ]
            for row in bundle.candidate_rows
        ],
    )
    top_antibiotics = sorted(
        bundle.antibiotic_rows, key=lambda row: row["binary_total"], reverse=True
    )[:15]
    coverage_table = _markdown_table(
        ["Antibiotic", "Unique samples", "R", "S", "I", "Excluded"],
        [
            [
                row["antibiotic"],
                row["unique_biosamples"],
                row["resistant"],
                row["susceptible"],
                row["intermediate"],
                row["excluded"],
            ]
            for row in top_antibiotics
        ],
    )
    retained = ", ".join(summary["candidate_antibiotics_retained"]) or "None"
    deferred = ", ".join(summary["candidate_antibiotics_deferred"]) or "None"
    if manifest["query"] == "file import":
        sample_scope = "Unique BioSamples after the configured species filter"
        reproducibility_note = (
            "The manifest records the source filename and SHA-256 fingerprint, configuration hash, "
            "retrieval time, and processed counts. The raw dataset remains outside version control."
        )
    else:
        sample_scope = "NCBI matching BioSamples at retrieval"
        reproducibility_note = (
            "The manifest records source URLs, query translation, configuration hash, "
            "retrieval time, and processed counts. Raw XML was parsed in batches and was "
            "not retained."
        )
    return f"""# Phase 1 NCBI AST Data Audit

Generated: {summary["generated_at_utc"]}

## Scope

- Species: *{summary["species"]}*
- Query: `{manifest["query"]}`
- {sample_scope}: {manifest["ncbi_reported_biosamples"]:,}
- BioSamples processed: {manifest["processed_biosamples"]:,}
- AST rows parsed: {summary["ast_records"]:,}
- Unique antibiotics observed: {summary["unique_antibiotics"]:,}

This is an availability and label-quality audit. It does not train a model and does not establish
clinical performance.

## Binary-label usability

- Resistant or susceptible rows: {summary["binary_records"]:,}
- Rows excluded from the binary endpoint: {summary["excluded_from_binary"]:,}
- Duplicate isolate-antibiotic pairs: {summary["duplicate_isolate_antibiotic_pairs"]:,}
- Conflicting binary-label pairs: {summary["conflicting_binary_label_pairs"]:,}

## Candidate antibiotic decision

{candidate_table}

**Retain for Phase 2 cohort construction:** {retained}.

**Defer pending more or different data:** {deferred}.

These decisions use only the configured count and class-balance thresholds. Phase 2 must still
verify assembly availability, duplicate genomes, lineage leakage, testing-standard consistency,
and genome quality before a final training cohort is frozen.

## Highest-coverage antibiotics

{coverage_table}

## Important limitations

- BioSample antibiograms are submitted by data providers and are not fully harmonized by NCBI.
- Recent or highly represented projects can dominate the records.
- This audit does not yet join BioSamples to GenBank assemblies.
- Counts do not guarantee that genomes will pass Phase 2 sequence-quality checks.
- `intermediate`, nonsusceptible, undefined, and other non-binary outcomes are excluded from the
  initial R/S endpoint.
- Laboratory method, platform, testing standard, source, time, and geography require cohort-level
  review before modeling.

## Reproducibility

{reproducibility_note}
"""


def write_audit_outputs(
    output_dir: Path,
    bundle: AuditBundle,
    manifest: dict[str, Any],
    *,
    records: list[ASTRecord] | None = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "audit_summary.json").write_text(
        json.dumps(bundle.summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "data_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_csv(output_dir / "antibiotic_summary.csv", bundle.antibiotic_rows)
    _write_csv(output_dir / "phenotype_summary.csv", bundle.phenotype_rows)
    _write_csv(output_dir / "missingness_summary.csv", bundle.missingness_rows)
    _write_csv(output_dir / "candidate_recommendations.csv", bundle.candidate_rows)
    (output_dir / "audit_report.md").write_text(
        render_markdown_report(bundle, manifest), encoding="utf-8"
    )
    if records is not None:
        _write_csv(output_dir / "records.csv", [record.as_dict() for record in records])
