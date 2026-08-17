from __future__ import annotations

import csv
from pathlib import Path

from genome_amr_predictor.data.labels import normalize_antibiotic, normalize_phenotype
from genome_amr_predictor.data.schema import ASTRecord

ALIASES = {
    "biosample_accession": ("biosample_accession", "biosample", "biosample_acc"),
    "organism": ("organism", "scientific_name", "scientific name"),
    "taxonomy_id": ("taxonomy_id", "taxid"),
    "bioproject_accession": ("bioproject_accession", "bioproject"),
    "sra_accession": ("sra_accession", "sra"),
    "isolate": ("isolate", "isolate_name"),
    "collection_date": ("collection_date", "collection date"),
    "geo_loc_name": ("geo_loc_name", "location", "geographic location"),
    "isolation_source": ("isolation_source", "isolation source"),
    "host": ("host",),
    "antibiotic": ("antibiotic",),
    "phenotype_raw": ("phenotype_raw", "resistance_phenotype", "resistance phenotype"),
    "measurement_sign": ("measurement_sign", "measurement sign"),
    "measurement": ("measurement", "mic", "mic (mg/l)", "disk diffusion (mm)"),
    "measurement_units": ("measurement_units", "measurement units", "units"),
    "laboratory_method": ("laboratory_method", "laboratory typing method"),
    "laboratory_platform": ("laboratory_platform", "laboratory typing platform"),
    "vendor": ("vendor",),
    "method_version_or_reagent": (
        "method_version_or_reagent",
        "laboratory typing method version or reagent",
    ),
    "testing_standard": ("testing_standard", "testing standard"),
}


def _normalized_key(value: str) -> str:
    return value.strip().lower().replace("-", "_")


def _get(row: dict[str, str], field: str) -> str:
    normalized_row = {_normalized_key(key): (value or "").strip() for key, value in row.items()}
    for alias in ALIASES[field]:
        value = normalized_row.get(_normalized_key(alias), "")
        if value:
            return value
    return ""


def _detect_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    sample = path.read_text(encoding="utf-8-sig", errors="replace")[:8192]
    return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter


def read_ast_file(path: Path) -> list[ASTRecord]:
    delimiter = _detect_delimiter(path)
    records: list[ASTRecord] = []
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError("AST input file has no header row")
        for row in reader:
            antibiotic = normalize_antibiotic(_get(row, "antibiotic"))
            phenotype_raw = _get(row, "phenotype_raw")
            if not antibiotic or not phenotype_raw:
                continue
            phenotype_normalized, binary_label = normalize_phenotype(phenotype_raw)
            records.append(
                ASTRecord(
                    biosample_accession=_get(row, "biosample_accession"),
                    organism=_get(row, "organism"),
                    taxonomy_id=_get(row, "taxonomy_id"),
                    bioproject_accession=_get(row, "bioproject_accession"),
                    sra_accession=_get(row, "sra_accession"),
                    isolate=_get(row, "isolate"),
                    collection_date=_get(row, "collection_date"),
                    geo_loc_name=_get(row, "geo_loc_name"),
                    isolation_source=_get(row, "isolation_source"),
                    host=_get(row, "host"),
                    antibiotic=antibiotic,
                    phenotype_raw=phenotype_raw,
                    phenotype_normalized=phenotype_normalized,
                    binary_label=binary_label,
                    measurement_sign=_get(row, "measurement_sign"),
                    measurement=_get(row, "measurement"),
                    measurement_units=_get(row, "measurement_units"),
                    laboratory_method=_get(row, "laboratory_method"),
                    laboratory_platform=_get(row, "laboratory_platform"),
                    vendor=_get(row, "vendor"),
                    method_version_or_reagent=_get(row, "method_version_or_reagent"),
                    testing_standard=_get(row, "testing_standard"),
                )
            )
    return records
