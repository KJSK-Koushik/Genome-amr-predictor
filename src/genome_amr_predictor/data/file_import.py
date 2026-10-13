from __future__ import annotations

import csv
import re
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

MALFORMED_QUOTED_FIELD = re.compile(r',""([^"\r\n]*)"",')


def _normalized_key(value: str) -> str:
    return value.strip().lower().lstrip("#").replace("-", "_")


NORMALIZED_ALIASES = {
    field: tuple(_normalized_key(alias) for alias in aliases) for field, aliases in ALIASES.items()
}


def _repair_export_quoting(line: str) -> str:
    """Repair doubled outer quotes observed in some NCBI AST Browser rows."""
    return MALFORMED_QUOTED_FIELD.sub(lambda match: f',"{match.group(1)}",', line)


def _matches_species(organism: str, species: str | None) -> bool:
    if not species:
        return True
    normalized_organism = " ".join(organism.casefold().split())
    normalized_species = " ".join(species.casefold().split())
    return normalized_organism == normalized_species or normalized_organism.startswith(
        f"{normalized_species} "
    )


def _get(row: dict[str, str], field: str) -> str:
    for alias in NORMALIZED_ALIASES[field]:
        value = row.get(alias, "")
        if value:
            return value
    return ""


def _detect_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".tab"}:
        return "\t"
    sample = path.read_text(encoding="utf-8-sig", errors="replace")[:8192]
    return csv.Sniffer().sniff(sample, delimiters=",\t").delimiter


def read_ast_file(path: Path, *, species: str | None = None) -> list[ASTRecord]:
    delimiter = _detect_delimiter(path)
    records: list[ASTRecord] = []
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(
            (_repair_export_quoting(line) for line in stream),
            delimiter=delimiter,
        )
        if not reader.fieldnames:
            raise ValueError("AST input file has no header row")
        for row in reader:
            if None in row:
                raise ValueError(f"Malformed AST row at line {reader.line_num}")
            normalized_row = {
                _normalized_key(key): (value or "").strip() for key, value in row.items()
            }
            organism = _get(normalized_row, "organism")
            if not _matches_species(organism, species):
                continue
            antibiotic = normalize_antibiotic(_get(normalized_row, "antibiotic"))
            phenotype_raw = _get(normalized_row, "phenotype_raw")
            if not antibiotic or not phenotype_raw:
                continue
            phenotype_normalized, binary_label = normalize_phenotype(phenotype_raw)
            records.append(
                ASTRecord(
                    biosample_accession=_get(normalized_row, "biosample_accession"),
                    organism=organism,
                    taxonomy_id=_get(normalized_row, "taxonomy_id"),
                    bioproject_accession=_get(normalized_row, "bioproject_accession"),
                    sra_accession=_get(normalized_row, "sra_accession"),
                    isolate=_get(normalized_row, "isolate"),
                    collection_date=_get(normalized_row, "collection_date"),
                    geo_loc_name=_get(normalized_row, "geo_loc_name"),
                    isolation_source=_get(normalized_row, "isolation_source"),
                    host=_get(normalized_row, "host"),
                    antibiotic=antibiotic,
                    phenotype_raw=phenotype_raw,
                    phenotype_normalized=phenotype_normalized,
                    binary_label=binary_label,
                    measurement_sign=_get(normalized_row, "measurement_sign"),
                    measurement=_get(normalized_row, "measurement"),
                    measurement_units=_get(normalized_row, "measurement_units"),
                    laboratory_method=_get(normalized_row, "laboratory_method"),
                    laboratory_platform=_get(normalized_row, "laboratory_platform"),
                    vendor=_get(normalized_row, "vendor"),
                    method_version_or_reagent=_get(normalized_row, "method_version_or_reagent"),
                    testing_standard=_get(normalized_row, "testing_standard"),
                )
            )
    return records
