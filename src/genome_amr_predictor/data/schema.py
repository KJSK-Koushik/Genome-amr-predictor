from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from typing import Any


@dataclass(frozen=True, slots=True)
class ASTRecord:
    biosample_accession: str
    organism: str
    taxonomy_id: str
    bioproject_accession: str
    sra_accession: str
    isolate: str
    collection_date: str
    geo_loc_name: str
    isolation_source: str
    host: str
    antibiotic: str
    phenotype_raw: str
    phenotype_normalized: str
    binary_label: int | None
    measurement_sign: str
    measurement: str
    measurement_units: str
    laboratory_method: str
    laboratory_platform: str
    vendor: str
    method_version_or_reagent: str
    testing_standard: str

    @classmethod
    def column_names(cls) -> tuple[str, ...]:
        return tuple(field.name for field in fields(cls))

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


REQUIRED_COLUMNS = {
    "biosample_accession",
    "organism",
    "antibiotic",
    "phenotype_raw",
}
