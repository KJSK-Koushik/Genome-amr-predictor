from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class SelectionThresholds:
    minimum_binary_records: int
    minimum_resistant_records: int
    minimum_susceptible_records: int
    minimum_minority_fraction: float

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> SelectionThresholds:
        thresholds = cls(
            minimum_binary_records=int(value["minimum_binary_records"]),
            minimum_resistant_records=int(value["minimum_resistant_records"]),
            minimum_susceptible_records=int(value["minimum_susceptible_records"]),
            minimum_minority_fraction=float(value["minimum_minority_fraction"]),
        )
        if thresholds.minimum_binary_records < 1:
            raise ValueError("minimum_binary_records must be positive")
        if thresholds.minimum_resistant_records < 1:
            raise ValueError("minimum_resistant_records must be positive")
        if thresholds.minimum_susceptible_records < 1:
            raise ValueError("minimum_susceptible_records must be positive")
        if not 0 < thresholds.minimum_minority_fraction <= 0.5:
            raise ValueError("minimum_minority_fraction must be in (0, 0.5]")
        return thresholds


@dataclass(frozen=True, slots=True)
class AuditConfig:
    species: str
    entrez_query: str
    candidate_antibiotics: tuple[str, ...]
    selection_thresholds: SelectionThresholds
    batch_size: int

    @classmethod
    def load(cls, path: Path) -> AuditConfig:
        raw = json.loads(path.read_text(encoding="utf-8"))
        config = cls(
            species=str(raw["species"]).strip(),
            entrez_query=str(raw["entrez_query"]).strip(),
            candidate_antibiotics=tuple(
                str(item).strip().lower() for item in raw["candidate_antibiotics"]
            ),
            selection_thresholds=SelectionThresholds.from_dict(raw["selection_thresholds"]),
            batch_size=int(raw.get("batch_size", 400)),
        )
        if not config.species:
            raise ValueError("species is required")
        if not config.entrez_query:
            raise ValueError("entrez_query is required")
        if not config.candidate_antibiotics:
            raise ValueError("candidate_antibiotics cannot be empty")
        if not 1 <= config.batch_size <= 500:
            raise ValueError("batch_size must be between 1 and 500")
        return config

    def as_dict(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "entrez_query": self.entrez_query,
            "candidate_antibiotics": list(self.candidate_antibiotics),
            "selection_thresholds": {
                "minimum_binary_records": self.selection_thresholds.minimum_binary_records,
                "minimum_resistant_records": self.selection_thresholds.minimum_resistant_records,
                "minimum_susceptible_records": (
                    self.selection_thresholds.minimum_susceptible_records
                ),
                "minimum_minority_fraction": (self.selection_thresholds.minimum_minority_fraction),
            },
            "batch_size": self.batch_size,
        }

    def sha256(self) -> str:
        canonical = json.dumps(self.as_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
