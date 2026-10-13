from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class AssemblyQCConfig:
    minimum_length: int
    maximum_length: int
    maximum_contigs: int
    minimum_contig_n50: int
    minimum_checkm_completeness: float
    maximum_checkm_contamination: float
    require_checkm: bool
    reject_ani_mismatch: bool


@dataclass(frozen=True, slots=True)
class SplitConfig:
    seed: str
    train_fraction: float
    validation_fraction: float
    test_fraction: float


@dataclass(frozen=True, slots=True)
class Phase2Config:
    species: str
    antibiotics: tuple[str, ...]
    assembly_request_batch_size: int
    assembly_qc: AssemblyQCConfig
    splits: SplitConfig

    @classmethod
    def load(cls, path: Path) -> Phase2Config:
        raw = json.loads(path.read_text(encoding="utf-8"))
        qc = AssemblyQCConfig(**raw["assembly_qc"])
        splits = SplitConfig(**raw["splits"])
        config = cls(
            species=str(raw["species"]).strip(),
            antibiotics=tuple(str(value).strip().lower() for value in raw["antibiotics"]),
            assembly_request_batch_size=int(raw.get("assembly_request_batch_size", 50)),
            assembly_qc=qc,
            splits=splits,
        )
        config.validate()
        return config

    def validate(self) -> None:
        if not self.species:
            raise ValueError("species is required")
        if not self.antibiotics:
            raise ValueError("antibiotics cannot be empty")
        if len(set(self.antibiotics)) != len(self.antibiotics):
            raise ValueError("antibiotics must be unique")
        if not 1 <= self.assembly_request_batch_size <= 100:
            raise ValueError("assembly_request_batch_size must be between 1 and 100")
        if self.assembly_qc.minimum_length >= self.assembly_qc.maximum_length:
            raise ValueError("minimum_length must be smaller than maximum_length")
        if self.assembly_qc.minimum_length < 1:
            raise ValueError("minimum_length must be positive")
        if self.assembly_qc.maximum_contigs < 1:
            raise ValueError("maximum_contigs must be positive")
        if self.assembly_qc.minimum_contig_n50 < 1:
            raise ValueError("minimum_contig_n50 must be positive")
        if not 0 <= self.assembly_qc.minimum_checkm_completeness <= 100:
            raise ValueError("minimum_checkm_completeness must be between 0 and 100")
        if not 0 <= self.assembly_qc.maximum_checkm_contamination <= 100:
            raise ValueError("maximum_checkm_contamination must be between 0 and 100")
        fractions = (
            self.splits.train_fraction,
            self.splits.validation_fraction,
            self.splits.test_fraction,
        )
        if any(value <= 0 for value in fractions):
            raise ValueError("all split fractions must be positive")
        if abs(sum(fractions) - 1.0) > 1e-9:
            raise ValueError("split fractions must sum to 1")
        if not self.splits.seed:
            raise ValueError("split seed is required")

    def as_dict(self) -> dict[str, Any]:
        return {
            "species": self.species,
            "antibiotics": list(self.antibiotics),
            "assembly_request_batch_size": self.assembly_request_batch_size,
            "assembly_qc": asdict(self.assembly_qc),
            "splits": asdict(self.splits),
        }

    def sha256(self) -> str:
        canonical = json.dumps(self.as_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()
