from __future__ import annotations

import re

_WHITESPACE = re.compile(r"\s+")

RESISTANT = {"r", "resistant"}
SUSCEPTIBLE = {"s", "susceptible", "sensitive"}
INTERMEDIATE = {"i", "intermediate"}


def normalize_text(value: str | None) -> str:
    if value is None:
        return ""
    return _WHITESPACE.sub(" ", value.strip().lower())


def normalize_antibiotic(value: str | None) -> str:
    return normalize_text(value).replace("\N{EN DASH}", "-").replace("\N{EM DASH}", "-")


def normalize_phenotype(value: str | None) -> tuple[str, int | None]:
    normalized = normalize_text(value)
    if normalized in RESISTANT:
        return "resistant", 1
    if normalized in SUSCEPTIBLE:
        return "susceptible", 0
    if normalized in INTERMEDIATE:
        return "intermediate", None
    return "excluded", None
