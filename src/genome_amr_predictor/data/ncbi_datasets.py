from __future__ import annotations

import http.client
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

DATASETS_API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
DATASETS_DOCS_URL = "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/"


def _integer(value: Any) -> int | None:
    if value in (None, ""):
        return None
    return int(value)


def _number(value: Any) -> float | None:
    if value in (None, ""):
        return None
    return float(value)


@dataclass(frozen=True, slots=True)
class AssemblyMetadata:
    biosample_accession: str
    assembly_accession: str
    paired_accession: str
    source_database: str
    organism: str
    tax_id: int | None
    assembly_level: str
    assembly_name: str
    assembly_status: str
    bioproject_accession: str
    release_date: str
    total_sequence_length: int | None
    total_ungapped_length: int | None
    number_of_contigs: int | None
    contig_n50: int | None
    gc_percent: float | None
    sequencing_technology: str
    assembly_method: str
    checkm_completeness: float | None
    checkm_contamination: float | None
    ani_taxonomy_check_status: str
    ani_match_status: str
    annotation_name: str

    @classmethod
    def from_report(cls, report: dict[str, Any]) -> AssemblyMetadata:
        organism = report.get("organism") or {}
        assembly_info = report.get("assembly_info") or {}
        biosample = assembly_info.get("biosample") or {}
        stats = report.get("assembly_stats") or {}
        checkm = report.get("checkm_info") or {}
        ani = report.get("average_nucleotide_identity") or {}
        annotation = report.get("annotation_info") or {}
        return cls(
            biosample_accession=str(biosample.get("accession") or ""),
            assembly_accession=str(report.get("accession") or ""),
            paired_accession=str(report.get("paired_accession") or ""),
            source_database=str(report.get("source_database") or ""),
            organism=str(organism.get("organism_name") or ""),
            tax_id=_integer(organism.get("tax_id")),
            assembly_level=str(assembly_info.get("assembly_level") or ""),
            assembly_name=str(assembly_info.get("assembly_name") or ""),
            assembly_status=str(assembly_info.get("assembly_status") or ""),
            bioproject_accession=str(assembly_info.get("bioproject_accession") or ""),
            release_date=str(assembly_info.get("release_date") or ""),
            total_sequence_length=_integer(stats.get("total_sequence_length")),
            total_ungapped_length=_integer(stats.get("total_ungapped_length")),
            number_of_contigs=_integer(stats.get("number_of_contigs")),
            contig_n50=_integer(stats.get("contig_n50")),
            gc_percent=_number(stats.get("gc_percent")),
            sequencing_technology=str(assembly_info.get("sequencing_tech") or ""),
            assembly_method=str(assembly_info.get("assembly_method") or ""),
            checkm_completeness=_number(checkm.get("completeness")),
            checkm_contamination=_number(checkm.get("contamination")),
            ani_taxonomy_check_status=str(ani.get("taxonomy_check_status") or ""),
            ani_match_status=str(ani.get("match_status") or ""),
            annotation_name=str(annotation.get("name") or ""),
        )

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> AssemblyMetadata:
        return cls(**value)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _chunks(values: Sequence[str], size: int) -> Iterator[Sequence[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


class NCBIDatasetsClient:
    def __init__(
        self,
        *,
        api_key: str | None = None,
        timeout: float = 60.0,
        retries: int = 5,
        retry_delay: float = 1.0,
    ) -> None:
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.timeout = timeout
        self.retries = retries
        self.retry_delay = retry_delay

    def _request_json(self, url: str) -> dict[str, Any]:
        headers = {
            "Accept": "application/json",
            "User-Agent": "genome-amr-predictor/0.1 (research software)",
        }
        if self.api_key:
            headers["api-key"] = self.api_key
        request = urllib.request.Request(url, headers=headers)
        retryable = (
            ConnectionError,
            TimeoutError,
            http.client.HTTPException,
            urllib.error.URLError,
        )
        for attempt in range(self.retries + 1):
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    return json.loads(response.read().decode("utf-8"))
            except urllib.error.HTTPError as exc:
                if exc.code not in {429, 500, 502, 503, 504} or attempt >= self.retries:
                    raise RuntimeError(
                        f"NCBI Datasets request failed with HTTP {exc.code}"
                    ) from exc
            except retryable as exc:
                if attempt >= self.retries:
                    raise RuntimeError("NCBI Datasets request failed after retries") from exc
            time.sleep(self.retry_delay * (2**attempt))
        raise RuntimeError("NCBI Datasets request failed")

    def assembly_reports(self, biosample_accessions: Sequence[str]) -> list[AssemblyMetadata]:
        if not biosample_accessions:
            return []
        encoded = urllib.parse.quote(",".join(biosample_accessions), safe=",")
        parameters = {
            "page_size": 1000,
            "filters.exclude_paired_reports": "true",
            "filters.exclude_atypical": "true",
            "filters.assembly_version": "current",
            "filters.is_metagenome_derived": "metagenome_derived_exclude",
        }
        reports: list[AssemblyMetadata] = []
        page_token = ""
        while True:
            if page_token:
                parameters["page_token"] = page_token
            query = urllib.parse.urlencode(parameters)
            url = f"{DATASETS_API_BASE}/genome/biosample/{encoded}/dataset_report?{query}"
            payload = self._request_json(url)
            reports.extend(
                AssemblyMetadata.from_report(report) for report in payload.get("reports", [])
            )
            page_token = str(payload.get("next_page_token") or "")
            if not page_token:
                return reports


_LEVEL_RANK = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1,
}


def select_preferred_assembly(candidates: Iterable[AssemblyMetadata]) -> AssemblyMetadata | None:
    values = [candidate for candidate in candidates if candidate.assembly_accession]
    if not values:
        return None

    def rank(value: AssemblyMetadata) -> tuple[int, int, int, str, str]:
        return (
            int("REFSEQ" in value.source_database.upper()),
            _LEVEL_RANK.get(value.assembly_level, 0),
            int(bool(value.annotation_name)),
            value.release_date,
            value.assembly_accession,
        )

    return max(values, key=rank)


def load_resolution_cache(path: Path) -> dict[str, AssemblyMetadata | None]:
    resolved: dict[str, AssemblyMetadata | None] = {}
    if not path.is_file():
        return resolved
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            item = json.loads(line)
            assembly = item.get("assembly")
            resolved[item["query_biosample"]] = (
                AssemblyMetadata.from_dict(assembly) if assembly else None
            )
    return resolved


def resolve_assemblies(
    biosample_accessions: Sequence[str],
    client: NCBIDatasetsClient,
    cache_path: Path,
    *,
    batch_size: int,
    progress: Callable[[int, int, int], None] | None = None,
) -> dict[str, AssemblyMetadata | None]:
    requested = tuple(dict.fromkeys(accession for accession in biosample_accessions if accession))
    resolved = load_resolution_cache(cache_path)
    pending = [accession for accession in requested if accession not in resolved]
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    total_batches = (len(pending) + batch_size - 1) // batch_size
    with cache_path.open("a", encoding="utf-8") as stream:
        for index, batch in enumerate(_chunks(pending, batch_size), start=1):
            reports = client.assembly_reports(batch)
            by_biosample: dict[str, list[AssemblyMetadata]] = {}
            for report in reports:
                by_biosample.setdefault(report.biosample_accession, []).append(report)
            for accession in batch:
                selected = select_preferred_assembly(by_biosample.get(accession, []))
                resolved[accession] = selected
                stream.write(
                    json.dumps(
                        {
                            "query_biosample": accession,
                            "candidate_count": len(by_biosample.get(accession, [])),
                            "assembly": selected.as_dict() if selected else None,
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )
            stream.flush()
            if progress:
                progress(index, total_batches, len(resolved))
    return {accession: resolved.get(accession) for accession in requested}
