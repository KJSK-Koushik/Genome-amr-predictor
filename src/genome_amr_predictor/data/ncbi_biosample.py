from __future__ import annotations

import http.client
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from typing import Any

from genome_amr_predictor.data.labels import normalize_antibiotic, normalize_phenotype
from genome_amr_predictor.data.schema import ASTRecord

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOSAMPLE_URL = "https://www.ncbi.nlm.nih.gov/biosample/"
AST_DOCS_URL = "https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/"


@dataclass(frozen=True, slots=True)
class SearchResult:
    count: int
    ids: tuple[str, ...]
    translated_query: str


class NCBIRequestError(RuntimeError):
    """Raised when an NCBI E-utilities request cannot be completed."""


class NCBIClient:
    def __init__(
        self,
        *,
        email: str | None = None,
        api_key: str | None = None,
        tool: str = "genome-amr-predictor",
        retries: int = 4,
        timeout_seconds: float = 60.0,
    ) -> None:
        self.email = email or os.getenv("NCBI_EMAIL", "")
        self.api_key = api_key or os.getenv("NCBI_API_KEY", "")
        self.tool = tool
        self.retries = retries
        self.timeout_seconds = timeout_seconds
        self._minimum_interval = 0.11 if self.api_key else 0.36
        self._last_request_at = 0.0

    def _wait_for_rate_limit(self) -> None:
        elapsed = time.monotonic() - self._last_request_at
        remaining = self._minimum_interval - elapsed
        if remaining > 0:
            time.sleep(remaining)

    def _request(self, endpoint: str, params: dict[str, Any]) -> bytes:
        common = {"tool": self.tool}
        if self.email:
            common["email"] = self.email
        if self.api_key:
            common["api_key"] = self.api_key
        encoded = urllib.parse.urlencode({**params, **common})
        base_url = f"{EUTILS_BASE}/{endpoint}"
        if len(encoded) > 1800:
            request = urllib.request.Request(
                base_url,
                data=encoded.encode("utf-8"),
                headers={
                    "Content-Type": "application/x-www-form-urlencoded",
                    "User-Agent": f"{self.tool}/0.1 (research data audit)",
                },
            )
        else:
            request = urllib.request.Request(
                f"{base_url}?{encoded}",
                headers={"User-Agent": f"{self.tool}/0.1 (research data audit)"},
            )
        last_error: Exception | None = None
        for attempt in range(self.retries):
            self._wait_for_rate_limit()
            try:
                with urllib.request.urlopen(request, timeout=self.timeout_seconds) as response:
                    payload = response.read()
                self._last_request_at = time.monotonic()
                return payload
            except (
                ConnectionError,
                TimeoutError,
                http.client.HTTPException,
                urllib.error.URLError,
            ) as exc:
                last_error = exc
                self._last_request_at = time.monotonic()
                if attempt + 1 < self.retries:
                    time.sleep(2**attempt)
        raise NCBIRequestError(f"NCBI request failed after {self.retries} attempts: {last_error}")

    def search(self, query: str, *, max_records: int = 0) -> SearchResult:
        initial = json.loads(
            self._request(
                "esearch.fcgi",
                {"db": "biosample", "term": query, "retmode": "json", "retmax": 0},
            )
        )
        result = initial["esearchresult"]
        count = int(result["count"])
        requested = count if max_records <= 0 else min(count, max_records)
        if requested == 0:
            return SearchResult(count=count, ids=(), translated_query=result["querytranslation"])
        payload = json.loads(
            self._request(
                "esearch.fcgi",
                {
                    "db": "biosample",
                    "term": query,
                    "retmode": "json",
                    "retmax": requested,
                },
            )
        )
        ids = tuple(payload["esearchresult"]["idlist"])
        return SearchResult(
            count=count,
            ids=ids,
            translated_query=payload["esearchresult"]["querytranslation"],
        )

    def fetch_batches(self, ids: Sequence[str], *, batch_size: int) -> Iterator[bytes]:
        for start in range(0, len(ids), batch_size):
            batch = ids[start : start + batch_size]
            yield self._request(
                "efetch.fcgi",
                {"db": "biosample", "id": ",".join(batch), "retmode": "xml"},
            )


def _attribute_map(sample: ET.Element) -> dict[str, str]:
    result: dict[str, str] = {}
    for attribute in sample.findall("./Attributes/Attribute"):
        names = (
            attribute.get("harmonized_name", ""),
            attribute.get("attribute_name", ""),
            attribute.get("display_name", ""),
        )
        value = (attribute.text or "").strip()
        for name in names:
            if name:
                result[name.strip().lower().replace(" ", "_")] = value
    return result


def _find_id(sample: ET.Element, *, db: str) -> str:
    for element in sample.findall("./Ids/Id"):
        if element.get("db", "").lower() == db.lower():
            return (element.text or "").strip()
    return ""


def _find_bioproject(sample: ET.Element) -> str:
    for link in sample.findall("./Links/Link"):
        if link.get("target", "").lower() == "bioproject":
            return link.get("label", "").strip()
    return ""


def _row_map(headers: list[str], cells: list[str]) -> dict[str, str]:
    return {headers[index]: value for index, value in enumerate(cells) if index < len(headers)}


def _first(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name, "")
        if value:
            return value
    return ""


def parse_biosample_xml(payload: bytes) -> list[ASTRecord]:
    root = ET.fromstring(payload)
    records: list[ASTRecord] = []
    for sample in root.findall("./BioSample"):
        accession = sample.get("accession", "").strip()
        organism_element = sample.find("./Description/Organism")
        organism = organism_element.get("taxonomy_name", "") if organism_element is not None else ""
        taxonomy_id = (
            organism_element.get("taxonomy_id", "") if organism_element is not None else ""
        )
        attributes = _attribute_map(sample)
        common = {
            "biosample_accession": accession,
            "organism": organism,
            "taxonomy_id": taxonomy_id,
            "bioproject_accession": _find_bioproject(sample),
            "sra_accession": _find_id(sample, db="SRA"),
            "isolate": attributes.get("isolate", ""),
            "collection_date": attributes.get("collection_date", ""),
            "geo_loc_name": attributes.get("geo_loc_name", ""),
            "isolation_source": attributes.get("isolation_source", ""),
            "host": attributes.get("host", ""),
        }
        for table in sample.findall("./Description/Comment/Table"):
            if not table.get("class", "").lower().startswith("antibiogram"):
                continue
            headers = [
                (cell.text or "").strip().lower().replace(" ", "_")
                for cell in table.findall("./Header/Cell")
            ]
            for row_element in table.findall("./Body/Row"):
                cells = [(cell.text or "").strip() for cell in row_element.findall("./Cell")]
                row = _row_map(headers, cells)
                antibiotic = normalize_antibiotic(_first(row, "antibiotic"))
                phenotype_raw = _first(row, "resistance_phenotype", "phenotype")
                phenotype_normalized, binary_label = normalize_phenotype(phenotype_raw)
                if not antibiotic or not phenotype_raw:
                    continue
                records.append(
                    ASTRecord(
                        **common,
                        antibiotic=antibiotic,
                        phenotype_raw=phenotype_raw,
                        phenotype_normalized=phenotype_normalized,
                        binary_label=binary_label,
                        measurement_sign=_first(row, "measurement_sign"),
                        measurement=_first(row, "measurement", "mic_(mg/l)", "mic"),
                        measurement_units=_first(row, "measurement_units", "units"),
                        laboratory_method=_first(row, "laboratory_typing_method"),
                        laboratory_platform=_first(row, "laboratory_typing_platform"),
                        vendor=_first(row, "vendor"),
                        method_version_or_reagent=_first(
                            row, "laboratory_typing_method_version_or_reagent"
                        ),
                        testing_standard=_first(row, "testing_standard"),
                    )
                )
    return records
