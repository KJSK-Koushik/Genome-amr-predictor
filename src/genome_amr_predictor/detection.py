"""AMRFinderPlus nucleotide adapter and lossless, versioned evidence normalization."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import tempfile
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path, PureWindowsPath
from typing import Any

from genome_amr_predictor import __version__
from genome_amr_predictor.data.fasta import inspect_fasta

SCHEMA_VERSION = "amrfinder-evidence-v1"
ONTOLOGY_VERSION = "ncbi-class-preservation-v1"
WARNING = "Detection evidence does not establish phenotypic resistance or susceptibility."
ALIASES = {
    "Gene symbol": "Element symbol",
    "Sequence name": "Element name",
    "Accession of closest sequence": "Closest reference accession",
    "Name of closest sequence": "Closest reference name",
}
REQUIRED = {
    "Element symbol",
    "Element name",
    "Scope",
    "Type",
    "Subtype",
    "Class",
    "Subclass",
    "Method",
    "Target length",
    "Reference sequence length",
    "% Coverage of reference",
    "% Identity to reference",
    "Alignment length",
    "Closest reference accession",
    "Closest reference name",
    "HMM accession",
    "HMM description",
}
KNOWN_METHODS = {
    "ALLELEP",
    "ALLELEX",
    "EXACTP",
    "EXACTX",
    "BLASTP",
    "BLASTX",
    "PARTIALP",
    "PARTIALX",
    "PARTIAL_CONTIG_ENDP",
    "PARTIAL_CONTIG_ENDX",
    "INTERNAL_STOP",
    "HMM",
    "POINTP",
    "POINTX",
    "POINTN",
    "ALLELE",
    "EXACT",
    "BLAST",
    "PARTIAL",
    "PARTIAL_CONTIG_END",
    "POINT",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _number(value: str, *, integer: bool = False, percent: bool = False) -> int | float | None:
    if value in {"", "NA"}:
        return None
    result = int(value) if integer else float(value)
    if not math.isfinite(result) or result < 0 or (percent and result > 100):
        raise ValueError(f"Invalid detector metric: {value}")
    return result


def parse_amrfinder(path: Path, sample_id: str) -> list[dict[str, Any]]:
    """Parse header-only no-hit results; reject truncation instead of treating it as no hits."""
    if not sample_id.strip():
        raise ValueError("sample_id must be nonempty")
    hits: list[dict[str, Any]] = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        original_header = next(reader, [])
        header = [ALIASES.get(name, name) for name in original_header]
        if len(set(header)) != len(header):
            raise ValueError("Duplicate AMRFinderPlus columns")
        missing = REQUIRED - set(header)
        if missing:
            raise ValueError(f"Missing AMRFinderPlus columns: {', '.join(sorted(missing))}")
        for line, values in enumerate(reader, start=2):
            if not values:
                continue
            if len(values) != len(header):
                raise ValueError(f"Malformed AMRFinderPlus row {line}: wrong column count")
            row = dict(zip(header, values, strict=True))
            if any(
                not row[key] or row[key] == "NA"
                for key in ("Element symbol", "Element name", "Scope", "Type", "Subtype", "Method")
            ):
                raise ValueError(f"Missing evidence identifier at row {line}")
            start = _number(row.get("Start", "NA"), integer=True)
            stop = _number(row.get("Stop", "NA"), integer=True)
            strand = row.get("Strand", "NA")
            if (start is None) != (stop is None) or (
                start is not None and (start < 1 or stop < start or strand not in {"+", "-"})
            ):
                raise ValueError(f"Invalid 1-based coordinates at row {line}")
            method = row["Method"]
            warnings = []
            if method.startswith("PARTIAL"):
                warnings.append("partial_hit")
            if method == "INTERNAL_STOP":
                warnings.append("internal_stop")
            if method not in KNOWN_METHODS:
                warnings.append("unrecognized_method")
            if row["Scope"] != "core":
                warnings.append("plus_or_unrecognized_scope")
            if row["Type"] not in {"AMR", "STRESS", "VIRULENCE"}:
                warnings.append("unrecognized_element_type")
            kind = "mutation" if row["Subtype"] in {"POINT", "POINT_DISRUPT"} else "gene"
            if row["Type"] != "AMR":
                kind = "other_element"
            # The native labels are the curated ontology; no drug-specific phenotype is inferred.
            hits.append(
                {
                    "sample_id": sample_id,
                    "element_symbol": row["Element symbol"],
                    "element_name": row["Element name"],
                    "evidence_kind": kind,
                    "protein_id": row.get("Protein id", "NA"),
                    "contig_id": row.get("Contig id", "NA"),
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                    "coordinate_system": "1-based-inclusive",
                    "scope": row["Scope"],
                    "element_type": row["Type"],
                    "element_subtype": row["Subtype"],
                    "drug_class": row["Class"],
                    "drug_subclass": row["Subclass"],
                    "hierarchy_nodes": row.get("Hierarchy node", "").split("::")
                    if row.get("Hierarchy node") not in {None, "", "NA"}
                    else [],
                    "mechanism": None,
                    "mechanism_status": "not_separately_reported_by_detector",
                    "ontology_source": "NCBI AMRFinderPlus",
                    "method": method,
                    "allele_symbol": row["Element symbol"] if method.startswith("ALLELE") else None,
                    "target_length": _number(row["Target length"], integer=True),
                    "reference_length": _number(row["Reference sequence length"], integer=True),
                    "coverage_percent": _number(row["% Coverage of reference"], percent=True),
                    "identity_percent": _number(row["% Identity to reference"], percent=True),
                    "alignment_length": _number(row["Alignment length"], integer=True),
                    "reference_accession": row["Closest reference accession"],
                    "reference_name": row["Closest reference name"],
                    "hmm_accession": row["HMM accession"],
                    "hmm_description": row["HMM description"],
                    "warnings": warnings,
                    "raw": dict(zip(original_header, values, strict=True)),
                }
            )
    return hits


def write_evidence(path: Path, hits: list[dict[str, Any]], provenance: dict[str, Any]) -> None:
    document = {
        "schema_version": SCHEMA_VERSION,
        "ontology_version": ONTOLOGY_VERSION,
        "package_version": __version__,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "status": "completed",
        "warning": WARNING,
        "no_hit_interpretation": "No supported determinant detected; susceptibility is unknown.",
        "summary": {
            "total_hits": len(hits),
            "amr_hits": sum(hit["element_type"] == "AMR" for hit in hits),
            "mutation_hits": sum(hit["evidence_kind"] == "mutation" for hit in hits),
            "methods": dict(sorted(Counter(hit["method"] for hit in hits).items())),
        },
        "provenance": provenance,
        "hits": hits,
    }
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


class AMRFinderRuntime:
    def __init__(self, executable: Path, database: Path, wsl_distribution: str | None = None):
        self.executable = executable.resolve(strict=True)
        self.database = database.resolve(strict=True)
        self.wsl_distribution = wsl_distribution
        if not self.executable.is_file() or not self.database.is_dir():
            raise ValueError("An executable file and explicit database directory are required")

    def runtime_path(self, path: Path) -> str:
        if not self.wsl_distribution:
            return str(path.resolve())
        windows_path = PureWindowsPath(str(path.resolve()))
        if windows_path.drive.lower().startswith(("\\\\wsl.localhost\\", "\\\\wsl$\\")):
            distribution = windows_path.drive.rsplit("\\", 1)[-1]
            if distribution.casefold() != self.wsl_distribution.casefold():
                raise ValueError("WSL filesystem path belongs to a different distribution")
            return "/" + "/".join(windows_path.parts[1:])
        try:
            result = subprocess.run(
                [
                    "wsl",
                    "-d",
                    self.wsl_distribution,
                    "--",
                    "wslpath",
                    "-a",
                    "-u",
                    str(path.resolve()),
                ],
                capture_output=True,
                text=True,
                check=True,
                timeout=30,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(f"WSL could not access path {path}: {exc}") from exc
        return result.stdout.strip()

    def command(self, arguments: list[str]) -> list[str]:
        command = [self.runtime_path(self.executable), *arguments]
        if self.wsl_distribution:
            command = ["wsl", "-d", self.wsl_distribution, "--", *command]
        return command

    def invoke(self, arguments: list[str], timeout: int) -> subprocess.CompletedProcess[str]:
        try:
            result = subprocess.run(
                self.command(arguments),
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise RuntimeError(f"AMRFinderPlus could not execute: {exc}") from exc
        if result.returncode:
            raise RuntimeError(f"AMRFinderPlus exited {result.returncode}: {result.stderr.strip()}")
        return result


def detect_genome(
    fasta: Path,
    output_dir: Path,
    sample_id: str,
    runtime: AMRFinderRuntime,
    *,
    threads: int = 4,
    timeout: int = 3600,
) -> dict[str, Any]:
    if threads < 1 or timeout < 1:
        raise ValueError("threads and timeout must be positive")
    metrics = inspect_fasta(fasta)
    # Duplicate contig identifiers make location evidence ambiguous.
    opener = gzip.open if fasta.suffix.lower() == ".gz" else open
    with opener(fasta, "rt", encoding="utf-8") as stream:
        ids = [
            line[1:].split()[0] if line[1:].split() else ""
            for line in stream
            if line.startswith(">")
        ]
    if "" in ids or len(ids) != len(set(ids)):
        raise ValueError("FASTA contig identifiers must be nonempty and unique")
    if not sample_id.strip():
        raise ValueError("sample_id must be nonempty")
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise ValueError("Detection output directory must be empty; use a new run directory")
    raw = output_dir / "amrfinder.tsv"
    manifest_path = output_dir / "run_manifest.json"
    manifest: dict[str, Any] = {
        "status": "running",
        "sample_id": sample_id,
        "started_at_utc": datetime.now(UTC).isoformat(),
        "package_version": __version__,
        "schema_version": SCHEMA_VERSION,
        "input": metrics.as_dict(),
        "executable": str(runtime.executable),
        "executable_sha256": sha256(runtime.executable),
        "database": str(runtime.database),
        "wsl_distribution": runtime.wsl_distribution,
        "organism": "Escherichia",
        "mode": "nucleotide",
        "plus": False,
        "threads": threads,
    }
    try:
        version = runtime.invoke(["--version"], 60)
        manifest["software_version"] = version.stdout.strip()
        db_version = runtime.invoke(["-d", runtime.runtime_path(runtime.database), "-V"], 60)
        manifest["database_version_report"] = (db_version.stdout + db_version.stderr).strip()
        if not manifest["software_version"] or not manifest["database_version_report"]:
            raise RuntimeError("AMRFinderPlus did not report software/database versions")
        database_hashes = {
            p.relative_to(runtime.database).as_posix(): sha256(p)
            for p in sorted(runtime.database.rglob("*"))
            if p.is_file()
        }
        if not database_hashes:
            raise ValueError("AMRFinderPlus database is empty")
        manifest["database_files_sha256"] = database_hashes
        manifest["adapter_sha256"] = sha256(Path(__file__))
        with tempfile.TemporaryDirectory(prefix="genome-amr-", dir=output_dir) as temporary:
            input_path = fasta.resolve()
            if fasta.suffix.lower() == ".gz":
                input_path = Path(temporary) / "genome.fna"
                with gzip.open(fasta, "rb") as source, input_path.open("wb") as target:
                    shutil.copyfileobj(source, target)
            arguments = [
                "-n",
                runtime.runtime_path(input_path),
                "-d",
                runtime.runtime_path(runtime.database),
                "-O",
                "Escherichia",
                "--threads",
                str(threads),
                "--print_node",
                "-o",
                runtime.runtime_path(raw),
            ]
            manifest["command"] = runtime.command(arguments)
            result = runtime.invoke(arguments, timeout)
        (output_dir / "amrfinder.stderr.log").write_text(result.stderr, encoding="utf-8")
        hits = parse_amrfinder(raw, sample_id)
        for hit in hits:
            if hit["contig_id"] not in ids:
                raise ValueError(f"Detector reported unknown contig: {hit['contig_id']}")
        if sha256(fasta) != metrics.file_sha256:
            raise ValueError("Input FASTA changed during detection")
        manifest["status"] = "completed"
        manifest["total_hits"] = len(hits)
        manifest["raw_output_sha256"] = sha256(raw)
        write_evidence(output_dir / "detections.json", hits, manifest)
        manifest["evidence_sha256"] = sha256(output_dir / "detections.json")
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
        raise
    finally:
        manifest["finished_at_utc"] = datetime.now(UTC).isoformat()
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    return manifest
