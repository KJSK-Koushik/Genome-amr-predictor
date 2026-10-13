from __future__ import annotations

import csv
import http.client
import json
import re
import shutil
import ssl
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections import defaultdict
from collections.abc import Callable, Sequence
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Any

from genome_amr_predictor.data.fasta import inspect_fasta

_FTP_SSL_CONTEXT = ssl.create_default_context()


def datasets_version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "version"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0])
    for row in rows[1:]:
        fieldnames.extend(field for field in row if field not in fieldnames)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _run(command: list[str], *, retries: int = 0) -> None:
    for attempt in range(retries + 1):
        result = subprocess.run(command, check=False)
        if result.returncode == 0:
            return
        if attempt < retries:
            time.sleep(min(2**attempt, 30))
    raise RuntimeError(
        f"Command failed after {retries + 1} attempts with exit code "
        f"{result.returncode}: {command[0]} {command[1]}"
    )


def _extract_safely(archive_path: Path, destination: Path) -> None:
    root = destination.resolve()
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if not target.is_relative_to(root):
                raise ValueError(f"Unsafe path in NCBI data package: {member.filename}")
        archive.extractall(destination)


def _shards(values: Sequence[str], size: int) -> list[Sequence[str]]:
    return [values[start : start + size] for start in range(0, len(values), size)]


def _ftp_url(accession: str, filename: str) -> str:
    match = re.fullmatch(r"(GC[AF])_(\d{9})\.\d+", accession)
    if match is None or not filename.startswith(accession + "_"):
        raise ValueError("Invalid genome FTP accession or filename")
    prefix, digits = match.groups()
    folder = filename.removesuffix("_genomic.fna")
    components = [prefix, digits[:3], digits[3:6], digits[6:], folder, filename + ".gz"]
    return "https://ftp.ncbi.nlm.nih.gov/genomes/all/" + "/".join(
        urllib.parse.quote(component, safe="") for component in components
    )


def _fetch_ftp_fasta(task: tuple[str, Path]) -> bool:
    url, target = task
    if target.is_file():
        try:
            inspect_fasta(target)
            return True
        except (EOFError, OSError, ValueError):
            target.replace(target.with_name(target.name + f".invalid-{time.time_ns()}"))
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".download.gz")
    for attempt in range(5):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "GenomeAMRPredictor/0.1"})
            with (
                urllib.request.urlopen(request, timeout=60, context=_FTP_SSL_CONTEXT) as response,
                temporary.open("wb") as out,
            ):
                shutil.copyfileobj(response, out)
            inspect_fasta(temporary)
            temporary.replace(target)
            return True
        except (EOFError, OSError, ValueError, urllib.error.URLError, http.client.HTTPException):
            temporary.unlink(missing_ok=True)
            if attempt < 4:
                time.sleep(2**attempt)
    return False


def _prefetch_ftp(package_dir: Path, accessions: Sequence[str], max_workers: int) -> bool:
    tasks = []
    allowed = set(accessions)
    dataset_dir = package_dir / "ncbi_dataset"
    for line in (dataset_dir / "fetch.txt").read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            continue
        relative = Path(fields[2])
        if relative.parent.name not in allowed or not relative.name.endswith("_genomic.fna"):
            continue
        target = (dataset_dir / relative).with_name(relative.name + ".gz").resolve()
        if not target.is_relative_to(dataset_dir.resolve()):
            raise ValueError("Unsafe path in NCBI fetch manifest")
        tasks.append((_ftp_url(relative.parent.name, relative.name), target))
    with ThreadPoolExecutor(max_workers=min(max_workers, 8)) as executor:
        results = list(executor.map(_fetch_ftp_fasta, tasks))
    return len(tasks) == len(allowed) and all(results)


def _inspect_downloaded_accession(
    task: tuple[str, Path, Path],
) -> dict[str, Any]:
    accession, package_dir, output_dir = task
    assembly_dir = package_dir / "ncbi_dataset" / "data" / accession
    fasta_paths = sorted(assembly_dir.glob("*_genomic.fna.gz"))
    if len(fasta_paths) != 1:
        return {
            "assembly_accession": accession,
            "download_status": "missing_or_ambiguous_fasta",
            "error": f"found_{len(fasta_paths)}_files",
        }
    try:
        metrics = inspect_fasta(fasta_paths[0]).as_dict()
    except (EOFError, OSError, UnicodeError, ValueError) as exc:
        return {
            "assembly_accession": accession,
            "download_status": "invalid_fasta",
            "error": str(exc),
        }
    metrics["path"] = fasta_paths[0].relative_to(output_dir).as_posix()
    return {
        "assembly_accession": accession,
        "download_status": "verified",
        "error": "",
        **metrics,
    }


def download_genome_shards(
    accessions: Sequence[str],
    output_dir: Path,
    datasets_executable: Path,
    *,
    shard_size: int = 250,
    max_shards: int = 0,
    max_workers: int = 10,
    prefer_ftp: bool = False,
    progress: Callable[[int, int, str], None] | None = None,
) -> list[dict[str, Any]]:
    if not datasets_executable.is_file():
        raise FileNotFoundError(f"NCBI Datasets executable not found: {datasets_executable}")
    if shard_size < 1:
        raise ValueError("shard_size must be positive")
    if not 1 <= max_workers <= 30:
        raise ValueError("max_workers must be between 1 and 30")
    unique_accessions = sorted(
        dict.fromkeys(value.strip() for value in accessions if value.strip())
    )
    shards = _shards(unique_accessions, shard_size)
    if max_shards:
        shards = shards[:max_shards]
    output_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, Any]] = []
    inspector: ProcessPoolExecutor | None = None

    for index, shard in enumerate(shards, start=1):
        name = f"shard-{index:04d}"
        shard_dir = output_dir / name
        complete_path = shard_dir / "complete.json"
        metrics_path = shard_dir / "genome_metrics.jsonl"
        if complete_path.is_file() and metrics_path.is_file():
            cached_rows = _read_jsonl(metrics_path)
            if {row["assembly_accession"] for row in cached_rows} != set(shard):
                raise ValueError(f"Cached accession list differs for {name}; use a new output dir")
            all_rows.extend(cached_rows)
            if progress:
                progress(index, len(shards), "cached")
            continue

        shard_dir.mkdir(parents=True, exist_ok=True)
        accession_path = shard_dir / "accessions.txt"
        package_path = shard_dir / "dehydrated.zip"
        package_dir = shard_dir / "package"
        extracted_marker = shard_dir / "extracted.json"
        accession_path.write_text("\n".join(shard) + "\n", encoding="utf-8")

        if not package_path.is_file():
            _run(
                [
                    str(datasets_executable),
                    "download",
                    "genome",
                    "accession",
                    "--inputfile",
                    str(accession_path),
                    "--include",
                    "genome",
                    "--dehydrated",
                    "--filename",
                    str(package_path),
                    "--no-progressbar",
                ],
                retries=3,
            )
        if not extracted_marker.is_file():
            _extract_safely(package_path, package_dir)
            extracted_marker.write_text('{"status":"extracted"}\n', encoding="utf-8")
        ftp_complete = prefer_ftp and _prefetch_ftp(package_dir, shard, max_workers)
        if not ftp_complete:
            _run(
                [
                    str(datasets_executable),
                    "rehydrate",
                    "--directory",
                    str(package_dir),
                    "--gzip",
                    "--max-workers",
                    str(max_workers),
                    "--no-progressbar",
                ],
                retries=5,
            )

        if inspector is None:
            inspector = ProcessPoolExecutor(max_workers=min(max_workers, 8))
        rows = list(
            inspector.map(
                _inspect_downloaded_accession,
                ((accession, package_dir, output_dir) for accession in shard),
            )
        )
        invalid = [row for row in rows if row["download_status"] == "invalid_fasta"]
        if invalid:
            for row in invalid:
                assembly_dir = package_dir / "ncbi_dataset" / "data" / row["assembly_accession"]
                for fasta in assembly_dir.glob("*_genomic.fna.gz"):
                    # Preserve damaged downloads for diagnosis; rehydrate replaces missing files.
                    fasta.replace(fasta.with_name(fasta.name + f".invalid-{time.time_ns()}"))
            _run(
                [
                    str(datasets_executable),
                    "rehydrate",
                    "--directory",
                    str(package_dir),
                    "--gzip",
                    "--max-workers",
                    str(max_workers),
                    "--no-progressbar",
                ],
                retries=5,
            )
            rows = list(
                inspector.map(
                    _inspect_downloaded_accession,
                    ((accession, package_dir, output_dir) for accession in shard),
                )
            )
        with metrics_path.open("w", encoding="utf-8") as stream:
            for row in rows:
                stream.write(json.dumps(row, sort_keys=True) + "\n")
        completion = {
            "accessions_requested": len(shard),
            "fastas_verified": sum(row["download_status"] == "verified" for row in rows),
        }
        completion_path = (
            complete_path
            if completion["fastas_verified"] == completion["accessions_requested"]
            else shard_dir / "incomplete.json"
        )
        completion_path.write_text(
            json.dumps(
                completion,
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        if completion_path == complete_path:
            (shard_dir / "incomplete.json").unlink(missing_ok=True)
        all_rows.extend(rows)
        if progress:
            progress(index, len(shards), "downloaded")

    if inspector is not None:
        inspector.shutdown()
    for row in all_rows:
        path_value = str(row.get("path") or "")
        if not path_value:
            continue
        candidate = Path(path_value)
        if candidate.is_absolute():
            try:
                row["path"] = candidate.relative_to(output_dir.resolve()).as_posix()
            except ValueError:
                pass
    _write_csv(output_dir / "genome_manifest.csv", all_rows)
    duplicate_rows = _duplicate_sequence_rows(all_rows)
    _write_csv(output_dir / "exact_sequence_duplicates.csv", duplicate_rows)
    return all_rows


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as stream:
        return [json.loads(line) for line in stream if line.strip()]


def _duplicate_sequence_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_digest: defaultdict[str, list[str]] = defaultdict(list)
    for row in rows:
        digest = str(row.get("sequence_sha256") or "")
        if digest:
            by_digest[digest].append(str(row["assembly_accession"]))
    return [
        {
            "sequence_sha256": digest,
            "assembly_count": len(accessions),
            "assembly_accessions": ";".join(sorted(accessions)),
        }
        for digest, accessions in sorted(by_digest.items())
        if len(accessions) > 1
    ]
