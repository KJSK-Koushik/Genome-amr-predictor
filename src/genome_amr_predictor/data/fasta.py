from __future__ import annotations

import gzip
import hashlib
from dataclasses import asdict, dataclass
from pathlib import Path

ALLOWED_BASE_BYTES = frozenset(b"ACGTNRYKMSWBDHV")


@dataclass(frozen=True, slots=True)
class FastaMetrics:
    path: str
    file_sha256: str
    sequence_sha256: str
    total_length: int
    contig_count: int
    contig_n50: int
    gc_percent: float
    ambiguous_bases: int
    ambiguous_percent: float

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _read_bytes(path: Path) -> bytes:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rb") as stream:
            return stream.read()
    return path.read_bytes()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _n50(lengths: list[int]) -> int:
    target = sum(lengths) / 2
    cumulative = 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= target:
            return length
    return 0


def inspect_fasta(path: Path) -> FastaMetrics:
    lengths: list[int] = []
    contig_hashes: list[str] = []
    total_gc = 0
    total_ambiguous = 0
    payload = _read_bytes(path).strip()
    if not payload:
        raise ValueError(f"FASTA contains no records: {path}")
    if not payload.startswith(b">"):
        raise ValueError("FASTA sequence appears before a header")

    for record_number, record in enumerate(payload.split(b"\n>"), start=1):
        if record_number > 1:
            record = b">" + record
        _header, separator, body = record.partition(b"\n")
        if not separator:
            raise ValueError(f"FASTA contains an empty sequence: {path}")
        sequence = b"".join(body.split()).upper()
        if not sequence:
            raise ValueError(f"FASTA contains an empty sequence: {path}")
        invalid = set(sequence) - ALLOWED_BASE_BYTES
        if invalid:
            rendered = "".join(chr(value) for value in sorted(invalid))
            raise ValueError(f"Invalid FASTA symbols in record {record_number}: {rendered}")
        length = len(sequence)
        adenine = sequence.count(b"A")
        cytosine = sequence.count(b"C")
        guanine = sequence.count(b"G")
        thymine = sequence.count(b"T")
        lengths.append(length)
        contig_hashes.append(hashlib.sha256(sequence).hexdigest())
        total_gc += cytosine + guanine
        total_ambiguous += length - adenine - cytosine - guanine - thymine

    total_length = sum(lengths)
    sequence_digest = hashlib.sha256("|".join(sorted(contig_hashes)).encode()).hexdigest()
    return FastaMetrics(
        path=str(path.resolve()),
        file_sha256=_file_sha256(path),
        sequence_sha256=sequence_digest,
        total_length=total_length,
        contig_count=len(lengths),
        contig_n50=_n50(lengths),
        gc_percent=round(100 * total_gc / total_length, 6),
        ambiguous_bases=total_ambiguous,
        ambiguous_percent=round(100 * total_ambiguous / total_length, 6),
    )
