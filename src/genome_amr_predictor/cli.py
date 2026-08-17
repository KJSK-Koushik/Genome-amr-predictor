from __future__ import annotations

import argparse
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

from genome_amr_predictor import __version__
from genome_amr_predictor.config import AuditConfig
from genome_amr_predictor.data.audit import audit_records, write_audit_outputs
from genome_amr_predictor.data.file_import import read_ast_file
from genome_amr_predictor.data.ncbi_biosample import (
    AST_DOCS_URL,
    BIOSAMPLE_URL,
    EUTILS_BASE,
    NCBIClient,
    parse_biosample_xml,
)
from genome_amr_predictor.data.schema import ASTRecord


def _positive_or_zero(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or positive")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genome-amr",
        description="Genome AMR Predictor research utilities.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)

    ncbi = commands.add_parser(
        "audit-ncbi",
        help="Audit structured NCBI BioSample antibiograms.",
    )
    ncbi.add_argument("--config", type=Path, required=True)
    ncbi.add_argument("--output-dir", type=Path, required=True)
    ncbi.add_argument(
        "--max-records",
        type=_positive_or_zero,
        default=1000,
        help="Maximum BioSamples to process; zero means all matching records.",
    )
    ncbi.add_argument("--write-records", action="store_true")
    ncbi.add_argument("--email", default=None, help="Optional NCBI contact email.")
    ncbi.add_argument("--api-key", default=None, help=argparse.SUPPRESS)

    file_command = commands.add_parser(
        "audit-file",
        help="Audit a CSV or TSV exported from the NCBI AST Browser.",
    )
    file_command.add_argument("--input", type=Path, required=True)
    file_command.add_argument("--config", type=Path, required=True)
    file_command.add_argument("--output-dir", type=Path, required=True)
    file_command.add_argument("--write-records", action="store_true")
    return parser


def _run_audit(
    records: list[ASTRecord],
    config: AuditConfig,
    output_dir: Path,
    manifest: dict[str, object],
    *,
    write_records: bool,
) -> int:
    if not records:
        print("No antibiogram rows were parsed; no audit was written.", file=sys.stderr)
        return 2
    bundle = audit_records(records, config)
    write_audit_outputs(
        output_dir,
        bundle,
        manifest,
        records=records if write_records else None,
    )
    print(json.dumps(bundle.summary, indent=2))
    print(f"Audit written to {output_dir.resolve()}")
    return 0


def audit_ncbi(args: argparse.Namespace) -> int:
    config = AuditConfig.load(args.config)
    client = NCBIClient(email=args.email, api_key=args.api_key)
    retrieved_at = datetime.now(UTC).isoformat()
    print(f"Searching NCBI BioSample: {config.entrez_query}")
    search = client.search(config.entrez_query, max_records=args.max_records)
    print(f"NCBI reports {search.count:,} matching BioSamples; processing {len(search.ids):,}.")

    records: list[ASTRecord] = []
    batches = (len(search.ids) + config.batch_size - 1) // config.batch_size
    for index, payload in enumerate(
        client.fetch_batches(search.ids, batch_size=config.batch_size), start=1
    ):
        batch_records = parse_biosample_xml(payload)
        records.extend(batch_records)
        print(
            f"Parsed batch {index}/{batches}: {len(batch_records):,} AST rows "
            f"({len(records):,} total)."
        )

    processed_accessions = {record.biosample_accession for record in records}
    manifest: dict[str, object] = {
        "source": "NCBI BioSample structured antibiograms",
        "source_urls": [BIOSAMPLE_URL, AST_DOCS_URL, EUTILS_BASE],
        "retrieved_at_utc": retrieved_at,
        "query": config.entrez_query,
        "translated_query": search.translated_query,
        "ncbi_reported_biosamples": search.count,
        "requested_biosamples": len(search.ids),
        "processed_biosamples": len(processed_accessions),
        "parsed_ast_rows": len(records),
        "raw_xml_retained": False,
        "config": config.as_dict(),
        "config_sha256": config.sha256(),
        "package_version": __version__,
    }
    return _run_audit(
        records,
        config,
        args.output_dir,
        manifest,
        write_records=args.write_records,
    )


def audit_file(args: argparse.Namespace) -> int:
    config = AuditConfig.load(args.config)
    input_path = args.input.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"AST input file not found: {input_path}")
    records = read_ast_file(input_path)
    manifest: dict[str, object] = {
        "source": "User-provided AST Browser-compatible file",
        "source_file": str(input_path),
        "retrieved_at_utc": datetime.now(UTC).isoformat(),
        "query": "file import",
        "translated_query": "file import",
        "ncbi_reported_biosamples": len(
            {record.biosample_accession for record in records if record.biosample_accession}
        ),
        "requested_biosamples": len(
            {record.biosample_accession for record in records if record.biosample_accession}
        ),
        "processed_biosamples": len(
            {record.biosample_accession for record in records if record.biosample_accession}
        ),
        "parsed_ast_rows": len(records),
        "raw_xml_retained": False,
        "config": config.as_dict(),
        "config_sha256": config.sha256(),
        "package_version": __version__,
    }
    return _run_audit(
        records,
        config,
        args.output_dir,
        manifest,
        write_records=args.write_records,
    )


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "audit-ncbi":
            return audit_ncbi(args)
        if args.command == "audit-file":
            return audit_file(args)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    parser.error(f"Unknown command: {args.command}")
    return 2
