from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

from genome_amr_predictor import __version__
from genome_amr_predictor.cohort_config import Phase2Config
from genome_amr_predictor.config import AuditConfig
from genome_amr_predictor.data.audit import audit_records, write_audit_outputs
from genome_amr_predictor.data.cohort import (
    NormalizedLabels,
    build_cohort,
    normalize_training_labels,
    write_cohort_outputs,
)
from genome_amr_predictor.data.file_import import read_ast_file
from genome_amr_predictor.data.finalize import finalize_sequence_cohort
from genome_amr_predictor.data.genome_download import datasets_version, download_genome_shards
from genome_amr_predictor.data.ncbi_biosample import (
    AST_DOCS_URL,
    BIOSAMPLE_URL,
    EUTILS_BASE,
    NCBIClient,
    parse_biosample_xml,
)
from genome_amr_predictor.data.ncbi_datasets import (
    DATASETS_API_BASE,
    DATASETS_DOCS_URL,
    NCBIDatasetsClient,
    resolve_assemblies,
)
from genome_amr_predictor.data.schema import ASTRecord
from genome_amr_predictor.detection import (
    AMRFinderRuntime,
    detect_genome,
    parse_amrfinder,
    sha256,
    write_evidence,
)


def _positive_or_zero(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or positive")
    return parsed


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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

    cohort = commands.add_parser(
        "construct-cohort",
        help="Join AST labels to NCBI assemblies, apply metadata QC, and freeze splits.",
    )
    cohort.add_argument("--input", type=Path, required=True)
    cohort.add_argument("--config", type=Path, required=True)
    cohort.add_argument("--output-dir", type=Path, required=True)
    cohort.add_argument("--assembly-cache", type=Path, required=True)
    cohort.add_argument(
        "--max-biosamples",
        type=_positive_or_zero,
        default=0,
        help="Maximum labeled BioSamples to resolve; zero means all.",
    )
    cohort.add_argument("--api-key", default=None, help=argparse.SUPPRESS)

    genomes = commands.add_parser(
        "download-genomes",
        help="Download and verify genome FASTA files in resumable NCBI Datasets shards.",
    )
    genomes.add_argument("--accessions", type=Path, required=True)
    genomes.add_argument("--output-dir", type=Path, required=True)
    genomes.add_argument("--datasets-executable", type=Path, required=True)
    genomes.add_argument("--shard-size", type=int, default=250)
    genomes.add_argument("--max-shards", type=_positive_or_zero, default=0)
    genomes.add_argument("--max-workers", type=int, default=10)
    genomes.add_argument(
        "--prefer-ftp",
        action="store_true",
        help="Fetch compressed FASTA from NCBI HTTPS FTP, with Datasets fallback.",
    )

    finalize = commands.add_parser(
        "finalize-cohort",
        help="Join verified FASTA hashes to the metadata cohort and remove exact sequences.",
    )
    finalize.add_argument("--cohort", type=Path, required=True)
    finalize.add_argument("--assemblies", type=Path, required=True)
    finalize.add_argument("--genome-manifest", type=Path, required=True)
    finalize.add_argument("--output-dir", type=Path, required=True)
    finalize.add_argument("--genome-dir", type=Path, required=True)
    finalize.add_argument("--config", type=Path, required=True)
    finalize.add_argument("--max-workers", type=int, default=8)
    detect = commands.add_parser("detect", help="Detect curated AMR evidence in an E. coli FASTA.")
    detect.add_argument("--fasta", type=Path, required=True)
    detect.add_argument("--sample-id", required=True)
    detect.add_argument("--output-dir", type=Path, required=True)
    detect.add_argument("--executable", type=Path, required=True)
    detect.add_argument("--database", type=Path, required=True)
    detect.add_argument("--wsl-distribution", default=None)
    detect.add_argument("--threads", type=int, default=4)
    detect.add_argument("--timeout", type=int, default=3600)
    parse = commands.add_parser("parse-detection", help="Normalize an existing AMRFinderPlus TSV.")
    parse.add_argument("--input", type=Path, required=True)
    parse.add_argument("--sample-id", required=True)
    parse.add_argument("--software-version", required=True)
    parse.add_argument("--database-version", required=True)
    parse.add_argument("--output", type=Path, required=True)
    feature = commands.add_parser(
        "build-features", help="Verify cohort calls and build ML features."
    )
    feature.add_argument("--assemblies", type=Path, required=True)
    feature.add_argument("--detection-dir", type=Path, required=True)
    feature.add_argument("--output-dir", type=Path, required=True)
    train = commands.add_parser(
        "train-models", help="Fit and evaluate frozen-split research models."
    )
    train.add_argument("--features", type=Path, required=True)
    train.add_argument("--cohort", type=Path, required=True)
    train.add_argument("--similarity-groups", type=Path, required=True)
    train.add_argument("--config", type=Path, required=True)
    train.add_argument("--model-dir", type=Path, required=True)
    train.add_argument("--report-dir", type=Path, required=True)
    analyze = commands.add_parser(
        "analyze", help="Run integrated genome detection and phenotype prediction."
    )
    analyze.add_argument("--fasta", type=Path, required=True)
    analyze.add_argument("--sample-id", required=True)
    analyze.add_argument("--output-dir", type=Path, required=True)
    analyze.add_argument("--config", type=Path, default=Path("configs/phase5_engine.json"))
    analyze.add_argument("--detection-cache", type=Path, default=None)
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
    records = read_ast_file(input_path, species=config.species)
    manifest: dict[str, object] = {
        "source": "User-provided AST Browser-compatible file",
        "source_file": input_path.name,
        "source_file_sha256": _sha256(input_path),
        "species_filter": config.species,
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


def construct_cohort(args: argparse.Namespace) -> int:
    config = Phase2Config.load(args.config)
    input_path = args.input.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"AST input file not found: {input_path}")
    records = read_ast_file(input_path, species=config.species)
    normalized = normalize_training_labels(records, config.antibiotics)
    biosamples = sorted({record.biosample_accession for record in normalized.records})
    if args.max_biosamples:
        biosamples = biosamples[: args.max_biosamples]
        selected = set(biosamples)
        normalized = NormalizedLabels(
            records=tuple(
                record for record in normalized.records if record.biosample_accession in selected
            ),
            duplicate_rows_removed=normalized.duplicate_rows_removed,
            conflicting_pairs_removed=normalized.conflicting_pairs_removed,
            rows_without_biosample_removed=normalized.rows_without_biosample_removed,
        )
    print(f"Resolving assemblies for {len(biosamples):,} BioSamples.")
    client = NCBIDatasetsClient(api_key=args.api_key)

    def progress(index: int, total: int, resolved: int) -> None:
        print(f"Assembly batch {index}/{total}: {resolved:,} cache entries available.")

    assemblies = resolve_assemblies(
        biosamples,
        client,
        args.assembly_cache.resolve(),
        batch_size=config.assembly_request_batch_size,
        progress=progress,
    )
    bundle = build_cohort(normalized, assemblies, config)
    manifest: dict[str, object] = {
        "source": "User-provided NCBI AST Browser export joined to NCBI Datasets assemblies",
        "source_file": input_path.name,
        "source_file_sha256": _sha256(input_path),
        "assembly_cache_file": args.assembly_cache.name,
        "assembly_cache_sha256": _sha256(args.assembly_cache.resolve()),
        "ncbi_datasets_api": DATASETS_API_BASE,
        "source_urls": [DATASETS_DOCS_URL, DATASETS_API_BASE],
        "retrieved_at_utc": datetime.now(UTC).isoformat(),
        "config": config.as_dict(),
        "config_sha256": config.sha256(),
        "package_version": __version__,
        "requested_biosamples": len(biosamples),
        "cohort_summary": bundle.summary,
    }
    write_cohort_outputs(args.output_dir, bundle, manifest)
    print(json.dumps(bundle.summary, indent=2))
    print(f"Cohort written to {args.output_dir.resolve()}")
    return 0


def download_genomes(args: argparse.Namespace) -> int:
    if not args.accessions.is_file():
        raise FileNotFoundError(f"Assembly accession file not found: {args.accessions}")
    if not 1 <= args.max_workers <= 30:
        raise ValueError("max-workers must be between 1 and 30")
    accessions = args.accessions.read_text(encoding="utf-8").splitlines()
    version = datasets_version(args.datasets_executable)
    print(f"Using {version}; {len(accessions):,} accessions listed.")

    def progress(index: int, total: int, status: str) -> None:
        print(f"Genome shard {index}/{total}: {status}.")

    rows = download_genome_shards(
        accessions,
        args.output_dir,
        args.datasets_executable,
        shard_size=args.shard_size,
        max_shards=args.max_shards,
        max_workers=args.max_workers,
        prefer_ftp=args.prefer_ftp,
        progress=progress,
    )
    verified = sum(row["download_status"] == "verified" for row in rows)
    download_manifest = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "package_version": __version__,
        "ncbi_datasets_version": version,
        "accession_file": args.accessions.name,
        "accession_file_sha256": _sha256(args.accessions),
        "accessions_listed": len(accessions),
        "genomes_reported": len(rows),
        "genomes_verified": verified,
        "shard_size": args.shard_size,
        "max_workers": args.max_workers,
        "sequence_transport": "ncbi-https-ftp-with-datasets-fallback"
        if args.prefer_ftp
        else "ncbi-datasets-rehydrate",
    }
    (args.output_dir / "download_manifest.json").write_text(
        json.dumps(download_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"genomes_listed": len(rows), "genomes_verified": verified}, indent=2))
    return 0 if verified == len(rows) else 2


def finalize_cohort(args: argparse.Namespace) -> int:
    config = Phase2Config.load(args.config)
    for path in (args.cohort, args.assemblies, args.genome_manifest):
        if not path.is_file():
            raise FileNotFoundError(f"Required cohort input not found: {path}")
    summary = finalize_sequence_cohort(
        args.cohort,
        args.assemblies,
        args.genome_manifest,
        args.output_dir,
        genome_dir=args.genome_dir,
        qc=config.assembly_qc,
        max_workers=args.max_workers,
        progress=lambda index, total: print(f"Sequence verification {index:,}/{total:,}."),
    )
    final_manifest = {
        "python_version": sys.version,
        "source_code_sha256": {
            path.relative_to(Path(__file__).parent).as_posix(): _sha256(path)
            for path in sorted(Path(__file__).parent.rglob("*.py"))
        },
        "config": config.as_dict(),
        "config_sha256": config.sha256(),
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "package_version": __version__,
        "inputs": {
            "cohort": {"name": args.cohort.name, "sha256": _sha256(args.cohort)},
            "assemblies": {
                "name": args.assemblies.name,
                "sha256": _sha256(args.assemblies),
            },
            "genome_manifest": {
                "name": args.genome_manifest.name,
                "sha256": _sha256(args.genome_manifest),
            },
        },
        "summary": summary,
        "outputs": {
            path.name: {"sha256": _sha256(path)}
            for path in sorted(args.output_dir.iterdir())
            if path.is_file() and path.name not in {"final_manifest.json", "verification.json"}
        },
    }
    for key, path in (
        ("cohort_provenance", args.cohort.parent / "data_manifest.json"),
        ("download_provenance", args.genome_dir / "download_manifest.json"),
    ):
        if path.is_file():
            final_manifest[key] = json.loads(path.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "final_manifest.json").write_text(
        json.dumps(final_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(f"Final cohort written to {args.output_dir.resolve()}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "analyze":
            from genome_amr_predictor.engine import analyze_genome

            result = analyze_genome(
                args.fasta,
                args.output_dir,
                args.sample_id,
                args.config,
                detection_cache=args.detection_cache,
            )
            print(
                json.dumps(
                    {
                        "status": result["status"],
                        "output_dir": str(args.output_dir),
                        "errors": result["errors"],
                    }
                )
            )
            return (
                0 if result["status"] == "completed" else 1 if result["status"] == "partial" else 2
            )
        if args.command == "build-features":
            from genome_amr_predictor.features import build_features

            summary = build_features(args.assemblies, args.detection_dir, args.output_dir)
            print(
                json.dumps({key: summary[key] for key in ("status", "genomes", "schema_version")})
            )
            return 0
        if args.command == "train-models":
            from genome_amr_predictor.modeling import train_models

            summary = train_models(
                args.features,
                args.cohort,
                args.similarity_groups,
                args.config,
                args.model_dir,
                args.report_dir,
            )
            print(json.dumps({"status": summary["status"], "models": list(summary["models"])}))
            return 0
        if args.command == "detect":
            runtime = AMRFinderRuntime(args.executable, args.database, args.wsl_distribution)
            manifest = detect_genome(
                args.fasta,
                args.output_dir,
                args.sample_id,
                runtime,
                threads=args.threads,
                timeout=args.timeout,
            )
            print(json.dumps({"status": manifest["status"], "total_hits": manifest["total_hits"]}))
            return 0
        if args.command == "parse-detection":
            hits = parse_amrfinder(args.input, args.sample_id)
            if args.output.exists():
                raise ValueError("Output already exists; choose a new evidence file")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            write_evidence(
                args.output,
                hits,
                {
                    "source": "imported_tsv",
                    "input_sha256": sha256(args.input),
                    "software_version": args.software_version,
                    "database_version": args.database_version,
                    "version_source": "user_supplied_not_runtime_verified",
                },
            )
            print(json.dumps({"status": "completed", "total_hits": len(hits)}))
            return 0
        if args.command == "audit-ncbi":
            return audit_ncbi(args)
        if args.command == "audit-file":
            return audit_file(args)
        if args.command == "construct-cohort":
            return construct_cohort(args)
        if args.command == "download-genomes":
            return download_genomes(args)
        if args.command == "finalize-cohort":
            return finalize_cohort(args)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))
    parser.error(f"Unknown command: {args.command}")
    return 2
