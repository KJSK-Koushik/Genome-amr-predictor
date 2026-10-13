"""Resume real cohort detection; failed or missing calls are never treated as absent genes."""

from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path

from genome_amr_predictor.detection import AMRFinderRuntime, detect_genome, sha256


def run_one(row, args, runtime):
    sample = row["biosample_accession"]
    if not sample.isalnum() or not sample.startswith(("SAMN", "SAME", "SAMD")):
        raise ValueError("Invalid BioSample path identifier")
    root = args.genome_dir.resolve()
    fasta = (root / row["fasta_path"]).resolve()
    if not fasta.is_relative_to(root) or sha256(fasta) != row["fasta_file_sha256"]:
        raise ValueError(f"Frozen FASTA changed: {row['biosample_accession']}")
    output = args.output_dir / row["biosample_accession"]
    saved = output / "run_manifest.json"
    if saved.exists():
        manifest = json.loads(saved.read_text())
        if manifest["status"] == "completed":
            if (
                manifest["input"]["file_sha256"] != row["fasta_file_sha256"]
                or manifest["input"]["sequence_sha256"] != row["sequence_sha256"]
                or manifest["executable_sha256"] != sha256(runtime.executable)
                or manifest["database_files_sha256"] != args.database_hashes
                or manifest["raw_output_sha256"] != sha256(output / "amrfinder.tsv")
                or manifest["evidence_sha256"] != sha256(output / "detections.json")
            ):
                raise ValueError("Cached detector provenance mismatch")
            return "cached"
        # Keep failed runs for diagnosis; retry into a fresh directory.
        backup = output.with_name(
            output.name + ".failed-" + datetime.now(UTC).strftime("%Y%m%d%H%M%S%f")
        )
        output.rename(backup)
    elif output.exists() and any(output.iterdir()):
        if not output.resolve().is_relative_to(args.output_dir.resolve()):
            raise ValueError("Unsafe interrupted-run path")
        backup = output.with_name(
            output.name + ".interrupted-" + datetime.now(UTC).strftime("%Y%m%d%H%M%S%f")
        )
        output.rename(backup)
    manifest = detect_genome(
        fasta,
        output,
        row["biosample_accession"],
        runtime,
        threads=args.threads,
        timeout=args.timeout,
    )
    if manifest["input"]["sequence_sha256"] != row["sequence_sha256"]:
        raise ValueError("Frozen sequence fingerprint mismatch")
    return "completed"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--genome-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 16 or args.limit < 0:
        parser.error("workers must be 1..16 and limit must be nonnegative")
    runtime = AMRFinderRuntime(args.executable, args.database)
    args.database_hashes = {
        p.relative_to(runtime.database).as_posix(): sha256(p)
        for p in sorted(runtime.database.rglob("*"))
        if p.is_file()
    }
    with args.assemblies.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if args.limit:
        rows = rows[: args.limit]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    progress = {
        "status": "running",
        "requested": len(rows),
        "completed": 0,
        "cached": 0,
        "failures": [],
        "assembly_table_sha256": sha256(args.assemblies),
        "started_at_utc": datetime.now(UTC).isoformat(),
        "workers": args.workers,
        "threads_per_genome": args.threads,
    }
    path = args.output_dir / "batch_progress.json"
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = {pool.submit(run_one, row, args, runtime): row for row in rows}
        for future in as_completed(jobs):
            row = jobs[future]
            try:
                progress[future.result()] += 1
            except Exception as exc:
                progress["failures"].append(
                    {"biosample_accession": row["biosample_accession"], "error": str(exc)}
                )
            progress["updated_at_utc"] = datetime.now(UTC).isoformat()
            temporary = path.with_suffix(".tmp")
            temporary.write_text(json.dumps(progress, indent=2) + "\n")
            temporary.replace(path)
            done = progress["completed"] + progress["cached"]
            print(f"Verified {done}/{len(rows)}; failures {len(progress['failures'])}", flush=True)
    progress["status"] = "failed" if progress["failures"] else "completed"
    path.write_text(json.dumps(progress, indent=2) + "\n")
    return 1 if progress["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
