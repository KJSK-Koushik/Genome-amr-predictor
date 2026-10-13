"""Re-run every delivery check and write one summary; never refits models or edits evidence."""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import platform
import re
import subprocess
import sys
import time
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "src")]


def _tests() -> dict:
    command = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"]
    done = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
    tail = done.stdout.strip().splitlines()[-1] if done.stdout.strip() else done.stderr[-500:]
    if done.returncode:
        raise RuntimeError(f"pytest failed: {tail}")
    return {"summary": tail, "passed": int(re.search(r"(\d+) passed", tail).group(1))}


def _lint() -> dict:
    done = subprocess.run(
        [sys.executable, "-m", "ruff", "check", "."],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if done.returncode:
        raise RuntimeError(done.stdout[-1000:])
    return {"summary": done.stdout.strip()}


def _phase2() -> dict:
    from verify_phase2 import verify

    return verify(ROOT / "reports/phase2/2026-08-18")


def _phase3() -> dict:
    from verify_phase3 import verify

    return verify(ROOT / "reports/phase3/2026-10-13/smoke-linux-runtime")


def _pretrained() -> dict:
    # The verifier uses project-relative paths and rewrites its deterministic result file.
    import os

    import verify_pretrained_hzi

    previous = Path.cwd()
    os.chdir(ROOT)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            verify_pretrained_hzi.main()
    finally:
        os.chdir(previous)
    return json.loads((ROOT / "reports/pretrained/2026-10-13/verification.json").read_text())


def _phase5() -> dict:
    from verify_phase5 import verify

    return verify(ROOT / "reports/phase5/2026-10-13")


def _saved_runs() -> dict:
    from genome_amr_predictor.reporting import find_runs, load_run, render_pdf

    runs = find_runs(ROOT / "reports")
    rows = []
    for directory in runs:
        run = load_run(directory)
        if not render_pdf(run).startswith(b"%PDF"):
            raise ValueError(f"PDF rendering failed for {directory}")
        rows.append(
            {
                "run": directory.relative_to(ROOT).as_posix(),
                "sample_id": run["sample_id"],
                "status": run["status"],
                "detection": run["detection"]["status"],
                "prediction": run["phenotype_prediction"]["status"],
            }
        )
    if not rows:
        raise ValueError("No committed analysis runs found")
    return {"runs_verified": len(rows), "runs": rows}


CHECKS = {
    "lint": _lint,
    "tests": _tests,
    "phase2_cohort": _phase2,
    "phase3_detection": _phase3,
    "pretrained_models": _pretrained,
    "phase5_engine": _phase5,
    "saved_runs_and_pdf": _saved_runs,
}
NEEDS_LOCAL_DATA = {"phase3_detection", "pretrained_models", "phase5_engine"}


def _environment() -> dict:
    def version(name: str) -> str | None:
        try:
            return metadata.version(name)
        except metadata.PackageNotFoundError:
            return None

    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=False
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain", "--", "src", "tests", "scripts", "configs"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()
    config = json.loads((ROOT / "configs/phase5_engine.json").read_text())
    return {
        "git_commit": commit or None,
        "uncommitted_code_changes": bool(dirty),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": {
            name: version(name)
            for name in ("genome-amr-predictor", "scikit-learn", "numpy", "streamlit", "fpdf2")
        },
        "amrfinderplus": config["detector"]["software_version"],
        "amrfinder_database": config["detector"]["database_version"],
    }


def run(skip: set[str]) -> dict:
    results = {}
    for name, check in CHECKS.items():
        if name in skip:
            results[name] = {"status": "skipped"}
            continue
        start = time.perf_counter()
        try:
            detail = check()
            results[name] = {"status": "passed", "detail": detail}
        except Exception as exc:
            results[name] = {"status": "failed", "error": f"{type(exc).__name__}: {exc}"}
        results[name]["seconds"] = round(time.perf_counter() - start, 1)
    statuses = {value["status"] for value in results.values()}
    return {
        "schema": "genome-amr-delivery-verification-v1",
        "generated_at_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "status": "failed" if "failed" in statuses else "passed",
        "environment": _environment(),
        "checks": results,
        "not_covered": [
            "Own Phase 4 cohort models: full-cohort training was suspended; no metrics exist.",
            "Population-level detection or phenotype accuracy: integration checks only.",
            "Browser file selection in the dashboard upload widget (manual step).",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--portable",
        action="store_true",
        help="Skip checks needing the local AMRFinderPlus install, genomes or model downloads.",
    )
    args = parser.parse_args()
    result = run(NEEDS_LOCAL_DATA if args.portable else set())
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text)
    sys.exit(0 if result["status"] == "passed" else 1)
