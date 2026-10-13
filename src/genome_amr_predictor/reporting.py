"""Load verified analysis runs and render them as display tables and PDF reports."""

from __future__ import annotations

import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from genome_amr_predictor.detection import sha256

RUN_SCHEMA = "genome-amr-analysis-v1"
ARTIFACTS = ("analysis.json", "predictions.csv", "detections.csv", "report.md")
DECISION_THRESHOLD = 0.5
BORDERLINE_MARGIN = 0.05


def load_run(directory: Path) -> dict:
    """Return a finished run only if every exported artifact matches its recorded hash."""
    manifest = json.loads((directory / "run_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema_version") != RUN_SCHEMA:
        raise ValueError("Not an integrated analysis run")
    if manifest.get("stage") != "finished" or "artifacts_sha256" not in manifest:
        raise ValueError("Analysis run has not finished")
    for name in ARTIFACTS:
        if sha256(directory / name) != manifest["artifacts_sha256"][name]:
            raise ValueError(f"{name} changed after the run finished")
    manifest["run_directory"] = str(directory.resolve())
    return manifest


def find_runs(*roots: Path) -> list[Path]:
    """Finished analysis run directories under the roots, newest first."""
    runs = set()
    for root in roots:
        if not root.is_dir():
            continue
        for path in root.rglob("run_manifest.json"):
            try:
                value = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if value.get("schema_version") == RUN_SCHEMA and value.get("stage") == "finished":
                runs.add(path.parent)
    return sorted(runs, key=lambda path: path.stat().st_mtime, reverse=True)


def detection_summary(run: dict) -> str:
    """One sentence that never turns an unavailable detector into gene absence."""
    detection = run["detection"]
    if detection["status"] == "completed":
        total = detection["evidence"]["summary"]["total_hits"]
        if total == 0:
            return "No known determinants detected. This does not establish susceptibility."
        source = " (verified reuse of saved detection)" if detection.get("reused") else ""
        return f"{total} known determinant elements detected{source}."
    if detection["status"] == "failed":
        return f"Detection failed: gene presence/absence is unknown. {detection['error']}"
    return "Detection was not available for this run: gene presence/absence is unknown."


def prediction_rows(run: dict) -> list[dict]:
    """Predictions with an explanatory note for anything other than a scored call."""
    result = run["phenotype_prediction"]["result"]
    if result is None:
        return []
    rows = []
    for row in result["predictions"]:
        if row["status"] == "unsupported":
            note = "No validated model for this antibiotic."
        elif row["status"] == "withheld_qc":
            note = "Withheld: assembly failed prediction quality checks."
        elif row["features_present"] == 0:
            note = "No model features present; score is the model baseline."
        elif abs(row["resistance_score"] - DECISION_THRESHOLD) < BORDERLINE_MARGIN:
            note = "Close to the decision threshold; treat the call as uncertain."
        else:
            note = ""
        rows.append({**row, "note": note})
    return rows


def determinant_classes(run: dict) -> dict[str, int]:
    """Count detected elements per curated drug class, most frequent first."""
    evidence = run["detection"]["evidence"]
    if evidence is None:
        return {}
    counts = Counter(hit["drug_class"] or "Unclassified" for hit in evidence["hits"])
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))


def _latin(text: object) -> str:
    # Core PDF fonts are Latin-1 only: dashes, comparison signs and curly quotes.
    replacements = dict(zip("—–≥≤’“”", ["-", "-", ">=", "<=", "'", '"', '"'], strict=True))  # noqa: RUF001
    value = str(text)
    for old, new in replacements.items():
        value = value.replace(old, new)
    return value.encode("latin-1", "replace").decode("latin-1")


def _database_version(report: str | None) -> str | None:
    for line in (report or "").splitlines():
        if line.startswith("Database version:"):
            return line.split(":", 1)[1].strip()
    return None


def render_pdf(run: dict, *, compress: bool = True) -> bytes:
    """A printable report of one run; scientific outputs and limitations stay separate."""
    from fpdf import FPDF
    from fpdf.fonts import FontFace

    pdf = FPDF(format="A4")
    pdf.set_compression(compress)
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.set_title(_latin(f"Genome AMR analysis: {run['sample_id']}"))
    pdf.add_page()
    width = pdf.epw
    header = FontFace(emphasis="BOLD", fill_color=(235, 235, 232))

    def heading(text: str) -> None:
        pdf.ln(3)
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(width, 7, _latin(text), new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", size=9)

    def paragraph(text: str) -> None:
        pdf.multi_cell(width, 5, _latin(text), new_x="LMARGIN", new_y="NEXT")

    def table(rows: list[list[object]], widths: tuple[int, ...]) -> None:
        pdf.set_font("Helvetica", size=8)
        pdf.set_fill_color(255, 255, 255)
        with pdf.table(col_widths=widths, headings_style=header, line_height=5) as grid:
            for values in rows:
                line = grid.row()
                for value in values:
                    line.cell(_latin("-" if value is None or value == "" else value))
        pdf.set_font("Helvetica", size=9)

    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(
        width, 9, _latin(f"Genome AMR analysis: {run['sample_id']}"), new_x="LMARGIN", new_y="NEXT"
    )
    pdf.set_font("Helvetica", size=9)
    paragraph(
        f"Run status: {run['status']}    Finished: {run.get('finished_at_utc', '-')}    "
        f"Report generated: {datetime.now(UTC).isoformat(timespec='seconds')}"
    )
    pdf.set_fill_color(253, 243, 224)
    pdf.set_font("Helvetica", "B", 9)
    pdf.multi_cell(
        width,
        5,
        "Research and educational use only. Not a clinical result; laboratory confirmation "
        "is required. Detected determinants and predicted phenotypes are separate outputs.",
        fill=True,
        new_x="LMARGIN",
        new_y="NEXT",
    )
    pdf.set_font("Helvetica", size=9)

    heading("Input genome")
    metrics = run["input"]
    if metrics:
        table(
            [
                ["Species", "Contigs", "Total length", "N50", "GC %", "Ambiguous %"],
                [
                    metrics["species"] + " (user-supplied)",
                    metrics["contig_count"],
                    f"{metrics['total_length']:,}",
                    f"{metrics['contig_n50']:,}",
                    f"{metrics['gc_percent']:.2f}",
                    f"{metrics['ambiguous_percent']:.3f}",
                ],
            ],
            (40, 15, 25, 20, 15, 20),
        )
        paragraph(f"File SHA-256: {metrics['file_sha256']}")
    else:
        paragraph("Input was rejected before analysis; no scientific outputs are available.")

    heading("Predicted phenotypes (pretrained models)")
    rows = prediction_rows(run)
    if rows:
        scored = [row for row in rows if row["resistance_score"] is not None]
        if scored:
            label_width, bar_width = 50.0, width - 70
            for row in scored:
                y = pdf.get_y()
                pdf.cell(label_width, 5, _latin(row["antibiotic"]))
                pdf.set_fill_color(240, 239, 236)
                pdf.rect(pdf.l_margin + label_width, y + 1, bar_width, 3, style="F")
                color = (235, 104, 52) if row["prediction"] == "R" else (42, 120, 214)
                pdf.set_fill_color(*color)
                pdf.rect(
                    pdf.l_margin + label_width,
                    y + 1,
                    bar_width * row["resistance_score"],
                    3,
                    style="F",
                )
                pdf.set_draw_color(80, 80, 80)
                threshold_x = pdf.l_margin + label_width + bar_width * DECISION_THRESHOLD
                pdf.line(threshold_x, y, threshold_x, y + 5)
                pdf.set_x(pdf.l_margin + label_width + bar_width + 2)
                pdf.cell(
                    18,
                    5,
                    f"{row['prediction']} {row['resistance_score']:.2f}",
                    new_x="LMARGIN",
                    new_y="NEXT",
                )
            paragraph(
                f"Bars show resistance score; the line marks the {DECISION_THRESHOLD} "
                "decision threshold."
            )
        table(
            [["Antibiotic", "Prediction", "Score", "Status", "Features", "Note"]]
            + [
                [
                    row["antibiotic"],
                    row["prediction"],
                    None if row["resistance_score"] is None else f"{row['resistance_score']:.4f}",
                    row["status"],
                    None
                    if row["features_total"] is None
                    else f"{row['features_present']}/{row['features_total']}",
                    row["note"],
                ]
                for row in rows
            ],
            (32, 16, 14, 18, 16, 60),
        )
        paragraph("Scores are external model outputs, not locally calibrated probabilities.")
    else:
        paragraph(
            f"Phenotype prediction {run['phenotype_prediction']['status']}: "
            + run["phenotype_prediction"].get("error", "no predictions are available.")
        )

    heading("Known resistance determinants (AMRFinderPlus)")
    paragraph(detection_summary(run))
    evidence = run["detection"]["evidence"]
    if evidence and evidence["hits"]:
        table(
            [
                [
                    "Element",
                    "Type",
                    "Evidence",
                    "Class",
                    "Method",
                    "Identity %",
                    "Coverage %",
                    "Contig",
                ]
            ]
            + [
                [
                    hit["element_symbol"],
                    hit["element_type"],
                    hit["evidence_kind"],
                    hit["drug_class"],
                    hit["method"],
                    hit["identity_percent"],
                    hit["coverage_percent"],
                    hit["contig_id"],
                ]
                for hit in evidence["hits"]
            ],
            (26, 13, 18, 38, 20, 14, 15, 36),
        )

    heading("Limitations")
    for warning in run["warnings"]:
        paragraph(f"- {warning}")
    if run["errors"]:
        heading("Errors")
        for error in run["errors"]:
            paragraph(f"- {error['stage']}: {error['message']}")

    heading("Provenance")
    manifest = (evidence or {}).get("verified_manifest", {})
    provenance = [
        ("Package version", run.get("package_version")),
        ("Engine SHA-256", run.get("engine_sha256")),
        ("Configuration SHA-256", run.get("config_sha256")),
        ("Detector version", manifest.get("software_version")),
        ("Detector database", _database_version(manifest.get("database_version_report"))),
        (
            "Model registry SHA-256",
            (run["phenotype_prediction"]["result"] or {}).get("registry_sha256"),
        ),
        ("analysis.json SHA-256", run["artifacts_sha256"]["analysis.json"]),
    ]
    for name, value in provenance:
        paragraph(f"{name}: {value if value else '-'}")
    return bytes(pdf.output())
