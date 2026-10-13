"""Streamlit dashboard for integrated genome AMR analysis runs.

Launch with ``./scripts/run_dashboard.ps1`` or
``streamlit run src/genome_amr_predictor/dashboard.py``.
"""

from __future__ import annotations

import json
import os
import re
import uuid
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUN_ROOTS = (
    ROOT / "data" / "processed" / "phase6" / "runs",
    ROOT / "data" / "processed" / "phase5",
    ROOT / "reports" / "phase5",
)
RUN_ROOTS_ENV = "GENOME_AMR_RUN_ROOTS"
UPLOAD_ROOT = ROOT / "data" / "processed" / "phase6" / "uploads"
CONFIG = ROOT / "configs" / "phase5_engine.json"
FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fas")
RESISTANT, SUSCEPTIBLE = "#eb6834", "#2a78d6"
NOT_CLINICAL = (
    "Research and educational use only. These results are not a clinical diagnosis and do not "
    "recommend treatment; laboratory confirmation is required."
)


def run_roots() -> tuple[Path, ...]:
    """Run folders to list; the first receives new analyses."""
    configured = os.environ.get(RUN_ROOTS_ENV)
    if configured:
        return tuple(Path(part) for part in configured.split(os.pathsep) if part)
    return DEFAULT_RUN_ROOTS


def upload_path(filename: str, base: Path) -> Path:
    """A fresh location for an uploaded FASTA; the user's file name is never used as a path."""
    name = filename.lower()
    compressed = name.endswith(".gz")
    stem = name[:-3] if compressed else name
    if not stem.endswith(FASTA_SUFFIXES):
        raise ValueError("Upload an assembled FASTA file (.fna, .fa, .fasta, optionally .gz)")
    return base / uuid.uuid4().hex / ("input.fna.gz" if compressed else "input.fna")


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "sample"


def _label(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def run_title(path: Path) -> str:
    """Sample, status and finish date of a run, for the run picker."""
    try:
        run = json.loads((path / "run_manifest.json").read_text(encoding="utf-8"))
        finished = run.get("finished_at_utc", "")[:16].replace("T", " ")
        reused = " · reused detection" if run["detection"].get("reused") else ""
        return f"{run['sample_id']} · {run['status']} · {finished} UTC{reused}"
    except (OSError, ValueError, KeyError):
        return _label(path)


def _analyze_new(st) -> None:
    from genome_amr_predictor.engine import analyze_genome

    with st.sidebar.form("analyze"):
        upload = st.file_uploader(
            "Assembled E. coli genome (FASTA)", type=["fna", "fa", "fasta", "fas", "gz"]
        )
        sample_id = st.text_input("Sample ID", placeholder="e.g. SAMN03075588")
        st.caption("Species is assumed to be Escherichia coli; it is not verified automatically.")
        submitted = st.form_submit_button("Run analysis", type="primary")
    if not submitted:
        return
    if upload is None or not sample_id.strip():
        st.sidebar.error("Provide both a FASTA file and a sample ID.")
        return
    try:
        target = upload_path(upload.name, UPLOAD_ROOT)
    except ValueError as exc:
        st.sidebar.error(str(exc))
        return
    target.parent.mkdir(parents=True)
    target.write_bytes(upload.getvalue())
    output = run_roots()[0] / target.parent.name
    with st.spinner("Running quality checks, pretrained models and AMRFinderPlus (1-2 minutes)"):
        try:
            analyze_genome(target, output, sample_id.strip(), CONFIG)
        except Exception as exc:
            st.sidebar.error(f"Analysis could not start: {exc}")
            return
    st.session_state["run_dir"] = str(output)


def _choose_saved(st, find_runs) -> None:
    runs = find_runs(*run_roots())
    if not runs:
        st.sidebar.info("No saved runs yet. Analyze a genome to create one.")
        return
    keys = [_label(path) for path in runs]
    current = st.session_state.get("run_dir")
    current_key = _label(Path(current)) if current else None
    choice = st.sidebar.selectbox(
        "Saved run",
        range(len(runs)),
        index=keys.index(current_key) if current_key in keys else 0,
        format_func=lambda i: run_title(runs[i]),
    )
    st.sidebar.caption(keys[choice])
    st.session_state["run_dir"] = str(runs[choice])


def _status(st, run: dict) -> None:
    messages = {
        "completed": (st.success, "Both detection and phenotype prediction completed."),
        "partial": (st.warning, "One stage failed; see the affected tab and Limitations."),
        "failed": (st.error, "The analysis failed; no scientific results are available."),
    }
    show, text = messages[run["status"]]
    show(f"Run status: **{run['status']}**. {text}")


def _metrics(st, run: dict) -> None:
    metrics = run["input"]
    if not metrics:
        return
    tiles = [
        ("Contigs", f"{metrics['contig_count']:,}"),
        ("Length", f"{metrics['total_length'] / 1e6:.1f} Mb"),
        ("N50", f"{metrics['contig_n50'] / 1e3:.0f} kb"),
        ("GC", f"{metrics['gc_percent']:.1f}%"),
        ("Ambiguous", f"{metrics['ambiguous_percent']:.2f}%"),
    ]
    for column, (label, value) in zip(st.columns(len(tiles)), tiles, strict=True):
        column.metric(label, value)


def _predictions(st, run: dict, reporting) -> None:
    import altair as alt
    import pandas as pd

    rows = reporting.prediction_rows(run)
    if not rows:
        stage = run["phenotype_prediction"]
        st.error(
            f"Phenotype prediction {stage['status']}: "
            + stage.get("error", "no predictions are available.")
        )
        return
    scored = [row for row in rows if row["resistance_score"] is not None]
    if scored:
        frame = pd.DataFrame(scored)
        frame["label"] = [f"{row['prediction']} {row['resistance_score']:.2f}" for row in scored]
        base = alt.Chart(frame).encode(
            y=alt.Y("antibiotic:N", sort=None, title=None),
            x=alt.X(
                "resistance_score:Q",
                scale=alt.Scale(domain=[0, 1]),
                title="Resistance score (not a calibrated probability)",
            ),
        )
        bars = base.mark_bar(cornerRadiusEnd=4, size=16).encode(
            color=alt.Color(
                "prediction:N",
                scale=alt.Scale(domain=["R", "S"], range=[RESISTANT, SUSCEPTIBLE]),
                legend=alt.Legend(title="Predicted", orient="top"),
            ),
            tooltip=[
                alt.Tooltip("antibiotic:N", title="Antibiotic"),
                alt.Tooltip("prediction:N", title="Prediction"),
                alt.Tooltip("resistance_score:Q", title="Score", format=".3f"),
                alt.Tooltip("features_present:Q", title="Features present"),
                alt.Tooltip("note:N", title="Note"),
            ],
        )
        labels = base.mark_text(align="left", dx=6).encode(text="label:N")
        threshold = (
            alt.Chart(pd.DataFrame({"x": [reporting.DECISION_THRESHOLD]}))
            .mark_rule(strokeDash=[4, 3], color="gray")
            .encode(x="x:Q")
        )
        st.altair_chart(
            (bars + labels + threshold).properties(height=alt.Step(36)),
            width="stretch",
        )
        st.caption(
            f"Dashed line: {reporting.DECISION_THRESHOLD} decision threshold. Scores come from "
            "externally trained models that have not been validated on this project's cohort."
        )
    for row in rows:
        if row["note"] and row["status"] == "predicted":
            st.warning(f"**{row['antibiotic']}**: {row['note']}")
    table = [
        {
            "Antibiotic": row["antibiotic"],
            "Prediction": row["prediction"] or "—",
            "Score": "—" if row["resistance_score"] is None else f"{row['resistance_score']:.4f}",
            "Status": row["status"],
            "Model features present": "—"
            if row["features_total"] is None
            else f"{row['features_present']} / {row['features_total']}",
            "Note": row["note"],
        }
        for row in rows
    ]
    st.dataframe(pd.DataFrame(table), hide_index=True, width="stretch")


def _determinants(st, run: dict, reporting) -> None:
    import altair as alt
    import pandas as pd

    summary = reporting.detection_summary(run)
    evidence = run["detection"]["evidence"]
    if evidence is None:
        st.error(summary)
        st.caption("An unavailable detector result is never interpreted as absent genes.")
        return
    st.markdown(summary)
    if not evidence["hits"]:
        return
    classes = reporting.determinant_classes(run)
    frame = pd.DataFrame({"drug_class": list(classes), "elements": list(classes.values())})
    base = alt.Chart(frame).encode(
        y=alt.Y("drug_class:N", sort=None, title=None, axis=alt.Axis(labelLimit=240)),
        x=alt.X("elements:Q", title="Detected elements", axis=alt.Axis(tickMinStep=1)),
    )
    chart = base.mark_bar(cornerRadiusEnd=4, size=16, color=SUSCEPTIBLE).encode(
        tooltip=[
            alt.Tooltip("drug_class:N", title="Drug class"),
            alt.Tooltip("elements:Q", title="Elements"),
        ]
    ) + base.mark_text(align="left", dx=6).encode(text="elements:Q")
    st.altair_chart(chart.properties(height=alt.Step(30)), width="stretch")
    st.caption(
        "Curated NCBI drug classes of detected genes, alleles and point mutations. A class "
        "listed here is evidence of a known determinant, not a predicted phenotype."
    )
    kinds = sorted({hit["evidence_kind"] for hit in evidence["hits"]})
    chosen = st.multiselect("Evidence type", kinds, default=kinds)
    hits = [hit for hit in evidence["hits"] if hit["evidence_kind"] in chosen]
    columns = [
        "element_symbol",
        "element_name",
        "element_type",
        "evidence_kind",
        "drug_class",
        "drug_subclass",
        "method",
        "identity_percent",
        "coverage_percent",
        "contig_id",
        "start",
        "stop",
        "strand",
    ]
    st.dataframe(
        pd.DataFrame(hits, columns=columns),
        hide_index=True,
        width="stretch",
        column_config={
            "element_symbol": "Element",
            "element_name": "Name",
            "element_type": "Type",
            "evidence_kind": "Evidence",
            "drug_class": "Class",
            "drug_subclass": "Subclass",
            "method": "Method",
            "identity_percent": st.column_config.NumberColumn("Identity %", format="%.2f"),
            "coverage_percent": st.column_config.NumberColumn("Coverage %", format="%.2f"),
            "contig_id": "Contig",
            "start": "Start",
            "stop": "Stop",
            "strand": "Strand",
        },
    )


def _limitations(st, run: dict) -> None:
    for warning in run["warnings"]:
        st.markdown(f"- {warning}")
    if run["errors"]:
        st.subheader("Errors")
        for error in run["errors"]:
            st.error(f"{error['stage']}: {error['message']}")


def _provenance(st, run: dict) -> None:
    import pandas as pd

    evidence = run["detection"]["evidence"] or {}
    manifest = evidence.get("verified_manifest", {})
    prediction = run["phenotype_prediction"]["result"] or {}
    rows = [
        ("Run directory", run["run_directory"]),
        ("Started (UTC)", run["started_at_utc"]),
        ("Finished (UTC)", run.get("finished_at_utc")),
        ("Package version", run["package_version"]),
        ("Engine SHA-256", run.get("engine_sha256")),
        ("Configuration SHA-256", run.get("config_sha256")),
        ("Input file SHA-256", (run["input"] or {}).get("file_sha256")),
        ("AMRFinderPlus version", manifest.get("software_version")),
        ("Detection reused", run["detection"].get("reused")),
        ("Model registry SHA-256", prediction.get("registry_sha256")),
    ]
    rows += [(f"{name} SHA-256", digest) for name, digest in run["artifacts_sha256"].items()]
    st.dataframe(
        pd.DataFrame(
            [(name, "-" if value is None else str(value)) for name, value in rows],
            columns=["Item", "Value"],
        ),
        hide_index=True,
        width="stretch",
    )


def _downloads(st, run: dict, directory: Path, reporting) -> None:
    prefix = safe_name(run["sample_id"])
    if run["status"] != "completed":
        st.warning(
            "Empty CSV files from failed stages do not indicate absent genes or susceptibility."
        )
    st.download_button(
        "PDF report",
        reporting.render_pdf(run),
        f"{prefix}_report.pdf",
        "application/pdf",
        type="primary",
    )
    files = (
        ("Full results (JSON)", "analysis.json", "application/json"),
        ("Predictions (CSV)", "predictions.csv", "text/csv"),
        ("Detected determinants (CSV)", "detections.csv", "text/csv"),
        ("Summary (Markdown)", "report.md", "text/markdown"),
    )
    for label, name, mime in files:
        st.download_button(label, (directory / name).read_bytes(), f"{prefix}_{name}", mime)
    st.caption("Exported files are byte-identical to the verified run artifacts.")


def main() -> None:
    import streamlit as st

    from genome_amr_predictor import reporting

    st.set_page_config(page_title="Genome AMR Analyzer", layout="wide")
    st.sidebar.title("Genome AMR Analyzer")
    mode = st.sidebar.radio("Source", ["Open a saved run", "Analyze a new genome"])
    if mode == "Analyze a new genome":
        _analyze_new(st)
    else:
        _choose_saved(st, reporting.find_runs)
    st.sidebar.caption(NOT_CLINICAL)

    if "run_dir" not in st.session_state:
        st.title("Genome AMR Analyzer")
        st.info(NOT_CLINICAL)
        st.write("Open a saved run or analyze a new assembled E. coli genome from the sidebar.")
        return
    directory = Path(st.session_state["run_dir"])
    try:
        run = reporting.load_run(directory)
    except (OSError, ValueError, KeyError) as exc:
        st.error(f"This run cannot be displayed because it failed an integrity check: {exc}")
        return

    st.title(f"Sample {run['sample_id']}")
    st.info(NOT_CLINICAL)
    _status(st, run)
    _metrics(st, run)
    tabs = st.tabs(
        ["Phenotype predictions", "Known determinants", "Limitations", "Provenance", "Downloads"]
    )
    with tabs[0]:
        _predictions(st, run, reporting)
    with tabs[1]:
        _determinants(st, run, reporting)
    with tabs[2]:
        _limitations(st, run)
    with tabs[3]:
        _provenance(st, run)
    with tabs[4]:
        _downloads(st, run, directory, reporting)


if __name__ == "__main__":
    main()
