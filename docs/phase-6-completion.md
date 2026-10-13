# Phase 6 Completion Note

Completed on 2026-10-14 (Asia/Calcutta).

## Outcome

A Streamlit dashboard now sits on top of the Phase 5 run format. It can open any saved run or
analyse a newly uploaded assembled *E. coli* FASTA. It shows input QC, phenotype predictions, and
known determinants on separate tabs, with charts, limitations, provenance, and downloads.

```powershell
python -m pip install -e ".[dev,ml,ui]"
./scripts/run_dashboard.ps1            # http://localhost:8501
```

A PDF can also be produced without the dashboard:

```powershell
$env:PYTHONPATH = "src"
.venv/Scripts/python.exe -m genome_amr_predictor export-pdf --run-dir <run dir> --output report.pdf
```

## What the interface shows

| Area | Content |
| --- | --- |
| Header | Research-use notice; run status (completed / partial / failed); contigs, length, N50, GC, ambiguous bases |
| Phenotype predictions | Resistance-score bar chart with the 0.5 threshold; R/S coloured and labelled; table with feature coverage and notes |
| Known determinants | Elements per curated drug class; filterable table with method, identity, coverage, contig, and coordinates |
| Limitations | All run warnings and stage errors |
| Provenance | Package, engine, configuration, input, detector, registry, and per-artifact SHA-256 values |
| Downloads | PDF report, `analysis.json`, `predictions.csv`, `detections.csv`, `report.md` |

## Safety behaviour

- Before display or export, every run's artifacts are checked against the hashes recorded at run
  end. An altered run is refused with an integrity error rather than shown.
- A failed or rejected detection stage reads "gene presence/absence is unknown". The dashboard
  never shows an empty gene table or "no genes found" in that case. A completed run with zero hits
  states that this does not establish susceptibility.
- Prediction notes flag unsupported antibiotics, QC-withheld predictions, and calls with no model
  features present. They also flag scores within 0.05 of the threshold as uncertain. In the
  `mismatched-cache` run, for example, ceftriaxone scored 0.491.
- An upload is saved under a generated name (`input.fna[.gz]`) in a fresh folder. The browser
  file name is never used as a path, and non-FASTA extensions are rejected.
- Charts use two colours validated for colour-vision deficiency, and the R/S letter is printed
  beside every bar, so colour is never the only cue. Scores are labelled as not calibrated
  probabilities.

## Validation

All **114 package tests passed** and **Ruff passed**. The 23 new tests cover:

- integrity refusal of edited runs and rejection of unfinished or foreign runs;
- run discovery that ignores nested detector manifests;
- failed-detection wording and no-hit wording;
- every prediction note;
- PDF content for completed, failed-detection, and rejected-input runs;
- the `export-pdf` CLI, including refusal for an edited run;
- upload path sanitisation and run-picker labels; and
- headless Streamlit `AppTest` renders of the landing page, a completed run (with all five
  download buttons), a failed-detection run, and a tampered run.

The dashboard was also driven in a browser on the Phase 5 `live` and `mismatched-cache` runs and a
new run; the `reuse` run appeared in the run picker. The new run went through
the dashboard's own upload code path on test-split genome SAMN11638310 (GCA_014625455.1):

| Stage | Result |
| --- | --- |
| Detection (live AMRFinderPlus) | 37 elements, 35 AMR, 7 mutations, matching the Phase 3 smoke run |
| Predictions | R for ampicillin 0.996, ceftriaxone 0.772, ciprofloxacin 0.998, gentamicin 0.962, tetracycline 0.867; 2 unsupported |

Browser file selection itself could not be automated in this environment, so the upload widget
was exercised up to file selection only. These are integration checks, not accuracy claims.

## Evidence

- [Example run](../reports/phase6/2026-10-14/SAMN11638310/report.md) and
  [PDF report](../reports/phase6/2026-10-14/SAMN11638310_report.pdf)
- [Dashboard](../src/genome_amr_predictor/dashboard.py) and
  [report rendering](../src/genome_amr_predictor/reporting.py)

Uploaded genomes and dashboard runs are written under `data/processed/phase6/` (ignored by git).
Set `GENOME_AMR_RUN_ROOTS` (path-separator list) to browse other run folders.

## Limitations and next phase

The dashboard is a local, single-user research tool with no authentication. Run it only on a
trusted machine. Species is still user-supplied. Analyses run synchronously, so the page waits
1-2 minutes per genome. Phase 7 covers final verification, examples, documentation, and the
demonstration guide.
