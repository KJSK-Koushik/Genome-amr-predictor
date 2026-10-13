# Demonstration Guide

A 10–12 minute live demonstration of the full pipeline. Every step uses committed evidence or
real public genomes, so it can be repeated.

## Before the session (15 minutes, once)

1. Complete [setup](user-guide.md#1-setup-windows-11) and fetch the examples:
   `.venv/Scripts/python.exe examples/fetch_examples.py`.
2. Run `.venv/Scripts/python.exe scripts/verify_delivery.py`. All checks should report `passed`.
3. Start WSL (`wsl -d Ubuntu -e true`) so the first detection does not wait for it to boot.
4. Start the dashboard with `./scripts/run_dashboard.ps1` and open http://localhost:8501.
5. Keep [the example PDF](../reports/phase6/2026-10-14/SAMN11638310_report.pdf) open as a
   fallback if a live run is slow.

## Script

| Time | Step | What to show | Point to make |
| ---: | --- | --- | --- |
| 0:00 | Problem and scope | README safety boundary | Research tool for *E. coli* assemblies; not clinical |
| 1:00 | Data foundation | [Phase 2 final report](../reports/phase2/2026-08-18/final/final_report.md) | 7,911 verified genomes, 47,485 R/S labels, project-grouped frozen splits |
| 2:30 | Open a saved run | Dashboard, *Open a saved run*, `SAMN03075588 · completed` | QC tiles; status banner |
| 3:30 | Predictions | *Phenotype predictions* tab | Scores vs 0.5 threshold; gentamicin note shows 0 features present; two unsupported drugs |
| 5:00 | Determinants | *Known determinants* tab, filter to `mutation` | gyrA/parC mutations agree with the ciprofloxacin R call; evidence and prediction stay separate |
| 6:30 | Failure safety | Select `mismatch-check · partial` | Cached detection for a different genome is rejected: "presence/absence is unknown", not "no genes" |
| 7:30 | Live analysis | *Analyze a new genome*, upload `data/raw/examples/SAMN03177615.fna.gz`, sample ID `SAMN03177615` | Takes 1–2 minutes; meanwhile open *Provenance* on another run |
| 9:30 | Result | New run: 6 elements; R ampicillin and tetracycline | Validation-split genome; whether the external models saw it in training is unknown |
| 10:30 | Downloads | *Downloads* tab, then the PDF | PDF, JSON, CSV and Markdown; outputs are hash-verified |
| 11:30 | Limitations | *Limitations* tab | External models unvalidated on our cohort; scores uncalibrated; species not verified |

## Talking points for questions

- **Why two outputs?** A known gene is evidence of a mechanism. A model score is a statistical
  prediction. Combining them would hide disagreement
  ([decision 0002](decisions/0002-separate-detection-and-prediction.md)).
- **How accurate is it?** Not yet measured on this cohort. The pretrained models come from the
  HZI benchmark ([delivery note](pretrained-deadline-delivery.md)). Local full-cohort training was
  paused before test evaluation, so the project makes no accuracy claim.
- **What prevents misleading results?** Pinned detector and database versions, checksums on every
  artifact, refusal of edited runs, and explicit unknown/unsupported states. QC withholds
  predictions for implausible assemblies.
- **Reproducibility.** `scripts/verify_delivery.py` re-runs the tests and every phase verifier, and
  records the environment ([latest result](../reports/phase7/2026-10-14/delivery_verification.json)).

## If something goes wrong

| Problem | Recovery |
| --- | --- |
| Upload run takes too long | Show the saved `SAMN11638310` run and its PDF instead |
| Detection fails in the live run | Use it as a demonstration of the failure state; predictions still appear |
| Dashboard will not start | Present the [example PDF](../reports/phase6/2026-10-14/SAMN11638310_report.pdf) and [live run report](../reports/phase5/2026-10-13/live/report.md) |
