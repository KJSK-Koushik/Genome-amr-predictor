# Phase 5 Completion Note

Completed on 2026-10-13 (Asia/Calcutta).

## Outcome

A single `analyze` command now connects genome input validation, the five pretrained phenotype
models, and AMRFinderPlus detection. The output is one versioned run directory with JSON, CSV,
and Markdown outputs. The stages fail independently: a detector failure leaves predictions intact,
and a model failure leaves detection evidence intact. Neither failure is reported as an absent gene
or a susceptible phenotype. Input integrity failures invalidate both outputs. Saved detection
results can be reused only after full provenance and raw-to-structured verification. Existing
run directories are never overwritten. See the [integrated analysis contract](phase-5-engine-contract.md).

## Validation

All **91 package tests passed** and **Ruff passed**. Engine tests cover verified reuse and export
hashes. They also cover detection rejection for a changed input, altered evidence, altered raw
output, a changed database, and a failed manifest. The remaining tests cover model failure that
retains detection, invalid FASTA that skips all stages, overwrite refusal, and the CLI path.

Three real runs on the frozen-cohort training genome SAMN03075588 (GCF_000770275.1):

| Run | Detection | Prediction | Status | Wall time |
| --- | --- | --- | --- | ---: |
| `live`: native AMRFinderPlus | completed, 18 hits (17 AMR, 7 mutations) | completed | completed | 84 s |
| `reuse`: cache from `live` | completed, reused; CSVs byte-identical to `live` | completed | completed | 20 s |
| `mismatched-cache`: `live` cache on GCA_059054765.1 | failed: different input genome; no hits exported | completed | partial (exit 1) | — |

Predictions for SAMN03075588:

| Antibiotic | Prediction | Resistance score | Model features present |
| --- | --- | ---: | ---: |
| Ampicillin | R | 0.988 | 52 / 1000 |
| Ceftriaxone | R | 0.948 | 429 / 1000 |
| Ciprofloxacin | R | 0.998 | 547 / 1000 |
| Gentamicin | S | 0.214 | 0 / 1000 |
| Tetracycline | R | 0.949 | 1000 / 1000 |
| Cefepime | unsupported | — | — |
| Trimethoprim/sulfamethoxazole | unsupported | — | — |

The gentamicin call had none of the model's features present, so its score reflects the model's
baseline. This is consistent with the absence of resistance-associated k-mers, but it carries
little positive evidence. An independent verifier passed artifact hashes, input/snapshot hashes,
raw/structured detection equality, prediction input binding, reuse equality, and
mismatched-cache rejection. This is integration validation on one genome, not an accuracy claim.

## Evidence and reproduction

- [Live run report](../reports/phase5/2026-10-13/live/report.md) and
  [analysis](../reports/phase5/2026-10-13/live/analysis.json)
- [Reuse run analysis](../reports/phase5/2026-10-13/reuse/analysis.json)
- [Mismatched-cache analysis](../reports/phase5/2026-10-13/mismatched-cache/analysis.json)
- [Independent verification](../reports/phase5/2026-10-13/verification.json)

```powershell
./scripts/run_phase5.ps1 -OutputDir reports/phase5/new/live
./scripts/run_phase5.ps1 -OutputDir reports/phase5/new/reuse -DetectionCache reports/phase5/new/live/detection
$env:PYTHONPATH = "src"; .venv/Scripts/python.exe scripts/verify_phase5.py --report-dir reports/phase5/new
```

The verifier also expects a `mismatched-cache` run that applies the `live` cache to another genome.
Input snapshots are ignored by git.

## Limitations and next phase

Species is user-supplied. Pretrained scores are not locally calibrated, and their accuracy on our
cohort is unevaluated. The original full-cohort Phase 4 modeling remains suspended, with 228
detection calls preserved. Phase 6 adds the Streamlit dashboard, charts, and JSON/CSV/PDF downloads
on top of this run directory format.
