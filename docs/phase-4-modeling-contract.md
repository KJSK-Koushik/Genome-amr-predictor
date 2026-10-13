# Phase 4 modeling contract

Status: implemented pipeline; full cohort computation in progress. No real model metrics or
approved models are claimed until the full run and independent verification finish.

## Frozen inputs and features

The input cohort contains 7,911 sequence-verified *E. coli* genomes and 47,485 observed R/S labels
for seven antibiotics. Original BioProject-grouped train/validation/test assignments remain fixed.
`configs/phase4_modeling.json` records modeling choices before test evaluation.

`scripts/run_phase4_detections.py` runs every genome through AMRFinderPlus 4.2.7 and database
2026-08-07.1. The native Ubuntu installation is `/opt/genome-amr/amrfinderplus/4.2.7`.
The active cohort run uses 16 independent workers with one detector thread each, BLAST+ 2.12.0,
and HMMER 3.4. Complete cached calls are accepted only after input, sequence, executable,
database, raw-output, and structured-output fingerprints match. Failed/interrupted runs are
preserved before retry; a failure is never encoded as absence of genes.

The feature builder verifies all detector calls against the frozen assemblies. It uses AMR gene
and mutation presence, distinct detected loci, hierarchy families, curated class/subclass labels,
native match methods, scope, and warnings. Stress/virulence elements do not become AMR features.
Detected loci are not claimed to be functional gene copy counts. No mechanism feature is inferred
when the curated mechanism annotation is absent. Sample IDs, project IDs, laboratory metadata,
phenotypes, geography, and lineage assignments are not model features.

Two feature sets are compared: detector-only and detector plus measured length, contig count,
N50, GC, and ambiguity. Dictionary vocabulary and scaling are fit on training genomes only,
including inside calibration folds. Unknown validation/test gene features are ignored by the
frozen vocabulary rather than adding columns after fitting.

## Genomic similarity sensitivity

`scripts/run_phase4_similarity.py` verifies FASTA fingerprints, sketches all genomes with Mash 2.3
(k=21, sketch size=1,000, hash seed=42), and computes single-linkage components at Mash distances
0.005 and 0.01 with p-value at most 1e-5. These are approximate whole-genome similarity groups,
not definitive lineage assignments or species calls.

The frozen project split remains the primary evaluation. For each distance cutoff, a stricter
test subset excludes every test genome whose component contains any training or validation
genome, including connections through other genomes. Its counts and metrics are reported using
the same frozen model and threshold. Empty/single-class subsets are explicitly unestimable for
ROC/AP; limited surviving cohorts cannot support broad claims. No test-driven regrouping or
retuning is performed. See the accepted [split decision](decisions/0003-project-grouped-frozen-splits.md)
and [Mash tutorial](https://mash.readthedocs.io/en/latest/tutorials.html).

## Selection, calibration, and evaluation

Each antibiotic compares balanced logistic regression (C=0.1, 1, 10) with random forest
(200 trees; maximum depth 10 or unrestricted; minimum leaf size 2) on both feature sets.
Validation average precision selects the candidate; validation Brier score breaks ties.
The prevalence-only baseline uses the training resistance prevalence.

The selected classifier's project-separated training out-of-fold probabilities fit a sigmoid
calibrator. There are up to five stratified group folds; training folds must contain both classes.
The calibrator is retained only when it preserves monotonicity and improves validation Brier
score. Validation MCC selects the probability threshold, with ties resolved toward 0.5.
Training/validation data may therefore influence these decisions; test labels cannot.
See [scikit-learn calibration guidance](https://scikit-learn.org/stable/modules/calibration.html).

All seven artifacts and choices are saved before any test evaluation. Outputs include sensitivity,
specificity, precision, F1, MCC, balanced accuracy, ROC AUC, average precision, Brier score, log loss,
confusion matrices, reliability bins, and project-bootstrap intervals. They also include quality
ablation comparisons, genomic-similarity sensitivity, and testing-standard/platform strata.
Missing methods/standard years and mixed standards remain limitations, not harmonized labels.

Models are research candidates. A flag records whether the selected model beats the prevalence
baseline on both test AP and Brier; failure must prevent claims of useful confidence. No artifact
is automatically approved for a user-facing registry. Influential features show association,
not causal explanation. Random forest fitting/inference uses one worker for reproducibility.

## Commands and artifacts

Install modeling dependencies with `uv pip install --python .venv/Scripts/python.exe -e ".[dev,ml]"`.
Existing .venv Python is used because the machine's global Python launcher is not reliable.

Native detection and similarity jobs run under Ubuntu WSL with `PYTHONPATH=src` from the project.
Their saved state is in `data/processed/phase4/detections/batch_progress.json` and
`data/processed/phase4/similarity/similarity_manifest.json`; diagnostics are in
`data/processed/phase4-detection.log` and `data/processed/phase4-similarity.log`.

After both full jobs complete:

```powershell
./scripts/run_phase4_models.ps1
```

The wrapper verifies complete cohort coverage, builds verified features, fits the models, evaluates
the frozen test set, and verifies saved model inference against saved predictions. Features and
models remain ignored in `data/processed/phase4/features` and `models/phase4/2026-10-13`.
Compact model cards, metrics, test predictions, frozen decisions, and hashes are written under
`reports/phase4/2026-10-13/models`. A completion note is written only after these artifacts pass.

Software tests use constructed fixtures only to verify behavior, including invariance of model
choices and weights when test labels are flipped. Fixtures never enter the real evaluation cohort.
# Deadline scope change — 2026-10-13

The user requested pretrained models for tomorrow morning's deadline. The original contract
below remains the deferred full-cohort research plan. Its detection batch and heartbeat are
suspended. The alternative pretrained delivery, coverage and verification are documented in
[pretrained-deadline-delivery.md](pretrained-deadline-delivery.md). It does not satisfy the
original independent evaluation/calibration gates and must not be reported as completed
full-cohort Phase 4 research.

