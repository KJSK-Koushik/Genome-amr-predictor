# Genome AMR Predictor

Genome AMR Predictor is a research and educational system for analyzing assembled bacterial
genomes. The planned system has two deliberately separate outputs:

1. detection of known antimicrobial-resistance genes and supported mutations; and
2. machine-learning prediction of resistant or susceptible phenotypes for selected antibiotics.

The first release is scoped to *Escherichia coli*. It is not a clinical diagnostic system and
does not recommend treatment. Laboratory confirmation is required.

## Current status

**Deadline scope: five pretrained phenotype models are runnable and verified.**

Run `./scripts/run_pretrained.ps1` in PowerShell for a real genome demonstration, or provide
`-Fasta "C:/path/genome.fna.gz" -OutputDir "reports/my-genome"`. JSON and CSV exports include
source provenance and quality checks. Ampicillin, ceftriaxone, ciprofloxacin, gentamicin and
tetracycline are supported; cefepime and trimethoprim/sulfamethoxazole explicitly abstain.
See the [pretrained delivery note](docs/pretrained-deadline-delivery.md) for attribution,
verification and limitations. These externally trained models have not been evaluated for
independent accuracy on our cohort.

The completed 2026-08-17 audit processed 11,075 *E. coli* BioSamples and 171,471 AST rows.
See the [Phase 1 completion note](docs/phase-1-completion.md) and the
[generated audit report](reports/phase1/2026-08-17/audit_report.md).

Phase 1 provided:

- a reproducible Python 3.11 package;
- typed records and label-normalization rules;
- a public NCBI BioSample antibiogram ingestion path;
- support for later AST Browser CSV/TSV imports;
- compact data-quality and class-balance reports; and
- tests and continuous integration.

Phase 2 produced 7,911 verified genomes and 47,485 observed R/S labels across seven antibiotics.
The frozen splits contain 5,592 training, 749 validation, and 1,570 test genomes; all 46 project
groups stay within one split. Sequence QC and artifact verification passed, with no exact sequence
duplicates. See the [Phase 2 completion note](docs/phase-2-completion.md),
[final cohort report](reports/phase2/2026-08-18/final/final_report.md), and
[data contract](docs/phase-2-data-contract.md). Genomic-similarity sensitivity analysis remains
required before final model performance claims.

Phase 3 adds the AMRFinderPlus adapter, versioned structured result parsing, NCBI class/subclass
and hierarchy mapping, and detection provenance. AMRFinderPlus 4.2.7 with database 2026-08-07.1
passed all four upstream reference tests and live detection on three frozen cohort genomes.
The package suite passes 66 tests. See the [Phase 3 completion note](docs/phase-3-completion.md)
and [detection contract](docs/phase-3-detection-contract.md) for commands, evidence, and limitations.
The web interface remains scheduled for a later phase.

The original Phase 4 detection batch is suspended to prioritize the user-requested pretrained
deadline scope, with 228 successful calls preserved. Genomic-similarity groups are complete
for all 7,911 genomes.
Its modeling pipeline compares logistic regression and random forest, fits calibration on
project-separated training folds, freezes model choices before test evaluation, and writes
research model cards. See the [Phase 4 modeling contract](docs/phase-4-modeling-contract.md).
Real cohort model metrics remain pending full computation and verification.

Phase 5 adds the integrated `analyze` command (`./scripts/run_phase5.ps1`). It validates the input
genome, runs the pretrained phenotype models and AMRFinderPlus, and writes one versioned run
directory with JSON, CSV and Markdown outputs. Stage failures stay isolated, and a failed or
unverifiable detector is never reported as "no genes found". See the
[Phase 5 completion note](docs/phase-5-completion.md) and
[integrated analysis contract](docs/phase-5-engine-contract.md).

## Quick start

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev,ml]"
python -m pytest
```

Run a small live NCBI audit:

```powershell
python -m genome_amr_predictor audit-ncbi `
  --config configs/phase1_audit.json `
  --output-dir reports/local-audit `
  --max-records 100
```

Run the full current BioSample antibiogram audit by setting `--max-records 0`. Raw XML is parsed
in batches and is not retained. The command writes only summaries unless `--write-records` is
provided.

Audit a CSV or TSV exported from the NCBI AST Browser:

```powershell
python -m genome_amr_predictor audit-file `
  --input path\to\ast-export.csv `
  --config configs/phase1_audit.json `
  --output-dir reports/local-file-audit
```

File audits apply the species configured in `phase1_audit.json`, recognize the AST Browser's
`#BioSample` header, repair a known doubled-quote export defect, and record the source file's
SHA-256 fingerprint. Raw downloaded datasets remain ignored by Git.

## Repository layout

```text
configs/                         Versioned audit configuration
data/                            Ignored raw/interim/processed data directories
docs/                            Architecture, data contract, decisions, overview report
models/                          Ignored future model artifacts
reports/phase1/                  Committed compact audit summaries
reports/phase2/                  Frozen cohort manifests and quality reports
scripts/                         Reproducible command wrappers
src/genome_amr_predictor/        Application package
tests/                           Unit and integration-style fixture tests
```

## Data sources

- [NCBI Pathogen Detection AST documentation](https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/)
- [NCBI BioSample](https://www.ncbi.nlm.nih.gov/biosample/)
- [NCBI Datasets genome download documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)
- [NCBI genome assembly report fields](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-assembly/)
- [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)

## Safety boundary

Outputs from this project are computational research results. A detected gene does not always
guarantee phenotypic resistance, and an ML prediction is not a treatment recommendation.
