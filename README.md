# Genome AMR Predictor

Genome AMR Predictor is a research and educational system for analyzing assembled bacterial
genomes. The planned system has two deliberately separate outputs:

1. detection of known antimicrobial-resistance genes and supported mutations; and
2. machine-learning prediction of resistant or susceptible phenotypes for selected antibiotics.

The first release is scoped to *Escherichia coli*. It is not a clinical diagnostic system and
does not recommend treatment. Laboratory confirmation is required.

## Current status

**Phase 1 - project foundation and data audit.**

The completed 2026-08-17 audit processed 11,075 *E. coli* BioSamples and 171,471 AST rows.
See the [Phase 1 completion note](docs/phase-1-completion.md) and the
[generated audit report](reports/phase1/2026-08-17/audit_report.md).

This phase provides:

- a reproducible Python 3.11 package;
- typed records and label-normalization rules;
- a public NCBI BioSample antibiogram ingestion path;
- support for later AST Browser CSV/TSV imports;
- compact data-quality and class-balance reports; and
- tests and continuous integration.

Genome downloads, AMRFinderPlus integration, model training, and the web interface are outside
Phase 1 and require later approval gates.

## Quick start

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
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

## Repository layout

```text
configs/                         Versioned audit configuration
data/                            Ignored raw/interim/processed data directories
docs/                            Architecture, data contract, decisions, overview report
models/                          Ignored future model artifacts
reports/phase1/                  Committed compact audit summaries
scripts/                         Reproducible command wrappers
src/genome_amr_predictor/        Application package
tests/                           Unit and integration-style fixture tests
```

## Data sources

- [NCBI Pathogen Detection AST documentation](https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/)
- [NCBI BioSample](https://www.ncbi.nlm.nih.gov/biosample/)
- [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)

## Safety boundary

Outputs from this project are computational research results. A detected gene does not always
guarantee phenotypic resistance, and an ML prediction is not a treatment recommendation.
