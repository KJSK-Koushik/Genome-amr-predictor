# Genome AMR Predictor

Genome AMR Predictor is a research and educational system for analysing assembled
*Escherichia coli* genomes. It produces two deliberately separate outputs:

1. **Known resistance determinants**: genes, alleles and point mutations detected by NCBI
   AMRFinderPlus 4.2.7 (database 2026-08-07.1); and
2. **Predicted phenotypes**: resistant/susceptible calls for ampicillin, ceftriaxone,
   ciprofloxacin, gentamicin and tetracycline from published pretrained models. Cefepime and
   trimethoprim/sulfamethoxazole are reported as unsupported.

It is not a clinical diagnostic system and does not recommend treatment. Laboratory confirmation
is required.

## Quick start

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev,ml,ui]"
python scripts/import_pretrained_hzi.py          # pinned pretrained models
python examples/fetch_examples.py                # three verified example genomes
python scripts/verify_delivery.py                # tests + every phase verifier
./scripts/run_dashboard.ps1                      # http://localhost:8501
```

Clone to a short path such as `D:\amr`, because of Windows path limits. AMRFinderPlus runs in
WSL Ubuntu. See the [user guide](docs/user-guide.md) for detector setup, inputs, the command line,
how to read results, and troubleshooting.

## Documentation

| Document | Purpose |
| --- | --- |
| [User guide](docs/user-guide.md) | Setup, running analyses, interpreting outputs |
| [Demonstration guide](docs/demonstration-guide.md) | 12-minute live demo script with fallbacks |
| [Examples](examples/README.md) | Public example genomes and saved outputs |
| [Phase 7 completion](docs/phase-7-completion.md) | Final verification results and known limitations |
| [Pretrained models](docs/pretrained-deadline-delivery.md) | Model source, attribution, verification |
| [Integrated analysis contract](docs/phase-5-engine-contract.md) | Run stages, failure rules, output files |
| [Detection contract](docs/phase-3-detection-contract.md) | AMRFinderPlus adapter and provenance |
| [Data contract](docs/phase-2-data-contract.md) | Cohort construction and frozen splits |

## Project status

All seven planned phases are delivered. The one exception is Phase 4's locally trained models:
full-cohort training was paused in favour of pretrained models, so the project makes **no
accuracy claim on its own cohort**.

| Phase | Outcome | Evidence |
| --- | --- | --- |
| 1 Foundation and audit | 11,075 BioSamples and 171,471 AST rows audited | [note](docs/phase-1-completion.md) |
| 2 Dataset | 7,911 verified genomes, 47,485 R/S labels; 5,592 / 749 / 1,570 project-grouped splits | [note](docs/phase-2-completion.md) |
| 3 Detection | AMRFinderPlus adapter; 4/4 NCBI reference tests; 3 live genomes | [note](docs/phase-3-completion.md) |
| 4 Modelling | Pipeline and similarity groups (7,911 genomes) built; cohort detection paused at 228 genomes; pretrained HZI models imported and verified | [contract](docs/phase-4-modeling-contract.md), [pretrained](docs/pretrained-deadline-delivery.md) |
| 5 Integrated engine | `analyze` command, isolated stage failures, verified detection reuse | [note](docs/phase-5-completion.md) |
| 6 Interface | Streamlit dashboard; PDF/JSON/CSV/Markdown downloads | [note](docs/phase-6-completion.md) |
| 7 Delivery | Clean-clone tests, all-phase verification, guides, examples | [note](docs/phase-7-completion.md) |

## Commands

| Command | Purpose |
| --- | --- |
| `./scripts/run_dashboard.ps1` | Dashboard: upload a genome or open a saved run |
| `genome-amr analyze --fasta … --sample-id … --output-dir …` | One integrated analysis (exit 0/1/2 = completed/partial/failed) |
| `genome-amr export-pdf --run-dir … --output …` | PDF from a verified run |
| `genome-amr detect` / `parse-detection` | AMRFinderPlus only |
| `python scripts/verify_delivery.py [--portable]` | Re-run all checks; `--portable` skips detector/genome checks |
| `genome-amr audit-ncbi` / `audit-file` | Phase 1 AST label audits ([note](docs/phase-1-completion.md)) |

Run `genome-amr` as `python -m genome_amr_predictor` if the console script is not on PATH.

## Repository layout

```text
configs/                         Frozen audit, cohort, modelling and engine configuration
data/                            Ignored raw genomes, interim files and analysis runs
docs/                            Guides, phase notes, contracts and decisions
examples/                        Example genome fetcher and guide
models/                          Ignored imported model artifacts
reports/                         Committed evidence per phase (genome copies are ignored)
scripts/                         Reproducible run and verification wrappers
src/genome_amr_predictor/        Application package, including dashboard.py
tests/                           Unit, integration and dashboard tests
```

## Data and model sources

- [NCBI Pathogen Detection AST documentation](https://www.ncbi.nlm.nih.gov/pathogens/docs/ast/)
- [NCBI BioSample](https://www.ncbi.nlm.nih.gov/biosample/)
- [NCBI Datasets genome download documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)
- [NCBI genome assembly report fields](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-assembly/)
- [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/)
- [HZI AMR prediction pipeline](https://github.com/hzi-bifo/AMR_prediction_pipeline) (MIT), source of the pretrained models

## Safety boundary

Outputs from this project are computational research results. A detected gene does not always
mean phenotypic resistance, and the absence of detected genes does not establish susceptibility.
An ML prediction is not a treatment recommendation. Prediction scores are not calibrated
probabilities, and the species is assumed rather than verified.
