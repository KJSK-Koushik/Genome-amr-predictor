# Pretrained model delivery — 2026-10-13

The user changed the modeling scope to pretrained models because the project is due tomorrow
morning. Five published E. coli phenotype models are downloaded, imported and runnable locally.
This completes the pretrained inference step, not the original full-cohort Phase 4 evaluation.

## Source and model coverage

Models come from [HZI's AMR prediction pipeline](https://github.com/hzi-bifo/AMR_prediction_pipeline),
commit `717fa43d002060a3759efaf616f9dd887f6c0db1`, under the MIT license. The associated
[2024 benchmark paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11070729/) describes the approach.
We fixed the upstream `kma` (low-similarity) scenario before any local phenotype evaluation.
No local test labels were used to select models, parameters, thresholds or this scenario.

| Antibiotic | Imported algorithm | Features | Availability |
| --- | --- | --- | --- |
| Ampicillin | Random forest, 60 trees | 1,000 binary 13-mers | Ready |
| Ceftriaxone | Random forest, 10 trees | 1,000 binary 13-mers | Ready |
| Ciprofloxacin | Random forest, 60 trees | 1,000 binary 13-mers | Ready |
| Gentamicin | Logistic regression | 1,000 binary 13-mers | Ready |
| Tetracycline | Random forest, 40 trees | 1,000 binary 13-mers | Ready |
| Cefepime | — | — | Explicit unsupported result |
| Trimethoprim/sulfamethoxazole | — | — | Explicit unsupported result |

Trimethoprim alone is a different endpoint and is not substituted for the combination drug.
RADAR was also reviewed, but its promoter/expression/copy-number features are unavailable in
our adapter; it was not integrated simply to fill the two missing endpoints.

## Run the demo

From the project folder in PowerShell:

```powershell
./scripts/run_pretrained.ps1
```

For another assembled E. coli genome:

```powershell
./scripts/run_pretrained.ps1 -Fasta "C:/path/genome.fna.gz" -OutputDir "reports/my-genome"
```

Outputs are `prediction.json` (including quality, input hash, model hashes, source commit,
limitations and unsupported endpoints) and `prediction.csv`. Inference reads portable JSON,
requires only the standard library, and handles FASTA or gzipped FASTA. The importer requires
the existing `.venv` ML dependencies and network access if the pinned files are absent.
The downloaded files and model artifacts are ignored by Git; the importer reconstructs them.
Upstream pickle classes are replaced with inert state containers under an explicit allowlist;
inference never loads pickle. The copied MIT license is in `models/pretrained-hzi/LICENSE`.

## Verification

- All 4,999 distinct selected vocabulary strings matched native GenomeTester4 13-mer queries
  on the real Phase 3 train demonstration genome: zero binary feature mismatches.
- All five portable models matched current scikit-learn execution of the original fitted
  tree/linear parameters across 35 vectors each, including the real genome: 175 comparisons.
  Largest score difference was 2.22e-16; all binary predictions matched.
- Original scikit-learn 0.24.1 tree state was adapted for the current tree executor by adding
  the unused missing-value flag and normalizing class counts. This verifies the conversion;
  it is not an execution of the entire legacy Snakemake environment.
- The package suite passed 81 tests. New checks cover strand invariance, FASTA line wrapping,
  contig boundaries, ambiguous bases, malformed input, model checksum failure and abstention.
- Three real Phase 3 demonstration genomes completed inference in 1.910, 1.970 and 2.256
  seconds each, with input hashes verified against the frozen assembly manifest. All produced
  five predictions and two unsupported results. Evidence is under
  `reports/pretrained/2026-10-13/`; timings are in `smoke_summary.json` and exact parity checks
  are in `verification.json`.

The native feature verification uses `glistmaker -w 13 -c 1` and `glistquery -f words.txt`.
The preserved inputs/logs are in ignored `data/external/pretrained-hzi/verification`.
`scripts/verify_pretrained_hzi.py` checks these outputs and the model conversion.

## Scientific limits and retained work

These are externally trained research models, not models trained by this project. Published
performance is not our cohort performance. Overlap between their PATRIC training genomes and
our NCBI cohort is unknown; do not present local demonstrations as independent accuracy tests.
Resistance scores are model outputs, not locally calibrated confidence or treatment advice.
Species is supplied by the user, not automatically identified. Assemblies outside 3–7 Mb or
with more than 5% ambiguous bases receive withheld predictions. These pragmatic gates are
not a validated clinical QC protocol. Known-gene detection remains separate from phenotype
prediction; the verified AMRFinderPlus adapter and three detection examples remain available.

Phase 2 and Phase 3 remain complete. Full-cohort similarity computation is complete for 7,911
genomes, but the small purged holdouts reveal a major generalization limitation. The original
AMRFinder cohort batch is suspended with 228 successful calls preserved; interrupted outputs
remain incomplete and are not interpreted as gene absence. The old training heartbeat is
paused. Full-cohort training, calibration, independent accuracy evaluation, the integrated
engine and web UI remain separate unfinished work. No Phase 4 completion claim is made.
