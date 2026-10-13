# Examples

Three public *E. coli* assemblies from the frozen Phase 2 cohort, one per split, are the project's
reference examples. Genomes are not committed; fetch them from NCBI (about 4.8 MB in total):

```powershell
.venv/Scripts/python.exe examples/fetch_examples.py
```

Files are saved to `data/raw/examples/<BioSample>.fna.gz` (ignored by git). Each download is checked
against the SHA-256 of its uncompressed sequence, so a changed or corrupted file is rejected.

| BioSample | Assembly | Split | Detected elements (AMR / mutations) | Pretrained predictions |
| --- | --- | --- | --- | --- |
| SAMN03075588 | GCF_000770275.1 | train | 18 (17 / 7) | R: AMP, CRO, CIP, TET; S: GEN* |
| SAMN03177615 | GCF_000819375.1 | validation | 6 (6 / 0) | R: AMP, TET; S: CRO, CIP, GEN |
| SAMN11638310 | GCA_014625455.1 | test | 37 (35 / 7) | R: AMP, CRO, CIP, GEN, TET |

\* Gentamicin had no model features present; the score is the model baseline.
Cefepime and trimethoprim/sulfamethoxazole are always reported as unsupported.

## Run an example

```powershell
./scripts/run_phase5.ps1 -Fasta data/raw/examples/SAMN11638310.fna.gz -SampleId SAMN11638310 `
    -OutputDir data/processed/examples/SAMN11638310
$env:PYTHONPATH = "src"
.venv/Scripts/python.exe -m genome_amr_predictor export-pdf `
    --run-dir data/processed/examples/SAMN11638310 --output SAMN11638310_report.pdf
```

Or start the dashboard (`./scripts/run_dashboard.ps1`) and upload the file.

## Saved example outputs

| Output | What it shows |
| --- | --- |
| [Phase 5 live run](../reports/phase5/2026-10-13/live/report.md) | SAMN03075588, fresh detection and predictions |
| [Phase 5 reuse run](../reports/phase5/2026-10-13/reuse/analysis.json) | Same genome with verified reuse of saved detection |
| [Mismatched cache](../reports/phase5/2026-10-13/mismatched-cache/report.md) | Detection rejected for a different genome; status `partial` |
| [Phase 6 PDF report](../reports/phase6/2026-10-14/SAMN11638310_report.pdf) | SAMN11638310, full printable report |

These are integration examples, not an assessment of predictive accuracy. Results are for research
and education only and require laboratory confirmation.
