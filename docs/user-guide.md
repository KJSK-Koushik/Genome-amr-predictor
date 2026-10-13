# User Guide

Genome AMR Analyzer takes one assembled *Escherichia coli* genome. It reports known antimicrobial
resistance determinants (genes, alleles and point mutations from NCBI AMRFinderPlus). It also
reports predicted resistant/susceptible (R/S) phenotypes for five antibiotics from published
pretrained models. The two outputs are kept separate. Results are for research and education only,
are not a clinical diagnosis, and require laboratory confirmation.

## 1. Setup (Windows 11)

Clone to a short path (for example `D:\amr`). Some committed report paths are 79 characters
deep, and Windows limits full paths to 260 characters.

```powershell
py -3.11 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev,ml,ui]"
python scripts/import_pretrained_hzi.py      # pinned HZI models -> models/pretrained-hzi
```

**AMRFinderPlus** is needed for determinant detection. It runs inside Ubuntu on WSL:

1. Install WSL Ubuntu, then BLAST+, HMMER and Perl inside it.
2. Install AMRFinderPlus 4.2.7 under `/opt/genome-amr/amrfinderplus/4.2.7` (no spaces in the path),
   following [NCBI's instructions](https://github.com/ncbi/amr/wiki/Installing-AMRFinder).
3. Install database release `2026-08-07.1` into `data/2026-08-07.1` beside the executable.

[`configs/phase5_engine.json`](../configs/phase5_engine.json) pins the executable SHA-256 and
database version. A different installation is refused rather than silently used. Without the
detector, phenotype predictions still run and detection is reported as failed (unknown).

Check the installation (takes about 2 minutes):

```powershell
.venv/Scripts/python.exe scripts/verify_delivery.py
```

## 2. Input requirements

- One assembled genome in FASTA (`.fna`, `.fa`, `.fasta`, `.fas`), optionally gzip-compressed.
  Raw reads (FASTQ) are not accepted.
- Contig identifiers must be nonempty and unique.
- The species must be *E. coli*; this is assumed, not verified.
- Predictions are withheld if total length is outside 3–7 Mb or more than 5% of bases are
  ambiguous.

## 3. Running an analysis

**Dashboard** (recommended):

```powershell
./scripts/run_dashboard.ps1        # then open http://localhost:8501
```

Choose *Analyze a new genome*, upload the FASTA file, enter a sample ID, and select
*Run analysis* (1–2 minutes). Choose *Open a saved run* to revisit earlier results.

**Command line:**

```powershell
$env:PYTHONPATH = "src"
.venv/Scripts/python.exe -m genome_amr_predictor analyze --fasta genome.fna.gz `
    --sample-id MY_SAMPLE --output-dir data/processed/my-run
.venv/Scripts/python.exe -m genome_amr_predictor export-pdf --run-dir data/processed/my-run `
    --output my-run.pdf
```

Exit codes: `0` completed, `1` partial (one stage failed), `2` failed. Output folders must be new
or empty; earlier results are never overwritten.

## 4. Reading the results

| Item | Meaning |
| --- | --- |
| Run status | `completed`: both stages ran. `partial`: one failed. `failed`: input rejected |
| Resistance score | Model output from 0 to 1; R when above 0.5. **Not a calibrated probability** |
| Model features present | How many of the model's 1,000 k-mer features occur in the genome |
| Note: *no model features present* | Score is the model's baseline; little positive evidence |
| Note: *close to the decision threshold* | Score within 0.05 of 0.5; treat the call as uncertain |
| `unsupported` | No validated model (cefepime, trimethoprim/sulfamethoxazole) |
| `withheld_qc` | Assembly failed prediction QC; no call is made |
| Detected element | Known determinant matched by AMRFinderPlus, with method, identity and coverage |
| *Detection failed: presence/absence is unknown* | The detector did not produce trusted evidence; this is **not** "no genes" |
| *No known determinants detected* | Does not establish susceptibility |

A detected determinant and a predicted phenotype can disagree. Neither replaces laboratory
antimicrobial susceptibility testing.

## 5. Outputs and integrity

Each run folder contains `analysis.json` (everything), `predictions.csv`, `detections.csv`,
`report.md`, a copy of the input, detector output under `detection/`, and `run_manifest.json`.
The manifest records hashes of every output plus software, database, model and configuration
versions. The dashboard and `export-pdf` verify these hashes and refuse a run that was edited.

## 6. Troubleshooting

| Symptom | Fix |
| --- | --- |
| Detection failed: checksum or version differs | AMRFinderPlus or its database differs from the pinned configuration; reinstall the recorded versions |
| Detection failed after a long wait | Start WSL (`wsl -d Ubuntu -e true`) and retry; the timeout is 30 minutes |
| Prediction failed: registry differs | Re-run `scripts/import_pretrained_hzi.py`; do not edit `models/` |
| *Filename too long* when cloning or testing | Clone to a shorter path, or run `git config core.longpaths true` |
| Port 8501 in use | `./scripts/run_dashboard.ps1 -Port 8502` |
| *Output directory must be empty* | Choose a new run folder |
