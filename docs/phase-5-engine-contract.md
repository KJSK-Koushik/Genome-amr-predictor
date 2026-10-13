# Phase 5 Integrated Analysis Contract

The `analyze` command runs one assembled *E. coli* genome through input validation, pretrained
phenotype prediction, and AMRFinderPlus known-determinant detection, then writes one run
directory. Detection evidence and phenotype predictions stay separate outputs throughout.

## Command

```powershell
./scripts/run_phase5.ps1 -Fasta "C:/path/genome.fna.gz" -SampleId "my-sample" -OutputDir "reports/my-run"
```

or directly:

```powershell
$env:PYTHONPATH = "src"
.venv/Scripts/python.exe -m genome_amr_predictor analyze --fasta genome.fna.gz `
    --sample-id my-sample --output-dir reports/my-run [--detection-cache <detection dir>]
```

The frozen configuration is [`configs/phase5_engine.json`](../configs/phase5_engine.json).
Exit codes: `0` completed, `1` partial (one stage failed), `2` failed (input rejected or both
stages failed). The output directory must be new or empty; an existing run is never overwritten.

## Stages

1. **Input validation**: the configuration schema and species must match. The FASTA (plain or gzip)
   must parse, and contig identifiers must be nonempty and unique. The input is copied to a
   snapshot whose SHA-256 is checked immediately and again at the end. Any failure here
   makes both scientific outputs `not_available`.
2. **Phenotype prediction**: the pretrained registry hash must match the configuration. Prediction
   runs on the snapshot, and its recorded input hash must equal the snapshot hash.
3. **Detection**: the executable SHA-256 and database `version.txt` must match the configuration.
   Without `--detection-cache`, AMRFinderPlus runs live into `detection/`. With a cache, existing
   evidence is reused only after verification (below).

A failure in one scientific stage is recorded and the other stage still runs.

## Detection reuse verification

Cached or freshly produced detection evidence is accepted only if all of these hold:

- manifest status `completed` and evidence schema `amrfinder-evidence-v1`;
- input file SHA-256 equals the current snapshot;
- software version, database version, organism `Escherichia`, nucleotide mode, and no `--plus`;
- adapter code, executable, and every reference-database file hash match;
- `detections.json` and `amrfinder.tsv` hashes match the manifest; and
- re-parsing the raw TSV reproduces the normalized hits exactly.

Otherwise detection is `failed` with the reason, evidence is `null`, `detections.csv` is header-only,
and the warning "gene presence/absence is unknown" is added. A rejected or failed detector is never
reported as "no genes found". Reused hits keep the `sample_id` of the run that produced them.

## Run directory

| File | Content |
| --- | --- |
| `run_manifest.json` | Progress record updated at each stage; final copy adds artifact SHA-256 values |
| `analysis.json` | Status, input QC metrics, both stage results, errors, warnings, code/config hashes |
| `predictions.csv` | One row per antibiotic: prediction, resistance score, status, feature coverage |
| `detections.csv` | One row per detected gene, allele, or point mutation |
| `report.md` | Human-readable summary with limitations and errors |
| `input.fna[.gz]` | Input snapshot (ignored by git) |
| `detection/` | Live AMRFinderPlus TSV, stderr log, normalized evidence, and manifest |

Run status is `completed` when both stages succeed, `partial` when one succeeds, and `failed`
otherwise. Empty CSVs from failed stages carry no biological meaning.

## Versioning

Every run records the package version, configuration hash, engine hash, prediction-code hash,
pretrained registry hash (via configuration), detector executable/database fingerprints, and
detection adapter hash. Changing any of these invalidates cached detection reuse.

## Scientific boundary

- Species is user-supplied; automatic species identification is not performed.
- Pretrained scores are external model outputs, not locally calibrated probabilities. Their
  accuracy on this cohort and their training overlap with it are unknown.
- Cefepime and trimethoprim/sulfamethoxazole return explicit `unsupported` rows.
- Absence of a detected determinant does not establish susceptibility.
- Research and educational use only; laboratory confirmation is required.
