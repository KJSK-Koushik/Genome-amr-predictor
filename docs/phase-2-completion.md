# Phase 2 Completion Note

Completed on 2026-10-13 (Asia/Calcutta), using the cohort frozen on 2026-08-18.

## Outcome

All 7,911 selected *Escherichia coli* assemblies were downloaded and verified. Finalization reread
every FASTA, checked its file and normalized sequence SHA-256 fingerprints, and applied the
configured sequence QC thresholds: 4.0–6.5 Mb, at most 500 contigs, and N50 of at least 20 kb.
No additional genomes were excluded. Exact sequence deduplication found no duplicate groups or
duplicate-sequence label conflicts.

| Measure | Result |
|---|---:|
| Verified and retained genomes | 7,911 |
| Observed isolate-antibiotic R/S rows | 47,485 |
| Antibiotics | 7 |
| Indivisible project/BioSample groups | 46 |
| Training genomes | 5,592 (70.7%) |
| Validation genomes | 749 (9.5%) |
| Test genomes | 1,570 (19.8%) |
| Exact sequence duplicate groups | 0 |
| Project overlap between splits | 0 |
| Changed frozen endpoints | 0 |

The seven endpoints are ampicillin, ceftriaxone, ciprofloxacin, gentamicin, tetracycline,
trimethoprim-sulfamethoxazole, and cefepime. Both susceptible and resistant observations are
present for every antibiotic in each split. Whole-project grouping accounts for the departure
from the target 70/15/15 proportions; the original frozen assignments remain unchanged.

Before sequence download, cohort construction removed 29 repeated label rows and one conflicting
BioSample-antibiotic pair. Of 8,273 BioSamples with usable labels, 342 lacked assemblies and
20 selected assemblies failed metadata QC. Exclusion reasons may overlap for an assembly.

## Reproducible evidence

- [Final cohort report](../reports/phase2/2026-08-18/final/final_report.md)
- [Final training rows](../reports/phase2/2026-08-18/final/final_cohort.csv)
- [Final assembly metrics](../reports/phase2/2026-08-18/final/final_assemblies.csv)
- [Per-antibiotic split/class counts](../reports/phase2/2026-08-18/final/split_summary.csv)
- [Final provenance and output hashes](../reports/phase2/2026-08-18/final/final_manifest.json)
- [Independent integrity check](../reports/phase2/2026-08-18/final/verification.json)
- [Data contract](phase-2-data-contract.md) and [split decision](decisions/0003-project-grouped-frozen-splits.md)

The source AST CSV SHA-256 remains
`f9fa051f3fda97420c75d56dfac23784abcb4a8858b63d7f3ece25716b90326e`.
The final manifest records source/configuration hashes, NCBI Datasets 18.35.0, Python 3.12.14,
package source hashes, download provenance, and output fingerprints. FASTA paths are relative to
`data/raw/phase2-genomes`; consumers join by assembly accession. Raw AST, genome downloads,
and the local virtual environment remain ignored by Git.

Run the independent artifact checks with:

```powershell
.venv/Scripts/python.exe scripts/verify_phase2.py --report-dir reports/phase2/2026-08-18
```

Verification passed for output fingerprints, frozen endpoint preservation, assembly joins,
sequence QC metrics, unique BioSamples/accessions/sequence digests, split separation, and class
coverage. The package test suite passed all 43 tests, and Ruff passed.

## Limitations and next phase

Project grouping and exact sequence hashes do not establish that closely related strains are
absent across splits. The accepted split decision requires a genomic-similarity sensitivity
analysis before final model performance claims.

Laboratory methods and testing-standard years are not separately reported in the source export.
Testing standards are mixed: CLSI 44,804 rows, EUCAST 2,561, SFM 96, and unspecified 24.
Reported R/S labels are preserved; breakpoint harmonization is not claimed. The frozen table's
measurement-unit fields are unpopulated, so measurement-based analysis must consult the raw MIC
(mg/L) and disk-diffusion (mm) columns. Missing CheckM values remain documented QC warnings.

Phase 3 is AMRFinderPlus detection integration: executable/database provenance, detector output
parsing, ontology mapping, and fixture-based validation. Phase 2 does not produce detector results
or train phenotype models.
