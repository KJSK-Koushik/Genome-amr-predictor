# Phase 1 Completion Note

## Outcome

Phase 1 establishes a reproducible repository, a documented data contract, automated tests,
continuous integration, and an NCBI BioSample antibiogram audit pipeline. It does not download
genomes, detect resistance genes, or train a prediction model.

The full audit was retrieved on 2026-08-17 using the query
`Escherichia coli[Organism] AND antibiogram[filter]`.

| Measure | Result |
|---|---:|
| Matching and processed BioSamples | 11,075 |
| Parsed AST rows | 171,471 |
| Unique antibiotics | 96 |
| Resistant or susceptible rows | 131,831 |
| Rows excluded from the binary endpoint | 39,640 |
| Duplicate isolate-antibiotic pairs | 52 |
| Conflicting binary-label pairs | 6 |

## Antibiotic decision for Phase 2

The following antibiotics pass the configured minimum count and class-balance thresholds:

- ampicillin;
- ceftriaxone;
- ciprofloxacin;
- gentamicin;
- tetracycline;
- trimethoprim-sulfamethoxazole; and
- cefepime.

Meropenem is deferred because resistant records form only 3.4% of its binary labels, below the
10% Phase 1 threshold. A highly imbalanced endpoint can produce misleading accuracy and unstable
resistance estimates.

These seven antibiotics are a provisional shortlist, not a frozen modeling cohort. Phase 2 must
join labels to downloadable assemblies, remove conflicts and duplicate genomes, check genome
quality and lineage leakage, and review laboratory-method and testing-standard consistency.

## Label policy refinement

Synthetic genome sequences may later be useful for software tests or controlled augmentation
experiments, but synthetic phenotype labels must not be mixed into the primary evaluation data.
Fabricated clinical labels would make model performance unreliable. The main labels therefore
remain observed resistant or susceptible results; intermediate and ambiguous outcomes are kept
for auditing but excluded from the first binary endpoint.

## Reproducible evidence

- `reports/phase1/2026-08-17/audit_report.md` contains the readable audit.
- `reports/phase1/2026-08-17/data_manifest.json` records provenance and configuration.
- `reports/phase1/2026-08-17/candidate_recommendations.csv` contains the selection evidence.
- `configs/phase1_audit.json` contains the versioned thresholds and candidate list.

Phase 2 requires explicit approval before implementation begins.
