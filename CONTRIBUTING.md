# Contributing

## Development workflow

1. Create a focused branch from `main`.
2. Keep data provenance and version information with every generated report.
3. Do not commit raw genomes, raw BioSample exports, credentials, or trained model binaries.
4. Add or update tests for behavior changes.
5. Run `python -m pytest` and `python -m ruff check .` before opening a pull request.

## Scientific rules

- Keep ARG detection separate from phenotype prediction.
- Never fabricate laboratory phenotype labels.
- Never place synthetic examples in validation or test sets.
- Prevent duplicate and lineage leakage before reporting model performance.
- Record dataset, reference-database, code, and model versions.
- Describe the project as research and educational software, not clinical decision support.
