# Phase 3 detection contract

The detection layer accepts assembled *Escherichia coli* DNA FASTA and produces curated gene,
allele, and supported mutation evidence using NCBI AMRFinderPlus. The adapter uses nucleotide
mode with `--organism Escherichia`, `--print_node`, an explicit database directory, and native
curated thresholds. It leaves `--plus` disabled for production calls. The species setting is a
scope assumption, not an independent species classifier.

## Runtime

This workspace uses AMRFinderPlus **4.2.7**, database **2026-08-07.1**, Ubuntu through WSL,
BLAST+ 2.12.0, and HMMER 3.4. The release archive is cached in ignored
`data/external/amrfinderplus`; the active Linux installation is
`/opt/genome-amr/amrfinderplus/4.2.7` inside Ubuntu. The AMRFinderPlus executable and database are
accessible from Windows through `\\wsl.localhost\Ubuntu\opt\genome-amr\amrfinderplus\4.2.7`.
The AMRFinderPlus executable and complete database file set are
fingerprinted with each run. Database updates are explicit setup operations; detection never
updates the database or follows a moving `latest` path implicitly.

For setup on another machine, use NCBI's [installation instructions](https://github.com/ncbi/amr/wiki/Installing-AMRFinder).
The official binary archive used here is
[amrfinder_binaries_v4.2.7.tar.gz](https://github.com/ncbi/amr/releases/download/amrfinder_v4.2.7/amrfinder_binaries_v4.2.7.tar.gz).
It requires Linux with BLAST+, HMMER, and Perl. Install the Linux binaries in a directory without
spaces or shell metacharacters: the upstream StxTyper subprocess does not safely quote all internal
paths. The genome input and result directories may remain under the original workspace path.
On Windows the adapter translates mounted-drive paths with
`wslpath` and invokes the chosen WSL distribution without constructing a shell command.
Install the database using `amrfinder -u`, then select a dated database directory explicitly.
Future reproduction must use the recorded database release, not assume that `-u` selects this version.

## Running detection

```powershell
.venv/Scripts/python.exe -m genome_amr_predictor detect `
  --fasta path/to/genome.fna.gz `
  --sample-id SAMN_example `
  --executable '\\wsl.localhost\Ubuntu\opt\genome-amr\amrfinderplus\4.2.7\amrfinder' `
  --database '\\wsl.localhost\Ubuntu\opt\genome-amr\amrfinderplus\4.2.7\data\2026-08-07.1' `
  --wsl-distribution Ubuntu `
  --output-dir reports/local-detection/new-run
```

Omit `--wsl-distribution` on Linux when using native paths. Output directories must be empty.
FASTA validation checks syntax, DNA alphabet, sequence metrics, and unique contig identifiers.
Compressed FASTA is decompressed temporarily inside the workspace and removed after execution.
The adapter records failures explicitly and rejects malformed or missing detector output.

## Outputs and ontology

- `amrfinder.tsv`: untouched native detector output.
- `amrfinder.stderr.log`: native successful execution diagnostics.
- `detections.json`: `amrfinder-evidence-v1` structured evidence and provenance.
- `run_manifest.json`: completion/failure, command, input/sequence hashes, adapter hash,
  executable hash, database file hashes, version report, and output hashes.

Each hit retains symbol/name, native method, scope/type/subtype, 1-based inclusive coordinates,
strand, identity/coverage/alignment fields, reference accession/name, HMM fields, hierarchy nodes,
and its original row. An allele symbol is exposed only for the detector's ALLELE method.
POINT/POINT_DISRUPT records are mutations; stress and virulence calls imported from plus output
remain separate element types. Partial, internal-stop, unfamiliar-method, and plus-scope hits
carry warnings without being silently discarded. Missing numeric fields remain null.

The versioned `ncbi-class-preservation-v1` ontology mapping preserves NCBI class/subclass and
hierarchy identifiers. It does not turn a broad class into a claim about an individual antibiotic.
AMRFinderPlus TSV does not provide a separate curated mechanism column: mechanism is explicitly
null with status `not_separately_reported_by_detector`. Mechanism enrichment would require a
separately sourced, versioned mapping; mechanism names are not inferred from gene-name substrings.
See [NCBI output fields](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#fields) and
[reference hierarchy documentation](https://www.ncbi.nlm.nih.gov/pathogens/docs/gene_hierarchy/).

A valid header-only table is a completed no-hit result. Empty files, truncated rows, invalid
numeric metrics, and failed commands are errors. No hits means that no supported determinant was
detected; susceptibility remains unknown. Detection does not establish phenotypic resistance.
See [NCBI interpretation guidance](https://github.com/ncbi/amr/wiki/Interpreting-results).

## Offline imports and validation

`parse-detection` accepts a TSV plus explicit software/database version strings. Its provenance
marks those versions as user supplied and not verified against an executable. Current headers and
older Gene symbol/Sequence name/closest-sequence headers are supported; unfamiliar extra columns
are retained in each raw row.

The fixture tests exercise the subprocess boundary without inventing biological validation data.
`scripts/run_phase3_smoke.py --output-dir reports/phase3/new-smoke` selects one real genome per
frozen split, checks its frozen file/sequence fingerprints, executes detection, and writes a smoke
summary. This verifies integration, not population-scale detection accuracy or phenotype prediction.
