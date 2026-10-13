# AMRFinderPlus fixture provenance

`amrfinder_v4.2.7_dna.expected.tsv` is the unchanged `test_dna.expected` file extracted from
NCBI's [4.2.7 Linux binary release](https://github.com/ncbi/amr/releases/download/amrfinder_v4.2.7/amrfinder_binaries_v4.2.7.tar.gz).
It is an upstream expected detector output for software testing, not a labeled cohort member or
biological validation set. NCBI's [source and license](https://github.com/ncbi/amr) apply.

The small constructed rows in `test_detection.py` exercise parsing and failure handling only.
