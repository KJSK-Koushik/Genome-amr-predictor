# Phase 7 Completion Note

Completed on 2026-10-14 (Asia/Calcutta).

## Outcome

The project is packaged for delivery:

- one verification command,
- a clean-clone check of the committed tree,
- verified public example genomes,
- a user guide and a rehearsed demonstration guide, and
- a README that describes the delivered system rather than the plan.

| Deliverable | Location |
| --- | --- |
| Delivery verification | [`scripts/verify_delivery.py`](../scripts/verify_delivery.py), [result](../reports/phase7/2026-10-14/delivery_verification.json) |
| Example genomes | [`examples/`](../examples/README.md): fetcher with sequence-hash checks |
| User guide | [user-guide.md](user-guide.md) |
| Demonstration guide | [demonstration-guide.md](demonstration-guide.md) |

## Verification

`scripts/verify_delivery.py` runs these checks and records the git commit, Python, package and
detector versions. It never refits models or edits evidence.

| Check | Result |
| --- | --- |
| Ruff lint | passed |
| Test suite | 114 passed |
| Phase 2 cohort | frozen sources, outputs, splits and label rules match the manifest |
| Phase 3 detection | 3 live genomes: raw/structured equality, installed executable and database hashes |
| Pretrained models | 4,999 k-mer features match GenomeTester4; scores match scikit-learn within 2.3e-16 for 5 models |
| Phase 5 engine | live, reuse and mismatched-cache runs |
| Saved runs and PDF | all 4 committed runs pass integrity checks and render to PDF |

**Clean clone.** A fresh clone of commit `cd3fdda` passed lint, all 114 tests, and the portable
checks, using the project's interpreter with the clone's own `src/` first on the import path.
The Phase 2 check could not run in the clone, because one file's full path there reached
Windows' 260-character limit. That came from the deep scratch folder, not the repository, and
the same check passes in the project folder. The deepest committed path is 79 characters, so the guides advise
a short clone location. A second virtual environment could not be built in that scratch folder
for the same path-length reason; CI installs from scratch on every push.

**Demo rehearsal.** The live demo step was run on validation genome SAMN03177615. It completed in
80 seconds with 6 detected elements and the predictions given in the guide: R for ampicillin and
tetracycline, S for ceftriaxone, ciprofloxacin and gentamicin. The example fetcher accepted the
three local genomes by sequence hash and rejected a corrupted file. All three NCBI download URLs
responded.

## Known limitations

- **No local accuracy claim.** Full-cohort detection stopped at 228 of 7,911 genomes, and the
  project's own models were not trained or test-evaluated. The pretrained models are external,
  their training overlap with this cohort is unknown, and their scores are not calibrated.
- **Species is assumed**, not identified.
- **Two antibiotics unsupported**: cefepime and trimethoprim/sulfamethoxazole.
- **Single-user local dashboard** with no authentication; analyses run synchronously
  (1–2 minutes).
- **Upload widget**: browser file selection was not automated; the code path after selection was
  tested with real genomes.
- **Windows path length**: clone to a short path.

## Possible next steps

1. Resume Phase 4 cohort detection and train, calibrate and test-evaluate the project's own models
   on the frozen splits, using the similarity groups for leakage-aware evaluation.
2. Evaluate the pretrained models on the frozen test split to measure their accuracy locally.
3. Add species identification (for example, a Mash or average-nucleotide-identity check) before
   analysis.
