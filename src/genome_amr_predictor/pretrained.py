"""Portable inference for the published HZI PhenotypeSeeker models."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from datetime import UTC, datetime
from pathlib import Path

from genome_amr_predictor.detection import sha256

UNSUPPORTED = ["cefepime", "trimethoprim/sulfamethoxazole"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical(kmer: str) -> str:
    return min(kmer, kmer.translate(COMPLEMENT)[::-1])


def extract_presence(path: Path, vocabulary: set[str], k: int = 13) -> tuple[set[str], dict]:
    """Keep FASTA contigs separate; ambiguous bases break windows."""
    wanted = {canonical(word) for word in vocabulary}
    found = set()
    bases = ambiguous = contigs = 0
    ids = set()
    tail = ""
    digest = hashlib.sha256()
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        for line in stream:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                identifier = line[1:].split()
                if not identifier or identifier[0] in ids:
                    raise ValueError("Missing or duplicate FASTA contig identifier")
                ids.add(identifier[0])
                contigs += 1
                tail = ""
                digest.update(b"\n>")
                continue
            if not contigs:
                raise ValueError("FASTA must start with a contig header")
            sequence = line.upper()
            if set(sequence) - set("ACGTRYSWKMBDHVN"):
                raise ValueError("Invalid nucleotide characters in FASTA")
            digest.update(sequence.encode("ascii"))
            bases += len(sequence)
            ambiguous += sum(base not in "ACGT" for base in sequence)
            sequence = tail + sequence
            for index in range(len(sequence) - k + 1):
                word = sequence[index : index + k]
                # Membership checks avoid computing reverse complements for every window.
                if word in vocabulary:
                    found.add(canonical(word))
            tail = sequence[-(k - 1) :]
    if bases < k:
        raise ValueError("Genome is empty or shorter than one feature")
    # Query lists identify either strand of each canonical k-mer.
    return found & wanted, {
        "bases": bases,
        "contigs": contigs,
        "ambiguous_bases": ambiguous,
        "sequence_sha256": digest.hexdigest(),
    }


def model_score(model: dict, features: list[int]) -> float:
    if model["kind"] == "logistic_regression":
        value = model["intercept"] + sum(
            w * x for w, x in zip(model["weights"], features, strict=True)
        )
        return 1 / (1 + math.exp(-value)) if value >= 0 else math.exp(value) / (1 + math.exp(value))
    scores = []
    for tree in model["trees"]:
        node = 0
        while tree["left"][node] != -1:
            branch = (
                "left" if features[tree["feature"][node]] <= tree["threshold"][node] else "right"
            )
            node = tree[branch][node]
        scores.append(tree["p_r"][node])
    return sum(scores) / len(scores)


def predict_pretrained(fasta: Path, model_dir: Path) -> dict:
    manifest_path = model_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["schema"] != "hzi-pretrained-v1":
        raise ValueError("Unsupported pretrained registry")
    models = []
    vocabulary = set()
    for filename, expected in manifest["artifacts"].items():
        path = model_dir / filename
        if sha256(path) != expected:
            raise ValueError(f"Model checksum mismatch: {filename}")
        model = json.loads(path.read_text(encoding="utf-8"))
        models.append(model)
        for word in model["vocabulary"]:
            vocabulary.update((word, word.translate(COMPLEMENT)[::-1]))
    found, quality = extract_presence(fasta, vocabulary)
    warnings = list(manifest["limitations"])
    if quality["bases"] < 3_000_000 or quality["bases"] > 7_000_000:
        warnings.append("Assembly length outside the expected E. coli range; predictions withheld.")
    if quality["ambiguous_bases"] / quality["bases"] > 0.05:
        warnings.append("More than 5% ambiguous bases; predictions withheld.")
    withhold = len(warnings) > len(manifest["limitations"])
    rows = []
    for model in models:
        features = [int(canonical(word) in found) for word in model["vocabulary"]]
        score = model_score(model, features) if not withhold else None
        # sklearn binary LR uses >0 decision; RF uses class-order argmax (tie -> S).
        rows.append(
            {
                "antibiotic": model["antibiotic"],
                "prediction": ("R" if score > 0.5 else "S") if score is not None else None,
                "resistance_score": score,
                "status": "withheld_qc" if withhold else "predicted",
                "features_present": sum(features),
                "features_total": len(features),
            }
        )
    rows.extend(
        {
            "antibiotic": drug,
            "prediction": None,
            "resistance_score": None,
            "status": "unsupported",
            "features_present": None,
            "features_total": None,
        }
        for drug in UNSUPPORTED
    )
    return {
        "schema": "pretrained-prediction-v1",
        "created_at": datetime.now(UTC).isoformat(),
        "input": str(fasta),
        "input_sha256": sha256(fasta),
        "quality": quality,
        "species": "Escherichia coli (user-supplied; not automatically verified)",
        "predictions": rows,
        "warnings": warnings,
        "model_provenance": manifest,
        "registry_sha256": sha256(manifest_path),
        "use": "Research/education only; requires laboratory confirmation.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, default=Path("models/pretrained-hzi"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = predict_pretrained(args.fasta, args.model_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "prediction.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    with (args.output_dir / "prediction.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result["predictions"][0]))
        writer.writeheader()
        writer.writerows(result["predictions"])
    print(json.dumps(result["predictions"], indent=2))


if __name__ == "__main__":
    main()
