"""Check portable inference against sklearn tree/LR execution and GenomeTester features."""

from __future__ import annotations

import io
import json
from pathlib import Path

import numpy as np
from import_pretrained_hzi import DRUGS, RestrictedReader
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier, _tree

from genome_amr_predictor.pretrained import canonical, extract_presence, model_score


def main():
    cache = Path("data/external/pretrained-hzi")
    model_dir = Path("models/pretrained-hzi")
    work = cache / "verification"
    models = {drug: json.loads((model_dir / f"{drug}.json").read_text()) for drug in DRUGS}
    vocab = {word for model in models.values() for word in model["vocabulary"]}
    vocab |= {word.translate(str.maketrans("ACGT", "TGCA"))[::-1] for word in vocab}
    found, _ = extract_presence(work / "genome.fna", vocab)
    native = {}
    for line in (work / "native_counts.tsv").read_text().splitlines():
        word, count = line.split()
        native[word] = int(count) > 0
    expected_words = {word for model in models.values() for word in model["vocabulary"]}
    if set(native) != expected_words:
        raise ValueError("Native query did not return the full vocabulary")
    mismatches = [word for word in native if native[word] != (canonical(word) in found)]
    if mismatches:
        raise ValueError(f"Feature mismatch: {mismatches[:5]}")
    rng = np.random.default_rng(42)
    checks = []
    node_dtype = (
        DecisionTreeClassifier().fit([[0], [1]], [0, 1]).tree_.__getstate__()["nodes"].dtype
    )
    for drug, portable in models.items():
        original = RestrictedReader(io.BytesIO((cache / f"{drug}_model.sav").read_bytes())).load()
        model = original.steps[0][1]
        real = [int(canonical(word) in found) for word in portable["vocabulary"]]
        x = np.vstack([np.zeros(1000), np.ones(1000), rng.integers(0, 2, (32, 1000)), real])
        if hasattr(model, "coef_"):
            reference = LogisticRegression()
            reference.__dict__.update(model.__dict__)
            expected = reference.predict_proba(x)[:, 1]
        else:
            outputs = []
            for estimator in model.estimators_:
                old = estimator.tree_
                nodes = np.zeros(old.node_count, dtype=node_dtype)
                for name in old.nodes.dtype.names:
                    nodes[name] = old.nodes[name]
                values = old.values / old.values.sum(axis=2, keepdims=True)
                tree = _tree.Tree(1000, np.array([2], dtype=np.intp), 1)
                tree.__setstate__(
                    {
                        "max_depth": old.max_depth,
                        "node_count": old.node_count,
                        "nodes": nodes,
                        "values": values,
                    }
                )
                outputs.append(tree.predict(x.astype(np.float32))[:, 1])
            expected = np.mean(outputs, axis=0)
        actual = np.array([model_score(portable, row.tolist()) for row in x])
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
        np.testing.assert_array_equal(actual > 0.5, expected > 0.5)
        checks.append(
            {
                "antibiotic": drug,
                "vectors_checked": len(x),
                "max_score_difference": float(np.max(np.abs(actual - expected))),
            }
        )
    output = Path("reports/pretrained/2026-10-13/verification.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "status": "passed",
        "native_feature_words_checked": len(native),
        "feature_mismatches": len(mismatches),
        "models": checks,
        "reference": "Current sklearn executes original weights/trees; legacy tree state "
        "adapted by adding missing-value flag and normalizing stored class counts. "
        "No fitting or parameter selection performed.",
        "native_features": "GenomeTester4 glistmaker -w 13 -c 1; glistquery -f words.txt",
    }
    output.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
