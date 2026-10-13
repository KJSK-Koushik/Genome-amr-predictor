"""Import the pinned MIT HZI models as portable JSON; never execute upstream classes."""

from __future__ import annotations

import io
import json
import pickle
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

from genome_amr_predictor.detection import sha256

COMMIT = "717fa43d002060a3759efaf616f9dd887f6c0db1"
BASE = f"https://raw.githubusercontent.com/hzi-bifo/AMR_prediction_pipeline/{COMMIT}/"
DRUGS = ["ampicillin", "ceftriaxone", "ciprofloxacin", "gentamicin", "tetracycline"]


class State:
    def __init__(self, *args):
        pass

    def __setstate__(self, state):
        self.__dict__.update(state)


class RestrictedReader(pickle.Unpickler):
    def find_class(self, module, name):
        states = {
            ("sklearn.pipeline", "Pipeline"),
            ("sklearn.ensemble._forest", "RandomForestClassifier"),
            ("sklearn.tree._classes", "DecisionTreeClassifier"),
            ("sklearn.tree._tree", "Tree"),
            ("sklearn.linear_model._logistic", "LogisticRegression"),
        }
        if (module, name) in states:
            return State
        arrays = {
            ("numpy.core.multiarray", "_reconstruct"): np._core.multiarray._reconstruct,
            ("numpy.core.multiarray", "scalar"): np._core.multiarray.scalar,
            ("numpy", "ndarray"): np.ndarray,
            ("numpy", "dtype"): np.dtype,
        }
        if (module, name) not in arrays:
            raise ValueError(f"Unsupported pickle class: {module}.{name}")
        return arrays[module, name]


def import_models(cache=Path("data/external/pretrained-hzi"), output=Path("models/pretrained-hzi")):
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    files = {}
    for drug in DRUGS:
        paths = {
            f"{drug}_model.sav": "AMR_software/PhenotypeSeeker/bin/model/kma/"
            f"Escherichia_coli/Escherichia_coli_{drug}_f1_macro_finalized_model.sav",
            f"{drug}_vocab.zip": "AMR_software/PhenotypeSeeker/bin/feature/"
            f"Escherichia_coli/{drug}_vocab.npy.zip",
        }
        for local, remote in paths.items():
            path = cache / local
            if not path.exists():
                with urllib.request.urlopen(BASE + remote, timeout=60) as response:
                    path.write_bytes(response.read())
            files[local] = {"url": BASE + remote, "sha256": sha256(path)}
        with zipfile.ZipFile(cache / f"{drug}_vocab.zip") as archive:
            data = archive.read(f"Escherichia_coli/{drug}_vocab.npy")
        stream = io.BytesIO(data)
        if np.lib.format.read_magic(stream) != (1, 0):
            raise ValueError("Unexpected vocabulary NPY format")
        shape, _, _ = np.lib.format.read_array_header_1_0(stream)
        vocab = RestrictedReader(stream).load().tolist()
        pipeline = RestrictedReader(io.BytesIO((cache / f"{drug}_model.sav").read_bytes())).load()
        if len(pipeline.steps) != 1 or pipeline.steps[0][0] != "cl":
            raise ValueError("Unexpected preprocessing pipeline")
        model = pipeline.steps[0][1]
        if model.classes_.tolist() != [0, 1] or shape != (model.n_features_in_,):
            raise ValueError("Model/vocabulary mismatch")
        if any(not isinstance(k, str) or len(k) != 13 or set(k) - set("ACGT") for k in vocab):
            raise ValueError("Invalid vocabulary")
        result = {"antibiotic": drug, "vocabulary": vocab, "classes": ["S", "R"]}
        if hasattr(model, "coef_"):
            result.update(
                kind="logistic_regression",
                weights=model.coef_[0].tolist(),
                intercept=float(model.intercept_[0]),
            )
        else:
            trees = []
            for estimator in model.estimators_:
                state = estimator.tree_
                values = state.values[:, 0, :]
                trees.append(
                    {
                        "left": state.nodes["left_child"].tolist(),
                        "right": state.nodes["right_child"].tolist(),
                        "feature": state.nodes["feature"].tolist(),
                        "threshold": state.nodes["threshold"].tolist(),
                        "p_r": (values[:, 1] / values.sum(axis=1)).tolist(),
                    }
                )
            result.update(kind="random_forest", trees=trees)
        (output / f"{drug}.json").write_text(json.dumps(result), encoding="utf-8")
    license_path = cache / "LICENSE"
    if not license_path.exists():
        with urllib.request.urlopen(BASE + "LICENSE", timeout=60) as response:
            license_path.write_bytes(response.read())
    (output / "LICENSE").write_bytes(license_path.read_bytes())
    manifest = {
        "schema": "hzi-pretrained-v1",
        "species": "Escherichia coli",
        "repository": "https://github.com/hzi-bifo/AMR_prediction_pipeline",
        "commit": COMMIT,
        "license": "MIT",
        "scenario": "kma",
        "k": 13,
        "feature_transform": "binary canonical k-mer presence in upstream vocabulary order",
        "source_files": files,
        "artifacts": {f"{d}.json": sha256(output / f"{d}.json") for d in DRUGS},
        "limitations": [
            "External models; not validated on our cohort.",
            "Training overlap with our NCBI/PATRIC cohort is unknown.",
            "Scores are not locally calibrated probabilities.",
            "Cefepime and trimethoprim/sulfamethoxazole unsupported.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"imported": DRUGS, "output": str(output)}, indent=2))


if __name__ == "__main__":
    import_models()
