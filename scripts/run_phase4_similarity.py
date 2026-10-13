"""Build Mash sketches and two conservative genomic-similarity component assignments."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

from genome_amr_predictor.detection import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemblies", type=Path, required=True)
    parser.add_argument("--genome-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    with args.assemblies.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    root = args.genome_dir.resolve()
    paths = []
    by_path = {}
    for i, row in enumerate(rows):
        path = (root / row["fasta_path"]).resolve()
        if not path.is_relative_to(root) or sha256(path) != row["fasta_file_sha256"]:
            raise ValueError("Frozen genome fingerprint mismatch")
        paths.append(str(path))
        by_path[str(path)] = i
    file_list = args.output_dir / "genome_files.txt"
    file_list.write_text("\n".join(paths) + "\n")
    sketch = args.output_dir / "cohort.msh"
    progress = {
        "status": "running",
        "genomes": len(rows),
        "assembly_table_sha256": sha256(args.assemblies),
        "mash_version": subprocess.check_output(["mash", "--version"], text=True).strip(),
        "k": 21,
        "sketch_size": 1000,
        "seed": 42,
        "distance_thresholds": [0.005, 0.01],
        "maximum_p_value": 1e-5,
    }
    progress_path = args.output_dir / "similarity_manifest.json"
    progress_path.write_text(json.dumps(progress, indent=2) + "\n")
    with (args.output_dir / "mash.stderr.log").open("w") as log:
        sketch_command = [
            "mash",
            "sketch",
            "-l",
            "-p",
            str(args.threads),
            "-k",
            "21",
            "-s",
            "1000",
            "-S",
            "42",
            "-o",
            str(args.output_dir / "cohort"),
            str(file_list),
        ]
        subprocess.run(sketch_command, check=True, stderr=log)
        progress["sketch_command"] = sketch_command
        progress["sketch_sha256"] = sha256(sketch)
        progress["stage"] = "distances"
        progress_path.write_text(json.dumps(progress, indent=2) + "\n")
        parents = {threshold: list(range(len(rows))) for threshold in (0.005, 0.01)}

        def find(parent, index):
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        command = [
            "mash",
            "dist",
            "-p",
            str(args.threads),
            "-d",
            "0.01",
            "-v",
            "0.00001",
            str(sketch),
            str(sketch),
        ]
        count = 0
        with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True) as process:
            for line in process.stdout:
                reference, query, distance, pvalue, _shared = line.rstrip("\n").split("\t")
                i, j = by_path[reference], by_path[query]
                if i >= j:
                    continue
                distance = float(distance)
                if float(pvalue) > 1e-5:
                    raise ValueError("Unexpected Mash p-value")
                count += 1
                for threshold, parent in parents.items():
                    if distance <= threshold:
                        left, right = find(parent, i), find(parent, j)
                        parent[max(left, right)] = min(left, right)
            if process.wait() != 0:
                raise RuntimeError("Mash distance calculation failed")
        components = args.output_dir / "similarity_groups.csv"
        groups = []
        for i, row in enumerate(rows):
            result = {"biosample_accession": row["biosample_accession"], "split": row["split"]}
            for threshold, parent in parents.items():
                result[f"mash_{threshold}"] = rows[find(parent, i)]["biosample_accession"]
            groups.append(result)
        with components.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(groups[0]))
            writer.writeheader()
            writer.writerows(groups)
        progress.update(
            status="completed",
            stage="completed",
            close_pairs=count,
            distance_command=command,
            groups_sha256=sha256(components),
        )
        progress_path.write_text(json.dumps(progress, indent=2) + "\n")
        print(json.dumps(progress, indent=2))


if __name__ == "__main__":
    main()
