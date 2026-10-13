"""Download the three public example genomes from NCBI and verify their sequence content."""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from genome_amr_predictor.data.fasta import inspect_fasta  # noqa: E402

BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
# BioSample, frozen split, NCBI path, and SHA-256 of the uncompressed sequence content.
EXAMPLES = {
    "SAMN03075588": (
        "train",
        "GCF/000/770/275/GCF_000770275.1_ASM77027v1",
        "8ab44ed7bc638925e60012bf7986ea1c002dc61dee0eeb9919559cb8b991eebc",
    ),
    "SAMN03177615": (
        "validation",
        "GCF/000/819/375/GCF_000819375.1_Escherichia_coli_CVM_N33429PS_v1.0",
        "b3dc50717420810e66052f4233466f6ba700a569349a694ba08a355519c7eacd",
    ),
    "SAMN11638310": (
        "test",
        "GCA/014/625/455/GCA_014625455.1_PDT000838265.1",
        "a480dc19d6c98786838f9803d8c70cdc9269e47c218b098d01560a3f2aca6c38",
    ),
}


def fetch(output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for sample, (split, folder, expected) in EXAMPLES.items():
        name = folder.rsplit("/", 1)[1] + "_genomic.fna.gz"
        target = output_dir / f"{sample}.fna.gz"
        if not target.is_file():
            temporary = target.with_suffix(".part")
            print(f"Downloading {sample} ({split}) from NCBI ...")
            urllib.request.urlretrieve(f"{BASE}{folder}/{name}", temporary)
            temporary.replace(target)
        actual = inspect_fasta(target).sequence_sha256
        if actual != expected:
            target.unlink()
            raise ValueError(f"{sample}: sequence differs from the verified genome; removed")
        print(f"{sample}: verified -> {target}")
        paths.append(target)
    return paths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data/raw/examples")
    fetch(parser.parse_args().output_dir)
