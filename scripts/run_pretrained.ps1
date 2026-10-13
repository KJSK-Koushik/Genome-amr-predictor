param(
    [string]$Fasta = "data/raw/phase2-genomes/shard-0028/package/ncbi_dataset/data/GCF_000770275.1/GCF_000770275.1_ASM77027v1_genomic.fna.gz",
    [string]$OutputDir = "reports/pretrained/demo"
)

$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "src"
if (-not (Test-Path -LiteralPath "models/pretrained-hzi/manifest.json")) {
    & .venv/Scripts/python.exe scripts/import_pretrained_hzi.py
    if ($LASTEXITCODE -ne 0) { throw "Pretrained model import failed." }
}
& .venv/Scripts/python.exe -m genome_amr_predictor.pretrained --fasta $Fasta --output-dir $OutputDir
if ($LASTEXITCODE -ne 0) { throw "Pretrained prediction failed." }
Write-Output "Saved predictions and provenance to $OutputDir/prediction.json and prediction.csv"
