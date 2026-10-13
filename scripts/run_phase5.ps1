param(
    [string]$Fasta = "data/raw/phase2-genomes/shard-0028/package/ncbi_dataset/data/GCF_000770275.1/GCF_000770275.1_ASM77027v1_genomic.fna.gz",
    [string]$SampleId = "SAMN03075588",
    [string]$OutputDir = "",
    [string]$DetectionCache = ""
)

$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "src"
if (-not $OutputDir) { $OutputDir = "data/processed/phase5/run-" + [guid]::NewGuid().ToString("N") }
$taskArgs = @("-m", "genome_amr_predictor", "analyze", "--fasta", $Fasta,
    "--sample-id", $SampleId, "--output-dir", $OutputDir)
if ($DetectionCache) { $taskArgs += @("--detection-cache", $DetectionCache) }
& .venv/Scripts/python.exe @taskArgs
if ($LASTEXITCODE -ne 0) { throw "Analysis is incomplete. Inspect $OutputDir/analysis.json for stage errors." }
Write-Output "Analysis saved to $OutputDir/analysis.json"
