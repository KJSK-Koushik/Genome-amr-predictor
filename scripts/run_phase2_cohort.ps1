param(
    [string]$AstFile,
    [string]$PythonExecutable = ".venv/Scripts/python.exe",
    [string]$DatasetsExecutable = "data/external/ncbi-datasets/datasets.exe",
    [ValidateRange(1, 30)]
    [int]$MaxWorkers = 30,
    [string]$OutputDate = (Get-Date -Format "yyyy-MM-dd"),
    [switch]$ResumeExistingCohort
)

$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "src"
$reportDir = "reports/phase2/$OutputDate"
$genomeDir = "data/raw/phase2-genomes"

if (-not $ResumeExistingCohort) {
    if (-not $AstFile) { throw "AstFile is required unless ResumeExistingCohort is set." }
    & $PythonExecutable -m genome_amr_predictor construct-cohort `
    --input $AstFile `
    --config configs/phase2_cohort.json `
    --output-dir $reportDir `
    --assembly-cache data/interim/phase2/full_assemblies.jsonl
    if ($LASTEXITCODE -ne 0) { throw "Cohort construction failed: $LASTEXITCODE" }
}

& $PythonExecutable -m genome_amr_predictor download-genomes `
    --accessions "$reportDir/assembly_accessions.txt" `
    --output-dir $genomeDir `
    --datasets-executable $DatasetsExecutable `
    --shard-size 250 `
    --max-workers $MaxWorkers
if ($LASTEXITCODE -ne 0) { throw "Genome download/verification failed: $LASTEXITCODE" }

& $PythonExecutable -m genome_amr_predictor finalize-cohort `
    --cohort "$reportDir/cohort.csv" `
    --assemblies "$reportDir/assemblies.csv" `
    --genome-manifest "$genomeDir/genome_manifest.csv" `
    --genome-dir $genomeDir `
    --config configs/phase2_cohort.json `
    --output-dir "$reportDir/final"
if ($LASTEXITCODE -ne 0) { throw "Sequence-cohort finalization failed: $LASTEXITCODE" }

& $PythonExecutable scripts/verify_phase2.py --report-dir $reportDir
if ($LASTEXITCODE -ne 0) { throw "Final cohort verification failed: $LASTEXITCODE" }
