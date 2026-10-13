param(
    [string]$DetectionDir = "data/processed/phase4/detections",
    [string]$FeatureDir = "data/processed/phase4/features",
    [string]$SimilarityDir = "data/processed/phase4/similarity",
    [string]$OutputDate = "2026-10-13"
)

$ErrorActionPreference = "Stop"
$taskPython = ".venv/Scripts/python.exe"
$taskAssemblies = "reports/phase2/2026-08-18/final/final_assemblies.csv"
$taskCohort = "reports/phase2/2026-08-18/final/final_cohort.csv"
$taskModelDir = "models/phase4/$OutputDate"
$taskReportDir = "reports/phase4/$OutputDate/models"
$taskBatch = Get-Content -LiteralPath "$DetectionDir/batch_progress.json" -Raw | ConvertFrom-Json
$taskGenomeCount = @(Import-Csv -LiteralPath $taskAssemblies).Count
if ($taskBatch.status -ne "completed" -or $taskBatch.requested -ne $taskGenomeCount -or
    ($taskBatch.completed + $taskBatch.cached) -ne $taskGenomeCount -or $taskBatch.failures.Count -ne 0) {
    throw "The full frozen detection cohort is not complete."
}
$taskSimilarity = Get-Content -LiteralPath "$SimilarityDir/similarity_manifest.json" -Raw | ConvertFrom-Json
if ($taskSimilarity.status -ne "completed" -or $taskSimilarity.genomes -ne $taskGenomeCount) {
    throw "The full genomic similarity analysis is not complete."
}
if (-not (Test-Path -LiteralPath "$FeatureDir/feature_manifest.json")) {
    & $taskPython -m genome_amr_predictor build-features --assemblies $taskAssemblies `
        --detection-dir $DetectionDir --output-dir $FeatureDir
    if ($LASTEXITCODE -ne 0) { throw "Feature construction failed." }
}
if (-not (Test-Path -LiteralPath "$taskReportDir/training_manifest.json")) {
    & $taskPython -m genome_amr_predictor train-models --features "$FeatureDir/features.jsonl" `
        --cohort $taskCohort --similarity-groups "$SimilarityDir/similarity_groups.csv" `
        --config configs/phase4_modeling.json --model-dir $taskModelDir --report-dir $taskReportDir
    if ($LASTEXITCODE -ne 0) { throw "Model fitting/evaluation failed; inspect the saved run." }
}
& $taskPython scripts/verify_phase4.py --features "$FeatureDir/features.jsonl" `
    --cohort $taskCohort --similarity-groups "$SimilarityDir/similarity_groups.csv" `
    --model-dir $taskModelDir --report-dir $taskReportDir
if ($LASTEXITCODE -ne 0) { throw "Independent model verification failed." }
