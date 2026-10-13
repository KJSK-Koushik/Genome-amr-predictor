param(
    [int]$MaxRecords = 0,
    [string]$OutputDirectory = "reports/phase1/current",
    [switch]$WriteRecords
)

$arguments = @(
    "-m", "genome_amr_predictor",
    "audit-ncbi",
    "--config", "configs/phase1_audit.json",
    "--output-dir", $OutputDirectory,
    "--max-records", $MaxRecords
)

if ($WriteRecords) {
    $arguments += "--write-records"
}

python @arguments
