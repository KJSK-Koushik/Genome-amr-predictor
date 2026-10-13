param([int]$Port = 8501)

$ErrorActionPreference = "Stop"
$env:PYTHONPATH = "src"
& .venv/Scripts/python.exe -m streamlit run src/genome_amr_predictor/dashboard.py `
    --server.port $Port --server.headless true --browser.gatherUsageStats false
