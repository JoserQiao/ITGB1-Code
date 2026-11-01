# Use 'Continue' and check $LASTEXITCODE per step so benign stderr output from
# native Python modules (for example, an openpyxl warning while reading xlsx)
# does not abort the pipeline as if it were a real failure.
$ErrorActionPreference = 'Continue'
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$env:PYTHONNOUSERSITE = '1'
$WindowsPythonExecutable = Join-Path $ProjectRoot '.venv\Scripts\python.exe'
$UnixPythonExecutable = Join-Path $ProjectRoot '.venv/bin/python'
$PythonExecutable = $null

if (Test-Path -LiteralPath $WindowsPythonExecutable -PathType Leaf) {
    $PythonExecutable = $WindowsPythonExecutable
}
elseif (Test-Path -LiteralPath $UnixPythonExecutable -PathType Leaf) {
    $PythonExecutable = $UnixPythonExecutable
}

if (-not $PythonExecutable -or -not (Test-Path -LiteralPath $PythonExecutable -PathType Leaf)) {
    throw "Virtual-environment Python was not found under .venv. Follow the README quick-start instructions first."
}

function Invoke-PipelineStep {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Module,
        [string[]]$ExtraArguments = @()
    )

    Write-Host "`n==> $Module" -ForegroundColor Cyan
    & $PythonExecutable -m $Module --project-root $ProjectRoot @ExtraArguments
    if ($LASTEXITCODE -ne 0) {
        throw "Pipeline step '$Module' failed with exit code $LASTEXITCODE."
    }
}

Push-Location $ProjectRoot
try {
    # 1. Public-data acquisition
    Invoke-PipelineStep 'scripts.acquisition.download_data'
    Invoke-PipelineStep 'scripts.acquisition.download_cmp'
    Invoke-PipelineStep 'scripts.acquisition.download_external_solid_tumor'
    Invoke-PipelineStep 'scripts.acquisition.download_imvigor210'

    # 2. Cohort preprocessing
    Invoke-PipelineStep 'scripts.preprocessing.prepare_geo'
    Invoke-PipelineStep 'scripts.preprocessing.prepare_external_ici'
    Invoke-PipelineStep 'scripts.preprocessing.prepare_external_solid_tumor'
    $Rscript = Get-Command Rscript -ErrorAction SilentlyContinue
    if (-not $Rscript) {
        throw 'Rscript was not found. IMvigor210 preprocessing requires R with Bioconductor DESeq2.'
    }
    & $Rscript.Source (Join-Path $ProjectRoot 'scripts/preprocessing/export_imvigor210.R') `
        (Join-Path $ProjectRoot 'data/raw/IMvigor210/cds.RData') `
        (Join-Path $ProjectRoot 'config/signatures.tsv') `
        (Join-Path $ProjectRoot 'data/processed/IMvigor210')
    if ($LASTEXITCODE -ne 0) { throw 'IMvigor210 R preprocessing failed.' }

    # 3. Primary analyses
    Invoke-PipelineStep 'scripts.analysis.analyze_response'
    Invoke-PipelineStep 'scripts.analysis.analyze_rectal_response'
    Invoke-PipelineStep 'scripts.analysis.analyze_single_cell'
    Invoke-PipelineStep 'scripts.analysis.analyze_imvigor210'
    Invoke-PipelineStep 'scripts.analysis.analyze_external_ici'
    Invoke-PipelineStep 'scripts.analysis.analyze_power_precision'
    Invoke-PipelineStep 'scripts.analysis.analyze_gdsc_lineage_adjusted'
    Invoke-PipelineStep 'scripts.analysis.analyze_external_solid_tumor'
    Invoke-PipelineStep 'scripts.analysis.analyze_external_marker_controls'

    # 4. Robustness and external validation
    Invoke-PipelineStep 'scripts.validation.robustness_analysis' @('--bootstrap', '5000', '--cv-repeats', '200')
    Invoke-PipelineStep 'scripts.validation.extended_auc_bootstrap' @('--bootstrap', '5000')
    Invoke-PipelineStep 'scripts.validation.alternative_mcpcounter'
    Invoke-PipelineStep 'scripts.validation.validate_integrated_extension'

    # 5. Figures
    Invoke-PipelineStep 'scripts.visualization.plot_all_figures_ordered'
    Invoke-PipelineStep 'scripts.visualization.plot_integrated_extension'
}
finally {
    Pop-Location
}
