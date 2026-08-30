# ITGB1 figures 4-6 reproducibility package

This repository regenerates manuscript Figures 4, 5 and 6 and Supplementary Figures S4, S5 and S6 from the included plot-ready source data. It contains no machine-specific paths and no superseded exploratory plotting scripts.

## Contents

- `data/figure4/`: source tables for the GSE96058 survival, GDSC pharmacogenomic, and corrected three-cohort pathway panels.
- `data/figure5/`: source tables for TCGA immune correlations, TCGA project-adjusted regressions, GEO immune-exclusion analyses, and the Supplementary Figure S5 specificity analysis.
- `data/figure6/`: processed expression TPM matrices and clinical matching files for GSE71340, GSE240671 and GSE212810.
- `scripts/`: one script for Figures 4-5, one for Supplementary Figures S4-S5, and one for Figure 6/S6.
- `outputs/`: regenerated PNG/PDF figures and Figure 6 source tables, created by the pipeline.

## Reproduce all figures

From the repository root:

```bash
Rscript scripts/install_r_dependencies.R
python3 -m pip install -r requirements.txt
bash scripts/run_all.sh
```

The pipeline writes `Figure4`, `Figure5`, `Figure6`, `Supplementary_FigureS4`, `Supplementary_FigureS5`, and `Supplementary_FigureS6` as both PNG and PDF files in `outputs/`.

## Software

The plotting scripts were tested with R 4.4.2 and Python 3.11. The R workflow uses only CRAN packages. Python dependencies are pinned in `requirements.txt`.

## Analysis scope

The repository contains plot-ready data rather than all raw sequencing files. The included source tables are exactly the inputs used for plotting, and data provenance is listed in `data/README.md`. Figure 4G is calculated from independently ranked gastric, pancreatic and liver tumour gene lists. Figure 5F is a descriptive running-enrichment curve and does not report a permutation-derived enrichment P value.

## Data provenance

Public datasets include GSE96058, GSE113255, GSE172356, GSE202069, GSE71340, GSE240671, GSE212810, TCGA and Cell Model Passports/GDSC. Please cite the original dataset resources when reusing the source data.
