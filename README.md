# ITGB1 cellular-origin and composition-aware analysis

This repository reproduces the manuscript:

> Cellular sources, composition dependence and context-specific drug-sensitivity associations of bulk ITGB1 in public cancer datasets

## Questions addressed

1. Which tumour-microenvironment compartments contribute to bulk `ITGB1` expression?
2. How much bulk `ITGB1` variation is associated with transcriptomic composition?
3. Does `ITGB1` predict objective response, progression-free survival, or pathological complete response after composition adjustment?
4. Does `ITGB1` add discrimination beyond established clinical variables in a larger ICI cohort?
5. Do direct pathology estimates and CPTAC protein measurements support the transcript-based interpretation?
6. Which external pathway and GDSC associations remain after composition or cancer-lineage adjustment?

## Evidence layers

- GSE195832: paired single-cell source attribution.
- GSE103322: independent patient-level single-cell source validation (5,902 annotated cells, 18 patients).
- GSE91061, GSE78220, GSE176307, GSE126044 and IMvigor210: pretreatment ICI objective response (total n=477).
- IMvigor210 additionally supports clinical-baseline and TMB-adjusted incremental-value analyses.
- GSE135222: ICI progression-free survival.
- GSE190826: neoadjuvant chemoradiotherapy pathological complete response.
- GSE172356, GSE113255 and GSE76427: continuous composition-aware pathway validation.
- Ten TCGA cancers: pathologist-estimated stromal, tumour-cell and tumour-nuclei percentages matched to ITGB1 expression.
- Seven CPTAC studies: paired ITGB1 RNA and mass-spectrometry protein abundance.
- Cell Model Passports/GDSC: 14 prespecified cytotoxic drugs with duplicate-assay control, lineage fixed effects, within-lineage random effects and leave-one-lineage-out sensitivity analysis.

## Reproduce the workflow

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
# R 4.6.1, Bioconductor 3.23 and DESeq2 1.52.0 are used for IMvigor210 export.
.\scripts\run_all.ps1
```

See `R_DEPENDENCIES.md` for the exact R setup and version check.

`config/datasets.tsv` records every public URL and expected local path. Raw downloads and large processed matrices are intentionally not committed. The acquisition and preprocessing stages recreate them.

## Repository layout

- `config/`: dataset manifest, cohort decisions, signatures and the prespecified GDSC drug family.
- `scripts/acquisition/`: public-data downloads.
- `scripts/preprocessing/`: cohort-specific preparation and patient-level filtering.
- `scripts/analysis/`: cellular-origin, composition, clinical-endpoint, negative-control, external-pathway and lineage-aware GDSC analyses.
- `scripts/validation/`: bootstrap, cross-validation, alternative-signature and invariant checks.
- `scripts/visualization/`: final manuscript and supplementary figures.
- `results/`: machine-readable final estimates used in the manuscript.
- `figures/`: final PNG and PDF figures.
- `MANUSCRIPT_CODE_MAP.md`: figure-by-figure mapping.

## Statistical guardrails

- Pretreatment samples only for response prediction.
- One analysis record per patient.
- Cells are nested within patients; clinical inference is patient level.
- Continuous standardized `ITGB1` in primary models; no median split for clinical inference.
- Higher `ITGB1` is prespecified as the nonresponse AUC direction; AUCs are never flipped by cohort.
- IMvigor210 incremental models use 100 repeated stratified five-fold validations and missing-indicator sensitivity analyses.
- Minimum detectable effects are design-based approximations; observed post-hoc power is not reported.
- `ITGB1` and each control marker are excluded from their own signatures.
- Transcriptomic composition scores are proxies, not DNA-derived purity estimates.
- REML tau-squared with modified Hartung-Knapp inference; DerSimonian-Laird and leave-one-cohort-out estimates are sensitivity analyses.
- Duplicate GDSC assays are selected by sample coverage before association testing.
- GDSC lineage labels come from Cell Model Passports sample/tissue metadata; DepMap is not used.
- FDR is controlled within prespecified analysis families.
- Association is not interpreted as causal resistance or therapeutic utility.

## Main result

Bulk `ITGB1` is a composition-sensitive, multicellular abundance marker and does not show stable stand-alone clinical-response discrimination. Direct TCGA pathology estimates show small stromal and tumour-content associations, while seven CPTAC studies show moderate RNA-protein correspondence. Five-cohort simulation indicates that 80% power requires an OR magnitude of about 1.55 without heterogeneity and 1.68-1.90 under increasing heterogeneity. Oxaliplatin is the only GDSC drug with FDR-significant within-dataset cross-lineage support and is not interpreted as independent replication or causal resistance.
