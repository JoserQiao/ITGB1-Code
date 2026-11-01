# ITGB1 cellular-origin and composition-aware analysis

This repository reproduces the manuscript:

> Cellular origins, stromal confounding and context-dependent drug-sensitivity associations of bulk ITGB1 across cancers

## Questions addressed

1. Which tumour-microenvironment compartments contribute to bulk `ITGB1` expression?
2. How much bulk `ITGB1` variation is associated with transcriptomic composition?
3. Does `ITGB1` predict objective response, progression-free survival, or pathological complete response after composition adjustment?
4. Which external pathway and GDSC associations remain after composition or cancer-lineage adjustment?

## Evidence layers

- GSE195832: paired single-cell source attribution.
- GSE91061, GSE78220, GSE176307 and GSE126044: pretreatment ICI objective response.
- GSE135222: ICI progression-free survival.
- GSE190826: neoadjuvant chemoradiotherapy pathological complete response.
- GSE172356, GSE113255 and GSE76427: continuous composition-aware pathway validation.
- Cell Model Passports/GDSC: 14 prespecified cytotoxic drugs with duplicate-assay control, lineage fixed effects, within-lineage random effects and leave-one-lineage-out sensitivity analysis.

## Reproduce the workflow

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.\scripts\run_all.ps1
```

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
- `ITGB1` and each control marker are excluded from their own signatures.
- Transcriptomic composition scores are proxies, not DNA-derived purity estimates.
- DerSimonian-Laird tau-squared with modified Hartung-Knapp confidence intervals for small-k meta-analysis.
- Duplicate GDSC assays are selected by sample coverage before association testing.
- FDR is controlled within prespecified analysis families.
- Association is not interpreted as causal resistance or therapeutic utility.

## Main result

Bulk `ITGB1` is composition sensitive and does not show stable stand-alone clinical-response discrimination. Composition and cancer-lineage adjustment attenuate pathway and drug-sensitivity associations, supporting a context-dependent rather than universal biomarker interpretation.
