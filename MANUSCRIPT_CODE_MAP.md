# Manuscript-to-code map

| Manuscript item | Analysis | Reproduction entry points |
| --- | --- | --- |
| Figure 1 | GSE195832 source attribution and expanded study framework | `scripts.analysis.analyze_single_cell`; `scripts.visualization.plot_all_figures_ordered` |
| Figure 2 | Six-cohort composition modelling, bootstrap, cross-validation and MCP-counter sensitivity | `scripts.analysis.analyze_response`; `scripts.validation.robustness_analysis`; `scripts.validation.alternative_mcpcounter`; `scripts.visualization.plot_all_figures_ordered` |
| Figure 3 | Five-cohort REML ICI response meta-analysis, prespecified-direction AUC, GSE135222 PFS and GSE190826 pCR | `scripts.analysis.analyze_response`; `scripts.analysis.analyze_imvigor210`; `scripts.analysis.analyze_external_ici`; `scripts.analysis.analyze_rectal_response`; `scripts.visualization.plot_all_figures_ordered` |
| Figure 4 | Continuous external pathway analysis before and after full-composition adjustment | `scripts.preprocessing.prepare_external_solid_tumor`; `scripts.analysis.analyze_external_solid_tumor`; `scripts.analysis.analyze_external_marker_controls`; `scripts.visualization.plot_integrated_extension` |
| Figure 5 | Global, lineage-fixed-effect and within-lineage GDSC analysis of 14 drugs | `scripts.analysis.analyze_gdsc_lineage_adjusted`; `scripts.visualization.plot_integrated_extension` |
| Supplementary Figure S1 | Detailed composition bootstrap, CV and alternative-method results | `scripts.validation.robustness_analysis`; `scripts.validation.alternative_mcpcounter`; `scripts.visualization.plot_all_figures_ordered` |
| Supplementary Figure S2 | Direction-preserving bootstrap AUC sensitivity | `scripts.validation.extended_auc_bootstrap`; `scripts.visualization.plot_all_figures_ordered` |
| Supplementary Figure S3 | Design-based minimum detectable effects and power curves | `scripts.analysis.analyze_power_precision` |
| Supplementary Figure S4 | Independent GSE103322 patient-level single-cell validation | `scripts.acquisition.download_gse103322_xena`; `scripts.analysis.analyze_gse103322_single_cell`; `scripts.visualization.plot_strengthening_analyses` |
| Supplementary Figure S5 | Direct TCGA pathology-composition validation | `scripts.acquisition.download_tcga_pathology_purity`; `scripts.analysis.analyze_tcga_pathology_purity`; `scripts.visualization.plot_strengthening_analyses` |
| Supplementary Figure S6 | IMvigor210 repeated cross-validation and missingness | `scripts.analysis.analyze_imvigor210_validation`; `scripts.visualization.plot_strengthening_analyses` |
| Supplementary Figure S7 | CPTAC ITGB1 RNA-protein validation | `scripts.acquisition.download_cptac_itgb1`; `scripts.analysis.analyze_cptac_itgb1`; `scripts.visualization.plot_strengthening_analyses` |

All figures have corresponding machine-readable estimates in `results/`. Legacy GEPIA/HPA portal panels, GSE96058 summary models, the previous signature-counting Fisher test and the discovery-only coexpression panel are not part of the revised manuscript or repository workflow.
