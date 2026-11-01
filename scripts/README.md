# Pipeline order

Run the complete workflow with `scripts/run_all.ps1`.

1. `scripts.acquisition.download_data`
2. `scripts.acquisition.download_cmp`
3. `scripts.acquisition.download_external_solid_tumor`
4. `scripts.preprocessing.prepare_geo`
5. `scripts.preprocessing.prepare_external_ici`
6. `scripts.preprocessing.prepare_external_solid_tumor`
7. `scripts.analysis.analyze_response`
8. `scripts.analysis.analyze_rectal_response`
9. `scripts.analysis.analyze_single_cell`
10. `scripts.analysis.analyze_external_ici`
11. `scripts.analysis.analyze_gdsc_lineage_adjusted`
12. `scripts.analysis.analyze_external_solid_tumor`
13. `scripts.analysis.analyze_external_marker_controls`
14. `scripts.validation.robustness_analysis`
15. `scripts.validation.extended_auc_bootstrap`
16. `scripts.validation.alternative_mcpcounter`
17. `scripts.validation.validate_integrated_extension`
18. `scripts.visualization.plot_all_figures_ordered`
19. `scripts.visualization.plot_integrated_extension`

Each command accepts `--project-root`. Analysis outputs are written to `results/`; figures are written to `figures/`.
