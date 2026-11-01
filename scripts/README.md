# Pipeline order

Run the complete workflow with `scripts/run_all.ps1`.

1. `scripts.acquisition.download_data`
2. `scripts.acquisition.download_cmp`
3. `scripts.acquisition.download_external_solid_tumor`
4. `scripts.acquisition.download_imvigor210`
5. `scripts.preprocessing.prepare_geo`
6. `scripts.preprocessing.prepare_external_ici`
7. `scripts.preprocessing.prepare_external_solid_tumor`
8. `scripts.preprocessing.export_imvigor210.R`
9. `scripts.analysis.analyze_response`
10. `scripts.analysis.analyze_rectal_response`
11. `scripts.analysis.analyze_single_cell`
12. `scripts.analysis.analyze_imvigor210`
13. `scripts.analysis.analyze_external_ici`
14. `scripts.analysis.analyze_power_precision`
15. `scripts.analysis.analyze_gdsc_lineage_adjusted`
16. `scripts.analysis.analyze_external_solid_tumor`
17. `scripts.analysis.analyze_external_marker_controls`
18. `scripts.validation.robustness_analysis`
19. `scripts.validation.extended_auc_bootstrap`
20. `scripts.validation.alternative_mcpcounter`
21. `scripts.validation.validate_integrated_extension`
22. `scripts.visualization.plot_all_figures_ordered`
23. `scripts.visualization.plot_integrated_extension`

Each command accepts `--project-root`. Analysis outputs are written to `results/`; figures are written to `figures/`.
