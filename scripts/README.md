# Pipeline order

Run the complete workflow with `scripts/run_all.ps1`.

1. `scripts.acquisition.download_data`
2. `scripts.acquisition.download_cmp`
3. `scripts.acquisition.download_external_solid_tumor`
4. `scripts.acquisition.download_imvigor210`
5. `scripts.acquisition.download_gse103322_xena`
6. `scripts.acquisition.download_tcga_pathology_purity`
7. `scripts.acquisition.download_cptac_itgb1`
8. `scripts.preprocessing.prepare_geo`
9. `scripts.preprocessing.prepare_external_ici`
10. `scripts.preprocessing.prepare_external_solid_tumor`
11. `scripts.preprocessing.export_imvigor210.R`
12. `scripts.analysis.analyze_response`
13. `scripts.analysis.analyze_rectal_response`
14. `scripts.analysis.analyze_single_cell`
15. `scripts.analysis.analyze_imvigor210`
16. `scripts.analysis.analyze_imvigor210_validation`
17. `scripts.analysis.analyze_imvigor210_multiple_imputation`
18. `scripts.analysis.analyze_external_ici`
19. `scripts.analysis.analyze_power_precision`
20. `scripts.analysis.analyze_gse103322_single_cell`
21. `scripts.analysis.analyze_tcga_pathology_purity`
22. `scripts.analysis.analyze_cptac_itgb1`
23. `scripts.analysis.analyze_gdsc_lineage_adjusted`
24. `scripts.analysis.analyze_external_solid_tumor`
25. `scripts.analysis.analyze_external_marker_controls`
26. `scripts.reporting.build_supplementary_tables`
27. `scripts.validation.verify_data_provenance --strict`
28. `scripts.validation.robustness_analysis`
29. `scripts.validation.extended_auc_bootstrap`
30. `scripts.validation.alternative_mcpcounter`
31. `scripts.validation.validate_integrated_extension`
32. `scripts.visualization.plot_all_figures_ordered`
33. `scripts.visualization.plot_integrated_extension`
34. `scripts.visualization.plot_strengthening_analyses`

Each command accepts `--project-root`. Analysis outputs are written to `results/`; figures are written to `figures/`.
