"""Build the machine-readable source tables used in Supplementary Tables S1-S10."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.core.cli import parse_project_root


ITGB1_FEATURES = [
    "ITGB1",
    "ITGB1_residual_CAF_ECM",
    "ITGB1_residual_Endothelial",
    "ITGB1_residual_StromalVascular",
    "ITGB1_residual_FullComposition",
]

CPTAC_LABELS = {
    "paad_cptac_2021": "Pancreatic cancer (2021)",
    "luad_cptac_2020": "Lung adenocarcinoma (2020)",
    "lusc_cptac_2021": "Lung squamous carcinoma (2021)",
    "coad_cptac_2019": "Colon cancer (2019)",
    "ucec_cptac_2020": "Endometrial cancer (2020)",
    "gbm_cptac_2021": "Glioblastoma (2021)",
    "brca_cptac_2020": "Breast cancer (2020)",
}


def read(root: Path, relative: str) -> pd.DataFrame:
    return pd.read_csv(root / relative, sep="\t")


def write(frame: pd.DataFrame, output: Path, name: str) -> None:
    frame.to_csv(output / name, sep="\t", index=False, na_rep="NA")


def table_s1(root: Path) -> pd.DataFrame:
    five = read(root, "results/clinical_response/five_ici_cohort_response_results.tsv")
    five = five.loc[five.feature.eq("ITGB1")].set_index("cohort")
    pfs = read(root, "results/clinical_response/GSE135222_pfs_results.tsv").iloc[0]
    pcr = read(root, "results/clinical_response/rectal_chemoradiotherapy_results_v3.tsv").iloc[0]
    sc = read(root, "results/single_cell/GSE103322_celltype_summary.tsv")
    external = read(root, "results/external_solid_tumor/preprocessing_qc.tsv").set_index("cohort")
    pathology = read(root, "results/pathology_purity/TCGA_ITGB1_pathology_correlations.tsv")
    pathology_n = pathology.groupby("variable").n.sum()
    protein = read(root, "results/protein_validation/CPTAC_ITGB1_RNA_protein_correlations.tsv")
    gdsc = read(root, "results/pharmacogenomics/gdsc_ITGB1_lineage_adjusted_summary.tsv")

    def response_summary(cohort: str) -> str:
        row = five.loc[cohort]
        return f"{int(row.n)}; {int(row.responders)} responders"

    rows = [
        ["GSE195832", "HNSCC", "4 paired patients", "Single-cell source attribution", "10x scRNA-seq; paired pre/post nivolumab"],
        ["GSE103322", "HNSCC", f"{int(sc.n_cells.sum()):,} cells; 18 patients", "Independent single-cell validation", "Original cell-type and patient annotations; Xena ITGB1 query"],
        ["GSE91061", "Melanoma", response_summary("GSE91061"), "ICI response", "Pretreatment FPKM"],
        ["GSE78220", "Melanoma", response_summary("GSE78220"), "ICI response", "Pretreatment FPKM"],
        ["GSE176307", "Urothelial cancer", response_summary("GSE176307"), "ICI response", "Pretreatment TPM"],
        ["GSE126044", "NSCLC", response_summary("GSE126044"), "ICI response", "Pretreatment counts"],
        ["IMvigor210", "Urothelial cancer", response_summary("IMvigor210"), "ICI response and incremental value", "Pretreatment DESeq2-normalized RNA-seq; atezolizumab"],
        ["GSE135222", "NSCLC", f"{int(pfs.n)}; {int(pfs.events)} PFS events", "PFS", "Pretreatment RNA-seq"],
        ["GSE190826", "Rectal cancer", f"{int(pcr.n)}; {int(pcr.responders)} pCR", "Chemoradiotherapy pCR", "Pretreatment bulk RNA"],
        ["GSE172356", "PDAC", str(int(external.loc["GSE172356", "n_samples"])), "External pathway validation", "HiSeq X Ten RNA-seq"],
        ["GSE113255", "Gastric cancer", f"{int(external.loc['GSE113255', 'n_samples'])} tumours", "External pathway validation", "NextSeq 500 RNA-seq"],
        ["GSE76427", "HCC", f"{int(external.loc['GSE76427', 'n_samples'])} tumours", "External pathway validation", "GPL10558 expression array; normalized series matrix; no TPM conversion"],
        ["TCGA/GDC pathology", "10 cancers", f"{int(pathology_n.min()):,}-{int(pathology_n.max()):,} samples", "Direct composition validation", "Pathologist-estimated stromal, tumour-cell and tumour-nuclei percentages"],
        ["CPTAC/cBioPortal", "7 cancers", f"{int(protein.n.sum())} paired primary-tumour samples", "RNA-protein validation", "Mass-spectrometry protein and matched RNA abundance"],
        ["Cell Model Passports/GDSC", "Cancer cell lines", f"{int(gdsc.n_total.min()):,}-{int(gdsc.n_total.max()):,} per drug", "Lineage-aware ln(IC50)", "14 prespecified cytotoxic drugs"],
    ]
    return pd.DataFrame(rows, columns=["Dataset", "Cancer", "Analysis n", "Role", "Processing/endpoint"])


def table_s3(root: Path) -> pd.DataFrame:
    rows = []
    meta = read(root, "results/clinical_response/random_effects_meta_analysis.tsv")
    for row in meta.itertuples():
        rows.append([row.feature, row.k, "OR", row.pooled_or, row.ci_low, row.ci_high, row.p_value, "ICI objective response"])
    pfs = read(root, "results/clinical_response/GSE135222_pfs_results.tsv")
    for row in pfs.itertuples():
        rows.append([row.feature, 1, "HR", row.hazard_ratio_per_sd, row.ci_low, row.ci_high, row.p_value, "GSE135222 PFS"])
    pcr = read(root, "results/clinical_response/rectal_chemoradiotherapy_results_v3.tsv")
    for row in pcr.loc[pcr.feature.isin(ITGB1_FEATURES)].itertuples():
        rows.append([row.feature, 1, "OR", row.odds_ratio_per_sd, row.ci_low, row.ci_high, row.p_logistic, "GSE190826 pCR"])
    return pd.DataFrame(rows, columns=["Feature", "k", "Effect measure", "Effect per SD", "CI low", "CI high", "P", "Endpoint"])


def table_s4(root: Path) -> pd.DataFrame:
    rows = []
    primary = read(root, "results/external_solid_tumor/external_pathway_random_effects_meta.tsv")
    for row in primary.itertuples():
        rows.append(["Primary pathway family", row.feature, row.pathway.replace("HALLMARK_", ""), row.meta_rho, row.meta_ci_low, row.meta_ci_high, row.meta_fdr, row.i_squared])
    controls = read(root, "results/external_solid_tumor/external_marker_control_meta.tsv")
    for row in controls.itertuples():
        rows.append(["Marker-control family", f"{row.marker} residual", row.pathway, row.meta_rho, row.meta_ci_low, row.meta_ci_high, row.meta_fdr, row.i_squared])
    return pd.DataFrame(rows, columns=["Multiplicity family", "Feature/marker", "Pathway", "Pooled rho", "CI low", "CI high", "FDR", "I2 (%)"])


def table_s5(root: Path) -> pd.DataFrame:
    gdsc = read(root, "results/pharmacogenomics/gdsc_ITGB1_lineage_adjusted_summary.tsv")
    return gdsc[[
        "drug_name", "n_total", "lineage_adjusted_beta_per_sd", "lineage_adjusted_ci_low",
        "lineage_adjusted_ci_high", "lineage_adjusted_fdr", "meta_rho", "meta_fdr",
        "attenuation_fraction",
    ]].rename(columns={
        "drug_name": "Drug", "n_total": "n", "lineage_adjusted_beta_per_sd": "Lineage-adjusted beta",
        "lineage_adjusted_ci_low": "CI low", "lineage_adjusted_ci_high": "CI high",
        "lineage_adjusted_fdr": "FDR", "meta_rho": "Within-lineage rho",
        "meta_fdr": "Meta FDR", "attenuation_fraction": "Attenuation fraction",
    })


def table_s6(root: Path) -> pd.DataFrame:
    models = read(root, "results/clinical_response/IMvigor210_clinical_models.tsv")
    increments = read(root, "results/clinical_response/IMvigor210_incremental_value.tsv")
    return models.merge(
        increments[["extended", "delta_auc", "bootstrap_ci_low", "bootstrap_ci_high", "lr_p_value"]],
        left_on="model", right_on="extended", how="left",
    ).drop(columns="extended")


def table_s9(root: Path) -> pd.DataFrame:
    rows = []
    sensitivity = read(root, "results/clinical_response/random_effects_meta_sensitivity.tsv")
    for row in sensitivity.loc[sensitivity.feature.eq("ITGB1")].itertuples():
        rows.append(["Meta-analysis", row.feature, f"Meta {row.tau_method}", row.pooled_or, row.ci_low, row.ci_high, row.p_value, "Primary unadjusted ITGB1 only"])
    leave_one_out = read(root, "results/clinical_response/random_effects_meta_leave_one_cohort_out.tsv")
    for row in leave_one_out.loc[leave_one_out.feature.eq("ITGB1")].itertuples():
        rows.append(["Leave-one-cohort-out", row.feature, f"Omit {row.omitted_cohort}", row.pooled_or, row.ci_low, row.ci_high, row.p_value, "Primary unadjusted ITGB1 only"])
    simulation = read(root, "results/power_precision/five_cohort_meta_simulation_mde.tsv")
    for row in simulation.itertuples():
        rows.append(["Power simulation", "ITGB1", f"tau(log OR)={row.tau_log_or:.2f}", row.simulated_minimum_detectable_or_magnitude, pd.NA, pd.NA, pd.NA, f"{int(row.replicates_per_effect):,} replicates; 80% MDE OR"])
    return pd.DataFrame(rows, columns=["Analysis family", "Feature", "Analysis", "Estimate", "CI low", "CI high", "P", "Detail"])


def table_s10(root: Path) -> pd.DataFrame:
    rows = []
    comparisons = read(root, "results/single_cell/GSE103322_paired_celltype_comparisons.tsv")
    for row in comparisons.itertuples():
        rows.append(["GSE103322", f"{row.cell_type_1} vs {row.cell_type_2}", row.paired_patients, "Mean difference", row.mean_difference, pd.NA, pd.NA, row.fdr, "FDR"])
    pathology = read(root, "results/pathology_purity/TCGA_ITGB1_pathology_meta.tsv")
    for row in pathology.itertuples():
        rows.append(["TCGA pathology", row.variable, row.k_lineages, "Spearman rho", row.meta_rho, row.meta_ci_low, row.meta_ci_high, row.meta_fdr, "FDR"])
    protein = read(root, "results/protein_validation/CPTAC_ITGB1_RNA_protein_correlations.tsv")
    for row in protein.itertuples():
        rows.append(["CPTAC", CPTAC_LABELS.get(row.study, row.study), row.n, "RNA-protein Spearman rho", row.spearman_rho, pd.NA, pd.NA, row.fdr, "FDR"])
    meta = read(root, "results/protein_validation/CPTAC_ITGB1_RNA_protein_meta.tsv").iloc[0]
    rows.append(["CPTAC pooled", "7 studies", int(protein.n.sum()), "RNA-protein Spearman rho", meta.meta_rho, meta.meta_ci_low, meta.meta_ci_high, meta.meta_p_value, "P"])
    return pd.DataFrame(rows, columns=["Source", "Comparison/cohort", "n/k", "Effect measure", "Estimate", "CI low", "CI high", "P/FDR value", "P/FDR type"])


def main() -> None:
    root, _ = parse_project_root()
    output = root / "results/supplementary_tables"
    output.mkdir(parents=True, exist_ok=True)

    write(table_s1(root), output, "TableS1_datasets.tsv")
    write(read(root, "config/signatures.tsv"), output, "TableS2_signatures.tsv")
    write(table_s3(root), output, "TableS3_clinical_results.tsv")
    write(table_s4(root), output, "TableS4_external_meta.tsv")
    write(table_s5(root), output, "TableS5_GDSC.tsv")
    write(table_s6(root), output, "TableS6_IMvigor_apparent_models.tsv")

    mde = read(root, "results/power_precision/minimum_detectable_effects.tsv")
    mde = mde.loc[~mde.cohort_or_pool.isin(["Four-cohort total", "Five-cohort total"])]
    write(mde, output, "TableS7_individual_MDE.tsv")

    write(read(root, "results/clinical_response/IMvigor210_repeated_cv_difference_summary.tsv"), output, "TableS8A_cross_validation.tsv")
    missing = read(root, "results/clinical_response/IMvigor210_missingness_comparison.tsv")
    write(missing.loc[missing.variable.eq("response")], output, "TableS8B_missingness.tsv")
    write(read(root, "results/clinical_response/IMvigor210_multiple_imputation_pooled.tsv"), output, "TableS8C_multiple_imputation.tsv")
    write(table_s9(root), output, "TableS9_primary_ITGB1_sensitivity.tsv")
    write(table_s10(root), output, "TableS10_public_validation.tsv")

    source_rows = [
        ["Table S1", "TableS1_datasets.tsv", "multiple tracked result tables", "Counts are derived from tracked results; labels and roles are declared in this reporting script."],
        ["Table S2", "TableS2_signatures.tsv", "config/signatures.tsv", "Exact copy."],
        ["Table S3", "TableS3_clinical_results.tsv", "random_effects_meta_analysis.tsv; GSE135222_pfs_results.tsv; rectal_chemoradiotherapy_results_v3.tsv", "Restricted to ITGB1 and the four prespecified composition-adjusted ITGB1 features; composition-score-only pCR rows are outside this table's scope."],
        ["Table S4", "TableS4_external_meta.tsv", "external_pathway_random_effects_meta.tsv; external_marker_control_meta.tsv", "All rows retained; multiplicity family is explicit."],
        ["Table S5", "TableS5_GDSC.tsv", "gdsc_ITGB1_lineage_adjusted_summary.tsv", "All 14 prespecified drugs retained."],
        ["Table S6", "TableS6_IMvigor_apparent_models.tsv", "IMvigor210_clinical_models.tsv; IMvigor210_incremental_value.tsv", "All apparent models retained."],
        ["Table S7", "TableS7_individual_MDE.tsv", "minimum_detectable_effects.tsv", "Aggregate four- and five-cohort rows are excluded because individual-cohort approximations and meta-analytic simulations are reported separately."],
        ["Table S8A", "TableS8A_cross_validation.tsv", "IMvigor210_repeated_cv_difference_summary.tsv", "All prespecified comparisons retained."],
        ["Table S8B", "TableS8B_missingness.tsv", "IMvigor210_missingness_comparison.tsv", "Response-missingness rows retained."],
        ["Table S8C", "TableS8C_multiple_imputation.tsv", "IMvigor210_multiple_imputation_pooled.tsv", "All three prespecified models retained."],
        ["Table S9", "TableS9_primary_ITGB1_sensitivity.tsv", "random_effects_meta_sensitivity.tsv; random_effects_meta_leave_one_cohort_out.tsv; five_cohort_meta_simulation_mde.tsv", "Meta-method and leave-one-out rows are restricted to the primary unadjusted ITGB1 feature; adjusted-feature sensitivity outputs remain in their complete tracked result tables."],
        ["Table S10", "TableS10_public_validation.tsv", "GSE103322_paired_celltype_comparisons.tsv; TCGA_ITGB1_pathology_meta.tsv; CPTAC_ITGB1_RNA_protein_correlations.tsv; CPTAC_ITGB1_RNA_protein_meta.tsv", "All prespecified public-validation rows retained."],
    ]
    write(pd.DataFrame(source_rows, columns=["Supplementary item", "Generated output", "Tracked source files", "Selection/assembly rule"]), output, "supplementary_table_sources.tsv")


if __name__ == "__main__":
    main()
