"""Validate data and result invariants for the integrated extension analyses."""

from __future__ import annotations

import numpy as np
import pandas as pd

from scripts.core.cli import parse_project_root
from scripts.core.common import COMPOSITION_MODELS, read_signatures


EXPECTED_SAMPLES = {"GSE172356": 62, "GSE113255": 130, "GSE76427": 115}


def record(rows: list[dict], check: str, passed: bool, detail: str) -> None:
    rows.append({"check": check, "passed": bool(passed), "detail": detail})
    if not passed:
        raise AssertionError(f"{check}: {detail}")


def main() -> None:
    root, _ = parse_project_root()
    rows: list[dict] = []
    signatures = read_signatures(root / "config/signatures.tsv")
    signature_genes = {gene.upper() for genes in signatures.values() for gene in genes}
    record(rows, "ITGB1 excluded from composition signatures", "ITGB1" not in signature_genes,
           f"ITGB1_present={('ITGB1' in signature_genes)}")

    for cohort, expected_n in EXPECTED_SAMPLES.items():
        expression = pd.read_csv(
            root / f"data/processed/{cohort}/expression_log2.tsv.gz", sep="\t", index_col=0
        )
        metadata = pd.read_csv(root / f"data/processed/{cohort}/metadata.tsv", sep="\t")
        features = pd.read_csv(
            root / f"data/processed/{cohort}/analysis_features.tsv", sep="\t", index_col=0
        )
        record(rows, f"{cohort} sample count", expression.shape[1] == expected_n,
               f"observed={expression.shape[1]}; expected={expected_n}")
        record(rows, f"{cohort} sample identifiers unique",
               not expression.columns.duplicated().any() and not metadata.sample_id.duplicated().any(),
               f"expression_duplicates={expression.columns.duplicated().sum()}; metadata_duplicates={metadata.sample_id.duplicated().sum()}")
        record(rows, f"{cohort} complete expression matrix", not expression.isna().any().any(),
               f"missing_values={int(expression.isna().sum().sum())}")
        record(rows, f"{cohort} ITGB1 present and variable",
               "ITGB1" in expression.index and expression.loc["ITGB1"].std(ddof=0) > 0,
               f"ITGB1_sd={expression.loc['ITGB1'].std(ddof=0):.6g}")
        predictors = COMPOSITION_MODELS["FullComposition"]
        residual = features["ITGB1_residual_FullComposition"]
        correlations = [abs(residual.corr(features[predictor])) for predictor in predictors]
        record(rows, f"{cohort} full-composition residual orthogonality", max(correlations) < 1e-8,
               f"max_abs_pearson={max(correlations):.3g}")

    gdsc = pd.read_csv(
        root / "results/pharmacogenomics/gdsc_ITGB1_lineage_adjusted_summary.tsv", sep="\t"
    )
    record(rows, "GDSC one primary assay per drug", not gdsc.drug_name.duplicated().any(),
           f"drugs={len(gdsc)}; duplicate_names={gdsc.drug_name.duplicated().sum()}")
    record(rows, "GDSC minimum total sample size", gdsc.n_total.ge(100).all(),
           f"minimum_n={gdsc.n_total.min()}")
    finite_columns = ["global_beta_per_sd", "lineage_adjusted_beta_per_sd", "lineage_adjusted_fdr"]
    record(rows, "GDSC estimates finite", np.isfinite(gdsc[finite_columns].to_numpy(float)).all(),
           f"nonfinite={int((~np.isfinite(gdsc[finite_columns].to_numpy(float))).sum())}")
    fdr_columns = [column for column in gdsc if column.endswith("fdr")]
    valid_fdr = gdsc[fdr_columns].apply(lambda x: x.between(0, 1) | x.isna()).all().all()
    record(rows, "GDSC FDR range", valid_fdr, f"columns={','.join(fdr_columns)}")

    clinical = pd.read_csv(
        root / "results/clinical_response/five_ici_cohort_response_results.tsv", sep="\t"
    )
    record(rows, "Clinical AUC uses prespecified resistance direction",
           "auc_discrimination" not in clinical.columns
           and "auc_resistance_direction" in clinical.columns
           and np.allclose(clinical.auc_response_direction + clinical.auc_resistance_direction, 1),
           "response and resistance AUCs must sum to one; no direction-free maximum is allowed")
    meta = pd.read_csv(
        root / "results/clinical_response/random_effects_meta_analysis.tsv", sep="\t"
    )
    record(rows, "Primary clinical meta-analysis uses REML",
           meta.method.str.startswith("REML").all(),
           f"methods={','.join(meta.method.unique())}")

    gse103322 = pd.read_csv(
        root / "results/single_cell/GSE103322_patient_celltype_ITGB1.tsv", sep="\t"
    )
    record(rows, "GSE103322 independent single-cell patients",
           gse103322.patient_id.nunique() == 18,
           f"patients={gse103322.patient_id.nunique()}")
    gse103322_tests = pd.read_csv(
        root / "results/single_cell/GSE103322_paired_celltype_comparisons.tsv", sep="\t"
    )
    record(rows, "GSE103322 prespecified comparison FDR range",
           len(gse103322_tests) == 2 and gse103322_tests.fdr.between(0, 1).all(),
           f"comparisons={len(gse103322_tests)}")
    pathology = pd.read_csv(
        root / "results/pathology_purity/TCGA_ITGB1_pathology_meta.tsv", sep="\t"
    )
    record(rows, "TCGA pathology validation covers ten cancers",
           pathology.k_lineages.eq(10).all(),
           f"k={','.join(pathology.k_lineages.astype(str))}")
    protein = pd.read_csv(
        root / "results/protein_validation/CPTAC_ITGB1_RNA_protein_correlations.tsv", sep="\t"
    )
    record(rows, "CPTAC RNA-protein validation covers seven studies",
           len(protein) == 7 and protein.n.sum() == 728,
           f"studies={len(protein)}; paired_samples={protein.n.sum()}")
    cptac_pairs = pd.read_csv(
        root / "data/processed/CPTAC_ITGB1/CPTAC_ITGB1_RNA_protein.tsv", sep="\t"
    )
    record(rows, "CPTAC one paired sample per patient",
           not cptac_pairs.duplicated(["study", "patientId"]).any(),
           f"duplicate_study_patients={cptac_pairs.duplicated(['study', 'patientId']).sum()}")
    record(rows, "CPTAC study-level FDR range", protein.fdr.between(0, 1).all(),
           f"studies={len(protein)}")
    imputed = pd.read_csv(
        root / "results/clinical_response/IMvigor210_multiple_imputation_pooled.tsv", sep="\t"
    )
    record(rows, "IMvigor210 multiple-imputation result grid",
           len(imputed) == 3 and imputed.imputations.eq(50).all()
           and imputed[["odds_ratio", "ci_low", "ci_high", "p_value"]].notna().all().all(),
           f"models={len(imputed)}; imputations={','.join(imputed.imputations.astype(str).unique())}")

    external = pd.read_csv(
        root / "results/external_solid_tumor/external_pathway_associations.tsv", sep="\t"
    )
    record(rows, "External pathway result grid complete", len(external) == 30,
           f"observed_rows={len(external)}; expected_rows=30")
    record(rows, "External pathway FDR range",
           external[["fdr", "spearman_fdr"]].apply(lambda x: x.between(0, 1)).all().all(),
           "all adjusted P values must be within [0,1]")

    out = root / "results/validation"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "integrated_extension_checks.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
