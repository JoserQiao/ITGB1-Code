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
