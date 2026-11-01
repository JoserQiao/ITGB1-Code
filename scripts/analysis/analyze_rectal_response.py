"""Model chemoradiotherapy response in the rectal cancer cohort."""

from __future__ import annotations

import pandas as pd

from scripts.core.cli import parse_project_root
from scripts.core.common import build_composition_features, bh_fdr, fit_binary_response_feature, read_signatures
from scripts.core.data import align_expression_metadata, load_cohort_expression


def main():
    root, _ = parse_project_root()
    d = root / "data/processed/GSE190826"
    expr = load_cohort_expression(root, "GSE190826")
    meta = pd.read_csv(d / "metadata_v2.tsv", sep="\t")
    signatures = read_signatures(root / "config/signatures.tsv")
    expr, meta = align_expression_metadata(expr, meta, "gsm")
    frame, diagnostics, coefficients = build_composition_features(expr, signatures, "GSE190826")
    frame.insert(0, "response", meta.response_binary.astype(int))
    rows = []
    features = ["ITGB1", "ITGB1_residual_CAF_ECM", "ITGB1_residual_Endothelial",
                "ITGB1_residual_StromalVascular", "ITGB1_residual_FullComposition",
                "CAF_core", "ECM", "Endothelial", "iCAF", "myCAF"]
    for feature in features:
        rows.append({"feature": feature, **fit_binary_response_feature(frame, feature)})
    out = root / "results/clinical_response"
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(d / "analysis_features_v3.tsv", sep="\t", index_label="gsm")
    result = pd.DataFrame(rows)
    result["fdr_logistic"] = bh_fdr(result.p_logistic)
    result["fdr_mann_whitney"] = bh_fdr(result.p_mann_whitney)
    result.to_csv(out / "rectal_chemoradiotherapy_results_v3.tsv", sep="\t", index=False)
    diagnostics.to_csv(out / "rectal_itgb1_composition_explained_variance_v3.tsv", sep="\t", index=False)
    coefficients.to_csv(out / "rectal_itgb1_composition_coefficients_v3.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
