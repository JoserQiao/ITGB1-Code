"""Model anti-PD-1 treatment response across melanoma cohorts."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from scripts.core.cli import parse_project_root
from scripts.core.common import (
    build_composition_features,
    bh_fdr,
    fit_binary_response_feature,
    read_signatures,
)
from scripts.core.data import align_expression_metadata


def bootstrap_auc(y, score, n_boot=2000, seed=20260813):
    rng = np.random.default_rng(seed)
    vals = []
    y, score = np.asarray(y), np.asarray(score)
    for _ in range(n_boot):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) == 2:
            vals.append(roc_auc_score(y[idx], score[idx]))
    return np.quantile(vals, [0.025, 0.975]).tolist()


def fit_models(expr, meta, signatures, cohort):
    expr, meta = align_expression_metadata(expr, meta, "patient_id")
    frame, diagnostics, coefficients = build_composition_features(expr, signatures, cohort)
    frame.insert(0, "response", meta.response_binary.astype(int))
    rows = []
    features = ["ITGB1", "ITGB1_residual_CAF_ECM", "ITGB1_residual_Endothelial",
                "ITGB1_residual_StromalVascular", "ITGB1_residual_FullComposition",
                "CAF_core", "ECM", "Endothelial", "Cytotoxic", "IFNG_response"]
    for feature in features:
        sub = frame[["response", feature]].dropna()
        rows.append({
            "cohort": cohort,
            "feature": feature,
            **fit_binary_response_feature(frame, feature),
            "auc_ci": bootstrap_auc(sub.response, sub[feature]),
        })
    result = pd.DataFrame(rows)
    result["fdr_logistic"] = bh_fdr(result.p_logistic)
    result["fdr_mann_whitney"] = bh_fdr(result.p_mann_whitney)
    correlations = []
    for stromal in ["CAF_core", "ECM", "Endothelial", "T_cell", "Myeloid", "Epithelial", "iCAF", "myCAF"]:
        if stromal in frame:
            sub = frame[["ITGB1", stromal]].dropna()
            rho, p = spearmanr(sub.ITGB1, sub[stromal])
            correlations.append({"cohort": cohort, "stromal_feature": stromal,
                                 "n": len(sub), "spearman_rho": rho, "p_value": p})
    corr = pd.DataFrame(correlations)
    if not corr.empty:
        corr["fdr"] = bh_fdr(corr.p_value)
    return frame, result, corr, diagnostics, coefficients


def main():
    root, _ = parse_project_root()
    signatures = read_signatures(root / "config/signatures.tsv")
    all_results = []
    all_correlations = []
    all_diagnostics = []
    all_coefficients = []
    for cohort in ["GSE91061", "GSE78220"]:
        d = root / f"data/processed/{cohort}"
        expr = pd.read_csv(d / "expression_log2_fpkm.tsv.gz", sep="\t", index_col=0)
        meta = pd.read_csv(d / "metadata_v2.tsv", sep="\t")
        frame, result, corr, diagnostics, coefficients = fit_models(expr, meta, signatures, cohort)
        frame.to_csv(d / "analysis_features_v3.tsv", sep="\t", index_label="patient_id")
        all_results.append(result)
        all_correlations.append(corr)
        all_diagnostics.append(diagnostics)
        all_coefficients.append(coefficients)
    out = root / "results/treatment_response"
    out.mkdir(parents=True, exist_ok=True)
    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(out / "anti_pd1_model_results_v3.tsv", sep="\t", index=False)
    pd.concat(all_correlations, ignore_index=True).to_csv(
        out / "anti_pd1_itgb1_composition_correlations_v3.tsv", sep="\t", index=False)
    pd.concat(all_diagnostics, ignore_index=True).to_csv(
        out / "anti_pd1_itgb1_composition_explained_variance_v3.tsv", sep="\t", index=False)
    pd.concat(all_coefficients, ignore_index=True).to_csv(
        out / "anti_pd1_itgb1_composition_coefficients_v3.tsv", sep="\t", index=False)
    summary = {"cohorts": combined.groupby("cohort").n.first().to_dict(),
               "note": "GSE91061 discovery; GSE78220 independent validation; pretreatment samples only."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
