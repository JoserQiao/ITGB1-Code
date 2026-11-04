"""Analyze ITGB1 in the IMvigor210 atezolizumab cohort.

The script expects files exported from the public IMvigor210CoreBiologies
package by scripts/preprocessing/export_imvigor210.R. Primary models use
continuous one-SD ITGB1 or full-composition residuals. Clinical incremental
models adjust for ECOG, PD-L1 immune-cell category, sex and metastatic status;
TMB is added in a complete-case sensitivity model.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2
from sklearn.metrics import brier_score_loss, roc_auc_score

from scripts.core.cli import parse_project_root
from scripts.core.common import ITGB1_FEATURES, add_itgb1_composition_residuals, fit_binary_response_feature
from scripts.analysis.analyze_external_ici import (
    leave_one_cohort_out_meta,
    meta_sensitivity,
    random_effects_meta,
)


def model_metrics(fit, data):
    prediction = fit.predict(data)
    return {
        "n": len(data),
        "responders": int(data.response.sum()),
        "auc": roc_auc_score(data.response, prediction),
        "brier": brier_score_loss(data.response, prediction),
        "log_likelihood": fit.llf,
        "parameters": int(fit.df_model + 1),
    }


def bootstrap_auc_difference(data, baseline_formula, extended_formula, bootstrap=2000, seed=20261101):
    rng = np.random.default_rng(seed)
    values = []
    indices = np.arange(len(data))
    for _ in range(bootstrap):
        draw = rng.choice(indices, size=len(indices), replace=True)
        sample = data.iloc[draw].copy()
        if sample.response.nunique() < 2:
            continue
        try:
            base = smf.glm(baseline_formula, data=sample, family=sm.families.Binomial()).fit()
            ext = smf.glm(extended_formula, data=sample, family=sm.families.Binomial()).fit()
            values.append(
                roc_auc_score(sample.response, ext.predict(sample))
                - roc_auc_score(sample.response, base.predict(sample))
            )
        except Exception:
            continue
    arr = np.asarray(values)
    return float(np.mean(arr)), float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975)), len(arr)


def main():
    root, _ = parse_project_root()
    processed = root / "data/processed/IMvigor210"
    features = pd.read_csv(processed / "analysis_features.tsv", sep="\t").set_index("sample_id")
    metadata = pd.read_csv(processed / "metadata.tsv", sep="\t").set_index("sample_id")
    features, diagnostics, coefficients = add_itgb1_composition_residuals(features, "IMvigor210")
    metadata["response"] = metadata.binaryResponse.map({"CR/PR": 1, "SD/PD": 0})
    frame = features.join(metadata, how="inner")
    response_frame = frame[frame.response.notna()].copy()

    rows = []
    for feature in ITGB1_FEATURES:
        result = fit_binary_response_feature(response_frame.rename(columns={"response": "response"}), feature)
        rows.append({"cohort": "IMvigor210", "feature": feature, **result})
    out = root / "results/clinical_response"
    out.mkdir(parents=True, exist_ok=True)
    response_results = pd.DataFrame(rows)
    response_results.to_csv(out / "IMvigor210_response_results.tsv", sep="\t", index=False)
    four_path = out / "four_ici_cohort_response_results.tsv"
    if four_path.is_file():
        four = pd.read_csv(four_path, sep="\t")
        five = pd.concat([four, response_results], ignore_index=True)
        five.to_csv(out / "five_ici_cohort_response_results.tsv", sep="\t", index=False, na_rep="NA")
        random_effects_meta(five, "REML").to_csv(
            out / "random_effects_meta_analysis.tsv", sep="\t", index=False)
        meta_sensitivity(five).to_csv(
            out / "random_effects_meta_sensitivity.tsv", sep="\t", index=False)
        leave_one_cohort_out_meta(five).to_csv(
            out / "random_effects_meta_leave_one_cohort_out.tsv", sep="\t", index=False)
    diagnostics.to_csv(out / "IMvigor210_composition_diagnostics.tsv", sep="\t", index=False, na_rep="NA")
    coefficients.to_csv(out / "IMvigor210_composition_coefficients.tsv", sep="\t", index=False)

    response_frame["ITGB1_z"] = (response_frame.ITGB1 - response_frame.ITGB1.mean()) / response_frame.ITGB1.std(ddof=0)
    response_frame["ITGB1_residual_z"] = (
        response_frame.ITGB1_residual_FullComposition
        - response_frame.ITGB1_residual_FullComposition.mean()
    ) / response_frame.ITGB1_residual_FullComposition.std(ddof=0)
    response_frame["log1p_TMB"] = np.log1p(pd.to_numeric(response_frame["FMOne.mutation.burden.per.MB"], errors="coerce"))
    response_frame = response_frame.rename(columns={
        "Baseline.ECOG.Score": "ECOG", "IC.Level": "IC_Level",
        "Met.Disease.Status": "Met_Status",
    })
    baseline = "response ~ C(ECOG) + C(IC_Level) + C(Sex) + C(Met_Status)"
    model_specs = [
        ("clinical_baseline", baseline),
        ("clinical_plus_ITGB1", baseline + " + ITGB1_z"),
        ("clinical_plus_ITGB1_residual", baseline + " + ITGB1_residual_z"),
        ("clinical_TMB_baseline", baseline + " + log1p_TMB"),
        ("clinical_TMB_plus_ITGB1", baseline + " + log1p_TMB + ITGB1_z"),
    ]
    model_rows, fits, data_by_model = [], {}, {}
    for name, formula in model_specs:
        columns = ["response", "ECOG", "IC_Level", "Sex", "Met_Status"]
        if "ITGB1_z" in formula: columns.append("ITGB1_z")
        if "ITGB1_residual_z" in formula: columns.append("ITGB1_residual_z")
        if "log1p_TMB" in formula: columns.append("log1p_TMB")
        data = response_frame[columns].dropna().copy()
        fit = smf.glm(formula, data=data, family=sm.families.Binomial()).fit()
        fits[name], data_by_model[name] = fit, data
        row = {"model": name, "formula": formula, **model_metrics(fit, data)}
        for term in ["ITGB1_z", "ITGB1_residual_z"]:
            if term in fit.params:
                ci = fit.conf_int().loc[term]
                row.update({
                    "itgb1_term": term,
                    "itgb1_or": np.exp(fit.params[term]),
                    "itgb1_ci_low": np.exp(ci.iloc[0]),
                    "itgb1_ci_high": np.exp(ci.iloc[1]),
                    "itgb1_p_value": fit.pvalues[term],
                })
        model_rows.append(row)

    comparisons = []
    for baseline_name, extended_name in [
        ("clinical_baseline", "clinical_plus_ITGB1"),
        ("clinical_baseline", "clinical_plus_ITGB1_residual"),
        ("clinical_TMB_baseline", "clinical_TMB_plus_ITGB1"),
    ]:
        # Refit both models on the extended model's complete-case set.
        data = data_by_model[extended_name]
        base_formula = next(x[1] for x in model_specs if x[0] == baseline_name)
        ext_formula = next(x[1] for x in model_specs if x[0] == extended_name)
        base = smf.glm(base_formula, data=data, family=sm.families.Binomial()).fit()
        ext = smf.glm(ext_formula, data=data, family=sm.families.Binomial()).fit()
        lr = 2 * (ext.llf - base.llf)
        delta_auc = roc_auc_score(data.response, ext.predict(data)) - roc_auc_score(data.response, base.predict(data))
        boot_mean, boot_low, boot_high, successful = bootstrap_auc_difference(data, base_formula, ext_formula)
        comparisons.append({
            "baseline": baseline_name, "extended": extended_name,
            "n": len(data), "lr_statistic": lr, "lr_p_value": chi2.sf(lr, 1),
            "delta_auc": delta_auc, "bootstrap_delta_auc_mean": boot_mean,
            "bootstrap_ci_low": boot_low, "bootstrap_ci_high": boot_high,
            "successful_bootstraps": successful,
        })
    pd.DataFrame(model_rows).to_csv(out / "IMvigor210_clinical_models.tsv", sep="\t", index=False, na_rep="NA")
    pd.DataFrame(comparisons).to_csv(out / "IMvigor210_incremental_value.tsv", sep="\t", index=False)
    features.to_csv(processed / "analysis_features_with_residuals.tsv", sep="\t", index_label="sample_id")


if __name__ == "__main__":
    main()
