"""Shared statistical and signature-scoring utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, rankdata
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.outliers_influence import variance_inflation_factor


def bh_fdr(pvalues):
    values = np.asarray(pvalues, dtype=float)
    out = np.full(values.shape, np.nan)
    keep = np.isfinite(values)
    if keep.any():
        out[keep] = multipletests(values[keep], method="fdr_bh")[1]
    return out


def signature_score(expression: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Rank-based sample score; expression must be genes x samples."""
    index_upper = {str(g).upper(): g for g in expression.index}
    selected = [index_upper[g.upper()] for g in genes if g.upper() in index_upper]
    if len(selected) < 2:
        return pd.Series(np.nan, index=expression.columns)
    ranked = pd.DataFrame(
        np.apply_along_axis(rankdata, 0, expression.to_numpy(dtype=float)) / expression.shape[0],
        index=expression.index,
        columns=expression.columns,
    )
    return ranked.loc[selected].mean(axis=0)


def read_signatures(path) -> dict[str, list[str]]:
    tab = pd.read_csv(path, sep="\t")
    return {row.signature: [x.strip() for x in row.genes.split(",") if x.strip()]
            for row in tab.itertuples(index=False)}


COMPOSITION_MODELS = {
    "CAF_ECM": ["CAF_core", "ECM"],
    "Endothelial": ["Endothelial"],
    "StromalVascular": ["CAF_core", "ECM", "Endothelial"],
    "FullComposition": ["CAF_core", "ECM", "Endothelial", "T_cell", "Myeloid", "Epithelial"],
}

ITGB1_FEATURES = ["ITGB1", *(f"ITGB1_residual_{name}" for name in COMPOSITION_MODELS)]


def build_composition_features(
    expression: pd.DataFrame,
    signatures: dict[str, list[str]],
    cohort: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Score signatures and add composition-adjusted ITGB1 residuals."""
    frame = pd.DataFrame({"ITGB1": expression.loc["ITGB1"]})
    for name, genes in signatures.items():
        frame[name] = signature_score(expression, genes)
    return add_itgb1_composition_residuals(frame, cohort)


def fit_binary_response_feature(frame: pd.DataFrame, feature: str) -> dict[str, float | int]:
    """Fit one standardized logistic model and its nonparametric summaries."""
    data = frame[["response", feature]].dropna()
    standardized = (data[feature] - data[feature].mean()) / data[feature].std(ddof=0)
    model = sm.GLM(
        data.response,
        sm.add_constant(standardized),
        family=sm.families.Binomial(),
    ).fit()
    confidence_interval = model.conf_int().loc[feature]
    auc = roc_auc_score(data.response, data[feature])
    p_mann_whitney = mannwhitneyu(
        data.loc[data.response.eq(1), feature],
        data.loc[data.response.eq(0), feature],
        alternative="two-sided",
    ).pvalue
    return {
        "n": len(data),
        "responders": int(data.response.sum()),
        "odds_ratio_per_sd": float(np.exp(model.params[feature])),
        "ci_low": float(np.exp(confidence_interval.iloc[0])),
        "ci_high": float(np.exp(confidence_interval.iloc[1])),
        "p_logistic": float(model.pvalues[feature]),
        "p_mann_whitney": float(p_mann_whitney),
        "auc_raw_direction": float(auc),
        "auc_discrimination": float(max(auc, 1 - auc)),
    }


def add_itgb1_composition_residuals(frame: pd.DataFrame, cohort: str):
    """Add nested composition-adjusted ITGB1 residuals and return diagnostics."""
    diagnostics, coefficients = [], []
    for model_name, predictors in COMPOSITION_MODELS.items():
        complete = frame[["ITGB1", *predictors]].dropna()
        x = sm.add_constant(complete[predictors])
        fit = sm.OLS(complete.ITGB1, x).fit()
        column = f"ITGB1_residual_{model_name}"
        frame[column] = np.nan
        frame.loc[complete.index, column] = fit.resid
        diagnostics.append({"cohort": cohort, "composition_model": model_name,
                            "n": len(complete), "predictors": ",".join(predictors),
                            "r_squared": fit.rsquared, "adjusted_r_squared": fit.rsquared_adj,
                            "rmse": float(np.sqrt(np.mean(fit.resid ** 2)))})
        # Center/scale before VIF: rank-based signature scores share a non-zero
        # baseline, which otherwise produces a spurious no-intercept VIF.
        design = complete[predictors].apply(
            lambda s: (s - s.mean()) / s.std(ddof=0)).replace([np.inf, -np.inf], np.nan)
        for predictor in predictors:
            vif = (variance_inflation_factor(design.to_numpy(float), design.columns.get_loc(predictor))
                   if len(predictors) > 1 else 1.0)
            coefficients.append({"cohort": cohort, "composition_model": model_name,
                                 "predictor": predictor, "coefficient": fit.params[predictor],
                                 "p_value": fit.pvalues[predictor], "vif": vif})
    diagnostic_table = pd.DataFrame(diagnostics)
    r2 = diagnostic_table.set_index("composition_model").r_squared
    diagnostic_table["incremental_r_squared"] = np.nan
    diagnostic_table.loc[diagnostic_table.composition_model.eq("StromalVascular"), "incremental_r_squared"] = (
        r2["StromalVascular"] - r2["CAF_ECM"])
    diagnostic_table.loc[diagnostic_table.composition_model.eq("FullComposition"), "incremental_r_squared"] = (
        r2["FullComposition"] - r2["StromalVascular"])
    return frame, diagnostic_table, pd.DataFrame(coefficients)
