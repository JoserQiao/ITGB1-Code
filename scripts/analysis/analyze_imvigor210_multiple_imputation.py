"""Multiple-imputation sensitivity models for missing IMvigor210 covariates.

Predictive mean matching is applied to numeric encodings of the two categorical
covariates with missing values and to log-transformed TMB. Imputed categorical
values are therefore always observed category codes. ITGB1, response and fully
observed clinical variables enter every imputation model. Estimates are pooled
with Rubin's rules and are a missing-at-random sensitivity analysis, not a
solution to possible outcome-related non-ignorable missingness.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm
from sklearn.metrics import roc_auc_score
from statsmodels.imputation.mice import MICEData

from scripts.analysis.analyze_imvigor210_validation import load_frame
from scripts.core.cli import parse_project_root, project_root_parser


MODELS = [
    ("clinical_plus_ITGB1", "response ~ C(ECOG) + C(IC_code) + C(Sex_code) + C(Met_code) + ITGB1_z", "ITGB1_z"),
    ("clinical_plus_residual", "response ~ C(ECOG) + C(IC_code) + C(Sex_code) + C(Met_code) + ITGB1_residual_z", "ITGB1_residual_z"),
    ("clinical_TMB_plus_ITGB1", "response ~ C(ECOG) + C(IC_code) + C(Sex_code) + C(Met_code) + log1p_TMB + ITGB1_z", "ITGB1_z"),
]


def pool_estimates(rows, imputations):
    table = pd.DataFrame(rows)
    pooled = []
    for (model, term), group in table.groupby(["model", "term"]):
        qbar = group.coef.mean()
        ubar = np.mean(group.se ** 2)
        between = group.coef.var(ddof=1)
        total = ubar + (1 + 1 / imputations) * between
        pooled_se = np.sqrt(total)
        pooled.append({
            "model": model, "term": term, "imputations": imputations,
            "pooled_log_or": qbar, "pooled_se": pooled_se,
            "odds_ratio": np.exp(qbar),
            "ci_low": np.exp(qbar - 1.96 * pooled_se),
            "ci_high": np.exp(qbar + 1.96 * pooled_se),
            "p_value": 2 * norm.sf(abs(qbar / pooled_se)),
            "within_variance": ubar, "between_variance": between,
            "mean_apparent_auc": group.auc.mean(),
            "auc_imputation_low": group.auc.quantile(.025),
            "auc_imputation_high": group.auc.quantile(.975),
            "assumption": "missing at random; predictive mean matching",
        })
    return pd.DataFrame(pooled), table


def main():
    parser = project_root_parser(__doc__)
    parser.add_argument("--imputations", type=int, default=50)
    parser.add_argument("--burn-in", type=int, default=20)
    root, args = parse_project_root(parser)
    frame = load_frame(root)
    encoded = pd.DataFrame({
        "response": frame.response.astype(float),
        "ITGB1_z": frame.ITGB1_z,
        "ITGB1_residual_z": frame.ITGB1_residual_z,
        "ECOG": pd.to_numeric(frame.ECOG, errors="coerce"),
        "IC_code": pd.Categorical(frame.IC_Level, categories=["IC0", "IC1", "IC2+"]).codes.astype(float),
        "Sex_code": pd.Categorical(frame.Sex, categories=["F", "M"]).codes.astype(float),
        "Met_code": pd.Categorical(frame.Met_Status, categories=["LN Only", "Visceral", "Liver"]).codes.astype(float),
        "log1p_TMB": frame.log1p_TMB,
    }, index=frame.index)
    # pandas uses -1 for missing categorical codes; MICE requires NaN.
    for column in ["IC_code", "Sex_code", "Met_code"]:
        encoded.loc[encoded[column].lt(0), column] = np.nan
    np.random.seed(20261105)
    mice = MICEData(encoded, k_pmm=20, perturbation_method="gaussian")
    mice.update_all(args.burn_in)
    rows = []
    for imputation in range(1, args.imputations + 1):
        mice.update_all(5)
        data = mice.data.copy()
        for model_name, formula, term in MODELS:
            fit = smf.glm(formula, data=data, family=sm.families.Binomial()).fit()
            rows.append({
                "imputation": imputation, "model": model_name, "term": term,
                "coef": fit.params[term], "se": fit.bse[term],
                "auc": roc_auc_score(data.response, fit.predict(data)),
            })
    pooled, estimates = pool_estimates(rows, args.imputations)
    out = root / "results/clinical_response"
    out.mkdir(parents=True, exist_ok=True)
    pooled.to_csv(out / "IMvigor210_multiple_imputation_pooled.tsv", sep="\t", index=False)
    estimates.to_csv(out / "IMvigor210_multiple_imputation_estimates.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
