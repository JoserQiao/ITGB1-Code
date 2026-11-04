"""Cross-validated incremental-value and missingness analyses for IMvigor210.

All model comparisons use identical patient sets and identical stratified
folds. Apparent AUCs from the primary script are retained for description;
the outputs here provide the internal-validation estimates used for inference.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
from scipy.special import logit
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from scripts.core.cli import parse_project_root, project_root_parser


CATEGORICAL = ["ECOG", "IC_Level", "Sex", "Met_Status"]
warnings.filterwarnings("ignore", message="Found unknown categories")


def load_frame(root):
    processed = root / "data/processed/IMvigor210"
    features_path = processed / "analysis_features_with_residuals.tsv"
    if not features_path.is_file():
        raise FileNotFoundError("Run scripts.analysis.analyze_imvigor210 first")
    features = pd.read_csv(features_path, sep="\t").set_index("sample_id")
    metadata = pd.read_csv(processed / "metadata.tsv", sep="\t").set_index("sample_id")
    metadata["response"] = metadata.binaryResponse.map({"CR/PR": 1, "SD/PD": 0})
    frame = features.join(metadata, how="inner")
    frame = frame[frame.response.notna()].copy()
    frame = frame.rename(columns={
        "Baseline.ECOG.Score": "ECOG",
        "IC.Level": "IC_Level",
        "Met.Disease.Status": "Met_Status",
    })
    for source, target in [
        ("ITGB1", "ITGB1_z"),
        ("ITGB1_residual_FullComposition", "ITGB1_residual_z"),
    ]:
        frame[target] = (frame[source] - frame[source].mean()) / frame[source].std(ddof=0)
    frame["log1p_TMB"] = np.log1p(
        pd.to_numeric(frame["FMOne.mutation.burden.per.MB"], errors="coerce")
    )
    return frame


def make_pipeline(continuous, allow_missing):
    categorical_steps = []
    if allow_missing:
        categorical_steps.append(("impute", SimpleImputer(strategy="constant", fill_value="Missing")))
    categorical_steps.append(("encode", OneHotEncoder(drop="first", handle_unknown="ignore")))
    continuous_steps = []
    if allow_missing:
        continuous_steps.append(("impute", SimpleImputer(strategy="median", add_indicator=True)))
    continuous_steps.append(("scale", StandardScaler()))
    preprocess = ColumnTransformer([
        ("categorical", Pipeline(categorical_steps), CATEGORICAL),
        ("continuous", Pipeline(continuous_steps), continuous),
    ])
    model = LogisticRegression(
        penalty="l2", C=1.0, solver="liblinear", max_iter=2000, random_state=0
    )
    return Pipeline([("preprocess", preprocess), ("model", model)])


def calibration_metrics(y, prediction):
    p = np.clip(np.asarray(prediction, dtype=float), 1e-6, 1 - 1e-6)
    x = logit(p)
    try:
        slope_fit = sm.GLM(y, sm.add_constant(x), family=sm.families.Binomial()).fit()
        intercept_fit = sm.GLM(y, np.ones((len(y), 1)), family=sm.families.Binomial(),
                               offset=x).fit()
        return float(intercept_fit.params[0]), float(slope_fit.params[1])
    except Exception:
        return np.nan, np.nan


def repeated_cv(frame, comparison, baseline_continuous, extended_continuous,
                allow_missing, repeats, seed):
    required = ["response", *CATEGORICAL, *set(baseline_continuous + extended_continuous)]
    data = frame[required].copy()
    if not allow_missing:
        data = data.dropna()
    y = data.response.astype(int).to_numpy()
    rows, comparison_rows = [], []
    for repeat in range(repeats):
        splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed + repeat)
        predictions = {
            "baseline": np.full(len(data), np.nan),
            "extended": np.full(len(data), np.nan),
        }
        for train, test in splitter.split(data, y):
            for model_name, continuous in [
                ("baseline", baseline_continuous),
                ("extended", extended_continuous),
            ]:
                columns = [*CATEGORICAL, *continuous]
                pipeline = make_pipeline(continuous, allow_missing)
                pipeline.fit(data.iloc[train][columns], y[train])
                predictions[model_name][test] = pipeline.predict_proba(
                    data.iloc[test][columns]
                )[:, 1]
        metrics = {}
        for model_name, prediction in predictions.items():
            auc = roc_auc_score(y, prediction)
            brier = brier_score_loss(y, prediction)
            intercept, slope = calibration_metrics(y, prediction)
            metrics[model_name] = {"auc": auc, "brier": brier}
            rows.append({
                "comparison": comparison,
                "model": model_name,
                "repeat": repeat + 1,
                "n": len(data),
                "responders": int(y.sum()),
                "auc": auc,
                "brier": brier,
                "calibration_intercept": intercept,
                "calibration_slope": slope,
                "missing_indicator_analysis": allow_missing,
            })
        comparison_rows.append({
            "comparison": comparison,
            "repeat": repeat + 1,
            "n": len(data),
            "responders": int(y.sum()),
            "delta_auc": metrics["extended"]["auc"] - metrics["baseline"]["auc"],
            "delta_brier": metrics["extended"]["brier"] - metrics["baseline"]["brier"],
            "missing_indicator_analysis": allow_missing,
        })
    return rows, comparison_rows


def summarize_repeats(table, group_columns, metrics):
    rows = []
    for keys, group in table.groupby(group_columns, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_columns, keys))
        row.update({"repeats": len(group), "n": int(group.n.iloc[0]),
                    "responders": int(group.responders.iloc[0])})
        for metric in metrics:
            values = group[metric].dropna().to_numpy(float)
            row[f"{metric}_mean"] = np.mean(values)
            row[f"{metric}_median"] = np.median(values)
            row[f"{metric}_stability_low"] = np.quantile(values, 0.025)
            row[f"{metric}_stability_high"] = np.quantile(values, 0.975)
        rows.append(row)
    return pd.DataFrame(rows)


def missingness_comparison(frame):
    definitions = {
        "clinical_complete_case": frame[CATEGORICAL].notna().all(axis=1),
        "TMB_complete_case": frame[[*CATEGORICAL, "log1p_TMB"]].notna().all(axis=1),
    }
    rows = []
    for definition, complete in definitions.items():
        response_tab = pd.crosstab(
            complete.map({True: "complete", False: "incomplete"}), frame.response.astype(int))
        response_p = fisher_exact(response_tab.to_numpy()).pvalue
        rows.append({
            "definition": definition, "variable": "response", "level": "binary",
            "complete_n": int(response_tab.loc["complete"].sum()),
            "incomplete_n": int(response_tab.loc["incomplete"].sum()),
            "complete_summary": (f"{int(response_tab.loc['complete', 1])}/"
                                 f"{int(response_tab.loc['complete'].sum())} "
                                 f"({100 * response_tab.loc['complete', 1] / response_tab.loc['complete'].sum():.1f}%)"),
            "incomplete_summary": (f"{int(response_tab.loc['incomplete', 1])}/"
                                   f"{int(response_tab.loc['incomplete'].sum())} "
                                   f"({100 * response_tab.loc['incomplete', 1] / response_tab.loc['incomplete'].sum():.1f}%)"),
            "test": "Fisher exact", "p_value": response_p,
        })
        for variable in ["ITGB1", "ITGB1_residual_FullComposition", "log1p_TMB"]:
            a = pd.to_numeric(frame.loc[complete, variable], errors="coerce").dropna()
            b = pd.to_numeric(frame.loc[~complete, variable], errors="coerce").dropna()
            p_value = mannwhitneyu(a, b, alternative="two-sided").pvalue if len(a) and len(b) else np.nan
            rows.append({
                "definition": definition, "variable": variable, "level": "continuous",
                "complete_n": len(a), "incomplete_n": len(b),
                "complete_summary": f"median={a.median():.4g}; IQR={a.quantile(.25):.4g}-{a.quantile(.75):.4g}",
                "incomplete_summary": (f"median={b.median():.4g}; IQR={b.quantile(.25):.4g}-{b.quantile(.75):.4g}"
                                       if len(b) else "not observed"),
                "test": "Mann-Whitney", "p_value": p_value,
            })
        for variable in CATEGORICAL:
            tab = pd.crosstab(complete.map({True: "complete", False: "incomplete"}),
                              frame[variable].fillna("Missing"))
            p_value = chi2_contingency(tab).pvalue if tab.shape[0] == 2 and tab.shape[1] > 1 else np.nan
            for level in tab.columns:
                rows.append({
                    "definition": definition, "variable": variable, "level": str(level),
                    "complete_n": int(tab.loc["complete", level]) if "complete" in tab.index else 0,
                    "incomplete_n": int(tab.loc["incomplete", level]) if "incomplete" in tab.index else 0,
                    "complete_summary": f"{100 * tab.loc['complete', level] / tab.loc['complete'].sum():.1f}%",
                    "incomplete_summary": f"{100 * tab.loc['incomplete', level] / tab.loc['incomplete'].sum():.1f}%",
                    "test": "Chi-square", "p_value": p_value,
                })
    return pd.DataFrame(rows)


def main():
    parser = project_root_parser(__doc__)
    parser.add_argument("--repeats", type=int, default=100)
    root, args = parse_project_root(parser)
    frame = load_frame(root)
    specifications = [
        ("clinical_plus_ITGB1_complete_case", [], ["ITGB1_z"], False),
        ("clinical_plus_residual_complete_case", [], ["ITGB1_residual_z"], False),
        ("TMB_plus_ITGB1_complete_case", ["log1p_TMB"], ["log1p_TMB", "ITGB1_z"], False),
        ("clinical_plus_ITGB1_missing_indicator", [], ["ITGB1_z"], True),
        ("TMB_plus_ITGB1_missing_indicator", ["log1p_TMB"], ["log1p_TMB", "ITGB1_z"], True),
    ]
    model_rows, comparison_rows = [], []
    for index, spec in enumerate(specifications):
        rows, comparisons = repeated_cv(
            frame, *spec, repeats=args.repeats, seed=20261102 + index * 1000
        )
        model_rows.extend(rows)
        comparison_rows.extend(comparisons)
    model_table = pd.DataFrame(model_rows)
    comparison_table = pd.DataFrame(comparison_rows)
    out = root / "results/clinical_response"
    out.mkdir(parents=True, exist_ok=True)
    model_table.to_csv(out / "IMvigor210_repeated_cv_repeats.tsv", sep="\t", index=False)
    comparison_table.to_csv(out / "IMvigor210_repeated_cv_differences.tsv", sep="\t", index=False)
    summarize_repeats(
        model_table, ["comparison", "model", "missing_indicator_analysis"],
        ["auc", "brier", "calibration_intercept", "calibration_slope"],
    ).to_csv(out / "IMvigor210_repeated_cv_summary.tsv", sep="\t", index=False)
    summarize_repeats(
        comparison_table, ["comparison", "missing_indicator_analysis"],
        ["delta_auc", "delta_brier"],
    ).to_csv(out / "IMvigor210_repeated_cv_difference_summary.tsv", sep="\t", index=False)
    missingness_comparison(frame).to_csv(
        out / "IMvigor210_missingness_comparison.tsv", sep="\t", index=False, na_rep="NA"
    )


if __name__ == "__main__":
    main()
