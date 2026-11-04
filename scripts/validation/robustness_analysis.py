"""Run bootstrap and cross-validation robustness analyses."""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import r2_score, roc_auc_score
from sklearn.model_selection import KFold, LeaveOneOut

from scripts.core.cli import parse_project_root, project_root_parser
from scripts.core.common import COMPOSITION_MODELS, ITGB1_FEATURES
from scripts.core.data import load_analysis_features


def fit_r2(frame: pd.DataFrame, predictors: list[str]) -> float:
    d = frame[["ITGB1", *predictors]].dropna()
    return float(sm.OLS(d.ITGB1, sm.add_constant(d[predictors])).fit().rsquared)


def bootstrap_composition(frame: pd.DataFrame, cohort: str, b: int, seed: int):
    rng = np.random.default_rng(seed)
    records = []
    distributions = {name: [] for name in COMPOSITION_MODELS}
    increments = {"Endothelial_beyond_CAF_ECM": [],
                  "OtherComposition_beyond_StromalVascular": []}
    for _ in range(b):
        sample = frame.iloc[rng.integers(0, len(frame), len(frame))]
        vals = {name: fit_r2(sample, predictors) for name, predictors in COMPOSITION_MODELS.items()}
        for name, value in vals.items():
            distributions[name].append(value)
        increments["Endothelial_beyond_CAF_ECM"].append(vals["StromalVascular"] - vals["CAF_ECM"])
        increments["OtherComposition_beyond_StromalVascular"].append(
            vals["FullComposition"] - vals["StromalVascular"])
    for name, values in distributions.items():
        point = fit_r2(frame, COMPOSITION_MODELS[name])
        lo, hi = np.quantile(values, [0.025, 0.975])
        records.append({"cohort": cohort, "estimand": "r_squared", "model": name,
                        "estimate": point, "ci_low": lo, "ci_high": hi,
                        "bootstrap_replicates": b})
    for name, values in increments.items():
        if name.startswith("Endothelial"):
            point = fit_r2(frame, COMPOSITION_MODELS["StromalVascular"]) - fit_r2(frame, COMPOSITION_MODELS["CAF_ECM"])
        else:
            point = fit_r2(frame, COMPOSITION_MODELS["FullComposition"]) - fit_r2(frame, COMPOSITION_MODELS["StromalVascular"])
        lo, hi = np.quantile(values, [0.025, 0.975])
        records.append({"cohort": cohort, "estimand": "incremental_r_squared", "model": name,
                        "estimate": point, "ci_low": lo, "ci_high": hi,
                        "bootstrap_replicates": b})
    return records


def bootstrap_auc(frame: pd.DataFrame, cohort: str, b: int, seed: int):
    rng = np.random.default_rng(seed)
    y = frame.response.to_numpy(int)
    # The direction is fixed before analysis: higher ITGB1 is evaluated as a
    # marker of non-response (resistance). Never flip an AUC by cohort or score.
    y_resistance = 1 - y
    point_auc = {f: roc_auc_score(y_resistance, frame[f]) for f in ITGB1_FEATURES}
    distributions = {f: [] for f in ITGB1_FEATURES}
    deltas = {f: [] for f in ITGB1_FEATURES[1:]}
    for _ in range(b):
        idx = rng.integers(0, len(frame), len(frame))
        if np.unique(y[idx]).size < 2:
            continue
        aucs = {f: roc_auc_score(y_resistance[idx], frame[f].to_numpy()[idx]) for f in ITGB1_FEATURES}
        for f, value in aucs.items():
            distributions[f].append(value)
        for f in ITGB1_FEATURES[1:]:
            deltas[f].append(aucs[f] - aucs["ITGB1"])
    rows = []
    for f in ITGB1_FEATURES:
        lo, hi = np.quantile(distributions[f], [0.025, 0.975])
        rows.append({"cohort": cohort, "estimand": "auc_resistance_direction", "feature": f,
                     "estimate": point_auc[f], "ci_low": lo, "ci_high": hi,
                     "bootstrap_replicates": len(distributions[f]),
                     "orientation": "higher_nonresponse"})
    for f in ITGB1_FEATURES[1:]:
        lo, hi = np.quantile(deltas[f], [0.025, 0.975])
        rows.append({"cohort": cohort, "estimand": "auc_difference_vs_unadjusted",
                     "feature": f, "estimate": point_auc[f] - point_auc["ITGB1"],
                     "ci_low": lo, "ci_high": hi, "bootstrap_replicates": len(deltas[f]),
                     "orientation": "higher_nonresponse"})
    return rows


def cross_validated_r2(frame: pd.DataFrame, cohort: str, repeats: int, seed: int):
    rows = []
    for name, predictors in COMPOSITION_MODELS.items():
        d = frame[["ITGB1", *predictors]].dropna()
        x, y = d[predictors].to_numpy(float), d.ITGB1.to_numpy(float)
        values = []
        for repeat in range(repeats):
            pred = np.full(len(y), np.nan)
            splitter = KFold(n_splits=5, shuffle=True, random_state=seed + repeat)
            for train, test in splitter.split(x):
                fit = sm.OLS(y[train], sm.add_constant(x[train], has_constant="add")).fit()
                pred[test] = fit.predict(sm.add_constant(x[test], has_constant="add"))
            values.append(r2_score(y, pred))
        loo_pred = np.full(len(y), np.nan)
        for train, test in LeaveOneOut().split(x):
            fit = sm.OLS(y[train], sm.add_constant(x[train], has_constant="add")).fit()
            loo_pred[test] = fit.predict(sm.add_constant(x[test], has_constant="add"))
        rows.append({"cohort": cohort, "model": name, "n": len(y),
                     "repeated_5fold_mean_r2": np.mean(values),
                     "repeated_5fold_median_r2": np.median(values),
                     "repeated_5fold_ci_low": np.quantile(values, 0.025),
                     "repeated_5fold_ci_high": np.quantile(values, 0.975),
                     "repeats": repeats, "loocv_r2": r2_score(y, loo_pred)})
    return rows


def main() -> None:
    parser = project_root_parser()
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--cv-repeats", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260814)
    root, args = parse_project_root(parser)
    out = root / "results/robustness"
    out.mkdir(parents=True, exist_ok=True)
    boot_r2, boot_auc, cv = [], [], []
    for i, cohort in enumerate(["GSE91061", "GSE78220", "GSE190826"]):
        frame = load_analysis_features(root, cohort)
        boot_r2.extend(bootstrap_composition(frame, cohort, args.bootstrap, args.seed + i))
        boot_auc.extend(bootstrap_auc(frame, cohort, args.bootstrap, args.seed + 100 + i))
        cv.extend(cross_validated_r2(frame, cohort, args.cv_repeats, args.seed + 200 + i))
    pd.DataFrame(boot_r2).to_csv(out / "bootstrap_composition_r2.tsv", sep="\t", index=False)
    pd.DataFrame(boot_auc).to_csv(out / "bootstrap_auc_and_differences.tsv", sep="\t", index=False)
    pd.DataFrame(cv).to_csv(out / "cross_validated_composition_r2.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
