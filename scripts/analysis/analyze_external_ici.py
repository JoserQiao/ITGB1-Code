"""Analyze independent immune-checkpoint inhibitor validation cohorts."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from statsmodels.duration.hazard_regression import PHReg

from scripts.analysis.analyze_response import fit_models
from scripts.core.cli import parse_project_root
from scripts.core.common import ITGB1_FEATURES, build_composition_features, read_signatures
from scripts.core.data import align_expression_metadata, load_cohort_expression
from scripts.validation.robustness_analysis import bootstrap_composition, cross_validated_r2


def survival_analysis(expr, meta, signatures, cohort):
    expr, meta = align_expression_metadata(expr, meta, "patient_id")
    frame, diagnostics, coefficients = build_composition_features(expr, signatures, cohort)
    frame["pfs_time"] = meta.pfs_time.astype(float)
    frame["pfs_event"] = meta.pfs_event.astype(int)
    rows = []
    for feature in ITGB1_FEATURES:
        d = frame[[feature, "pfs_time", "pfs_event"]].dropna()
        z = (d[feature] - d[feature].mean()) / d[feature].std(ddof=0)
        fit = PHReg(d.pfs_time, z.to_frame(feature), status=d.pfs_event).fit(disp=0)
        coef, se = float(fit.params[0]), float(fit.bse[0])
        rows.append({"cohort": cohort, "feature": feature, "n": len(d), "events": int(d.pfs_event.sum()),
                     "hazard_ratio_per_sd": np.exp(coef), "ci_low": np.exp(coef-1.96*se),
                     "ci_high": np.exp(coef+1.96*se), "p_value": float(fit.pvalues[0])})
    return frame, pd.DataFrame(rows), diagnostics, coefficients


def random_effects_meta(response: pd.DataFrame) -> pd.DataFrame:
    """DerSimonian-Laird tau2 with modified Hartung-Knapp inference."""
    rows = []
    for feature, d in response.groupby("feature"):
        d = d[(d.ci_low > 0) & (d.ci_high > 0)].copy()
        yi = np.log(d.odds_ratio_per_sd.to_numpy())
        sei = (np.log(d.ci_high.to_numpy()) - np.log(d.ci_low.to_numpy())) / (2 * 1.96)
        vi = sei ** 2
        wi = 1 / vi
        fixed = np.sum(wi * yi) / np.sum(wi)
        q = np.sum(wi * (yi-fixed)**2)
        df = len(d)-1
        c = np.sum(wi) - np.sum(wi**2)/np.sum(wi)
        tau2 = max(0, (q-df)/c) if c > 0 else 0
        wr = 1/(vi+tau2)
        pooled = np.sum(wr*yi)/np.sum(wr)
        hk_scale = max(1.0, np.sum(wr * (yi - pooled) ** 2) / df) if df > 0 else 1.0
        se = np.sqrt(hk_scale / np.sum(wr))
        critical = student_t.ppf(0.975, df) if df > 0 else 1.96
        p_value = 2 * student_t.sf(abs(pooled / se), df) if df > 0 else np.nan
        rows.append({"feature": feature, "k": len(d), "pooled_or": np.exp(pooled),
                     "ci_low": np.exp(pooled-critical*se), "ci_high": np.exp(pooled+critical*se),
                     "p_value": p_value, "tau_squared": tau2,
                     "i_squared": max(0, (q-df)/q)*100 if q > 0 else 0,
                     "hk_scale": hk_scale, "df": df,
                     "method": "DerSimonian-Laird tau2; modified Hartung-Knapp CI"})
    return pd.DataFrame(rows)


def main() -> None:
    root, _ = parse_project_root()
    signatures = read_signatures(root / "config/signatures.tsv")
    out = root / "results/external_validation"
    out.mkdir(parents=True, exist_ok=True)
    response_results, boot_rows, cv_rows = [], [], []
    for i, cohort in enumerate(["GSE176307", "GSE126044"]):
        d = root / f"data/processed/{cohort}"
        expr = load_cohort_expression(root, cohort)
        meta = pd.read_csv(d / "metadata_v2.tsv", sep="\t")
        frame, result, _, diagnostics, _ = fit_models(expr, meta, signatures, cohort)
        frame.to_csv(d / "analysis_features_v3.tsv", sep="\t", index_label="patient_id")
        result.to_csv(out / f"{cohort}_response_results.tsv", sep="\t", index=False)
        diagnostics.to_csv(out / f"{cohort}_composition_r2.tsv", sep="\t", index=False)
        response_results.append(result[result.feature.isin(ITGB1_FEATURES)])
        boot_rows.extend(bootstrap_composition(frame, cohort, 5000, 20260820+i))
        cv_rows.extend(cross_validated_r2(frame, cohort, 200, 20260920+i))

    d = root / "data/processed/GSE135222"
    expr = load_cohort_expression(root, "GSE135222")
    meta = pd.read_csv(d / "metadata_v2.tsv", sep="\t")
    frame, surv, diagnostics, _ = survival_analysis(expr, meta, signatures, "GSE135222")
    frame.to_csv(d / "analysis_features_v3.tsv", sep="\t", index_label="patient_id")
    surv.to_csv(out / "GSE135222_pfs_results.tsv", sep="\t", index=False)
    diagnostics.to_csv(out / "GSE135222_composition_r2.tsv", sep="\t", index=False)
    boot_rows.extend(bootstrap_composition(frame, "GSE135222", 5000, 20260822))
    cv_rows.extend(cross_validated_r2(frame, "GSE135222", 200, 20260922))

    new_response = pd.concat(response_results, ignore_index=True)
    new_response.to_csv(out / "new_ici_response_results.tsv", sep="\t", index=False)
    old = pd.read_csv(root / "results/treatment_response/anti_pd1_model_results_v3.tsv", sep="\t")
    combined = pd.concat([old[old.feature.isin(ITGB1_FEATURES)], new_response], ignore_index=True)
    combined.to_csv(out / "four_ici_cohort_response_results.tsv", sep="\t", index=False)
    random_effects_meta(combined).to_csv(out / "random_effects_meta_analysis.tsv", sep="\t", index=False)
    pd.DataFrame(boot_rows).to_csv(out / "external_bootstrap_composition_r2.tsv", sep="\t", index=False)
    pd.DataFrame(cv_rows).to_csv(out / "external_cross_validated_composition_r2.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
