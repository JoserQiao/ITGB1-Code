"""Test ITGB1--GDSC drug associations after cancer-lineage adjustment.

The primary assay for a duplicated drug name is selected by the largest
available sample size, never by its effect size or P value.  Two complementary
lineage-aware analyses are reported: an OLS model with lineage fixed effects
and a random-effects meta-analysis of within-lineage Spearman correlations.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr, t as student_t

from scripts.core.cli import parse_project_root, project_root_parser
from scripts.core.common import bh_fdr


MIN_TOTAL_N = 100
MIN_LINEAGE_N_FIXED = 5
MIN_LINEAGE_N_META = 20
MIN_META_LINEAGES = 3


def jsonapi_table(path: Path) -> pd.DataFrame:
    """Read Cell Model Passports JSON:API exports saved as JSON or JSONL."""
    if path.suffix == ".jsonl":
        items = [json.loads(line) for line in path.open(encoding="utf-8") if line.strip()]
    else:
        items = json.loads(path.read_text(encoding="utf-8"))["data"]
    rows = []
    for item in items:
        row = {"id": item.get("id"), **item.get("attributes", {})}
        for name, relationship in item.get("relationships", {}).items():
            data = relationship.get("data")
            if isinstance(data, dict):
                row[f"{name}_id"] = data.get("id")
        rows.append(row)
    return pd.DataFrame(rows)


def cell_model_lineages(path: Path) -> pd.DataFrame:
    """Map Cell Model Passports model IDs to portal-native tissue lineages."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    tissues = {
        str(item["id"]): item.get("attributes", {}).get("name")
        for item in payload.get("included", [])
        if item.get("type") == "tissue"
    }
    rows = []
    for sample in payload.get("data", []):
        relationships = sample.get("relationships", {})
        tissue_data = relationships.get("tissue", {}).get("data") or {}
        lineage = tissues.get(str(tissue_data.get("id")))
        if not lineage:
            lineage = sample.get("attributes", {}).get("sample_site")
        for model in relationships.get("models", {}).get("data") or []:
            rows.append({"model_id": str(model.get("id")), "lineage": lineage})
    table = pd.DataFrame(rows).dropna(subset=["model_id", "lineage"])
    return table.drop_duplicates("model_id")


def random_effects_fisher_z(table: pd.DataFrame) -> dict[str, float | int]:
    """DerSimonian--Laird random-effects pooling of Fisher-z correlations."""
    d = table[(table.n >= MIN_LINEAGE_N_META) & table.spearman_rho.notna()].copy()
    if len(d) < MIN_META_LINEAGES:
        return {
            "k_lineages": len(d), "meta_rho": np.nan, "meta_ci_low": np.nan,
            "meta_ci_high": np.nan, "meta_p_value": np.nan,
            "tau_squared": np.nan, "i_squared": np.nan,
        }
    rho = d.spearman_rho.clip(-0.999999, 0.999999).to_numpy(float)
    yi = np.arctanh(rho)
    vi = 1.0 / (d.n.to_numpy(float) - 3.0)
    wi = 1.0 / vi
    fixed = np.sum(wi * yi) / np.sum(wi)
    q = np.sum(wi * (yi - fixed) ** 2)
    df = len(d) - 1
    c = np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    wr = 1.0 / (vi + tau2)
    pooled = np.sum(wr * yi) / np.sum(wr)
    hk_scale = max(1.0, np.sum(wr * (yi - pooled) ** 2) / df) if df > 0 else 1.0
    se = np.sqrt(hk_scale / np.sum(wr))
    critical = student_t.ppf(0.975, df) if df > 0 else 1.96
    p_value = 2 * student_t.sf(abs(pooled / se), df) if df > 0 else np.nan
    return {
        "k_lineages": len(d),
        "meta_rho": float(np.tanh(pooled)),
        "meta_ci_low": float(np.tanh(pooled - critical * se)),
        "meta_ci_high": float(np.tanh(pooled + critical * se)),
        "meta_p_value": float(p_value),
        "tau_squared": float(tau2),
        "i_squared": float(max(0.0, (q - df) / q) * 100 if q > 0 else 0.0),
        "hk_scale": float(hk_scale),
        "df": int(df),
        "meta_method": "DerSimonian-Laird tau2; modified Hartung-Knapp CI",
    }


def load_analysis_table(root: Path, config_root: Path) -> pd.DataFrame:
    classes = pd.read_csv(config_root / "config/gdsc_drug_classes.tsv", sep="\t")
    classes = classes[classes.include_primary.eq(1)].drop_duplicates("drug_name")

    expression = jsonapi_table(root / "data/raw/CellModelPassports/ITGB1_rnaseq.json")
    expression = (
        expression[expression.dataset_id.astype(str).eq("15")]
        .drop_duplicates("model_id")
        .assign(model_id=lambda x: x.model_id.astype(str),
                ITGB1_log2_tpm=lambda x: np.log2(pd.to_numeric(x.tpm) + 1))
        .set_index("model_id")
    )
    curves = jsonapi_table(root / "data/raw/CellModelPassports/GDSC_fitted_curves.jsonl")
    curves["drug_id"] = curves.drug_id.astype(str)
    curves["model_id"] = curves.model_id.astype(str)

    drugs = jsonapi_table(root / "data/raw/CellModelPassports/drugs.jsonl")
    drugs["drug_id"] = drugs.id.astype(str)
    drug_name_column = "drug_name" if "drug_name" in drugs.columns else "name"
    drugs = drugs.rename(columns={drug_name_column: "drug_name"})
    curves = curves.merge(drugs[["drug_id", "drug_name"]], on="drug_id", how="left")
    curves = curves.merge(classes[["drug_name", "drug_class"]], on="drug_name", how="inner")
    curves = curves.join(expression[["ITGB1_log2_tpm"]], on="model_id", how="inner")

    model = cell_model_lineages(root / "data/raw/CellModelPassports/samples.json")
    curves = curves.merge(model, on="model_id", how="inner")
    curves["ln_ic50"] = pd.to_numeric(curves.ln_ic50, errors="coerce")
    curves = curves.dropna(subset=["ln_ic50", "ITGB1_log2_tpm", "lineage"])

    # Select the best-covered assay for duplicated drug names without looking
    # at an association estimate.  Retain all assay-level data in a QC table.
    assay_n = (
        curves.groupby(["drug_name", "drug_id"], as_index=False)
        .agg(n=("model_id", "nunique"))
        .sort_values(["drug_name", "n", "drug_id"], ascending=[True, False, True])
    )
    primary = assay_n.drop_duplicates("drug_name")[["drug_name", "drug_id"]]
    curves = curves.merge(primary.assign(primary_assay=True), on=["drug_name", "drug_id"], how="left")
    curves["primary_assay"] = curves.primary_assay.notna()
    return curves


def analyze_drug(group: pd.DataFrame) -> tuple[dict, list[dict], list[dict]]:
    d = group.copy()
    lineage_counts = d.lineage.value_counts()
    keep_lineages = lineage_counts[lineage_counts >= MIN_LINEAGE_N_FIXED].index
    fixed = d[d.lineage.isin(keep_lineages)].copy()
    fixed["ITGB1_z"] = (
        (fixed.ITGB1_log2_tpm - fixed.ITGB1_log2_tpm.mean())
        / fixed.ITGB1_log2_tpm.std(ddof=0)
    )

    global_fit = smf.ols("ln_ic50 ~ ITGB1_z", data=fixed).fit(cov_type="HC3")
    adjusted_fit = smf.ols("ln_ic50 ~ ITGB1_z + C(lineage)", data=fixed).fit(cov_type="HC3")
    global_ci = global_fit.conf_int().loc["ITGB1_z"]
    adjusted_ci = adjusted_fit.conf_int().loc["ITGB1_z"]

    lineage_rows = []
    for lineage, ld in d.groupby("lineage"):
        if len(ld) < MIN_LINEAGE_N_META or ld.ITGB1_log2_tpm.nunique() < 3 or ld.ln_ic50.nunique() < 3:
            continue
        rho, p_value = spearmanr(ld.ITGB1_log2_tpm, ld.ln_ic50)
        lineage_rows.append({
            "drug_id": d.drug_id.iloc[0], "drug_name": d.drug_name.iloc[0],
            "drug_class": d.drug_class.iloc[0], "lineage": lineage,
            "n": len(ld), "spearman_rho": rho, "p_value": p_value,
        })
    lineage_table = pd.DataFrame(lineage_rows)
    meta = random_effects_fisher_z(lineage_table) if not lineage_table.empty else random_effects_fisher_z(
        pd.DataFrame(columns=["n", "spearman_rho"])
    )

    loo_rows = []
    for omitted in sorted(fixed.lineage.unique()):
        ld = fixed[fixed.lineage.ne(omitted)]
        if ld.lineage.nunique() < 2:
            continue
        fit = smf.ols("ln_ic50 ~ ITGB1_z + C(lineage)", data=ld).fit(cov_type="HC3")
        ci = fit.conf_int().loc["ITGB1_z"]
        loo_rows.append({
            "drug_id": d.drug_id.iloc[0], "drug_name": d.drug_name.iloc[0],
            "omitted_lineage": omitted, "n": len(ld),
            "beta": fit.params["ITGB1_z"], "ci_low": ci.iloc[0],
            "ci_high": ci.iloc[1], "p_value": fit.pvalues["ITGB1_z"],
        })

    row = {
        "drug_id": d.drug_id.iloc[0], "drug_name": d.drug_name.iloc[0],
        "drug_class": d.drug_class.iloc[0], "n_total": len(d),
        "n_fixed_effect_model": len(fixed), "n_lineages_fixed": fixed.lineage.nunique(),
        "global_beta_per_sd": global_fit.params["ITGB1_z"],
        "global_ci_low": global_ci.iloc[0], "global_ci_high": global_ci.iloc[1],
        "global_p_value": global_fit.pvalues["ITGB1_z"],
        "lineage_adjusted_beta_per_sd": adjusted_fit.params["ITGB1_z"],
        "lineage_adjusted_ci_low": adjusted_ci.iloc[0],
        "lineage_adjusted_ci_high": adjusted_ci.iloc[1],
        "lineage_adjusted_p_value": adjusted_fit.pvalues["ITGB1_z"],
        **meta,
    }
    return row, lineage_rows, loo_rows


def main() -> None:
    parser = project_root_parser(__doc__)
    parser.add_argument(
        "--config-root", default=None,
        help="Code/config repository; defaults to project-root.",
    )
    parser.add_argument(
        "--output-root", default=None,
        help="Destination project; defaults to project-root.",
    )
    root, args = parse_project_root(parser)
    config_root = Path(args.config_root).resolve() if args.config_root else root
    output_root = Path(args.output_root).resolve() if args.output_root else root
    out = output_root / "results/pharmacogenomics"
    out.mkdir(parents=True, exist_ok=True)

    curves = load_analysis_table(root, config_root)
    curves.to_csv(out / "gdsc_ITGB1_lineage_analysis_input.tsv.gz", sep="\t", index=False, compression="gzip")
    assay_qc = (
        curves.groupby(["drug_name", "drug_class", "drug_id", "primary_assay"], as_index=False)
        .agg(n=("model_id", "nunique"), n_lineages=("lineage", "nunique"))
    )
    assay_qc.to_csv(out / "gdsc_chemotherapy_assay_selection.tsv", sep="\t", index=False)

    summary_rows, lineage_rows, loo_rows = [], [], []
    primary = curves[curves.primary_assay].copy()
    for _, d in primary.groupby(["drug_name", "drug_id"], sort=True):
        if len(d) < MIN_TOTAL_N:
            continue
        row, within, loo = analyze_drug(d)
        summary_rows.append(row)
        lineage_rows.extend(within)
        loo_rows.extend(loo)
    summary = pd.DataFrame(summary_rows)
    summary["global_fdr"] = bh_fdr(summary.global_p_value)
    summary["lineage_adjusted_fdr"] = bh_fdr(summary.lineage_adjusted_p_value)
    summary["meta_fdr"] = bh_fdr(summary.meta_p_value)
    summary["attenuation_fraction"] = 1 - (
        summary.lineage_adjusted_beta_per_sd.abs() / summary.global_beta_per_sd.abs()
    )
    summary["direction_concordant"] = np.sign(summary.lineage_adjusted_beta_per_sd).eq(np.sign(summary.meta_rho))
    summary.sort_values("lineage_adjusted_fdr").to_csv(
        out / "gdsc_ITGB1_lineage_adjusted_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(lineage_rows).to_csv(
        out / "gdsc_ITGB1_within_lineage_correlations.tsv", sep="\t", index=False
    )
    pd.DataFrame(loo_rows).to_csv(
        out / "gdsc_ITGB1_leave_one_lineage_out.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
