"""Validate composition estimates with MCP-counter signatures."""

from __future__ import annotations

import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

from scripts.core.cli import parse_project_root
from scripts.core.data import load_analysis_features, load_cohort_expression


def mcp_scores(expression: pd.DataFrame, signatures: dict[str, list[str]]) -> pd.DataFrame:
    upper = {str(g).upper(): g for g in expression.index}
    result = {}
    for population, genes in signatures.items():
        selected = [upper[g.upper()] for g in genes if g.upper() in upper]
        if selected:
            result[population] = expression.loc[selected].mean(axis=0)
    return pd.DataFrame(result)


def main() -> None:
    root, _ = parse_project_root()
    markers = pd.read_csv(root / "config/MCPcounter_genes.txt", sep="\t")
    signatures = markers.groupby("Cell population")["HUGO symbols"].apply(list).to_dict()
    out = root / "results/robustness"
    out.mkdir(parents=True, exist_ok=True)
    diagnostics, concordance = [], []
    mapping = {"Fibroblasts": "CAF_core", "Endothelial cells": "Endothelial",
               "T cells": "T_cell", "Monocytic lineage": "Myeloid",
               "Cytotoxic lymphocytes": "Cytotoxic"}
    predictors = list(mapping)
    for cohort in ["GSE91061", "GSE78220", "GSE190826", "GSE176307", "GSE126044", "GSE135222"]:
        expr = load_cohort_expression(root, cohort)
        feature = load_analysis_features(root, cohort)
        samples = [s for s in expr.columns if s in feature.index]
        score = mcp_scores(expr[samples], signatures)
        score.to_csv(out / f"{cohort}_MCPcounter_scores.tsv", sep="\t", index_label="sample")
        joined = feature.loc[samples].join(score)
        for mcp, internal in mapping.items():
            d = joined[[mcp, internal]].dropna()
            rho, p = spearmanr(d[mcp], d[internal])
            concordance.append({"cohort": cohort, "MCPcounter_population": mcp,
                                "internal_signature": internal, "n": len(d),
                                "spearman_rho": rho, "p_value": p})
        d = joined[["ITGB1", *predictors]].dropna()
        fit = sm.OLS(d.ITGB1, sm.add_constant(d[predictors])).fit()
        diagnostics.append({"cohort": cohort, "method": "MCPcounter",
                            "model": "stromal_immune", "n": len(d),
                            "predictors": ",".join(predictors), "r_squared": fit.rsquared,
                            "adjusted_r_squared": fit.rsquared_adj})
        internal_predictors = list(mapping.values())
        d2 = joined[["ITGB1", *internal_predictors]].dropna()
        fit2 = sm.OLS(d2.ITGB1, sm.add_constant(d2[internal_predictors])).fit()
        diagnostics.append({"cohort": cohort, "method": "InternalRankSignatures",
                            "model": "matched_stromal_immune", "n": len(d2),
                            "predictors": ",".join(internal_predictors), "r_squared": fit2.rsquared,
                            "adjusted_r_squared": fit2.rsquared_adj})
    pd.DataFrame(diagnostics).to_csv(out / "alternative_method_composition_r2.tsv", sep="\t", index=False)
    pd.DataFrame(concordance).to_csv(out / "MCPcounter_signature_concordance.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
