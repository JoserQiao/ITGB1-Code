"""Composition-aware pathway validation in external solid-tumor cohorts."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

from scripts.analysis.analyze_gdsc_lineage_adjusted import random_effects_fisher_z
from scripts.core.cli import parse_project_root, project_root_parser
from scripts.core.common import bh_fdr, build_composition_features, read_signatures, signature_score


COHORTS = ["GSE172356", "GSE113255", "GSE76427"]
FEATURES = ["ITGB1", "ITGB1_residual_FullComposition"]
PATHWAYS = [
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
]


def load_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            gene_sets[fields[0]] = fields[2:]
    return gene_sets


def standardized_association(frame: pd.DataFrame, feature: str, pathway: str) -> dict:
    data = frame[[feature, pathway]].dropna()
    x = (data[feature] - data[feature].mean()) / data[feature].std(ddof=0)
    y = (data[pathway] - data[pathway].mean()) / data[pathway].std(ddof=0)
    fit = sm.OLS(y, sm.add_constant(x.rename(feature))).fit(cov_type="HC3")
    ci = fit.conf_int().loc[feature]
    rho, p_spearman = spearmanr(data[feature], data[pathway])
    return {
        "n": len(data), "standardized_beta": fit.params[feature],
        "ci_low": ci.iloc[0], "ci_high": ci.iloc[1], "p_value": fit.pvalues[feature],
        "spearman_rho": rho, "p_spearman": p_spearman,
    }


def main() -> None:
    parser = project_root_parser(__doc__)
    parser.add_argument("--gene-set-root", default=None)
    parser.add_argument("--permutations", type=int, default=500)
    root, args = parse_project_root(parser)
    gene_set_root = Path(args.gene_set_root).resolve() if args.gene_set_root else root
    signatures = read_signatures(root / "config/signatures.tsv")
    gene_sets = load_gmt(gene_set_root / "data/raw/gene_sets/h.all.v2025.1.Hs.symbols.gmt")
    selected_sets = {name: [gene for gene in gene_sets[name] if gene.upper() != "ITGB1"] for name in PATHWAYS}
    out = root / "results/external_solid_tumor"
    out.mkdir(parents=True, exist_ok=True)

    association_rows, diagnostic_rows, coefficient_rows = [], [], []
    for cohort_index, cohort in enumerate(COHORTS):
        processed = root / f"data/processed/{cohort}"
        expression = pd.read_csv(processed / "expression_log2.tsv.gz", sep="\t", index_col=0)
        frame, diagnostics, coefficients = build_composition_features(expression, signatures, cohort)
        for pathway, genes in selected_sets.items():
            frame[pathway] = signature_score(expression, genes)
        frame.to_csv(processed / "analysis_features.tsv", sep="\t", index_label="sample_id")
        diagnostic_rows.append(diagnostics)
        coefficient_rows.append(coefficients)

        for feature in FEATURES:
            for pathway in PATHWAYS:
                association_rows.append({
                    "cohort": cohort, "feature": feature, "pathway": pathway,
                    **standardized_association(frame, feature, pathway),
                })

    associations = pd.DataFrame(association_rows)
    associations["fdr"] = np.nan
    associations["spearman_fdr"] = np.nan
    for feature, index in associations.groupby("feature").groups.items():
        associations.loc[index, "fdr"] = bh_fdr(associations.loc[index, "p_value"])
        associations.loc[index, "spearman_fdr"] = bh_fdr(associations.loc[index, "p_spearman"])
    associations.to_csv(out / "external_pathway_associations.tsv", sep="\t", index=False)

    meta_rows = []
    for (feature, pathway), group in associations.groupby(["feature", "pathway"]):
        meta_rows.append({"feature": feature, "pathway": pathway, **random_effects_fisher_z(group)})
    meta = pd.DataFrame(meta_rows)
    meta["meta_fdr"] = np.nan
    for feature, index in meta.groupby("feature").groups.items():
        meta.loc[index, "meta_fdr"] = bh_fdr(meta.loc[index, "meta_p_value"])
    meta.to_csv(out / "external_pathway_random_effects_meta.tsv", sep="\t", index=False)

    pd.concat(diagnostic_rows, ignore_index=True).to_csv(
        out / "external_composition_diagnostics.tsv", sep="\t", index=False
    )
    pd.concat(coefficient_rows, ignore_index=True).to_csv(
        out / "external_composition_coefficients.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
