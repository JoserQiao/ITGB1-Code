"""Direct histological composition validation using TCGA pathology slides."""

from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

from scripts.analysis.analyze_gdsc_lineage_adjusted import random_effects_fisher_z
from scripts.core.cli import parse_project_root
from scripts.core.common import bh_fdr


VARIABLES = ["percent_tumor_nuclei", "percent_tumor_cells", "percent_stromal_cells"]


def main():
    root, _ = parse_project_root()
    path = root / "data/processed/TCGA_pathology_purity/TCGA_ITGB1_pathology_combined.tsv"
    frame = pd.read_csv(path, sep="\t")
    rows = []
    for variable in VARIABLES:
        for cancer, group in frame.groupby("cancer"):
            data = group[["ITGB1", variable]].dropna()
            if len(data) < 20:
                continue
            rho, p_value = spearmanr(data.ITGB1, data[variable])
            rows.append({
                "variable": variable, "cancer": cancer, "n": len(data),
                "spearman_rho": rho, "p_value": p_value,
            })
    associations = pd.DataFrame(rows)
    associations["fdr"] = np.nan
    for variable, index in associations.groupby("variable").groups.items():
        associations.loc[index, "fdr"] = bh_fdr(associations.loc[index, "p_value"])

    meta_rows = []
    for variable, group in associations.groupby("variable"):
        meta_rows.append({"variable": variable, **random_effects_fisher_z(group)})
    meta = pd.DataFrame(meta_rows)
    if not meta.empty:
        meta["meta_fdr"] = bh_fdr(meta.meta_p_value)

    model_rows = []
    for variable in VARIABLES:
        data = frame[["cancer", "ITGB1", variable]].dropna().copy()
        data["ITGB1_z_within_cancer"] = data.groupby("cancer").ITGB1.transform(
            lambda x: (x - x.mean()) / x.std(ddof=0))
        data["pathology_z_within_cancer"] = data.groupby("cancer")[variable].transform(
            lambda x: (x - x.mean()) / x.std(ddof=0))
        data = data.replace([np.inf, -np.inf], np.nan).dropna()
        fit = smf.ols("ITGB1_z_within_cancer ~ pathology_z_within_cancer + C(cancer)",
                      data=data).fit(cov_type="HC3")
        ci = fit.conf_int().loc["pathology_z_within_cancer"]
        model_rows.append({
            "variable": variable, "n": len(data), "cancers": data.cancer.nunique(),
            "beta": fit.params["pathology_z_within_cancer"],
            "ci_low": ci.iloc[0], "ci_high": ci.iloc[1],
            "p_value": fit.pvalues["pathology_z_within_cancer"],
        })

    out = root / "results/pathology_purity"
    out.mkdir(parents=True, exist_ok=True)
    associations.to_csv(out / "TCGA_ITGB1_pathology_correlations.tsv", sep="\t", index=False)
    meta.to_csv(out / "TCGA_ITGB1_pathology_meta.tsv", sep="\t", index=False)
    pd.DataFrame(model_rows).to_csv(
        out / "TCGA_ITGB1_pathology_fixed_effect_models.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
