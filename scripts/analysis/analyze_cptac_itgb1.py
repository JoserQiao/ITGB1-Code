"""Validate RNA-to-protein correspondence for ITGB1 in public CPTAC data."""

from __future__ import annotations

import pandas as pd
from scipy.stats import spearmanr

from scripts.analysis.analyze_gdsc_lineage_adjusted import random_effects_fisher_z
from scripts.core.cli import parse_project_root


def main():
    root, _ = parse_project_root()
    frame = pd.read_csv(
        root / "data/processed/CPTAC_ITGB1/CPTAC_ITGB1_RNA_protein.tsv", sep="\t")
    rows = []
    for study, group in frame.groupby("study"):
        data = group[["ITGB1_RNA", "ITGB1_protein"]].dropna()
        if len(data) < 20:
            continue
        rho, p_value = spearmanr(data.ITGB1_RNA, data.ITGB1_protein)
        rows.append({"study": study, "n": len(data), "spearman_rho": rho, "p_value": p_value})
    associations = pd.DataFrame(rows)
    meta = random_effects_fisher_z(associations)
    out = root / "results/protein_validation"
    out.mkdir(parents=True, exist_ok=True)
    associations.to_csv(out / "CPTAC_ITGB1_RNA_protein_correlations.tsv", sep="\t", index=False)
    pd.DataFrame([meta]).to_csv(out / "CPTAC_ITGB1_RNA_protein_meta.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
