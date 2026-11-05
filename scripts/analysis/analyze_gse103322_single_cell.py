"""Patient-level independent single-cell source validation in GSE103322."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from scripts.core.cli import parse_project_root
from scripts.core.common import bh_fdr


COMPARISONS = [
    ("Fibroblast", "tumor cell"),
    ("Endothelial", "tumor cell"),
]


def main():
    root, _ = parse_project_root()
    directory = root / "data/processed/GSE103322"
    metadata = pd.read_csv(directory / "metadata.tsv", sep="\t")
    expression = pd.read_csv(directory / "itgb1_expression.tsv", sep="\t")
    frame = metadata.merge(expression, on="cell_id", how="inner")
    # The Xena sample-metadata table maps 21 source identifiers to 18 patients.
    frame["patient_id"] = frame["patient"].astype("Int64").astype(str)
    frame.loc[frame["patient"].isna(), "patient_id"] = frame.loc[
        frame["patient"].isna(), "patient_id_source"].astype(str)
    frame["cell_type"] = frame["cell type"].replace({"-Fibroblast": "Fibroblast"})
    frame["ITGB1_detected"] = pd.to_numeric(frame.ITGB1, errors="coerce").gt(0)
    frame = frame.dropna(subset=["patient_id", "cell_type", "ITGB1"])

    patient = (
        frame.groupby(["patient_id", "cell_type"], as_index=False)
        .agg(n_cells=("cell_id", "size"), mean_ITGB1=("ITGB1", "mean"),
             median_ITGB1=("ITGB1", "median"), detected_fraction=("ITGB1_detected", "mean"))
    )
    summary = (
        patient.groupby("cell_type", as_index=False)
        .agg(n_patients=("patient_id", "nunique"), n_cells=("n_cells", "sum"),
             patient_mean_ITGB1=("mean_ITGB1", "mean"),
             patient_median_ITGB1=("mean_ITGB1", "median"),
             patient_mean_detected_fraction=("detected_fraction", "mean"))
        .sort_values("patient_mean_ITGB1", ascending=False)
    )

    comparisons = []
    wide = patient.pivot(index="patient_id", columns="cell_type", values="mean_ITGB1")
    for first, second in COMPARISONS:
        paired = wide[[first, second]].dropna()
        if len(paired) >= 3 and np.any((paired[first] - paired[second]).ne(0)):
            test = wilcoxon(paired[first], paired[second], alternative="two-sided",
                            zero_method="wilcox")
            statistic, p_value = float(test.statistic), float(test.pvalue)
        else:
            statistic, p_value = np.nan, np.nan
        comparisons.append({
            "cell_type_1": first, "cell_type_2": second, "paired_patients": len(paired),
            "median_difference": float((paired[first] - paired[second]).median()) if len(paired) else np.nan,
            "mean_difference": float((paired[first] - paired[second]).mean()) if len(paired) else np.nan,
            "wilcoxon_statistic": statistic, "p_value": p_value,
        })

    out = root / "results/single_cell"
    out.mkdir(parents=True, exist_ok=True)
    patient.to_csv(out / "GSE103322_patient_celltype_ITGB1.tsv", sep="\t", index=False)
    summary.to_csv(out / "GSE103322_celltype_summary.tsv", sep="\t", index=False)
    comparison_table = pd.DataFrame(comparisons)
    comparison_table["fdr"] = bh_fdr(comparison_table.p_value)
    comparison_table.to_csv(
        out / "GSE103322_paired_celltype_comparisons.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
