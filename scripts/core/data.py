"""Shared cohort loading and expression-normalization utilities."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


EXPRESSION_FILES = {
    "GSE176307": "expression_log2_tpm.tsv.gz",
    "GSE126044": "expression_log2_cpm.tsv.gz",
    "GSE135222": "expression_log2.tsv.gz",
}


def cpm_log2(counts: pd.DataFrame) -> pd.DataFrame:
    """Filter low-count genes and return log2 counts per million."""
    minimum_samples = max(3, int(0.1 * counts.shape[1]))
    filtered = counts.loc[(counts >= 10).sum(axis=1) >= minimum_samples]
    return np.log2(filtered.div(filtered.sum(axis=0), axis=1) * 1e6 + 1)


def cohort_sample_column(cohort: str) -> str:
    return "gsm" if cohort == "GSE190826" else "patient_id"


def load_cohort_expression(root: Path, cohort: str) -> pd.DataFrame:
    """Load the analysis-scale expression matrix for a processed cohort."""
    directory = root / "data" / "processed" / cohort
    if cohort == "GSE190826":
        counts = pd.read_csv(directory / "counts.tsv.gz", sep="\t", index_col=0)
        return cpm_log2(counts)
    filename = EXPRESSION_FILES.get(cohort, "expression_log2_fpkm.tsv.gz")
    return pd.read_csv(directory / filename, sep="\t", index_col=0)


def load_analysis_features(root: Path, cohort: str) -> pd.DataFrame:
    path = root / "data" / "processed" / cohort / "analysis_features_v3.tsv"
    return pd.read_csv(path, sep="\t", index_col=0)


def align_expression_metadata(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Align expression columns and metadata rows without changing column order."""
    available = set(metadata[sample_column])
    samples = [sample for sample in expression.columns if sample in available]
    aligned_metadata = metadata.set_index(sample_column).loc[samples]
    return expression[samples], aligned_metadata
