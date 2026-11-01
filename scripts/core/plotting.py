"""Shared publication-figure constants and helpers."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


COHORT_ORDER = ["GSE91061", "GSE78220", "GSE190826", "GSE176307", "GSE126044", "GSE135222"]
COHORT_COLORS = {
    "GSE91061": "#0072B2",
    "GSE78220": "#D55E00",
    "GSE190826": "#009E73",
    "GSE176307": "#CC79A7",
    "GSE126044": "#E69F00",
    "GSE135222": "#56B4E9",
}
FEATURE_LABELS = {
    "ITGB1": "Unadjusted ITGB1",
    "ITGB1_residual_CAF_ECM": "CAF/ECM-adjusted",
    "ITGB1_residual_Endothelial": "Endothelial-adjusted",
    "ITGB1_residual_StromalVascular": "Stromal/vascular-adjusted",
    "ITGB1_residual_FullComposition": "Full-composition-adjusted",
}


def setup_plot_style(style: str = "whitegrid", *, bold_titles: bool = False) -> None:
    sns.set_theme(style=style, context="paper")
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.titlesize": 10,
        "axes.titleweight": "bold" if bold_titles else "normal",
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.8,
        "figure.dpi": 150,
        "savefig.dpi": 450,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def panel_label(ax, label: str, *, x: float = -0.12) -> None:
    ax.text(x, 1.08, label, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def clean_axis(ax, grid: str | None = "x") -> None:
    sns.despine(ax=ax)
    ax.grid(False)
    if grid:
        ax.grid(axis=grid, color="#dddddd", linewidth=0.6, zorder=0)


def save_figure(fig, output_base: Path) -> None:
    for extension in ("png", "pdf", "svg"):
        fig.savefig(output_base.with_suffix(f".{extension}"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load_composition_results(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine original and external composition robustness results."""
    bootstrap = pd.concat([
        pd.read_csv(root / "results/robustness/bootstrap_composition_r2.tsv", sep="\t"),
        pd.read_csv(root / "results/external_validation/external_bootstrap_composition_r2.tsv", sep="\t"),
    ], ignore_index=True)
    cross_validation = pd.concat([
        pd.read_csv(root / "results/robustness/cross_validated_composition_r2.tsv", sep="\t"),
        pd.read_csv(root / "results/external_validation/external_cross_validated_composition_r2.tsv", sep="\t"),
    ], ignore_index=True)
    return bootstrap, cross_validation
