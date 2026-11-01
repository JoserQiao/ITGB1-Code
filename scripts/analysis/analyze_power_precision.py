"""Minimum detectable effects for the clinical analyses.

These are design-based approximations, not observed post-hoc power. For a
standardized continuous predictor in logistic regression, information is
approximated by n*p*(1-p). For the Cox model, information is approximated by
the number of events. Two-sided alpha=0.05 and 80% power are used.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from scripts.core.cli import parse_project_root


ALPHA = 0.05
TARGET_POWER = 0.80


DESIGNS = [
    ("GSE91061", "ICI objective response", 49, 10, "OR"),
    ("GSE78220", "ICI objective response", 27, 14, "OR"),
    ("GSE176307", "ICI objective response", 87, 16, "OR"),
    ("GSE126044", "ICI objective response", 16, 5, "OR"),
    ("Four-cohort total", "ICI objective response", 179, 45, "OR"),
    ("IMvigor210", "ICI objective response", 298, 68, "OR"),
    ("Five-cohort total", "ICI objective response", 477, 113, "OR"),
    ("GSE190826", "Pathological complete response", 105, 26, "OR"),
    ("GSE135222", "Progression-free survival", 27, 21, "HR"),
]


def standard_error(total: int, events: int, measure: str) -> float:
    if measure == "HR":
        return 1 / np.sqrt(events)
    fraction = events / total
    return 1 / np.sqrt(total * fraction * (1 - fraction))


def power_for_effect(effect: np.ndarray, se: float) -> np.ndarray:
    noncentrality = np.log(effect) / se
    critical = norm.ppf(1 - ALPHA / 2)
    return norm.cdf(-critical - noncentrality) + 1 - norm.cdf(critical - noncentrality)


def main() -> None:
    root, _ = parse_project_root()
    z = norm.ppf(1 - ALPHA / 2) + norm.ppf(TARGET_POWER)
    rows = []
    for cohort, endpoint, total, events, measure in DESIGNS:
        se = standard_error(total, events, measure)
        rows.append({
            "cohort_or_pool": cohort,
            "endpoint": endpoint,
            "n": total,
            "events_or_responders": events,
            "event_fraction": events / total,
            "effect_measure": measure,
            "alpha_two_sided": ALPHA,
            "target_power": TARGET_POWER,
            "approx_log_effect_se": se,
            "minimum_detectable_effect": np.exp(z * se),
            "reciprocal_minimum_detectable_effect": np.exp(-z * se),
        })
    table = pd.DataFrame(rows)
    output = root / "results/power_precision"
    output.mkdir(parents=True, exist_ok=True)
    table.to_csv(output / "minimum_detectable_effects.tsv", sep="\t", index=False)

    figure_dir = root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    display = table[table.cohort_or_pool.isin([
        "GSE126044", "GSE78220", "GSE176307", "Four-cohort total",
        "IMvigor210", "Five-cohort total", "GSE190826", "GSE135222",
    ])].copy()
    y = np.arange(len(display))[::-1]
    colors = ["#4C78A8" if m == "OR" else "#F58518" for m in display.effect_measure]
    axes[0].scatter(display.minimum_detectable_effect, y, c=colors, s=38)
    axes[0].axvline(1, color="#555555", linestyle="--", linewidth=0.8)
    axes[0].set_xscale("log")
    axes[0].set_yticks(y, display.cohort_or_pool)
    axes[0].set_xlabel("Approximate minimum detectable OR or HR")
    axes[0].set_title("A  Two-sided 80% detectable effects")
    axes[0].grid(axis="x", color="#E5E5E5", linewidth=0.7)

    effects = np.linspace(1.0, 4.8, 250)
    for name, color in [
        ("GSE126044", "#B279A2"),
        ("GSE176307", "#59A14F"),
        ("IMvigor210", "#E15759"),
        ("Five-cohort total", "#4E79A7"),
    ]:
        row = table[table.cohort_or_pool.eq(name)].iloc[0]
        axes[1].plot(effects, power_for_effect(effects, row.approx_log_effect_se), label=name, color=color)
    axes[1].axhline(TARGET_POWER, color="#555555", linestyle="--", linewidth=0.8)
    axes[1].set_ylim(0, 1.02)
    axes[1].set_xlim(1, 4.8)
    axes[1].set_xlabel("True OR per 1-SD higher ITGB1")
    axes[1].set_ylabel("Approximate power")
    axes[1].set_title("B  Precision gained from the larger cohort")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(color="#E5E5E5", linewidth=0.7)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"FigureS3_power_precision.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
