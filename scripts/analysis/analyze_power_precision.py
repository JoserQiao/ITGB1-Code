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
from scipy.stats import norm, t as student_t

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


def meta_analytic_power_simulation(replicates=20000, seed=20261102):
    """Simulate a normal approximation to the five-cohort meta-analytic design.

    Cohort-specific log-OR standard errors use each observed sample size and
    response fraction. Between-cohort heterogeneity is prespecified on the
    log-OR scale. Modified Hartung-Knapp inference is then applied using the
    corresponding fixed tau-based random-effects weights. Tau is not
    re-estimated in each replicate. This avoids treating n=477 as one
    homogeneous cohort while remaining a design approximation rather than a
    patient-level outcome simulation.
    """
    rng = np.random.default_rng(seed)
    clinical = [row for row in DESIGNS if row[0] in {
        "GSE91061", "GSE78220", "GSE176307", "GSE126044", "IMvigor210"
    }]
    se = np.asarray([standard_error(total, events, measure)
                     for _, _, total, events, measure in clinical])
    k = len(se)
    df = k - 1
    critical = student_t.ppf(1 - ALPHA / 2, df)
    effects = np.round(np.arange(1.0, 2.505, 0.025), 3)
    rows = []
    for tau in [0.0, 0.10, 0.20]:
        weights = 1 / (se ** 2 + tau ** 2)
        weight_sum = weights.sum()
        for effect in effects:
            beta = np.log(effect)
            heterogeneity = rng.normal(0, tau, size=(replicates, k))
            sampling = rng.normal(0, se, size=(replicates, k))
            estimates = beta + heterogeneity + sampling
            pooled = np.sum(estimates * weights, axis=1) / weight_sum
            scale = np.sum(weights * (estimates - pooled[:, None]) ** 2, axis=1) / df
            scale = np.maximum(1.0, scale)
            pooled_se = np.sqrt(scale / weight_sum)
            statistic = np.abs(pooled / pooled_se)
            power = np.mean(statistic > critical)
            rows.append({
                "true_or_magnitude": effect,
                "tau_log_or": tau,
                "replicates": replicates,
                "power": power,
                "method": "five-cohort normal-approximation simulation; modified Hartung-Knapp",
            })
    curves = pd.DataFrame(rows)
    mde_rows = []
    for tau, group in curves.groupby("tau_log_or"):
        eligible = group[group.power.ge(TARGET_POWER)]
        mde_rows.append({
            "tau_log_or": tau,
            "target_power": TARGET_POWER,
            "simulated_minimum_detectable_or_magnitude": (
                eligible.true_or_magnitude.iloc[0] if len(eligible) else np.nan
            ),
            "replicates_per_effect": replicates,
        })
    return curves, pd.DataFrame(mde_rows)


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
    curves, simulated_mde = meta_analytic_power_simulation()
    curves.to_csv(output / "five_cohort_meta_simulation_power.tsv", sep="\t", index=False)
    simulated_mde.to_csv(output / "five_cohort_meta_simulation_mde.tsv", sep="\t", index=False)

    figure_dir = root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    display = table[table.cohort_or_pool.isin([
        "GSE91061", "GSE126044", "GSE78220", "GSE176307",
        "IMvigor210", "GSE190826", "GSE135222",
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

    for tau, color in [(0.0, "#4E79A7"), (0.10, "#59A14F"), (0.20, "#B279A2")]:
        group = curves[curves.tau_log_or.eq(tau)]
        axes[1].plot(group.true_or_magnitude, group.power,
                     label=f"tau(log OR) = {tau:.2f}", color=color)
    axes[1].axhline(TARGET_POWER, color="#555555", linestyle="--", linewidth=0.8)
    axes[1].set_ylim(0, 1.02)
    axes[1].set_xlim(1, 2.5)
    axes[1].set_xlabel("True OR magnitude per 1-SD higher ITGB1")
    axes[1].set_ylabel("Simulated meta-analytic power")
    axes[1].set_title("B  Five-cohort design under heterogeneity")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(color="#E5E5E5", linewidth=0.7)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
    for suffix in ("png", "pdf"):
        fig.savefig(figure_dir / f"FigureS3_power_precision.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
