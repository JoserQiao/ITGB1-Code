"""Plot lineage-adjusted pharmacogenomics and external pathway validation."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from scripts.core.cli import parse_project_root
from scripts.core.plotting import clean_axis, panel_label, save_figure, setup_plot_style


PATHWAY_LABELS = {
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
    "HALLMARK_TGF_BETA_SIGNALING": "TGF-β",
    "HALLMARK_ANGIOGENESIS": "Angiogenesis",
    "HALLMARK_HYPOXIA": "Hypoxia",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-γ response",
}
COHORT_LABELS = {
    "GSE172356": "Pancreatic cancer\nGSE172356 (n=62)",
    "GSE113255": "Gastric cancer\nGSE113255 (n=130)",
    "GSE76427": "HCC\nGSE76427 (n=115)",
}


def gdsc_figure(root) -> None:
    data = pd.read_csv(root / "results/pharmacogenomics/gdsc_ITGB1_lineage_adjusted_summary.tsv", sep="\t")
    data = data.sort_values("lineage_adjusted_beta_per_sd")
    y = np.arange(len(data))
    setup_plot_style("ticks")
    fig, axes = plt.subplots(1, 3, figsize=(13, 6.4), constrained_layout=True)

    ax = axes[0]
    for i, row in enumerate(data.itertuples()):
        ax.plot([row.lineage_adjusted_beta_per_sd, row.global_beta_per_sd], [i, i], color="#cccccc", lw=1)
    ax.scatter(data.global_beta_per_sd, y, facecolors="white", edgecolors="#666666", s=30, label="Unadjusted")
    ax.errorbar(
        data.lineage_adjusted_beta_per_sd, y,
        xerr=[data.lineage_adjusted_beta_per_sd - data.lineage_adjusted_ci_low,
              data.lineage_adjusted_ci_high - data.lineage_adjusted_beta_per_sd],
        fmt="o", color="#0072B2", capsize=2, ms=4, label="Lineage-adjusted",
    )
    ax.axvline(0, color="#555555", ls="--", lw=.8)
    ax.set_yticks(y, data.drug_name)
    ax.set_xlabel("Change in ln(IC50) per 1-SD ITGB1")
    ax.set_title("Lineage adjustment attenuates global associations")
    ax.legend(frameon=False, loc="lower right")
    clean_axis(ax, "x"); panel_label(ax, "A", x=-.18)

    ax = axes[1]
    significant = data.meta_fdr.lt(.05)
    colors = np.where(significant, "#D55E00", "#888888")
    for i, row in enumerate(data.itertuples()):
        ax.errorbar(
            row.meta_rho, i,
            xerr=[[row.meta_rho-row.meta_ci_low], [row.meta_ci_high-row.meta_rho]],
            fmt="o", color=colors[i], capsize=2, ms=4,
        )
    ax.axvline(0, color="#555555", ls="--", lw=.8)
    ax.set_yticks(y, data.drug_name)
    ax.set_xlabel("Random-effects within-lineage Spearman ρ")
    ax.set_title("Cross-lineage replication")
    clean_axis(ax, "x"); panel_label(ax, "B", x=-.18)

    ax = axes[2]
    attenuation = data.attenuation_fraction.clip(lower=0, upper=1) * 100
    ax.barh(y, attenuation, color="#56B4E9")
    ax.axvline(50, color="#555555", ls="--", lw=.8)
    ax.set_yticks(y, data.drug_name)
    ax.set_xlabel("Effect attenuation after lineage adjustment (%)")
    ax.set_xlim(0, 100)
    ax.set_title("Magnitude attributable to lineage structure")
    clean_axis(ax, "x"); panel_label(ax, "C", x=-.18)
    save_figure(fig, root / "figures/Figure_extension_GDSC_lineage_adjusted")


def pathway_figure(root) -> None:
    association = pd.read_csv(root / "results/external_solid_tumor/external_pathway_associations.tsv", sep="\t")
    diagnostics = pd.read_csv(root / "results/external_solid_tumor/external_composition_diagnostics.tsv", sep="\t")
    meta = pd.read_csv(root / "results/external_solid_tumor/external_pathway_random_effects_meta.tsv", sep="\t")
    pathway_order = list(PATHWAY_LABELS)
    cohort_order = list(COHORT_LABELS)
    setup_plot_style("ticks")
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), constrained_layout=True)

    ax = axes[0, 0]
    comp = diagnostics[diagnostics.composition_model.eq("FullComposition")].set_index("cohort").loc[cohort_order]
    ax.bar(np.arange(3), comp.r_squared, color=["#0072B2", "#D55E00", "#009E73"])
    ax.set_xticks(np.arange(3), [COHORT_LABELS[x] for x in cohort_order])
    ax.set_ylim(0, .6); ax.set_ylabel("Full-composition model $R^2$")
    ax.set_title("External-cohort composition dependence")
    clean_axis(ax, "y"); panel_label(ax, "A")

    for ax, feature, letter, title in [
        (axes[0, 1], "ITGB1", "B", "Unadjusted ITGB1"),
        (axes[1, 0], "ITGB1_residual_FullComposition", "C", "Composition-adjusted ITGB1"),
    ]:
        subset = association[association.feature.eq(feature)]
        heat = subset.pivot(index="cohort", columns="pathway", values="spearman_rho").loc[cohort_order, pathway_order]
        fdr = subset.pivot(index="cohort", columns="pathway", values="spearman_fdr").loc[cohort_order, pathway_order]
        labels = heat.copy().astype(object)
        for r in range(heat.shape[0]):
            for c in range(heat.shape[1]):
                labels.iloc[r, c] = f"{heat.iloc[r,c]:.2f}{'*' if fdr.iloc[r,c] < .05 else ''}"
        sns.heatmap(
            heat, ax=ax, cmap="RdBu_r", center=0, vmin=-.65, vmax=.65,
            annot=labels, fmt="", linewidths=.5, linecolor="white",
            cbar=feature.endswith("FullComposition"),
            cbar_kws={"label": "Spearman ρ", "shrink": .8},
        )
        ax.set_xticklabels([PATHWAY_LABELS[x] for x in pathway_order], rotation=30, ha="right")
        ax.set_yticklabels([COHORT_LABELS[x] for x in cohort_order], rotation=0)
        ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title(f"{title}–pathway associations (* FDR < 0.05)")
        panel_label(ax, letter, x=-.25 if letter == "C" else -.13)

    ax = axes[1, 1]
    offsets = {"ITGB1": .12, "ITGB1_residual_FullComposition": -.12}
    colors = {"ITGB1": "#777777", "ITGB1_residual_FullComposition": "#0072B2"}
    labels = {"ITGB1": "Unadjusted", "ITGB1_residual_FullComposition": "Full-composition-adjusted"}
    for feature in ["ITGB1", "ITGB1_residual_FullComposition"]:
        d = meta[meta.feature.eq(feature)].set_index("pathway").loc[pathway_order]
        ypos = np.arange(len(d)) + offsets[feature]
        ax.errorbar(
            d.meta_rho, ypos,
            xerr=[d.meta_rho-d.meta_ci_low, d.meta_ci_high-d.meta_rho],
            fmt="o", color=colors[feature], capsize=2, ms=4, label=labels[feature],
        )
    ax.axvline(0, color="#555555", ls="--", lw=.8)
    ax.set_yticks(np.arange(len(pathway_order)), [PATHWAY_LABELS[x] for x in pathway_order])
    ax.set_xlabel("Random-effects pooled Spearman ρ")
    ax.set_title("Pooled pathway associations across three cohorts")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, .5))
    clean_axis(ax, "x"); panel_label(ax, "D", x=-.17)
    save_figure(fig, root / "figures/Figure_extension_external_pathways")


def main() -> None:
    root, _ = parse_project_root()
    (root / "figures").mkdir(parents=True, exist_ok=True)
    gdsc_figure(root)
    pathway_figure(root)


if __name__ == "__main__":
    main()
