"""Plot independent source, pathology-composition and internal-validation analyses."""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import PercentFormatter

from scripts.core.cli import parse_project_root
from scripts.core.plotting import clean_axis, panel_label, save_figure, setup_plot_style


def figure_s4(root, out):
    summary = pd.read_csv(root / "results/single_cell/GSE103322_celltype_summary.tsv", sep="\t")
    patient = pd.read_csv(root / "results/single_cell/GSE103322_patient_celltype_ITGB1.tsv", sep="\t")
    summary = summary.sort_values("patient_mean_ITGB1")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)

    ax = axes[0]
    panel_label(ax, "A")
    y = np.arange(len(summary))
    points = ax.scatter(summary.patient_mean_ITGB1, y,
                        s=35 + 12 * summary.n_patients,
                        c=summary.patient_mean_detected_fraction,
                        cmap="viridis", vmin=0, vmax=1, edgecolor="white", linewidth=.7)
    ax.set_yticks(y, summary.cell_type)
    ax.set_xlabel("Patient-averaged ITGB1 expression")
    ax.set_title("Independent GSE103322 cell-type profile")
    fig.colorbar(points, ax=ax, label="Mean detected fraction", shrink=.78)
    clean_axis(ax, "x")

    ax = axes[1]
    panel_label(ax, "B")
    pivot = patient.pivot(index="patient_id", columns="cell_type", values="mean_ITGB1")
    comparisons = [("Fibroblast", "tumor cell"), ("Endothelial", "tumor cell")]
    colors = ["#CC79A7", "#56B4E9"]
    for index, ((cell_type, reference), color) in enumerate(zip(comparisons, colors)):
        paired = pivot[[cell_type, reference]].dropna()
        differences = paired[cell_type] - paired[reference]
        rng = np.random.default_rng(20261102 + index)
        x = np.full(len(differences), index) + rng.normal(0, .035, len(differences))
        ax.scatter(x, differences, color=color, alpha=.75, s=25)
        ax.plot([index - .18, index + .18], [differences.median()] * 2,
                color="#222222", linewidth=2)
    ax.axhline(0, color="#555555", linestyle="--", linewidth=.8)
    ax.set_xticks([0, 1], ["Fibroblast - tumor", "Endothelial - tumor"])
    ax.set_ylabel("Within-patient mean-expression difference")
    ax.set_title("Patient-level paired contrasts")
    clean_axis(ax, "y")
    save_figure(fig, out / "FigureS4_GSE103322_single_cell_validation")


def correlation_interval(rho, n):
    z = np.arctanh(np.clip(rho, -.999999, .999999))
    se = 1 / np.sqrt(n - 3)
    return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)


def figure_s5(root, out):
    associations = pd.read_csv(
        root / "results/pathology_purity/TCGA_ITGB1_pathology_correlations.tsv", sep="\t")
    meta = pd.read_csv(root / "results/pathology_purity/TCGA_ITGB1_pathology_meta.tsv", sep="\t")
    variables = [
        ("percent_stromal_cells", "Pathologist-estimated stromal cells"),
        ("percent_tumor_nuclei", "Pathologist-estimated tumour nuclei"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.4), constrained_layout=True)
    for letter, ax, (variable, title) in zip(["A", "B"], axes, variables):
        panel_label(ax, letter)
        data = associations[associations.variable.eq(variable)].sort_values("spearman_rho")
        pooled = meta[meta.variable.eq(variable)].iloc[0]
        labels = data.cancer.tolist() + ["Random-effects pooled"]
        estimates = data.spearman_rho.tolist() + [pooled.meta_rho]
        intervals = [correlation_interval(r.spearman_rho, r.n) for r in data.itertuples()]
        intervals.append((pooled.meta_ci_low, pooled.meta_ci_high))
        y = np.arange(len(labels))[::-1]
        for i, (estimate, interval) in enumerate(zip(estimates, intervals)):
            marker = "D" if i == len(labels) - 1 else "o"
            color = "#222222" if i == len(labels) - 1 else "#4E79A7"
            ax.errorbar(estimate, y[i], xerr=[[estimate - interval[0]], [interval[1] - estimate]],
                        fmt=marker, color=color, capsize=2, markersize=5)
        ax.axvline(0, color="#555555", linestyle="--", linewidth=.8)
        ax.set_yticks(y, labels)
        ax.set_xlabel("Spearman rho with ITGB1 (95% CI)")
        ax.set_title(title)
        clean_axis(ax, "x")
    save_figure(fig, out / "FigureS5_TCGA_pathology_composition")


def figure_s6(root, out):
    model = pd.read_csv(
        root / "results/clinical_response/IMvigor210_repeated_cv_summary.tsv", sep="\t")
    diff = pd.read_csv(
        root / "results/clinical_response/IMvigor210_repeated_cv_difference_summary.tsv", sep="\t")
    selected = [
        "clinical_plus_ITGB1_complete_case",
        "clinical_plus_residual_complete_case",
        "TMB_plus_ITGB1_complete_case",
        "clinical_plus_ITGB1_missing_indicator",
        "TMB_plus_ITGB1_missing_indicator",
    ]
    labels = [
        "Clinical + ITGB1\ncomplete-case",
        "Clinical + residual\ncomplete-case",
        "Clinical + TMB + ITGB1\ncomplete-case",
        "Clinical + ITGB1\nmissing-indicator",
        "Clinical + TMB + ITGB1\nmissing-indicator",
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), constrained_layout=True)

    ax = axes[0, 0]
    panel_label(ax, "A")
    for offset, model_name, marker in [(-.10, "baseline", "o"), (.10, "extended", "s")]:
        data = model[model.model.eq(model_name)].set_index("comparison").loc[selected]
        y = np.arange(len(selected)) + offset
        ax.errorbar(data.auc_mean, y,
                    xerr=[data.auc_mean - data.auc_stability_low,
                          data.auc_stability_high - data.auc_mean],
                    fmt=marker, capsize=2, label=model_name.capitalize())
    ax.set_yticks(np.arange(len(selected)), labels)
    ax.set_xlabel("Repeated five-fold CV AUC for response")
    ax.set_title("Internally validated discrimination")
    ax.legend(frameon=False)
    clean_axis(ax, "x")

    data = diff.set_index("comparison").loc[selected]
    ax = axes[0, 1]
    panel_label(ax, "B")
    y = np.arange(len(selected))
    ax.errorbar(data.delta_auc_mean, y,
                xerr=[data.delta_auc_mean - data.delta_auc_stability_low,
                      data.delta_auc_stability_high - data.delta_auc_mean],
                fmt="o", color="#4E79A7", capsize=2)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Extended - baseline CV AUC")
    ax.set_title("Incremental discrimination")
    clean_axis(ax, "x")

    ax = axes[1, 0]
    panel_label(ax, "C")
    ax.errorbar(data.delta_brier_mean, y,
                xerr=[data.delta_brier_mean - data.delta_brier_stability_low,
                      data.delta_brier_stability_high - data.delta_brier_mean],
                fmt="o", color="#E15759", capsize=2)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Extended - baseline Brier score")
    ax.set_title("Incremental prediction error")
    clean_axis(ax, "x")

    ax = axes[1, 1]
    panel_label(ax, "D")
    metadata = pd.read_csv(root / "data/processed/IMvigor210/metadata.tsv", sep="\t")
    metadata["response"] = metadata.binaryResponse.map({"CR/PR": 1, "SD/PD": 0})
    data = metadata[metadata.response.notna()].copy()
    clinical_complete = data[["IC.Level", "Met.Disease.Status"]].notna().all(axis=1)
    tmb_complete = clinical_complete & data["FMOne.mutation.burden.per.MB"].notna()
    rates = []
    for name, complete in [("Clinical variables", clinical_complete), ("Clinical + TMB", tmb_complete)]:
        for status, mask in [("Complete", complete), ("Incomplete", ~complete)]:
            rates.append({"analysis": name, "status": status,
                          "response_rate": data.loc[mask, "response"].mean(), "n": int(mask.sum())})
    rates = pd.DataFrame(rates)
    sns.barplot(data=rates, x="analysis", y="response_rate", hue="status", ax=ax,
                palette=["#4E79A7", "#F28E2B"])
    for patch in ax.patches:
        if np.isfinite(patch.get_height()) and patch.get_height() > 0:
            ax.text(patch.get_x() + patch.get_width() / 2, patch.get_height() + .01,
                    f"{patch.get_height():.1%}", ha="center", va="bottom", fontsize=8)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel("Objective response rate")
    ax.set_xlabel("")
    ax.set_title("Outcome-related missingness")
    ax.legend(frameon=False, title="Analysis set")
    clean_axis(ax, "y")
    save_figure(fig, out / "FigureS6_IMvigor210_internal_validation")


def figure_s7(root, out):
    associations = pd.read_csv(
        root / "results/protein_validation/CPTAC_ITGB1_RNA_protein_correlations.tsv", sep="\t")
    meta = pd.read_csv(
        root / "results/protein_validation/CPTAC_ITGB1_RNA_protein_meta.tsv", sep="\t").iloc[0]
    frame = pd.read_csv(
        root / "data/processed/CPTAC_ITGB1/CPTAC_ITGB1_RNA_protein.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)

    ax = axes[0]
    panel_label(ax, "A")
    data = associations.sort_values("spearman_rho")
    labels = [name.replace("_cptac_", " ").replace("_", " ").upper()
              for name in data.study] + ["RANDOM-EFFECTS POOLED"]
    estimates = data.spearman_rho.tolist() + [meta.meta_rho]
    intervals = [correlation_interval(r.spearman_rho, r.n) for r in data.itertuples()]
    intervals.append((meta.meta_ci_low, meta.meta_ci_high))
    y = np.arange(len(labels))[::-1]
    for i, (estimate, interval) in enumerate(zip(estimates, intervals)):
        marker = "D" if i == len(labels) - 1 else "o"
        color = "#222222" if i == len(labels) - 1 else "#59A14F"
        ax.errorbar(estimate, y[i], xerr=[[estimate - interval[0]], [interval[1] - estimate]],
                    fmt=marker, color=color, capsize=2)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("RNA-protein Spearman rho (95% CI)")
    ax.set_title("CPTAC cohort-specific correspondence")
    clean_axis(ax, "x")

    ax = axes[1]
    panel_label(ax, "B")
    frame["RNA_z"] = frame.groupby("study").ITGB1_RNA.transform(
        lambda x: (x - x.mean()) / x.std(ddof=0))
    frame["protein_z"] = frame.groupby("study").ITGB1_protein.transform(
        lambda x: (x - x.mean()) / x.std(ddof=0))
    sns.regplot(data=frame, x="RNA_z", y="protein_z", scatter_kws={"s": 10, "alpha": .35},
                line_kws={"color": "#222222", "linewidth": 1.2}, ax=ax)
    ax.set_xlabel("Within-cohort standardized ITGB1 RNA")
    ax.set_ylabel("Within-cohort standardized ITGB1 protein")
    ax.set_title("Pooled standardized samples")
    clean_axis(ax, None)
    save_figure(fig, out / "FigureS7_CPTAC_ITGB1_RNA_protein")


def main():
    root, _ = parse_project_root()
    out = root / "figures"
    out.mkdir(parents=True, exist_ok=True)
    setup_plot_style()
    figure_s4(root, out)
    figure_s5(root, out)
    figure_s6(root, out)
    figure_s7(root, out)


if __name__ == "__main__":
    main()
