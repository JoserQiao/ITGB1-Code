#!/usr/bin/env python3
"""Replot Figure 6 and Supplementary Figure S6 from included source tables."""

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/itgb1-matplotlib-cache")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "figure6"
OUT = ROOT / "outputs"
DATASETS = {"GSE71340": "Ovarian cancer", "GSE240671": "Breast cancer", "GSE212810": "Bladder cancer"}
ORDER = list(DATASETS)
COLORS = {"ITGB1-low": "#4C78A8", "ITGB1-high": "#D64F3A"}


def stars(q):
    if not np.isfinite(q): return ""
    if q < .001: return "***"
    if q < .01: return "**"
    if q < .05: return "*"
    return ""


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def load_tables():
    corr = pd.read_csv(DATA / "Figure6_selected_gene_correlations.tsv", sep="\t")
    panel = pd.read_csv(DATA / "Figure6_GEO_derived_top10_panel.tsv", sep="\t")
    scores = pd.read_csv(DATA / "Figure6_panel_sample_scores.tsv", sep="\t")
    statistics = pd.read_csv(DATA / "Figure6_high_low_statistics.tsv", sep="\t")
    return corr, panel, scores, statistics


def plot_main(corr, panel, scores, statistics):
    genes = panel["gene"].tolist()
    sns.set_theme(style="white", font="DejaVu Sans")
    fig = plt.figure(figsize=(12.4, 7.4))
    grid = fig.add_gridspec(2, 3, width_ratios=[1.05, 1.15, 1.0], height_ratios=[1, 1.05], wspace=.55, hspace=.55)
    ax = fig.add_subplot(grid[0, 0:2])
    display = panel.sort_values("pooled_spearman_rho", ascending=True)
    y = np.arange(len(display))
    for ds, color, offset in zip(ORDER, ["#4C78A8", "#8E8E8E", "#D64F3A"], [-.22, 0, .22]):
        ax.scatter(display[ds], y + offset, color=color, s=28, label=ds, zorder=3)
    ax.scatter(display["pooled_spearman_rho"], y, color="#202020", marker="D", s=42, label="Pooled", zorder=4)
    ax.set_yticks(y); ax.set_yticklabels(display["gene"], fontsize=8.5); ax.set_xlabel("Spearman correlation with ITGB1"); ax.axvline(0, color="#777777", lw=.8); ax.set_xlim(0,1); ax.legend(frameon=False, ncol=4, fontsize=7.5, loc="lower right"); style(ax)
    ax = fig.add_subplot(grid[0,2])
    matrix = corr.pivot(index="dataset", columns="gene", values="spearman_rho").reindex(index=ORDER, columns=genes)
    sns.heatmap(matrix, ax=ax, cmap="vlag", center=0, vmin=-.1, vmax=.9, linewidths=.4, linecolor="white", cbar_kws={"label":"Spearman rho", "shrink":.78})
    for yi, ds in enumerate(ORDER):
        for xi, gene in enumerate(genes):
            q = corr.loc[(corr.dataset == ds) & (corr.gene == gene), "fdr"].iloc[0]
            ax.text(xi+.5, yi+.5, stars(q), ha="center", va="center", fontsize=7.5)
    ax.set_xlabel(""); ax.set_ylabel(""); ax.tick_params(axis="x", rotation=60, labelsize=7); ax.tick_params(axis="y", rotation=0, labelsize=7)
    ax = fig.add_subplot(grid[1,0:2])
    sns.violinplot(data=scores, x="dataset", y="panel_score", hue="ITGB1_group", order=ORDER, hue_order=["ITGB1-low","ITGB1-high"], split=True, inner=None, palette=COLORS, linewidth=.8, ax=ax)
    sns.stripplot(data=scores, x="dataset", y="panel_score", hue="ITGB1_group", order=ORDER, hue_order=["ITGB1-low","ITGB1-high"], dodge=True, palette=COLORS, size=2, alpha=.38, ax=ax)
    handles, labels = ax.get_legend_handles_labels(); ax.legend(handles[:2], labels[:2], frameon=False, fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(.5,-.16))
    ax.axhline(0, color="#777777", lw=.7, ls="--"); ax.set_xlabel(""); ax.set_ylabel("Mean z-score of selected genes")
    ymax = scores.panel_score.max()
    for x, ds in enumerate(ORDER):
        row = statistics[(statistics.dataset == ds) & (statistics.gene == "__PANEL__")].iloc[0]
        count = scores[scores.dataset == ds].ITGB1_group.value_counts()
        ax.text(x, ymax+.18, f"n={count.get('ITGB1-low',0)}/{count.get('ITGB1-high',0)}\nq={row.fdr:.2g}", ha="center", va="bottom", fontsize=7)
    ax.set_ylim(scores.panel_score.min()-.25, ymax+.65); style(ax)
    ax = fig.add_subplot(grid[1,2])
    panel_rows = statistics[statistics.gene == "__PANEL__"].set_index("dataset").reindex(ORDER)
    values = panel_rows["median_high_minus_low"].to_numpy(); ypos = np.arange(len(ORDER))
    ax.hlines(ypos, 0, values, color="#B8B8B8", lw=1.2); ax.scatter(values, ypos, s=145, color="#D64F3A", edgecolor="white", zorder=3); ax.axvline(0, color="#333333", lw=.8)
    for y0, value in zip(ypos, values): ax.text(value+.035, y0, f"{value:.2f}", va="center", fontsize=8)
    ax.set_yticks(ypos); ax.set_yticklabels([f"{ds}\n{DATASETS[ds]}" for ds in ORDER], fontsize=8); ax.set_xlabel("Median score difference\nITGB1-high minus ITGB1-low"); ax.set_xlim(-.05, max(values)+.35); style(ax)
    fig.savefig(OUT / "Figure6.png", dpi=350, bbox_inches="tight"); fig.savefig(OUT / "Figure6.pdf", bbox_inches="tight"); plt.close(fig)


def plot_supplement(corr, panel, scores, statistics):
    genes = panel["gene"].tolist(); sns.set_theme(style="white", font="DejaVu Sans")
    fig = plt.figure(figsize=(12.8,9)); grid = fig.add_gridspec(2,2,wspace=.35,hspace=.5)
    ax = fig.add_subplot(grid[0,0]); display = panel.sort_values("pooled_spearman_rho", ascending=True); y = np.arange(len(display)); ax.barh(y, display.pooled_spearman_rho, color="#4C78A8", alpha=.88); ax.scatter(display.minimum_cohort_rho, y, color="#D64F3A", s=28, zorder=3, label="Minimum cohort rho"); ax.set_yticks(y); ax.set_yticklabels([f"{rank}. {gene}" for rank,gene in zip(display["rank"],display.gene)],fontsize=8); ax.set_xlabel("Pooled rho (bar) / minimum cohort rho (dot)"); ax.set_xlim(0,.8); ax.legend(frameon=False,fontsize=7.5,loc="lower right"); style(ax)
    ax = fig.add_subplot(grid[0,1]); matrix = corr.pivot(index="gene",columns="dataset",values="spearman_rho").reindex(index=genes,columns=ORDER); sns.heatmap(matrix,ax=ax,cmap="vlag",center=0,vmin=-.1,vmax=.9,linewidths=.4,linecolor="white",cbar_kws={"label":"Spearman rho","shrink":.76}); ax.set_xlabel("");ax.set_ylabel("");ax.tick_params(axis="x",rotation=30,labelsize=8);ax.tick_params(axis="y",rotation=0,labelsize=8)
    ax = fig.add_subplot(grid[1,0]); sns.boxplot(data=scores,x="dataset",y="panel_score",hue="ITGB1_group",order=ORDER,hue_order=["ITGB1-low","ITGB1-high"],palette=COLORS,fliersize=1.5,linewidth=.8,ax=ax); handles,labels=ax.get_legend_handles_labels();ax.legend(handles,labels,frameon=False,fontsize=8,ncol=2,loc="upper center",bbox_to_anchor=(.5,-.17));ax.axhline(0,color="#777777",lw=.7,ls="--");ax.set_xlabel("");ax.set_ylabel("Mean z-score of selected genes");style(ax)
    ax = fig.add_subplot(grid[1,1]); diff=statistics[statistics.gene!="__PANEL__"].pivot(index="gene",columns="dataset",values="median_high_minus_low").reindex(index=genes,columns=ORDER);sns.heatmap(diff,ax=ax,cmap="vlag",center=0,vmin=-1.2,vmax=1.2,linewidths=.4,linecolor="white",cbar_kws={"label":"Median log2(TPM+1) difference","shrink":.76});ax.set_xlabel("");ax.set_ylabel("");ax.tick_params(axis="x",rotation=30,labelsize=8);ax.tick_params(axis="y",rotation=0,labelsize=8)
    fig.savefig(OUT / "Supplementary_FigureS6.png",dpi=350,bbox_inches="tight");fig.savefig(OUT / "Supplementary_FigureS6.pdf",bbox_inches="tight");plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    plot_main(*load_tables())
    plot_supplement(*load_tables())
    print("Wrote Figure6 and Supplementary_FigureS6 to outputs/")
