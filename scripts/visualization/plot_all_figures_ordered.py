"""Generate all publication figures in manuscript order."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
import seaborn as sns

from scripts.core.cli import parse_project_root
from scripts.core.common import ITGB1_FEATURES
from scripts.core.plotting import (
    COHORT_COLORS,
    COHORT_ORDER,
    FEATURE_LABELS as FEATURE_LABEL_MAP,
    clean_axis as clean,
    load_composition_results as load_composition,
    panel_label,
    save_figure as save,
    setup_plot_style,
)


COLORS = COHORT_COLORS
FEATURES = ITGB1_FEATURES
FEATURE_LABELS = [FEATURE_LABEL_MAP[feature] for feature in FEATURES]
MODEL_ORDER = ["CAF_ECM", "Endothelial", "StromalVascular", "FullComposition"]
MODEL_LABELS = ["CAF/ECM", "Endothelial", "Stromal/vascular", "Full composition"]


def panel(ax, letter):
    panel_label(ax, letter, x=-0.13)


def figure1(root, out, source):
    sc = pd.read_csv(root / "results/single_cell/pretreatment_celltype_ranking.tsv", sep="\t")
    sc = sc[sc.cell_type != "Unclassified"].copy()
    paired = pd.read_csv(root / "results/single_cell/paired_pre_post_ITGB1.tsv", sep="\t")
    patient = pd.read_csv(root / "results/single_cell/patient_celltype_ITGB1.tsv", sep="\t")
    patient = patient[patient.cell_type != "Unclassified"]
    sc.to_csv(source / "Figure1_single_cell_source.tsv", sep="\t", index=False)
    paired.to_csv(source / "Figure1_paired_source.tsv", sep="\t", index=False)
    fig = plt.figure(figsize=(11, 7.3), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[.8, 1.25])

    ax = fig.add_subplot(gs[0, :]); ax.axis("off"); panel(ax, "A")
    boxes = [(0.01, "Single-cell source\nHNSCC, 4 paired patients", "#D9EAF7"),
             (0.205, "Bulk composition\n6 cohorts, n=311", "#DDF1E8"),
             (0.40, "Clinical response\n5 ICI cohorts, n=477", "#F5E3D3"),
             (0.595, "External pathways\n3 solid-tumour cohorts", "#F8E7A8"),
             (0.79, "Pharmacogenomics\nlineage-aware GDSC", "#EDE2F2")]
    for x, txt, color in boxes:
        patch = FancyBboxPatch((x, .27), .17, .45, boxstyle="round,pad=.012,rounding_size=.02",
                               facecolor=color, edgecolor="#555555", lw=.9, transform=ax.transAxes)
        ax.add_patch(patch); ax.text(x+.085, .495, txt, ha="center", va="center", fontsize=8.3,
                                    transform=ax.transAxes)
    for x in [.18, .375, .57, .765]:
        ax.add_patch(FancyArrowPatch((x, .495), (x+.025, .495), arrowstyle="-|>", mutation_scale=11,
                                    lw=1, color="#555555", transform=ax.transAxes))
    ax.text(.5, .90, "Analysis framework: from cellular origin to clinical relevance",
            ha="center", va="center", fontsize=11, fontweight="bold", transform=ax.transAxes)

    ax = fig.add_subplot(gs[1, 0]); panel(ax, "B")
    d = sc.sort_values("mean_ITGB1")
    y = np.arange(len(d))
    ax.scatter(d.mean_ITGB1, y, s=45+270*d.detected_fraction, c=d.detected_fraction,
               cmap="viridis", vmin=.35, vmax=.85, edgecolor="white", lw=.7, zorder=3)
    labels = d.cell_type.replace({"Fibroblast_CAF":"Fibroblast/CAF", "Endothelial":"Endothelial",
                                  "Epithelial_malignant_candidate":"Malignant epithelial",
                                  "T_NK":"T/NK", "B_Plasma":"B/plasma"})
    ax.set_yticks(y, labels); ax.set_xlabel("Mean normalized ITGB1 expression")
    ax.set_title("Pretreatment cellular sources of ITGB1")
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(.35, .85), cmap="viridis")
    cb = fig.colorbar(sm, ax=ax, fraction=.05, pad=.03); cb.set_label("Detected fraction", fontsize=8)
    clean(ax, "x")

    ax = fig.add_subplot(gs[1, 1]); panel(ax, "C")
    cell_order = ["Fibroblast_CAF", "Endothelial", "Epithelial_malignant_candidate", "T_NK", "Myeloid", "B_Plasma"]
    display = ["Fibroblast/CAF", "Endothelial", "Malignant epithelial", "T/NK", "Myeloid", "B/plasma"]
    for j, ct in enumerate(cell_order):
        d = patient[patient.cell_type.eq(ct)]
        wide = d.pivot(index="patient_id", columns="timepoint", values="mean_ITGB1").dropna()
        x0, x1 = j-.14, j+.14
        for _, r in wide.iterrows():
            ax.plot([x0, x1], [r.Pre, r.Post], color="#999999", lw=.8, alpha=.75)
            ax.scatter([x0], [r.Pre], color="#0072B2", s=18, zorder=3)
            ax.scatter([x1], [r.Post], color="#D55E00", s=18, zorder=3)
    ax.set_xticks(np.arange(len(cell_order)), display, rotation=25, ha="right")
    ax.set_ylabel("Patient-level mean ITGB1 expression")
    ax.set_title("Paired pre/post-treatment changes")
    ax.legend(handles=[Line2D([0],[0], marker='o', color='none', markerfacecolor='#0072B2', label='Pretreatment'),
                       Line2D([0],[0], marker='o', color='none', markerfacecolor='#D55E00', label='Post-treatment')],
              title="Timepoint", frameon=False, ncol=2, loc="upper right")
    clean(ax, "y")
    save(fig, out / "Figure1_study_design_and_single_cell_origin")


def figure2(root, out, source):
    boot, cv = load_composition(root)
    alt = pd.read_csv(root / "results/robustness/alternative_method_composition_r2.tsv", sep="\t")
    conc = pd.read_csv(root / "results/robustness/MCPcounter_signature_concordance.tsv", sep="\t")
    boot.to_csv(source / "Figure2_bootstrap_source.tsv", sep="\t", index=False)
    cv.to_csv(source / "Figure2_cv_source.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

    ax = axes[0,0]; panel(ax,"A")
    d = boot[(boot.estimand.eq("r_squared")) & boot.model.eq("FullComposition")].set_index("cohort").loc[COHORT_ORDER].reset_index()
    y = np.arange(len(d))[::-1]
    for i,r in d.iterrows():
        ax.errorbar(r.estimate,y[i],xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt='o',
                    color=COLORS[r.cohort],capsize=2,ms=5)
    ns = cv[cv.model.eq("FullComposition")].set_index("cohort").n
    ax.set_yticks(y,[f"{r.cohort} (n={int(ns[r.cohort])})" for r in d.itertuples()])
    ax.set_xlabel("Full-composition in-sample $R^2$ (bootstrap 95% interval)")
    ax.set_xlim(0,1.02); ax.set_title("Composition dependence across six cohorts"); clean(ax,"x")

    ax=axes[0,1]; panel(ax,"B")
    d=cv[cv.model.eq("FullComposition")].set_index("cohort").loc[COHORT_ORDER].reset_index(); y=np.arange(len(d))[::-1]
    for i,r in d.iterrows():
        ax.errorbar(r.repeated_5fold_mean_r2,y[i],xerr=[[r.repeated_5fold_mean_r2-r.repeated_5fold_ci_low],
                    [r.repeated_5fold_ci_high-r.repeated_5fold_mean_r2]],fmt='o',color=COLORS[r.cohort],capsize=2,ms=5)
    ax.axvline(0,color="#555555",ls="--",lw=.8); ax.set_yticks(y,[f"{r.cohort} (n={int(r.n)})" for r in d.itertuples()])
    ax.set_xlabel("Repeated 5-fold cross-validated $R^2$"); ax.set_title("Out-of-sample generalizability"); clean(ax,"x")

    ax=axes[1,0]; panel(ax,"C")
    d=alt[alt.method.eq("MCPcounter")].set_index("cohort").loc[COHORT_ORDER].reset_index()
    ax.bar(np.arange(len(d)),d.r_squared,color=[COLORS[x] for x in d.cohort],width=.68,zorder=3)
    ax.set_xticks(np.arange(len(d)),d.cohort,rotation=25,ha="right"); ax.set_ylabel("MCP-counter model $R^2$")
    ax.set_title("Independent deconvolution sensitivity analysis"); clean(ax,"y")

    ax=axes[1,1]; panel(ax,"D")
    pops=["Fibroblasts","Endothelial cells","T cells","Monocytic lineage","Cytotoxic lymphocytes"]
    mat=conc.pivot(index="MCPcounter_population",columns="cohort",values="spearman_rho").reindex(index=pops,columns=COHORT_ORDER)
    sns.heatmap(mat,ax=ax,cmap="Blues",vmin=0,vmax=1,annot=True,fmt=".2f",linewidths=.5,
                cbar_kws={"label":"Spearman ρ", "shrink":.8})
    ax.set_xlabel("Cohort"); ax.set_ylabel("MCP-counter population"); ax.set_title("Agreement with internal signatures")
    ax.set_xticklabels(ax.get_xticklabels(),rotation=25,ha="right")
    save(fig,out/"Figure2_bulk_composition_and_robustness")


def figure3(root,out,source):
    response=pd.read_csv(root/"results/clinical_response/five_ici_cohort_response_results.tsv",sep="\t")
    meta=pd.read_csv(root/"results/clinical_response/random_effects_meta_analysis.tsv",sep="\t")
    pfs=pd.read_csv(root/"results/clinical_response/GSE135222_pfs_results.tsv",sep="\t")
    rect=pd.read_csv(root/"results/clinical_response/rectal_chemoradiotherapy_results_v3.tsv",sep="\t")
    response.to_csv(source/"Figure3_response_source.tsv",sep="\t",index=False,na_rep="NA")
    meta.to_csv(source/"Figure3_meta_source.tsv",sep="\t",index=False)
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)

    ax=axes[0,0];panel(ax,"A")
    cohort_order=["GSE91061","GSE78220","GSE176307","GSE126044","IMvigor210"]
    d=response[response.feature.eq("ITGB1")].set_index("cohort").loc[cohort_order].reset_index()
    pooled=meta[meta.feature.eq("ITGB1")].iloc[0]; labels=[f"{r.cohort} (n={r.n})" for r in d.itertuples()]+["Random-effects pooled"]
    est=np.r_[d.odds_ratio_per_sd,pooled.pooled_or];lo=np.r_[d.ci_low,pooled.ci_low];hi=np.r_[d.ci_high,pooled.ci_high];y=np.arange(6)[::-1]
    for i,r in d.iterrows(): ax.errorbar(est[i],y[i],xerr=[[est[i]-lo[i]],[hi[i]-est[i]]],fmt='o',color=COLORS[r.cohort],capsize=2)
    ax.errorbar(est[-1],y[-1],xerr=[[est[-1]-lo[-1]],[hi[-1]-est[-1]]],fmt='D',color="#222222",capsize=3)
    ax.axvline(1,color="#555",ls="--",lw=.8);ax.set_xscale("log");ax.set_yticks(y,labels)
    ax.set_xlabel("Response OR per 1-SD ITGB1 (95% CI)");ax.set_title("Unadjusted ITGB1 and objective ICI response");clean(ax,"x")

    ax=axes[0,1];panel(ax,"B")
    d=meta.set_index("feature").loc[FEATURES].reset_index();y=np.arange(5)[::-1]
    ax.errorbar(d.pooled_or,y,xerr=[d.pooled_or-d.ci_low,d.ci_high-d.pooled_or],fmt='D',color="#333333",capsize=3)
    ax.axvline(1,color="#555",ls="--",lw=.8);ax.set_xscale("log");ax.set_yticks(y,FEATURE_LABELS)
    ax.set_xlabel("Random-effects pooled response OR (95% CI)");ax.set_title("Effect of composition adjustment");clean(ax,"x")

    ax=axes[1,0];panel(ax,"C")
    cohorts=["GSE91061","GSE78220","GSE176307","GSE126044","IMvigor210"]; offsets=np.linspace(-.28,.28,5)
    for off,c in zip(offsets,cohorts):
        d=response[response.cohort.eq(c)].set_index("feature").loc[FEATURES]
        ax.plot(np.arange(5)+off,d.auc_discrimination,'o-',color=COLORS[c],lw=1,ms=4,label=c)
    ax.axhline(.5,color="#555",ls="--",lw=.8);ax.set_xticks(np.arange(5),["Unadjusted","CAF/ECM","Endothelial","Stromal/vascular","Full"],rotation=25,ha="right")
    ax.set_ylabel("Direction-free AUC");ax.set_ylim(.45,.82);ax.set_title("Cohort-specific discrimination")
    ax.legend(title="ICI cohort",frameon=False,ncol=2,loc="upper left");clean(ax,"y")

    ax=axes[1,1];panel(ax,"D")
    d=pfs.set_index("feature").loc[FEATURES].reset_index();y=np.arange(5)[::-1]
    ax.errorbar(d.hazard_ratio_per_sd,y,xerr=[d.hazard_ratio_per_sd-d.ci_low,d.ci_high-d.hazard_ratio_per_sd],fmt='o',color=COLORS["GSE135222"],capsize=2,label="GSE135222 PFS")
    rd=rect[rect.feature.isin(FEATURES)].set_index("feature").loc[FEATURES].reset_index()
    # Rectal response is shown as OR, offset for readability; it is a distinct endpoint.
    ax.errorbar(rd.odds_ratio_per_sd,y-.16,xerr=[rd.odds_ratio_per_sd-rd.ci_low,rd.ci_high-rd.odds_ratio_per_sd],fmt='s',color=COLORS["GSE190826"],capsize=2,label="GSE190826 pCR")
    ax.axvline(1,color="#555",ls="--",lw=.8);ax.set_xscale("log");ax.set_yticks(y,FEATURE_LABELS)
    ax.set_xlabel("Effect estimate per 1-SD score (95% CI)");ax.set_title("Independent clinical endpoints")
    ax.legend(title="Endpoint",frameon=False,loc="center left",bbox_to_anchor=(1.02,0.5));clean(ax,"x")
    save(fig,out/"Figure3_treatment_response_and_meta_analysis")


def figure_s1(root,out,source):
    boot,cv=load_composition(root);alt=pd.read_csv(root/"results/robustness/alternative_method_composition_r2.tsv",sep="\t")
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True);offsets=np.linspace(-.22,.22,6)
    ax=axes[0,0];panel(ax,"A")
    for off,c in zip(offsets,COHORT_ORDER):
        d=boot[(boot.cohort.eq(c))&boot.estimand.eq("r_squared")].set_index("model").loc[MODEL_ORDER]
        x=np.arange(4)+off;ax.errorbar(x,d.estimate,yerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o',color=COLORS[c],capsize=1.5,ms=4,label=c)
    ax.set_xticks(np.arange(4),MODEL_LABELS,rotation=25,ha="right");ax.set_ylabel("In-sample $R^2$ (bootstrap 95% interval)");ax.set_title("All composition models");clean(ax,"y")
    ax.legend(title="Cohort",frameon=False,ncol=3,loc="upper left")
    ax=axes[0,1];panel(ax,"B")
    for off,c in zip(offsets,COHORT_ORDER):
        d=cv[cv.cohort.eq(c)].set_index("model").loc[MODEL_ORDER];x=np.arange(4)+off
        ax.errorbar(x,d.repeated_5fold_mean_r2,yerr=[d.repeated_5fold_mean_r2-d.repeated_5fold_ci_low,d.repeated_5fold_ci_high-d.repeated_5fold_mean_r2],fmt='o',color=COLORS[c],capsize=1.5,ms=4,label=c)
    ax.axhline(0,color="#555",ls="--",lw=.8);ax.set_xticks(np.arange(4),MODEL_LABELS,rotation=25,ha="right");ax.set_ylabel("Repeated 5-fold cross-validated $R^2$");ax.set_title("Model-specific generalizability");clean(ax,"y")
    ax.legend(title="Cohort",frameon=False,ncol=3,loc="lower left")
    ax=axes[1,0];panel(ax,"C")
    inc=boot[boot.estimand.eq("incremental_r_squared")];mods=["Endothelial_beyond_CAF_ECM","OtherComposition_beyond_StromalVascular"]
    for off,c in zip(offsets,COHORT_ORDER):
        d=inc[inc.cohort.eq(c)].set_index("model").loc[mods];x=np.arange(2)+off
        ax.errorbar(x,d.estimate,yerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o',color=COLORS[c],capsize=1.5,ms=4,label=c)
    ax.set_xticks(np.arange(2),["Endothelial beyond\nCAF/ECM","Other composition beyond\nstromal/vascular"]);ax.set_ylabel("Incremental $R^2$ (bootstrap 95% interval)");ax.set_title("Incremental composition contribution");clean(ax,"y")
    ax.legend(title="Cohort",frameon=False,ncol=3,loc="upper left")
    ax=axes[1,1];panel(ax,"D")
    methods=["MCPcounter","InternalRankSignatures"];w=.12
    for i,c in enumerate(COHORT_ORDER):
        d=alt[alt.cohort.eq(c)].set_index("method").loc[methods]
        ax.bar(np.arange(2)+(i-2.5)*w,d.r_squared,width=w,color=COLORS[c],label=c,zorder=3)
    ax.set_xticks(np.arange(2),["MCP-counter","Internal matched\nsignatures"]);ax.set_ylabel("Composition model $R^2$");ax.set_title("Alternative-method comparison");clean(ax,"y")
    ax.legend(title="Cohort",frameon=False,ncol=3,loc="upper left")
    boot.to_csv(source/"FigureS1_bootstrap_source.tsv",sep="\t",index=False);cv.to_csv(source/"FigureS1_CV_source.tsv",sep="\t",index=False)
    save(fig,out/"FigureS1_detailed_composition_robustness")


def figure_s2(root,out,source):
    auc=pd.read_csv(root/"results/robustness/extended_bootstrap_auc_and_differences.tsv",sep="\t")
    auc.to_csv(source/"FigureS2_AUC_bootstrap_source.tsv",sep="\t",index=False)
    fig,axes=plt.subplots(2,2,figsize=(10.5,8),constrained_layout=True);ici=["GSE91061","GSE78220","GSE176307","GSE126044"];offs=np.linspace(-.21,.21,4)
    ax=axes[0,0];panel(ax,"A")
    for off,c in zip(offs,ici):
        d=auc[(auc.cohort.eq(c))&auc.estimand.eq("auc_oriented_to_original_ITGB1")].set_index("feature").loc[FEATURES]
        x=np.arange(5)+off;ax.errorbar(x,d.estimate,yerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o-',color=COLORS[c],capsize=2,ms=4,label=c)
    ax.axhline(.5,color="#555",ls="--",lw=.8);ax.set_xticks(np.arange(5),["Unadjusted","CAF/ECM","Endothelial","Stromal/vascular","Full"],rotation=25,ha="right")
    ax.set_ylabel("AUC oriented to original ITGB1 direction");ax.set_title("Bootstrap AUC estimates");ax.legend(title="Cohort/endpoint",frameon=False,loc="lower left");clean(ax,"y")
    ax=axes[0,1];panel(ax,"B")
    deltas=auc[auc.estimand.eq("auc_difference_vs_unadjusted")];adj=FEATURES[1:]
    for off,c in zip(offs,ici):
        d=deltas[deltas.cohort.eq(c)].set_index("feature").loc[adj];y=np.arange(4)+off
        ax.errorbar(d.estimate,y,xerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o',color=COLORS[c],capsize=2,ms=4,label=c)
    ax.axvline(0,color="#555",ls="--",lw=.8);ax.set_yticks(np.arange(4),FEATURE_LABELS[1:]);ax.set_xlabel("AUC difference versus unadjusted ITGB1 (95% interval)")
    ax.set_title("Paired bootstrap AUC change");ax.legend(title="Cohort/endpoint",frameon=False,loc="lower right");clean(ax,"x")
    ax=axes[1,0];panel(ax,"C")
    d=auc[(auc.cohort.eq("GSE190826"))&auc.estimand.eq("auc_oriented_to_original_ITGB1")].set_index("feature").loc[FEATURES]
    ax.errorbar(np.arange(5),d.estimate,yerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o-',color=COLORS["GSE190826"],capsize=2,ms=4,label="GSE190826 pCR")
    ax.axhline(.5,color="#555",ls="--",lw=.8);ax.set_xticks(np.arange(5),["Unadjusted","CAF/ECM","Endothelial","Stromal/vascular","Full"],rotation=25,ha="right")
    ax.set_ylabel("AUC oriented to original ITGB1 direction");ax.set_title("Independent chemoradiotherapy endpoint");ax.legend(title="Endpoint",frameon=False);clean(ax,"y")
    ax=axes[1,1];panel(ax,"D")
    d=auc[(auc.cohort.eq("GSE190826"))&auc.estimand.eq("auc_difference_vs_unadjusted")].set_index("feature").loc[FEATURES[1:]];y=np.arange(4)
    ax.errorbar(d.estimate,y,xerr=[d.estimate-d.ci_low,d.ci_high-d.estimate],fmt='o',color=COLORS["GSE190826"],capsize=2,ms=4,label="GSE190826 pCR")
    ax.axvline(0,color="#555",ls="--",lw=.8);ax.set_yticks(y,FEATURE_LABELS[1:]);ax.set_xlabel("AUC difference versus unadjusted ITGB1 (95% interval)");ax.set_title("Chemoradiotherapy AUC change");ax.legend(title="Endpoint",frameon=False);clean(ax,"x")
    save(fig,out/"FigureS2_bootstrap_AUC_sensitivity")


def main():
    root, _ = parse_project_root()
    out = root / "figures"
    source = root / "results/figure_source_data"
    out.mkdir(parents=True, exist_ok=True)
    source.mkdir(parents=True, exist_ok=True)
    setup_plot_style(style="ticks", bold_titles=True)
    figure1(root, out, source)
    figure2(root, out, source)
    figure3(root, out, source)
    figure_s1(root, out, source)
    figure_s2(root, out, source)


if __name__ == "__main__":
    main()
