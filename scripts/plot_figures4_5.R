#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(forcats)
})

root <- normalizePath(".")
d4 <- file.path(root, "data", "figure4")
d5 <- file.path(root, "data", "figure5")
out <- file.path(root, "outputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
theme_pub <- function(base_size = 8) theme_classic(base_size = base_size) + theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black"), plot.margin = margin(4,4,4,4))
star <- function(p) case_when(is.na(p) ~ "", p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "")
plabel <- function(p) ifelse(p < .001, "P<0.001", paste0("P=", signif(p, 2)))
save_plot <- function(p, name, w, h) { ggsave(file.path(out, paste0(name, ".png")), p, width=w, height=h, dpi=450); ggsave(file.path(out, paste0(name, ".pdf")), p, width=w, height=h) }

# Figure 4 --------------------------------------------------------------------
model4 <- ggplot() +
  annotate("label",x=.10,y=.55,label="ITGB1-high",fontface="bold",size=3.7,fill="#f8d6d6") +
  annotate("label",x=.37,y=.55,label="Focal adhesion / ECM",fontface="bold",size=3.7,fill="#e8f2f8") +
  annotate("label",x=.64,y=.55,label="Lower chemo sensitivity",fontface="bold",size=3.7,fill="#fbefd9") +
  annotate("label",x=.86,y=.55,label="Poor prognosis",fontface="bold",size=3.7,fill="#e0f0df") +
  annotate("segment",x=.19,xend=.28,y=.55,yend=.55,arrow=arrow(length=grid::unit(.09,"in"))) +
  annotate("segment",x=.49,xend=.56,y=.55,yend=.55,arrow=arrow(length=grid::unit(.09,"in"))) +
  annotate("segment",x=.77,xend=.80,y=.55,yend=.55,arrow=arrow(length=grid::unit(.09,"in"))) +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+theme_void()

cox <- read_tsv(file.path(d4,"GSE96058_ITGB1_OS_Cox.tsv"),show_col_types=FALSE) %>% mutate(model=ifelse(adjusted,"Adjusted","Univariable"),row=paste0(cohort," (",model,")"),highlight=cohort=="Chemo-treated",row=fct_rev(factor(row,levels=row)),label=paste0("n=",n,", events=",events,", ",plabel(p)))
p4b <- ggplot(cox,aes(hr,row))+geom_vline(xintercept=1,linetype=2,color="grey55")+geom_errorbarh(aes(xmin=ci_low,xmax=ci_high,color=highlight),height=0,linewidth=.45)+geom_point(aes(color=highlight),size=1.9)+geom_text(aes(label=label),x=max(cox$ci_high)*1.08,hjust=0,size=2.2)+scale_color_manual(values=c("FALSE"="grey45","TRUE"="#C93F4A"),guide="none")+scale_x_log10(limits=c(min(cox$ci_low)*.85,max(cox$ci_high)*1.85))+labs(x="Hazard ratio per 1 SD ITGB1",y=NULL)+theme_pub(8)

chemo <- read_tsv(file.path(d4,"Figure4E_GDSC_chemotherapy_dotplot_source.tsv"),show_col_types=FALSE) %>% filter(fdr<.05,spearman_rho>0) %>% arrange(desc(spearman_rho)) %>% slice_head(n=12) %>% arrange(spearman_rho) %>% mutate(drug_drug_name=fct_inorder(drug_drug_name))
p4c <- ggplot(chemo,aes(spearman_rho,drug_drug_name))+geom_vline(xintercept=0,linetype=2,color="grey55")+geom_point(aes(size=-log10(fdr),color=spearman_rho))+scale_color_gradient(low="#F6B26B",high="#B2182B")+scale_size_continuous(range=c(2.4,6))+labs(x="Spearman rho: ITGB1 vs ln(IC50)",y=NULL,color="rho",size="-log10(FDR)")+theme_pub(8)

ratio <- read_tsv(file.path(d4,"GDSC_chemo_ITGB1_high_low_sensitivity_summary.tsv"),show_col_types=FALSE) %>% filter(mann_whitney_fdr<.05,IC50_ratio_high_vs_low>1) %>% arrange(desc(IC50_ratio_high_vs_low)) %>% slice_head(n=10) %>% arrange(IC50_ratio_high_vs_low) %>% mutate(drug=fct_inorder(drug))
p4d <- ggplot(ratio,aes(IC50_ratio_high_vs_low,drug,fill=IC50_ratio_high_vs_low))+geom_vline(xintercept=1,linetype=2,color="grey35")+geom_col(width=.68)+geom_text(aes(label=paste0(round(IC50_ratio_high_vs_low,2),"x")),hjust=-.08,size=2.4)+scale_fill_gradient(low="#F6B26B",high="#B2182B",guide="none")+coord_cartesian(xlim=c(0,max(ratio$IC50_ratio_high_vs_low)*1.18))+labs(x="Median IC50 ratio: ITGB1-high / ITGB1-low",y=NULL)+theme_pub(8)

scatter <- read_tsv(file.path(d4,"Figure4F_GDSC_representative_scatter_source.tsv"),show_col_types=FALSE) %>% filter(drug_drug_name %in% c("Oxaliplatin","Paclitaxel","Gemcitabine")) %>% mutate(drug_drug_name=factor(drug_drug_name,levels=c("Oxaliplatin","Paclitaxel","Gemcitabine")))
ann <- scatter %>% group_by(drug_drug_name) %>% summarise(x=min(ITGB1_log2_TPM_plus_1),y=max(ln_ic50),label=paste0("rho=",signif(first(spearman_rho),2),"\nFDR=",signif(first(fdr),2)),.groups="drop")
p4e <- ggplot(scatter,aes(ITGB1_log2_TPM_plus_1,ln_ic50))+geom_point(color="#4C78A8",alpha=.32,size=.5)+geom_smooth(method="lm",color="#D95F02",linewidth=.7)+geom_label(data=ann,aes(x=x,y=y,label=label),inherit.aes=FALSE,hjust=0,vjust=1,size=2.1,label.size=.15,fill="white")+facet_wrap(~drug_drug_name,scales="free_y",nrow=1)+labs(x="ITGB1 expression, log2(TPM + 1)",y="ln(IC50)")+theme_pub(8)

pathway <- read_tsv(file.path(d4,"GEO_pathway_high_low_statistics.tsv"),show_col_types=FALSE) %>% mutate(Cancer=factor(Cancer,c("Gastric cancer","Pancreatic cancer","Liver cancer")),Pathway=factor(Pathway,rev(c("ECM / collagen","Integrin / focal adhesion","EMT","TGF-beta","Hypoxia"))))
p4f <- ggplot(pathway,aes(Cancer,Pathway,fill=Delta))+geom_tile(color="white",linewidth=.45)+geom_text(aes(label=star(FDR)),size=3)+scale_fill_gradient2(low="#3b6fb6",mid="white",high="#c43b3b",midpoint=0,name="High-low\nscore")+labs(x=NULL,y=NULL)+theme_pub(8)+theme(axis.text.x=element_text(angle=25,hjust=1))

integrin_genes <- c("ITGB1","ITGA5","ITGAV","ITGA2","ITGA3","ITGA6","ITGB3","ITGB5","PTK2","PXN","TLN1","VCL","ACTN1","CAV1","LAMA4")
ranks <- read_tsv(file.path(d4,"GEO_integrin_focal_adhesion_ranked_genes.tsv"),show_col_types=FALSE)
gsea <- read_tsv(file.path(d4,"GEO_integrin_focal_adhesion_fgsea.tsv"),show_col_types=FALSE)
curve_one <- function(cohort,show_y=FALSE) {
  d <- ranks %>% filter(Cohort==cohort) %>% arrange(Rank); hits <- d$Gene %in% integrin_genes; w <- abs(d$RankMetric); running <- cumsum(ifelse(hits,w/sum(w[hits]),0)-ifelse(!hits,1/sum(!hits),0)); stat <- gsea %>% filter(Cohort==cohort)
  ggplot(tibble(Rank=d$Rank,ES=running,Hit=hits),aes(Rank,ES))+geom_hline(yintercept=0,linetype=2,linewidth=.25)+geom_line(color="#32c73e",linewidth=.55)+geom_segment(data=function(x)filter(x,Hit),aes(x=Rank,xend=Rank,y=min(running)-.05,yend=min(running)-.015),inherit.aes=FALSE,linewidth=.18)+annotate("label",x=Inf,y=Inf,label=sprintf("%s\nNES = %.2f\nFDR = %.3g",stat$Cancer,stat$NES,stat$padj),hjust=1.04,vjust=1.15,size=2.1,label.size=.12,fill="white")+labs(x=NULL,y=if(show_y)"Enrichment score" else NULL)+theme_pub(7)+theme(axis.title=element_text(size=6.5),axis.text.x=element_text(size=5.5))
}
p4g <- (curve_one("GSE113255",TRUE)|curve_one("GSE172356")|curve_one("GSE202069"))+plot_annotation(caption="Genes ranked by ITGB1-high minus ITGB1-low mean expression")&theme(plot.caption=element_text(hjust=.5,size=6.5,margin=margin(t=2)))

scores <- read_tsv(file.path(d4,"GSE96058_signature_scores.tsv"),show_col_types=FALSE) %>% filter(signature %in% c("ECM_Collagen","Integrin_Focal_Adhesion","EMT","TGF_beta","Drug_Efflux_ABC","Hypoxia")) %>% mutate(signature=recode(signature,ECM_Collagen="ECM/collagen",Integrin_Focal_Adhesion="Integrin/focal adhesion",TGF_beta="TGF-beta",Drug_Efflux_ABC="Drug efflux"),signature=factor(signature,levels=c("ECM/collagen","Integrin/focal adhesion","EMT","TGF-beta","Drug efflux","Hypoxia")),ITGB1_group=factor(ITGB1_group,levels=c("ITGB1-low","ITGB1-high")))
p4h <- ggplot(scores,aes(signature,score,fill=ITGB1_group))+geom_violin(position=position_dodge(.78),width=.72,alpha=.55,linewidth=.25)+geom_boxplot(position=position_dodge(.78),width=.18,outlier.shape=NA,linewidth=.25)+scale_fill_manual(values=c("ITGB1-low"="#4C78A8","ITGB1-high"="#D65F5F"))+labs(x=NULL,y="Signature score",fill=NULL)+theme_pub(8)+theme(axis.text.x=element_text(angle=35,hjust=1))
fig4 <- (model4/(p4b|p4c)/(p4d|p4e)/(p4f|wrap_elements(full=p4g)|p4h))+plot_layout(heights=c(.45,1,1,1),widths=c(.9,1.6,1))+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(face="bold",size=16))
save_plot(fig4,"Figure4",12,12)

# Figure 5 --------------------------------------------------------------------
model5 <- ggplot()+annotate("label",x=.14,y=.60,label="ITGB1-high",fontface="bold",size=3.7,fill="#f8d6d6")+annotate("label",x=.46,y=.60,label="Immune exclusion",fontface="bold",size=3.7,fill="#e8f2f8")+annotate("label",x=.77,y=.60,label="Suppressive state",fontface="bold",size=3.7,fill="#e0f0df")+annotate("segment",x=.25,xend=.35,y=.60,yend=.60,arrow=arrow(length=grid::unit(.09,"in")))+annotate("segment",x=.57,xend=.66,y=.60,yend=.60,arrow=arrow(length=grid::unit(.09,"in")))+coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand=FALSE)+theme_void()

tcga <- read_tsv(file.path(d5,"TCGA_pancancer_ITGB1_immune_gene_correlations.tsv"),show_col_types=FALSE)
genes <- c("TGFB1","CD47","HAVCR2","CD274","IDO1","PDCD1","CTLA4","LAG3","TIGIT","ARG1")
heat <- tcga %>% filter(Gene %in% genes) %>% mutate(Cancer=factor(Cancer,levels=unique(Cancer)),Gene=factor(Gene,levels=rev(genes)))
p5b <- ggplot(heat,aes(Cancer,Gene,fill=Rho))+geom_tile(color="white",linewidth=.25)+geom_text(aes(label=star(FDR)),size=2)+scale_fill_gradient2(low="#3b6fb6",mid="white",high="#c43b3b",midpoint=0,name="Spearman rho")+labs(x=NULL,y=NULL)+theme_pub(7)+theme(axis.text.x=element_text(angle=60,hjust=1,size=6))
counts <- tcga %>% filter(Rho>0,FDR<.05) %>% count(Cancer,name="Positive_genes") %>% right_join(tibble(Cancer=unique(tcga$Cancer)),by="Cancer") %>% mutate(Positive_genes=replace_na(Positive_genes,0),Cancer=fct_reorder(Cancer,Positive_genes))
p5c <- ggplot(counts,aes(Cancer,Positive_genes))+geom_col(fill="#315f7d")+geom_hline(yintercept=5,linetype=2,color="grey55")+labs(x=NULL,y="Positive genes, FDR < 0.05")+theme_pub(8)+theme(axis.text.x=element_text(angle=60,hjust=1))
expr <- read_tsv(file.path(d5,"TCGA_pancancer_ITGB1_immune_genes_expression.tsv"),show_col_types=FALSE) %>% pivot_wider(names_from=gene,values_from=expression_log2_tpm_plus_1) %>% filter(!is.na(ITGB1))
reg_genes <- c("TGFB1","CD47","HAVCR2","CD274","CTLA4","TIGIT","ARG1")
reg <- bind_rows(lapply(reg_genes,function(g) { d <- expr %>% filter(!is.na(.data[[g]])) %>% mutate(outcome_z=as.numeric(scale(.data[[g]])),itgb1_z=as.numeric(scale(ITGB1)),project_id=factor(project_id)); ct <- summary(lm(outcome_z~itgb1_z+project_id,data=d))$coefficients["itgb1_z",]; tibble(Gene=g,Beta=ct[["Estimate"]],SE=ct[["Std. Error"]],P=ct[["Pr(>|t|)"]]) })) %>% mutate(FDR=p.adjust(P,"BH"),CI_low=Beta-1.96*SE,CI_high=Beta+1.96*SE,Gene=factor(Gene,rev(reg_genes))) %>% filter(Beta>0,FDR<.05)
p5d <- ggplot(reg,aes(Beta,Gene))+geom_vline(xintercept=0,linetype=2,color="grey55",linewidth=.3)+geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.18,color="#c43b3b",linewidth=.45)+geom_point(color="#c43b3b",size=2)+geom_text(aes(label=ifelse(FDR<.001,"FDR<0.001",sprintf("FDR=%.3f",FDR))),hjust=-.08,size=2.5)+scale_x_continuous(expand=expansion(mult=c(.05,.22)))+labs(x="Standardized beta for ITGB1",y=NULL)+theme_pub(8)
immune <- read_tsv(file.path(d5,"GEO_immune_exclusion_ITGB1_regression.tsv"),show_col_types=FALSE) %>% mutate(Cancer=factor(Cancer,rev(c("Gastric cancer","Pancreatic cancer","Liver cancer"))),label=ifelse(FDR<.001,"FDR < 0.001",sprintf("FDR = %.3f",FDR)))
p5e <- ggplot(immune,aes(Beta,Cancer))+geom_vline(xintercept=0,linetype=2,color="grey55",linewidth=.3)+geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.18,color="#c43b3b",linewidth=.55)+geom_point(color="#c43b3b",size=2.4)+geom_text(aes(label=label),hjust=-.08,size=2.8)+scale_x_continuous(expand=expansion(mult=c(.05,.27)))+labs(x="Standardized beta for ITGB1",y=NULL)+theme_pub(8)
curve <- read_tsv(file.path(d5,"GEO_immune_exclusion_running_enrichment.tsv"),show_col_types=FALSE) %>% mutate(Cancer=factor(Cancer,c("Gastric cancer","Pancreatic cancer","Liver cancer")))
p5f <- ggplot(curve,aes(Rank,RunningES))+geom_hline(yintercept=0,linetype=2,linewidth=.25,color="grey50")+geom_line(color="#c43b3b",linewidth=.55)+geom_segment(data=function(d)filter(d,Hit),aes(x=Rank,xend=Rank,y=min(RunningES)-.05,yend=min(RunningES)-.015),inherit.aes=FALSE,linewidth=.20,color="black")+geom_label(data=function(d)distinct(d,Cancer,ES),aes(x=Inf,y=Inf,label=sprintf("ES = %.2f",ES)),hjust=1.04,vjust=1.12,size=2.4,label.size=.15,fill="white")+facet_wrap(~Cancer,nrow=1,scales="free_x")+labs(x="Ranked genes: ITGB1-high versus ITGB1-low",y="Running enrichment score")+theme_pub(8)
fig5 <- (model5/(p5b|p5c)/(p5d|p5e)/p5f)+plot_layout(heights=c(.42,1,.82,.95))+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(face="bold",size=16))
save_plot(fig5,"Figure5",12,10.8)
