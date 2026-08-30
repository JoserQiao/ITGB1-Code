#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork); library(forcats) })
root <- normalizePath("."); d4 <- file.path(root,"data","figure4"); d5 <- file.path(root,"data","figure5"); out <- file.path(root,"outputs"); dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_pub <- function(base_size=8) theme_classic(base_size=base_size)+theme(axis.text=element_text(color="black"),axis.title=element_text(color="black"),plot.margin=margin(4,4,4,4))
star <- function(p) case_when(is.na(p)~"",p<.001~"***",p<.01~"**",p<.05~"*",TRUE~"")
save_plot <- function(p,n,w,h) { ggsave(file.path(out,paste0(n,".png")),p,width=w,height=h,dpi=450);ggsave(file.path(out,paste0(n,".pdf")),p,width=w,height=h) }

# Supplementary Figure S4: corrected three-cohort pathway results.
pathway <- read_tsv(file.path(d4,"GEO_pathway_high_low_statistics.tsv"),show_col_types=FALSE) %>% mutate(Cancer=factor(Cancer,c("Gastric cancer","Pancreatic cancer","Liver cancer")),Pathway=factor(Pathway,rev(c("ECM / collagen","Integrin / focal adhesion","EMT","TGF-beta","Hypoxia"))))
gsea <- read_tsv(file.path(d4,"GEO_integrin_focal_adhesion_fgsea.tsv"),show_col_types=FALSE) %>% mutate(Cancer=factor(Cancer,c("Gastric cancer","Pancreatic cancer","Liver cancer")))
s4a <- ggplot(gsea,aes(Cancer,NES,color=NES,size=-log10(padj)))+geom_hline(yintercept=0,linetype=2,color="grey55")+geom_point()+geom_text(aes(label=sprintf("FDR=%.2g",padj)),vjust=-1,size=2.8)+scale_color_gradient(low="#f4a582",high="#b2182b")+labs(x=NULL,y="Integrin/focal-adhesion NES",size="-log10(FDR)")+theme_pub(8)
s4b <- ggplot(pathway,aes(Cancer,Pathway,fill=Delta))+geom_tile(color="white",linewidth=.45)+geom_text(aes(label=star(FDR)),size=3)+scale_fill_gradient2(low="#3b6fb6",mid="white",high="#c43b3b",midpoint=0,name="High-low\nscore")+labs(x=NULL,y=NULL)+theme_pub(8)+theme(axis.text.x=element_text(angle=25,hjust=1))
s4 <- (s4a|s4b)+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(face="bold",size=16));save_plot(s4,"Supplementary_FigureS4",9,3.6)

# Supplementary Figure S5: source-table-based specificity analysis.
tcga <- read_tsv(file.path(d5,"TCGA_pancancer_ITGB1_immune_gene_correlations.tsv"),show_col_types=FALSE) %>% filter(Rho>0,FDR<.05) %>% mutate(Cancer=factor(Cancer,levels=unique(Cancer)),Gene=factor(Gene,levels=rev(unique(Gene))))
s5a <- ggplot(tcga,aes(Cancer,Gene))+geom_point(aes(size=-log10(pmax(FDR,1e-10)),color=Rho))+scale_color_gradient(low="#f4a582",high="#b2182b",name="Spearman rho")+scale_size_continuous(range=c(.8,3),name="-log10(FDR)")+labs(x=NULL,y=NULL)+theme_pub(7)+theme(axis.text.x=element_text(angle=60,hjust=1,size=6))
consistency <- read_tsv(file.path(d5,"External_full_pathway_consistency.tsv"),show_col_types=FALSE) %>% filter(Delta>0) %>% mutate(Cohort=factor(Cohort,c("Gastric cancer","Pancreatic cancer","Liver cancer","NSCLC anti-PD-(L)1")),Pathway=fct_rev(factor(Pathway)))
s5b <- ggplot(consistency,aes(Cohort,Pathway))+geom_point(aes(size=-log10(pmax(FDR,1e-10)),color=Delta))+geom_text(aes(label=star(FDR)),nudge_y=.22,size=2.7)+scale_color_gradient(low="#f4a582",high="#b2182b",name="High-low")+scale_size_continuous(range=c(1.2,4),name="-log10(FDR)")+labs(x=NULL,y=NULL)+theme_pub(8)+theme(axis.text.x=element_text(angle=30,hjust=1))
meta <- read_tsv(file.path(d5,"Specificity_contrast_random_effects_meta.tsv"),show_col_types=FALSE) %>% mutate(Signature=fct_rev(factor(Signature)),Category=factor(Category,c("Stromal / exclusion","Effector immune")))
s5c <- ggplot(meta,aes(Pooled_effect,Signature,color=Category))+geom_vline(xintercept=0,linetype=2,linewidth=.25,color="grey50")+geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.18,linewidth=.45)+geom_point(size=2.4)+geom_text(aes(label=star(FDR)),nudge_x=.12,size=2.8,show.legend=FALSE)+scale_color_manual(values=c("Stromal / exclusion"="#b2182b","Effector immune"="#2166ac"))+labs(x="Pooled standardized effect size\n(ITGB1-high minus ITGB1-low)",y=NULL,color=NULL)+theme_pub(8)+theme(legend.position="top")
sumtab <- read_tsv(file.path(d5,"External_specificity_enrichment_summary.tsv"),show_col_types=FALSE)
s5d <- ggplot(sumtab,aes(Category,Fraction,fill=Category))+geom_col(width=.62,color="black",linewidth=.25)+geom_text(aes(label=paste0(Significant_positive,"/",Total_tests)),vjust=-.45,size=3.1)+annotate("text",x=1.5,y=min(1.02,max(sumtab$Fraction)+.18),label=paste0("Fisher P = ",signif(first(sumtab$Fisher_P),3)),size=3.1)+scale_fill_manual(values=c("Stromal / exclusion"="#b2182b","Effector immune"="#2166ac"))+scale_y_continuous(labels=scales::percent_format(accuracy=1),limits=c(0,1.05),expand=c(0,0))+labs(x=NULL,y="Fraction of significant positive tests")+theme_pub(8)+theme(legend.position="none")
s5 <- (s5a/(s5b|s5c)/s5d)+plot_layout(heights=c(1,1.05,.75))+plot_annotation(tag_levels="A")&theme(plot.tag=element_text(face="bold",size=16));save_plot(s5,"Supplementary_FigureS5",12,11.5)
