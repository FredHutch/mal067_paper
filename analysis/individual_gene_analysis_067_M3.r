library(here)
library(tidyverse)
library(Biobase)
library(GSEABase)
library(limma)
library(magrittr)
library(SummarizedExperiment)

fdr_cut <- 0.2

## load mal067 data
library(mal067data)
meta_dt <- as_tibble(pData(mal067_eset))

## significant genesets (down-selection)
sig_genesets <- read_csv(here("output/fig2_sig_genesets_downselection.csv")) %>%
  pull(geneset)

m3_disease <- read_csv(here("output/disease_rtss_dmso_M3.csv")) %>%
  filter(geneset %in% sig_genesets) %>%
  dplyr::select(geneset, NGenes, Direction, PValue, FDR)

m3_genesets <- m3_disease %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  filter(FDR <= fdr_cut) %>%
  pull(geneset)

## set up GSEA analysis
min_gene_set <- 5
btm_gtm <- getGmt(here("data/BTM_for_GSEA_20131008.gmt"))
btm_ind <- ids2indices(geneIds(btm_gtm), rownames(mal067_eset))
btm_ind <- btm_ind[sapply(btm_ind, length) > min_gene_set]


genesets_btm <- btm_gtm[m3_genesets]
genesets_ind <- btm_ind[m3_genesets] %>% unlist() %>% unique()

genesets_genes <- geneIds(genesets_btm) %>% unlist() %>% unique()

## Get voomed genes
sig_eset <- mal067_eset[rownames(mal067_eset) %in% genesets_genes,]
sig_voom <- mal067_voom[rownames(mal067_voom$E) %in% genesets_genes,]

## parameters set by parent: stim, time
stim="dmso"
time="M3"

################### DMSO M3 RTS,S vs Comparator ########################

form_main_dis_both <- "~plate + total_reads + age + case"
coeff <- "casecase"

## disease: select subset of expressionSet
meta_dt %>%
  filter(stimulation == stim,
         case == "case" | case == "control",
         visit == time) %$%
  col_id ->
  sample_dis

sig_meta <- meta_dt %>%
  filter(stimulation == stim,
         case == "case" | case == "control",
         visit == time)

sig_voom_analysis <- sig_voom[, meta_dt$col_id %in% sample_dis]

## drop extra case levels
sig_meta$case <- fct_drop(sig_meta$case)


design_dis <- model.matrix(~ plate + total_reads + age + case,
                           sig_meta)

# Gene level based analysis
dupcor_imm <- duplicateCorrelation(sig_voom_analysis$E,
                                   design_dis,
                                   block = sig_meta$pid)

lm_imm <- lmFit(sig_voom_analysis$E,
                design_dis,
                correlation = dupcor_imm)
eb_imm <- eBayes(lm_imm)


mal067_res <- topTable(eb_imm,
                       coef = "casecase",
                       number = "inf",
                       sort.by = "none") %>%
  arrange(P.Value)

mal067_res <- mal067_res %>%
  mutate(gene_name = rownames(mal067_res)) %>%
  as_tibble()

write_csv(mal067_res, here("output/single_gene/mal067_results.csv"))

