library(here)
library(tidyverse)
library(Biobase)
library(GSEABase)
library(limma)
library(magrittr)
library(SummarizedExperiment)

fdr_cut <- 0.2

## load Vahey data
vahey_eset <- readRDS(here("data", "gpl_eset.Rds"))
meta_dt <- as_tibble(pData(vahey_eset))

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
btm_ind <- ids2indices(geneIds(btm_gtm), rownames(vahey_eset))
btm_ind <- btm_ind[sapply(btm_ind, length) > min_gene_set]


genesets_btm <- btm_gtm[m3_genesets]
genesets_ind <- btm_ind[m3_genesets] %>% unlist() %>% unique()

genesets_genes <- geneIds(genesets_btm) %>% unlist() %>% unique()

## Get voomed genes
sig_eset <- vahey_eset[rownames(vahey_eset) %in% genesets_genes,]


sample_index <- pData(sig_eset) %>%
  filter(!is.na(visit)) %>%           ## remove controls
  filter(disease != "Unknown") %$%
    geo_id

gsea_eset <- sig_eset[,sample_index]

pData(gsea_eset) <- pData(gsea_eset) %>%
  mutate(disease=factor(disease, levels=c("Protected", "NonProtected"))) %>%
  set_rownames(pData(gsea_eset)$geo_id)

form <- "~disease"        ## formula used in camera
coef <- "diseaseNonProtected"

## get subset
pData(gsea_eset) %>%
  filter(visit == "T4") %$%
    geo_id ->
  visit_index

gsea_sub <- gsea_eset[,visit_index]

sig_meta <- as_tibble(pData(gsea_sub))

design_dis <- model.matrix(~disease,
                           sig_meta)

v1 <- voom(exprs(gsea_sub),
           design=design_dis)

# Gene level based analysis
dupcor_imm <- duplicateCorrelation(v1$E,
                                   design_dis,
                                   block = sig_meta$pid)

lm_imm <- lmFit(v1$E,
                design_dis,
                correlation = dupcor_imm)
eb_imm <- eBayes(lm_imm)


vahey_res <- topTable(eb_imm,
                       coef = coef,
                       number = "inf",
                       sort.by = "none") %>%
  arrange(P.Value)

vahey_res <- vahey_res %>%
  mutate(gene_name = rownames(vahey_res)) %>%
  as_tibble()

write_csv(vahey_res, here("output/single_gene/vahey_results.csv"))

