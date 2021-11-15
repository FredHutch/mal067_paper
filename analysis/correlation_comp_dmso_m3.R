  library(tidyverse)
  library(data.table)
  library(here)
  library(DT)
  library(magrittr)
  library(readr)
  library(glue)
  library(Biobase)
  library(GSEABase)

library(mal067data)
data(mal067_eset)
meta <- pData(mal067_eset) %>%
  as_tibble()
exprAll <- mal067_voom$E

load(here("data/m067_seattle_data.RData"))

names(l_df_rna) <- c("cytokines", "IAVI", "CEVAC", "antibody")

cd4_csp_mag <- read_rds(here("data/CD4_CSP_PFS.rds"))
cd4_hbs_mag <- read_rds(here("data/CD4_HBS_PFS.rds"))
cd8_csp_mag <- read_rds(here("data/CD8_CSP_PFS.rds"))
cd8_hbs_mag <- read_rds(here("data/CD8_HBS_PFS.rds"))

## significant genesets (down-selection) based on RTS,S vs Comparator, M3 (FDR <= 0.2, not TBA)
sig_genesets <- read_csv(here("output/fig2_sig_genesets_downselection.csv")) %>%
  pull(geneset)

minGeneSetCut <- 5

geneSet <- getGmt(here("data/BTM_for_GSEA_20131008.gmt"))
geneIds <- geneIds(geneSet)
setsIndices <- ids2indices(geneIds, rownames(exprAll))
setsIndices <- setsIndices[sapply(setsIndices, length) > minGeneSetCut]

breaks <- seq(-1,1,0.25)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

## RNAseq data, M3, DMSO, RTS,S
m3_dmso_comp_cols <- meta %>%
  filter(stimulation == "dmso",
         vaccine == "comparator",
         visit == "M3") %>%
  pull(col_id) %>%
  unique()

m3_dmso_comp_pids <- meta %>%
  filter(col_id %in% m3_dmso_comp_cols) %>%
  pull(pid) %>%
  unique()

m3_expr <- exprAll[,m3_dmso_comp_cols]

m3_geneset_scores <- NULL
for (gs in sig_genesets) {
  idxs <- unlist(setsIndices[gs])
  gs_scores <- apply(m3_expr[idxs,],2,mean)
  m3_geneset_scores <- bind_rows(
    m3_geneset_scores,
    tibble(col_id = names(gs_scores),
           geneset = gs,
           score = gs_scores)
  )
}

m3_geneset_scores <- m3_geneset_scores %>%
  left_join(meta %>% dplyr::select(col_id, pid),
            by = "col_id") %>%
  dplyr::select(col_id, pid, geneset, score)

## Other datasets

m3_iavi <- l_df_rna$IAVI %>%
  filter(visit == "M3",
         pid %in% m3_dmso_comp_pids) %>%
  dplyr::select(pid, eu_elisa_total.nanp, eu_elisa_total.cterm) %>%
  as_tibble()

m3_cevac <- l_df_rna$CEVAC %>%
  filter(visit == "M3",
         pid %in% m3_dmso_comp_pids) %>%
  dplyr::select(pid, `log_val.Anti-CS`, `log_val.HBV.S AB`) %>%
  as_tibble()

m3_antibody <- l_df_rna$antibody %>%
  filter(visit == "m3",
         pid %in% m3_dmso_comp_pids) %>%
  mutate(full_name = glue("{isotype}, {analyte}")) %>%
  dplyr::select(pid, full_name, log10_mfi_no_dil) %>%
  as_tibble()

m3_cytokines <- l_df_rna$cytokines %>%
  filter(visit == "m3",
         pid %in% m3_dmso_comp_pids) %>%
  mutate(full_name = glue("{stimulation}, {analyte}")) %>%
  dplyr::select(pid, full_name, Log10_ratio_conc_imp_trunc) %>%
  as_tibble()


##  correlations
iavi_cors <- m3_geneset_scores %>%
  inner_join(m3_iavi, by = "pid") %>%
  group_by(geneset) %>%
  summarize(cterm_n    = sum(!is.na(eu_elisa_total.cterm)),
            cterm_cor  = cor.test(score, eu_elisa_total.cterm, method = "spearman")$estimate,
            cterm_pval = cor.test(score, eu_elisa_total.cterm, method = "spearman")$p.value,
            nanp_n     = sum(!is.na(eu_elisa_total.nanp)),
            nanp_cor   = cor.test(score, eu_elisa_total.nanp, method = "spearman")$estimate,
            nanp_pval  = cor.test(score, eu_elisa_total.nanp, method = "spearman")$p.value,
            .groups = "drop")

cevac_cors <- m3_geneset_scores %>%
  inner_join(m3_cevac, by = "pid") %>%
  group_by(geneset) %>%
  summarize(anti_cs_n     = sum(!is.na(`log_val.Anti-CS`)),
            anti_cs_cor   = cor.test(score, `log_val.Anti-CS`, method = "spearman")$estimate,
            anti_cs_pval  = cor.test(score, `log_val.Anti-CS`, method = "spearman")$p.value,
            hbv_s_ab_n    = sum(!is.na(`log_val.HBV.S AB`)),
            hbv_s_ab_cor  = cor.test(score, `log_val.HBV.S AB`, method = "spearman")$estimate,
            hbv_s_ab_pval = cor.test(score, `log_val.HBV.S AB`, method = "spearman")$p.value,
            .groups = "drop")

antibody_cors <- m3_geneset_scores %>%
  inner_join(m3_antibody, by = "pid") %>%
  group_by(geneset, full_name) %>%
  summarize(n    = sum(!is.na(log10_mfi_no_dil)),
            cor  = cor.test(score, log10_mfi_no_dil, method = "spearman")$estimate,
            pval = cor.test(score, log10_mfi_no_dil, method = "spearman")$p.value,
            .groups = "drop")

cytokine_cors <- m3_geneset_scores %>%
  inner_join(m3_cytokines, by = "pid") %>%
  group_by(geneset, full_name) %>%
  summarize(n    = sum(!is.na(Log10_ratio_conc_imp_trunc)),
            cor  = cor.test(score, Log10_ratio_conc_imp_trunc, method = "spearman")$estimate,
            pval = cor.test(score, Log10_ratio_conc_imp_trunc, method = "spearman")$p.value,
            .groups = "drop")

write_csv(iavi_cors, here("output/m3_dmso_comp_iavi_cors.csv"))
write_csv(cevac_cors, here("output/m3_dmso_comp_cevac_cors.csv"))
write_csv(antibody_cors, here("output/m3_dmso_comp_antibody_cors.csv"))
write_csv(cytokine_cors, here("output/m3_dmso_comp_cytokine_cors.csv"))
