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

cd4_mag <- read_rds(here("data/PRIMARY_CD4_MAGNITUDE.rds"))
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


stims <- unique(meta$stimulation)
non_dmso_stims <- stims[stims != "dmso"]

m3_stim_geneset_scores <- NULL

for (stim in stims) {
  stim_cols <- meta %>%
    filter(stimulation == stim,
           vaccine == "comparator",
           visit == "M3") %>%
    pull(col_id) %>%
    unique()

  stim_pids <- meta %>%
    filter(col_id %in% stim_cols) %>%
    pull(pid) %>%
    unique()

  stim_expr <- exprAll[,stim_cols]

  stim_geneset_scores <- NULL
  for (gs in sig_genesets) {
    idxs <- unlist(setsIndices[gs])
    gs_scores <- apply(stim_expr[idxs,],2,mean)
    stim_geneset_scores <- bind_rows(
      stim_geneset_scores,
      tibble(col_id = names(gs_scores),
             geneset = gs,
             score = gs_scores)
    )
  }

  m3_stim_geneset_scores[[stim]] <- stim_geneset_scores %>%
    left_join(meta %>% dplyr::select(col_id, pid),
              by = "col_id") %>%
    dplyr::select(col_id, pid, geneset, score)
}

dmso_geneset_scores <- m3_stim_geneset_scores$dmso %>%
  dplyr::rename(dmso_score = score) %>%
  dplyr::select(-col_id)

m3_stim_geneset_scores$dmso <- NULL

stim_dmso_geneset_scores <- NULL

for (stim in non_dmso_stims) {
  stim_geneset_scores <- m3_stim_geneset_scores[[stim]]
  stim_dmso_geneset_scores[[stim]] <- inner_join(stim_geneset_scores,
                                 dmso_geneset_scores,
                                 by = c("pid", "geneset")) %>%
    mutate(diff_score = score - dmso_score)
}

## Other datasets

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

m3_iavi <- l_df_rna$IAVI %>%
  dplyr::filter(visit == "M3",
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


## correlations
for (s in non_dmso_stims) {
  stim_geneset_scores <- stim_dmso_geneset_scores[[s]]
  iavi_cors <- stim_geneset_scores %>%
    inner_join(m3_iavi, by = "pid") %>%
    group_by(geneset) %>%
    summarize(cterm_n    = sum(!is.na(eu_elisa_total.cterm)),
              cterm_cor  = cor.test(diff_score, eu_elisa_total.cterm, method = "spearman")$estimate,
              cterm_pval = cor.test(diff_score, eu_elisa_total.cterm, method = "spearman")$p.value,
              nanp_n     = sum(!is.na(eu_elisa_total.nanp)),
              nanp_cor   = cor.test(diff_score, eu_elisa_total.nanp, method = "spearman")$estimate,
              nanp_pval  = cor.test(diff_score, eu_elisa_total.nanp, method = "spearman")$p.value,
              .groups = "drop")

  cevac_cors <- stim_geneset_scores %>%
    inner_join(m3_cevac, by = "pid") %>%
    group_by(geneset) %>%
    summarize(anti_cs_n     = sum(!is.na(`log_val.Anti-CS`)),
              anti_cs_cor   = cor.test(diff_score, `log_val.Anti-CS`, method = "spearman")$estimate,
              anti_cs_pval  = cor.test(diff_score, `log_val.Anti-CS`, method = "spearman")$p.value,
              hbv_s_ab_n    = sum(!is.na(`log_val.HBV.S AB`)),
              hbv_s_ab_cor  = cor.test(diff_score, `log_val.HBV.S AB`, method = "spearman")$estimate,
              hbv_s_ab_pval = cor.test(diff_score, `log_val.HBV.S AB`, method = "spearman")$p.value,
              .groups = "drop")

  antibody_cors <- stim_geneset_scores %>%
    inner_join(m3_antibody, by = "pid") %>%
    group_by(geneset, full_name) %>%
    summarize(n    = sum(!is.na(log10_mfi_no_dil)),
              cor  = cor.test(diff_score, log10_mfi_no_dil, method = "spearman")$estimate,
              pval = cor.test(diff_score, log10_mfi_no_dil, method = "spearman")$p.value,
              .groups = "drop")

  cytokine_cors <- stim_geneset_scores %>%
    inner_join(m3_cytokines, by = "pid") %>%
    group_by(geneset, full_name) %>%
    summarize(n    = sum(!is.na(Log10_ratio_conc_imp_trunc)),
              cor  = cor.test(diff_score, Log10_ratio_conc_imp_trunc, method = "spearman")$estimate,
              pval = cor.test(diff_score, Log10_ratio_conc_imp_trunc, method = "spearman")$p.value,
              .groups = "drop")

  write_csv(iavi_cors, here(glue("output/m3_{s}_dmso_diff_comp_iavi_cors.csv")))
  write_csv(cevac_cors, here(glue("output/m3_{s}_dmso_diff_comp_cevac_cors.csv")))
  write_csv(antibody_cors, here(glue("output/m3_{s}_dmso_diff_comp_antibody_cors.csv")))
  write_csv(cytokine_cors, here(glue("output/m3_{s}_dmso_diff_comp_cytokine_cors.csv")))
}
