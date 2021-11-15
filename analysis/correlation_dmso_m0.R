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

## RNAseq data, M0, DMSO, RTS,S
m0_dmso_rtss_cols <- meta %>%
  filter(stimulation == "dmso",
         vaccine == "rtss",
         visit == "M0") %>%
  pull(col_id) %>%
  unique()

m0_dmso_rtss_pids <- meta %>%
  filter(col_id %in% m0_dmso_rtss_cols) %>%
  pull(pid) %>%
  unique()

m0_expr <- exprAll[,m0_dmso_rtss_cols]

m0_geneset_scores <- NULL
for (gs in sig_genesets) {
  idxs <- unlist(setsIndices[gs])
  gs_scores <- apply(m0_expr[idxs,],2,mean)
  m0_geneset_scores <- bind_rows(
    m0_geneset_scores,
    tibble(col_id = names(gs_scores),
           geneset = gs,
           score = gs_scores)
  )
}

m0_geneset_scores <- m0_geneset_scores %>%
  left_join(meta %>% dplyr::select(col_id, pid),
            by = "col_id") %>%
  dplyr::select(col_id, pid, geneset, score)

## Other datasets

m3_iavi <- l_df_rna$IAVI %>%
  filter(visit == "M3",
         pid %in% m0_dmso_rtss_pids) %>%
  dplyr::select(pid, eu_elisa_total.nanp, eu_elisa_total.cterm) %>%
  as_tibble()

m3_cevac <- l_df_rna$CEVAC %>%
  filter(visit == "M3",
         pid %in% m0_dmso_rtss_pids) %>%
  dplyr::select(pid, `log_val.Anti-CS`, `log_val.HBV.S AB`) %>%
  as_tibble()

m3_antibody <- l_df_rna$antibody %>%
  filter(visit == "m3",
         pid %in% m0_dmso_rtss_pids) %>%
  mutate(full_name = glue("{isotype}, {analyte}")) %>%
  dplyr::select(pid, full_name, log10_mfi_no_dil) %>%
  as_tibble()

m3_cytokines <- l_df_rna$cytokines %>%
  filter(visit == "m3",
         pid %in% m0_dmso_rtss_pids) %>%
  mutate(full_name = glue("{stimulation}, {analyte}")) %>%
  dplyr::select(pid, full_name, Log10_ratio_conc_imp_trunc) %>%
  as_tibble()

ptid_trans <- cd4_csp_mag %>%
  dplyr::select(PTID, pid) %>%
  distinct()

m3_cd4_csp_pfs <- cd4_csp_mag %>%
  filter(visit == "M3",
         pid %in% m0_dmso_rtss_pids) %>%
  dplyr::select(pid, PFS) %>%
  as_tibble()

m3_cd4_hbs_pfs <- cd4_hbs_mag %>%
  left_join(ptid_trans, by = "PTID") %>%
  filter(pid %in% m0_dmso_rtss_pids) %>%
  dplyr::select(pid, PFS) %>%
  as_tibble()

m3_cd8_csp_pfs <- cd8_csp_mag %>%
  left_join(ptid_trans, by = "PTID") %>%
  filter(pid %in% m0_dmso_rtss_pids) %>%
  dplyr::select(pid, PFS) %>%
  as_tibble()

m3_cd8_hbs_pfs <- cd8_hbs_mag %>%
  left_join(ptid_trans, by = "PTID") %>%
  filter(pid %in% m0_dmso_rtss_pids) %>%
  dplyr::select(pid, PFS) %>%
  as_tibble()

m3_pfs <- bind_rows(
  m3_cd4_csp_pfs %>%
    mutate(stim = "csp", tcell = "CD4"),
  m3_cd4_hbs_pfs %>%
    mutate(stim = "hbs", tcell = "CD4"),
  m3_cd8_csp_pfs %>%
    mutate(stim = "csp", tcell = "CD8"),
  m3_cd8_hbs_pfs %>%
    mutate(stim = "hbs", tcell = "CD8")
)

##  correlations
iavi_cors <- m0_geneset_scores %>%
  inner_join(m3_iavi, by = "pid") %>%
  group_by(geneset) %>%
  summarize(cterm_n    = sum(!is.na(eu_elisa_total.cterm)),
            cterm_cor  = cor.test(score, eu_elisa_total.cterm, method = "spearman")$estimate,
            cterm_pval = cor.test(score, eu_elisa_total.cterm, method = "spearman")$p.value,
            nanp_n     = sum(!is.na(eu_elisa_total.nanp)),
            nanp_cor   = cor.test(score, eu_elisa_total.nanp, method = "spearman")$estimate,
            nanp_pval  = cor.test(score, eu_elisa_total.nanp, method = "spearman")$p.value,
            .groups = "drop")

antibody_cors <- m0_geneset_scores %>%
  inner_join(m3_antibody, by = "pid") %>%
  group_by(geneset, full_name) %>%
  summarize(n    = sum(!is.na(log10_mfi_no_dil)),
            cor  = cor.test(score, log10_mfi_no_dil, method = "spearman")$estimate,
            pval = cor.test(score, log10_mfi_no_dil, method = "spearman")$p.value,
            .groups = "drop")

cytokine_cors <- m0_geneset_scores %>%
  inner_join(m3_cytokines, by = "pid") %>%
  group_by(geneset, full_name) %>%
  summarize(n    = sum(!is.na(Log10_ratio_conc_imp_trunc)),
            cor  = cor.test(score, Log10_ratio_conc_imp_trunc, method = "spearman")$estimate,
            pval = cor.test(score, Log10_ratio_conc_imp_trunc, method = "spearman")$p.value,
            .groups = "drop")


pfs_cors <- m0_geneset_scores %>%
  inner_join(m3_pfs, by = "pid") %>%
  group_by(geneset, stim, tcell) %>%
  summarize(n    = sum(!is.na(PFS)),
            cor  = cor.test(score, PFS, method = "spearman")$estimate,
            pval = cor.test(score, PFS, method = "spearman")$p.value,
            .groups = "drop")

write_csv(iavi_cors, here("output/m0_dmso_rtss_iavi_cors.csv"))
write_csv(antibody_cors, here("output/m0_dmso_rtss_antibody_cors.csv"))
write_csv(cytokine_cors, here("output/m0_dmso_rtss_cytokine_cors.csv"))
write_csv(pfs_cors, here("output/m0_dmso_rtss_pfs_cors.csv"))
