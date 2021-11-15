suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(here)
  library(DT)
  library(magrittr)
  library(readr)
  library(heatmap3)
  library(glue)
  library(ComplexHeatmap)
  library(cowplot)
  library(readxl)
  library(Biobase)
  library(GSEABase)
  library(RColorBrewer)
  library(circlize)
})

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
           vaccine == "rtss",
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

m3_dmso_rtss_cols <- meta %>%
  filter(stimulation == "dmso",
         vaccine == "rtss",
         visit == "M3") %>%
  pull(col_id) %>%
  unique()

m3_dmso_rtss_pids <- meta %>%
  filter(col_id %in% m3_dmso_rtss_cols) %>%
  pull(pid) %>%
  unique()

m3_iavi <- l_df_rna$IAVI %>%
  dplyr::filter(visit == "M3",
                pid %in% m3_dmso_rtss_pids) %>%
  dplyr::select(pid, eu_elisa_total.nanp, eu_elisa_total.cterm) %>%
  as_tibble()

m3_cevac <- l_df_rna$CEVAC %>%
  filter(visit == "M3",
         pid %in% m3_dmso_rtss_pids) %>%
  dplyr::select(pid, `log_val.Anti-CS`, `log_val.HBV.S AB`) %>%
  as_tibble()

m3_antibody <- l_df_rna$antibody %>%
  filter(visit == "m3",
         pid %in% m3_dmso_rtss_pids) %>%
  mutate(full_name = glue("{isotype}, {analyte}")) %>%
  dplyr::select(pid, full_name, log10_mfi_no_dil) %>%
  as_tibble()

m3_cytokines <- l_df_rna$cytokines %>%
  filter(visit == "m3",
         pid %in% m3_dmso_rtss_pids) %>%
  mutate(full_name = glue("{stimulation}, {analyte}")) %>%
  dplyr::select(pid, full_name, Log10_ratio_conc_imp_trunc) %>%
  as_tibble()

ptid_trans <- cd4_csp_mag %>%
  dplyr::select(PTID, pid) %>%
  distinct()

m3_cd4_csp_pfs <- cd4_csp_mag %>%
  filter(visit == "M3",
         pid %in% m3_dmso_rtss_pids) %>%
  dplyr::select(pid, PFS) %>%
  as_tibble()

m3_cd4_hbs_pfs <- cd4_hbs_mag %>%
  left_join(ptid_trans, by = "PTID") %>%
  filter(pid %in% m3_dmso_rtss_pids) %>%
  dplyr::select(pid, PFS) %>%
  as_tibble()

m3_cd8_csp_pfs <- cd8_csp_mag %>%
  left_join(ptid_trans, by = "PTID") %>%
  filter(pid %in% m3_dmso_rtss_pids) %>%
  dplyr::select(pid, PFS) %>%
  as_tibble()

m3_cd8_hbs_pfs <- cd8_hbs_mag %>%
  left_join(ptid_trans, by = "PTID") %>%
  filter(pid %in% m3_dmso_rtss_pids) %>%
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

stim_cors <- NULL

for (stim in non_dmso_stims) {
  iavi_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_rtss_iavi_cors.csv")))
  cevac_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_rtss_cevac_cors.csv")))
  antibody_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_rtss_antibody_cors.csv")))
  cytokine_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_rtss_cytokine_cors.csv")))

  iavi_cordata <- bind_rows(iavi_cors %>%
                              mutate(assay = "IAVI", variable = "C-term") %>%
                              dplyr::rename(n = cterm_n,
                                            cor = cterm_cor,
                                            pval = cterm_pval) %>%
                              dplyr::select(geneset, assay, variable, n, cor, pval),
                            iavi_cors %>%
                              mutate(assay = "IAVI", variable = "NANP") %>%
                              dplyr::rename(n = nanp_n,
                                            cor = nanp_cor,
                                            pval = nanp_pval) %>%
                              dplyr::select(geneset, assay, variable, n, cor, pval)
  )

  cevac_cordata <- bind_rows(cevac_cors %>%
                               mutate(assay = "CEVAC", variable = "Anti-CS") %>%
                               dplyr::rename(n = anti_cs_n,
                                             cor = anti_cs_cor,
                                             pval = anti_cs_pval) %>%
                               dplyr::select(geneset, assay, variable, n, cor, pval),
                             cevac_cors %>%
                               mutate(assay = "CEVAC", variable = "Hbv-s-Ab") %>%
                               dplyr::rename(n = hbv_s_ab_n,
                                             cor = hbv_s_ab_cor,
                                             pval = hbv_s_ab_pval) %>%
                               dplyr::select(geneset, assay, variable, n, cor, pval)
  )

  antibody_cordata <- antibody_cors %>%
    mutate(assay = "Antibody") %>%
    dplyr::rename(variable = full_name) %>%
    dplyr::select(geneset, assay, variable, n, cor, pval)

  cytokine_cordata <- cytokine_cors %>%
    mutate(assay = "Cytokine") %>%
    dplyr::rename(variable = full_name) %>%
    dplyr::select(geneset, assay, variable, n, cor, pval)

  pfs_cordata <- NULL
  pfs_file <- here(glue("output/m3_{stim}_dmso_diff_rtss_pfs_cors.csv"))
  if (file.exists(pfs_file)) {
    pfs_cordata <- read_csv(pfs_file) %>%
      mutate(assay = "PFS",
             variable = tcell) %>%
      dplyr::select(geneset, assay, variable, n, cor, pval)
  }

  stim_cors[[stim]] <- bind_rows(iavi_cordata,
                                 cevac_cordata,
                                 antibody_cordata,
                                 cytokine_cordata,
                                 pfs_cordata)
}

iavi_cors <- read_csv(here(glue("output/m3_dmso_rtss_iavi_cors.csv")))
cevac_cors <- read_csv(here(glue("output/m3_dmso_rtss_cevac_cors.csv")))
antibody_cors <- read_csv(here(glue("output/m3_dmso_rtss_antibody_cors.csv")))
cytokine_cors <- read_csv(here(glue("output/m3_dmso_rtss_cytokine_cors.csv")))
pfs_cors <- read_csv(here(glue("output/m3_dmso_rtss_pfs_cors.csv")))

iavi_cordata <- bind_rows(iavi_cors %>%
                            mutate(assay = "IAVI", variable = "C-term") %>%
                            dplyr::rename(n = cterm_n,
                                          cor = cterm_cor,
                                          pval = cterm_pval) %>%
                            dplyr::select(geneset, assay, variable, n, cor, pval),
                          iavi_cors %>%
                            mutate(assay = "IAVI", variable = "NANP") %>%
                            dplyr::rename(n = nanp_n,
                                          cor = nanp_cor,
                                          pval = nanp_pval) %>%
                            dplyr::select(geneset, assay, variable, n, cor, pval)
)

cevac_cordata <- bind_rows(cevac_cors %>%
                             mutate(assay = "CEVAC", variable = "Anti-CS") %>%
                             dplyr::rename(n = anti_cs_n,
                                           cor = anti_cs_cor,
                                           pval = anti_cs_pval) %>%
                             dplyr::select(geneset, assay, variable, n, cor, pval),
                           cevac_cors %>%
                             mutate(assay = "CEVAC", variable = "Hbv-s-Ab") %>%
                             dplyr::rename(n = hbv_s_ab_n,
                                           cor = hbv_s_ab_cor,
                                           pval = hbv_s_ab_pval) %>%
                             dplyr::select(geneset, assay, variable, n, cor, pval)
)

antibody_cordata <- antibody_cors %>%
  mutate(assay = "Antibody") %>%
  dplyr::rename(variable = full_name) %>%
  dplyr::select(geneset, assay, variable, n, cor, pval)

cytokine_cordata <- cytokine_cors %>%
  mutate(assay = "Cytokine") %>%
  dplyr::rename(variable = full_name) %>%
  dplyr::select(geneset, assay, variable, n, cor, pval)

pfs_cordata <- pfs_cors %>%
  mutate(assay = "PFS",
         variable = glue("{stim}, {tcell}")) %>%
  dplyr::select(geneset, assay, variable, n, cor, pval)

dmso_cors <- bind_rows(iavi_cordata,
                       cevac_cordata,
                       antibody_cordata,
                       cytokine_cordata,
                       pfs_cordata)


## heatmap color scheme
col_fun = colorRamp2(c(-1, 0, 1),  rev(RColorBrewer::brewer.pal(3, "RdBu")))

## Set colors for annotations
ann_colors <- groupings %>%
  dplyr::select(annotation) %>%
  distinct() %>%
  arrange(annotation) %>%
  mutate(col = c(brewer.pal(9, "Set1"),
                 brewer.pal(9, "Pastel1"),
                 brewer.pal(8, "Set2"))[1:n()]) %>%
  deframe()



all_cordata <- bind_rows(
  stim_cors$csp %>%
    mutate(stim = "CSP"),
  stim_cors$hbs %>%
    mutate(stim = "HBS"),
  stim_cors$ama1 %>%
    mutate(stim = "AMA1"),
  dmso_cors %>%
    mutate(stim = "DMSO")
  )

sig_cells <- all_cordata %>%
  filter(assay != "IAVI" | (assay == "IAVI" & stim == "CSP"),
         assay != "Cytokine" | (assay == "Cytokine" &
                                  str_to_lower(stim) == str_to_lower(str_sub(variable, 1, 3)))) %>%
  group_by(assay, variable) %>%
  mutate(pval_adj = p.adjust(pval, "fdr")) %>%
  ungroup() %>%
  filter(pval_adj <= 0.2) %>%
  dplyr::select(stim, geneset, variable) %>%
  mutate(stim_var = glue("{stim} {variable}"))

filtered_cordata <- all_cordata %>%
  filter(assay != "IAVI" | (assay == "IAVI" & stim == "CSP"),
         assay != "Cytokine" | (assay == "Cytokine" &
                                  str_to_lower(stim) == str_to_lower(str_sub(variable, 1, 3)))) %>%
  group_by(assay, variable) %>%
  mutate(pval_adj = p.adjust(pval, "fdr")) %>%
  ungroup() %>%
  filter(geneset %in% sig_cells$geneset & glue("{stim} {variable}") %in% sig_cells$stim_var)


filtered_cors <- filtered_cordata %>%
  mutate(var2 = glue("{stim} || {assay}][{variable}")) %>%
  pivot_wider(id_cols = "geneset", names_from = var2, values_from = cor) %>%
  left_join(groupings %>% dplyr::select(geneset, annotation), by = "geneset") %>%
  dplyr::filter(!is.na(annotation)) %>%
  arrange(annotation, geneset)

filtered_pvals <- filtered_cordata %>%
  mutate(var2 = glue("{stim} || {assay}][{variable}")) %>%
  pivot_wider(id_cols = "geneset", names_from = var2, values_from = pval_adj) %>%
  left_join(groupings %>% dplyr::select(geneset, annotation), by = "geneset") %>%
  dplyr::filter(!is.na(annotation)) %>%
  arrange(annotation, geneset)
pval_mat <- as.matrix(filtered_pvals %>%
                            dplyr::select(-geneset, -annotation))

## just the heatmap matrix
filtered_mat <- as.matrix(filtered_cors %>%
                            dplyr::select(-geneset, -annotation))
rownames(filtered_mat) <- filtered_cors$geneset

## row annotations
filtered_row_ha = rowAnnotation(f = filtered_cors$annotation,
                                show_legend = F,
                                annotation_label = "",
                                col = list(f = ann_colors[unique(filtered_cors$annotation)]))

## col annotations
filtered_col_ha = columnAnnotation(f = str_split(colnames(filtered_mat), " \\|\\| ") %>%
                                     lapply(`[[`, 2) %>%
                                     unlist() %>%
                                     str_split("\\]\\[") %>%
                                     lapply(`[[`, 1) %>%
                                     unlist(),
                                   show_legend = F,
                                   annotation_label = "",
                                   col = list(f = c(Antibody = "#F1A340", Cytokine = "#998EC3")))


filtered_col_groups <- str_split(colnames(filtered_mat), " \\|\\| ") %>% lapply(`[[`, 1) %>% unlist()
filtered_col_groups <- factor(filtered_col_groups,
                              levels = c("DMSO", "CSP", "HBS", "AMA1"))
colnames(filtered_mat) <- str_split(colnames(filtered_mat), " \\|\\| ") %>%
  lapply(`[[`, 2) %>%
  unlist()  %>%
  str_split("\\]\\[") %>%
  lapply(`[[`, 2) %>%
  unlist()

fig5 <- Heatmap(filtered_mat,
                cluster_columns = F,
                show_column_names = T,
                column_split = filtered_col_groups,
                column_names_rot = 45,
                bottom_annotation = filtered_col_ha,

                cluster_rows = F,
                show_row_names = F,
                row_split = filtered_cors$annotation,
                row_title_rot = 0,
                left_annotation = filtered_row_ha,

                heatmap_legend_param = list(title = "Correlation"),
                col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if (abs(pval_mat[i, j]) <= 0.2) {
                    grid.rect(x = x, y = y, width = width, height = height,
                              gp = gpar(col = "#333333", fill = "transparent"))
                  }
                })

pdf(here(glue("output/figures/Fig5.pdf")),
    width = 18, height = 12)
fig5
dev.off()

png(here(glue("output/figures/Fig5.png")),
    width = 1500, height = 1000)
fig5
dev.off()
