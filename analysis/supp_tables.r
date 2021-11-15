suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(glue)
})

out_folder <- "output"


# #Supplementary Tables
#
# Generate the various supplementary tables to go with figures
#
# ##Figure 2 Supplementary Tables  {.tabset .tabset-fade .tabset-pills}
#
# * RTS,S vs Comparator at M3 results for all stims
# * RTS,S M3 vs M0 for all stims
# * Significant genesets from RTS,S vs Comparator at M3 comparison (any stim, FDR <= 0.2)
#

## fig2_supp_tables
fdr_cut <- 0.2

## DMSO RTSS vs Comp at M3
dmso_rtss_comp <- read_csv(here(out_folder, "vaccine_both_dmso_M3.csv"))
## DMSO RTSS M3 vs M0
dmso_m3_m0 <- read_csv(here(out_folder, "M3-M0_rtss_dmso.csv"))
## AMA1 RTSS vs Comp at M3
ama1_rtss_comp <- read_csv(here(out_folder, "vaccine_both_ama1_M3.csv"))
## AMA1 RTSS M3 vs M0
ama1_m3_m0 <- read_csv(here(out_folder, "M3-M0_rtss_ama1.csv"))
## CSP RTSS vs Comp at M3
csp_rtss_comp <- read_csv(here(out_folder, "vaccine_both_csp_M3.csv"))
## CSP RTSS M3 vs M0
csp_m3_m0 <- read_csv(here(out_folder, "M3-M0_rtss_csp.csv"))
## HBS RTSS vs Comp at M3
hbs_rtss_comp <- read_csv(here(out_folder, "vaccine_both_hbs_M3.csv"))
## HBS RTSS M3 vs M0
hbs_m3_m0 <- read_csv(here(out_folder, "M3-M0_rtss_hbs.csv"))

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

rtss_vs_comp <- bind_rows(
  dmso_rtss_comp %>% mutate(Stimulation = "DMSO"),
  ama1_rtss_comp %>% mutate(Stimulation = "AMA1"),
  csp_rtss_comp %>% mutate(Stimulation = "CSP"),
  hbs_rtss_comp %>% mutate(Stimulation = "HBS")
) %>%
  filter(!str_detect(geneset, "TBA")) %>%
  mutate(is_significant = FDR <= fdr_cut) %>%
  left_join(groupings %>% select(annotation, geneset), by = "geneset") %>%
  arrange(FDR) %>%
  rename(Annotation = annotation, BTM = geneset) %>%
  select(Annotation, BTM, Stimulation, Direction, PValue, FDR, is_significant)

rtss_m3_vs_m0 <- bind_rows(
  dmso_m3_m0 %>% mutate(Stimulation = "DMSO"),
  ama1_m3_m0 %>% mutate(Stimulation = "AMA1"),
  csp_m3_m0 %>% mutate(Stimulation = "CSP"),
  hbs_m3_m0 %>% mutate(Stimulation = "HBS")
) %>%
  filter(!str_detect(geneset, "TBA")) %>%
  mutate(is_significant = FDR <= fdr_cut) %>%
  left_join(groupings %>% select(annotation, geneset), by = "geneset") %>%
  arrange(FDR) %>%
  rename(Annotation = annotation, BTM = geneset) %>%
  select(Annotation, BTM, Stimulation, Direction, PValue, FDR, is_significant)

sig_genesets <- rtss_vs_comp %>%
  filter(is_significant) %>%
  select(Annotation, BTM) %>%
  distinct() %>%
  arrange(Annotation, BTM)

write_csv(rtss_vs_comp, here(out_folder, "supp_tables/table_S1a_rtss_vs_comp_m3_results.csv"))
write_csv(rtss_m3_vs_m0, here(out_folder, "supp_tables/table_S1b_rtss_m3_vs_m0_results.csv"))
write_csv(sig_genesets, here(out_folder, "supp_tables/table_S2_sig_genesets.csv"))


## Figure 3 Supplementary Tables
## RTS,S M3 Cases vs Controls for all stims (after downselection)

##  fig3_supp_tables
## DMSO Cases vs Controls at M3
dmso_case_control <- read_csv(here(out_folder, "disease_rtss_dmso_M3.csv"))
## AMA1 Cases vs Controls at M3
ama1_case_control <- read_csv(here(out_folder, "disease_rtss_ama1_M3.csv"))
## CSP Cases vs Controls at M3
csp_case_control <- read_csv(here(out_folder, "disease_rtss_csp_M3.csv"))
## HBS Cases vs Controls at M3
hbs_case_control <- read_csv(here(out_folder, "disease_rtss_hbs_M3.csv"))

## significant genesets (down-selection)
sig_genesets <- read_csv(here(out_folder, "fig2_sig_genesets_downselection.csv")) %>%
  pull(geneset)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

## put the data together
rtss_case_vs_control <- bind_rows(
  dmso_case_control %>% mutate(Stimulation = "DMSO"),
  ama1_case_control %>% mutate(Stimulation = "AMA1"),
  csp_case_control %>% mutate(Stimulation = "CSP"),
  hbs_case_control %>% mutate(Stimulation = "HBS")
) %>%
  filter(geneset %in% sig_genesets) %>%
  group_by(Stimulation) %>%
  mutate(FDR = p.adjust(PValue, method = "fdr")) %>%
  ungroup() %>%
  mutate(is_significant = FDR <= fdr_cut) %>%
  left_join(groupings %>% select(annotation, geneset), by = "geneset") %>%
  arrange(FDR) %>%
  rename(Annotation = annotation, BTM = geneset) %>%
  select(Annotation, BTM, Stimulation, Direction, PValue, FDR, is_significant)

write_csv(rtss_case_vs_control, here(out_folder, "supp_tables/table_S3_rtss_case_vs_control_results.csv"))

##Figure 5 Supplementary Tables

## * M3 BTM and adaptive response correlations (significant only)

##  fig5_supp_tables
stim_cors <- NULL

for (stim in c("csp", "hbs", "ama1")) {
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

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))


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

filtered_cordata <- all_cordata %>%
  filter(assay != "IAVI" | (assay == "IAVI" & stim == "CSP"),
         assay != "Cytokine" | (assay == "Cytokine" &
                                  str_to_lower(stim) == str_to_lower(str_sub(variable, 1, 3)))) %>%
  group_by(assay, variable) %>%
  mutate(pval_adj = p.adjust(pval, "fdr")) %>%
  ungroup() %>%
  filter(pval_adj <= 0.2)

supp_tbl <- filtered_cordata %>%
  left_join(groupings %>%
              dplyr::select(annotation, geneset), by = "geneset") %>%
  dplyr::rename(Annotation = annotation,
         BTM = geneset,
         Stimulation = stim,
         Assay = assay,
         Variable = variable,
         Correlation = cor,
         PValue = pval,
         FDR = pval_adj) %>%
  arrange(FDR) %>%
  dplyr::select(Stimulation, Annotation, BTM, Assay, Variable, n, Correlation, PValue, FDR)

write_csv(supp_tbl,
          here(out_folder, "supp_tables/table_S4_M3_adaptive_correlation.csv"))

## fig5b_supp_tables
stim_cors <- NULL

for (stim in c("csp", "hbs", "ama1")) {
  iavi_cors <- read_csv(here(glue("output/m0_{stim}_dmso_diff_rtss_iavi_cors.csv")))
  antibody_cors <- read_csv(here(glue("output/m0_{stim}_dmso_diff_rtss_antibody_cors.csv")))
  cytokine_cors <- read_csv(here(glue("output/m0_{stim}_dmso_diff_rtss_cytokine_cors.csv")))

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

  antibody_cordata <- antibody_cors %>%
    mutate(assay = "Antibody") %>%
    dplyr::rename(variable = full_name) %>%
    dplyr::select(geneset, assay, variable, n, cor, pval)

  cytokine_cordata <- cytokine_cors %>%
    mutate(assay = "Cytokine") %>%
    dplyr::rename(variable = full_name) %>%
    dplyr::select(geneset, assay, variable, n, cor, pval)

  pfs_cordata <- NULL
  pfs_file <- here(glue("output/m0_{stim}_dmso_diff_rtss_pfs_cors.csv"))
  if (file.exists(pfs_file)) {
    pfs_cordata <- read_csv(pfs_file) %>%
      mutate(assay = "PFS",
             variable = glue("{tcell} PFS")) %>%
      dplyr::select(geneset, assay, variable, n, cor, pval)
  }

  stim_cors[[stim]] <- bind_rows(iavi_cordata,
                                 ## cevac_cordata,
                                 antibody_cordata,
                                 cytokine_cordata,
                                 pfs_cordata)
}

iavi_cors <- read_csv(here(glue("output/m0_dmso_rtss_iavi_cors.csv")))
antibody_cors <- read_csv(here(glue("output/m0_dmso_rtss_antibody_cors.csv")))
cytokine_cors <- read_csv(here(glue("output/m0_dmso_rtss_cytokine_cors.csv")))
pfs_cors <- read_csv(here(glue("output/m0_dmso_rtss_pfs_cors.csv")))

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
                       antibody_cordata,
                       cytokine_cordata,
                       pfs_cordata)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))


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

filtered_cordata <- all_cordata %>%
  filter(assay != "IAVI" | (assay == "IAVI" & stim == "CSP"),
         assay != "Cytokine" | (assay == "Cytokine" &
                                  str_to_lower(stim) == str_to_lower(str_sub(variable, 1, 3)))) %>%
  group_by(assay, variable) %>%
  mutate(pval_adj = p.adjust(pval, "fdr")) %>%
  ungroup() %>%
  filter(pval_adj <= 0.2)

supp_tbl <- filtered_cordata %>%
  left_join(groupings %>%
              dplyr::select(annotation, geneset), by = "geneset") %>%
  dplyr::rename(Annotation = annotation,
         BTM = geneset,
         Stimulation = stim,
         Assay = assay,
         Variable = variable,
         Correlation = cor,
         PValue = pval,
         FDR = pval_adj) %>%
  arrange(FDR) %>%
  dplyr::select(Stimulation, Annotation, BTM, Assay, Variable, n, Correlation, PValue, FDR)

write_csv(supp_tbl,
          here(out_folder, "supp_tables/table_S4b_M0_adaptive_correlation.csv"))

##Figure 6 Supplementary Tables

## M0 multi-study results
## M3 multi-study results

##  fig6_supp_tables_m0
fdr_cut <- 0.2

mal067 <- read_csv(here(out_folder, "disease_rtss_dmso_M0.csv")) %>%
  filter(FDR <= fdr_cut,
         !str_sub(geneset, 1, 3) == "TBA") %>%
  mutate(study = "mal067") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

sig_genesets <- unique(mal067$geneset)

sig_translation <- tibble(geneset = sig_genesets,
                          module = str_replace_all(sig_genesets, " ", "\\_") %>%
                            str_replace_all("\\(", "") %>%
                            str_replace_all("\\)", "") %>%
                            str_replace_all("\\,", "") %>%
                            str_replace_all(";", "") %>%
                            str_replace_all("\\_\\&\\_", "\\_") %>%
                            str_replace_all("\\_\\-\\_", "\\_") %>%
                            str_replace_all("\\-", "\\_") %>%
                            str_replace_all("\\.", "\\_") %>%
                            str_replace_all("\\/", "\\_") %>%
                            str_to_lower() %>%
                            str_replace_all("dcs", "d\\_cs") %>%
                            str_replace_all("gtpase", "gt\\_pase"))

vahey <- read_csv(here("data/Vahey_GSEA_disease.csv")) %>%
  filter(geneset %in% sig_genesets,
         time == "T0") %>%
  mutate(study = "vahey") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

mal068_rrr  <- read_csv(here(out_folder, "MAL068_RRR_baseline_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal068 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

mal071_rrr  <- read_csv(here(out_folder, "MAL071_RRR_baseline_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal071 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

supp_m0_multistudy <- bind_rows(vahey, mal068_rrr, mal071_rrr) %>%
  group_by(study) %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  ungroup() %>%
  bind_rows(mal067) %>%
  left_join(groupings, by = "geneset") %>%
  rename(Annotation = annotation,
         BTM = geneset,
         Study = study) %>%
  mutate(is_significant = FDR <= fdr_cut) %>%
  dplyr::select(Study, Annotation, BTM, Direction, PValue, FDR, is_significant) %>%
  arrange(FDR)

write_csv(supp_m0_multistudy,
          here(out_folder, "supp_tables/table_S5_m0_multistudy_comparison.csv"))

## fig6_supp_tables_m3
## significant genesets (down-selection)
sig_genesets <- read_csv(here(out_folder, "fig2_sig_genesets_downselection.csv")) %>%
  pull(geneset)

sig_translation <- tibble(geneset = sig_genesets,
                          module = str_replace_all(sig_genesets, " ", "\\_") %>%
                            str_replace_all("\\(", "") %>%
                            str_replace_all("\\)", "") %>%
                            str_replace_all("\\,", "") %>%
                            str_replace_all("\\_\\&\\_", "\\_") %>%
                            str_replace_all("\\_\\-\\_", "\\_") %>%
                            str_replace_all("\\-", "\\_") %>%
                            str_replace_all("\\.", "\\_") %>%
                            str_to_lower() %>%
                            str_replace_all("dcs", "d\\_cs") %>%
                            str_replace_all("gtpase", "gt\\_pase"))

mal067 <- read_csv(here(out_folder, "disease_rtss_dmso_M3.csv")) %>%
  filter(geneset %in% sig_genesets) %>%
  mutate(study = "mal067") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)
vahey <- read_csv(here("data/Vahey_GSEA_disease.csv")) %>%
  filter(geneset %in% sig_genesets,
         time == "T4") %>%
  mutate(study = "vahey") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)
mal068_rrr  <- read_csv(here(out_folder, "MAL068_RRR_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal068 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)
mal071_rrr  <- read_csv(here(out_folder, "MAL071_RRR_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal071 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))


mal067_sig <- mal067 %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  filter(FDR <= fdr_cut)

supp_m3_multistudy <- bind_rows(vahey, mal068_rrr,
                       mal071_rrr) %>%
  filter(geneset %in% mal067_sig$geneset) %>%
  group_by(study) %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  bind_rows(mal067_sig) %>%
  left_join(groupings, by = "geneset") %>%
  rename(Annotation = annotation,
         BTM = geneset,
         Study = study) %>%
  mutate(is_significant = FDR <= fdr_cut) %>%
  dplyr::select(Study, Annotation, BTM, Direction, PValue, FDR, is_significant) %>%
  arrange(FDR)

write_csv(supp_m3_multistudy,
          here(out_folder, "supp_tables/table_S6_m3_multistudy_comparison.csv"))




##Figure S2 Supplementary Tables

## Comparator M3 Cases vs Controls for all stims (after downselection)

 ## fig3_supp_tables
fdr_cut <- 0.2

## DMSO Cases vs Controls at M3
dmso_case_control <- read_csv(here(out_folder, "disease_comp_dmso_M3.csv"))
## AMA1 Cases vs Controls at M3
ama1_case_control <- read_csv(here("output/manuscript/GSEA_M3_ama1_dis_case_by_vac.csv")) %>%
  select(geneset, Direction_comp, PValue_comp, FDR_comp) %>%
  rename(Direction = Direction_comp,
         PValue = PValue_comp,
         FDR = FDR_comp)
## CSP Cases vs Controls at M3
csp_case_control <- read_csv(here("output/manuscript/GSEA_M3_csp_dis_case_by_vac.csv")) %>%
  select(geneset, Direction_comp, PValue_comp, FDR_comp) %>%
  rename(Direction = Direction_comp,
         PValue = PValue_comp,
         FDR = FDR_comp)
## HBS Cases vs Controls at M3
hbs_case_control <- read_csv(here("output/manuscript/GSEA_M3_hbs_dis_case_by_vac.csv")) %>%
  select(geneset, Direction_comp, PValue_comp, FDR_comp) %>%
  rename(Direction = Direction_comp,
         PValue = PValue_comp,
         FDR = FDR_comp)

## significant genesets (down-selection)
sig_genesets <- read_csv(here(out_folder, "fig2_sig_genesets_downselection.csv")) %>%
  pull(geneset)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

## put the data together
comp_case_vs_control <- bind_rows(
  dmso_case_control %>% mutate(Stimulation = "DMSO"),
  ama1_case_control %>% mutate(Stimulation = "AMA1"),
  csp_case_control %>% mutate(Stimulation = "CSP"),
  hbs_case_control %>% mutate(Stimulation = "HBS")
) %>%
  filter(geneset %in% sig_genesets) %>%
  group_by(Stimulation) %>%
  mutate(FDR = p.adjust(PValue, method = "fdr")) %>%
  ungroup() %>%
  mutate(is_significant = FDR <= fdr_cut) %>%
  left_join(groupings %>% select(annotation, geneset), by = "geneset") %>%
  arrange(FDR) %>%
  rename(Annotation = annotation, BTM = geneset) %>%
  select(Annotation, BTM, Stimulation, Direction, PValue, FDR, is_significant)

write_csv(comp_case_vs_control, here(out_folder, "supp_tables/table_S7_comp_case_vs_control_results.csv"))


##Figure S4 Supplementary Tables  {.tabset .tabset-fade .tabset-pills}

## M3 BTM and adaptive response correlations for Comparator recipients (significant only)

## fig5_supp_tables
stim_cors <- NULL

for (stim in c("csp", "hbs", "ama1")) {
  iavi_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_comp_iavi_cors.csv")))
  cevac_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_comp_cevac_cors.csv")))
  antibody_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_comp_antibody_cors.csv")))
  cytokine_cors <- read_csv(here(glue("output/m3_{stim}_dmso_diff_comp_cytokine_cors.csv")))

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


  stim_cors[[stim]] <- bind_rows(iavi_cordata,
                                 cevac_cordata,
                                 antibody_cordata,
                                 cytokine_cordata)
}

iavi_cors <- read_csv(here(glue("output/m3_dmso_comp_iavi_cors.csv")))
cevac_cors <- read_csv(here(glue("output/m3_dmso_comp_cevac_cors.csv")))
antibody_cors <- read_csv(here(glue("output/m3_dmso_comp_antibody_cors.csv")))
cytokine_cors <- read_csv(here(glue("output/m3_dmso_comp_cytokine_cors.csv")))

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


dmso_cors <- bind_rows(iavi_cordata,
                       cevac_cordata,
                       antibody_cordata,
                       cytokine_cordata)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))


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

filtered_cordata <- all_cordata %>%
  filter(assay != "IAVI" | (assay == "IAVI" & stim == "CSP"),
         assay != "Cytokine" | (assay == "Cytokine" &
                                  str_to_lower(stim) == str_to_lower(str_sub(variable, 1, 3)))) %>%
  group_by(assay, variable) %>%
  mutate(pval_adj = p.adjust(pval, "fdr")) %>%
  ungroup() %>%
  filter(pval_adj <= 0.2)

supp_tbl <- filtered_cordata %>%
  left_join(groupings %>%
              dplyr::select(annotation, geneset), by = "geneset") %>%
  dplyr::rename(Annotation = annotation,
         BTM = geneset,
         Stimulation = stim,
         Assay = assay,
         Variable = variable,
         Correlation = cor,
         PValue = pval,
         FDR = pval_adj) %>%
  arrange(FDR) %>%
  dplyr::select(Stimulation, Annotation, BTM, Assay, Variable, n, Correlation, PValue, FDR)

write_csv(supp_tbl,
          here(out_folder, "supp_tables/table_S8_M3_comparator_adaptive_correlation.csv"))


##Figure S5 Supplementary Tables  {.tabset .tabset-fade .tabset-pills}

## M0 multi-study results for comparator recipients in MAL067

##r fig6_supp_tables_m0
fdr_cut <- 0.2

mal067 <- read_csv(here("output/disease_comp_dmso_M0.csv")) %>%
  filter(FDR <= fdr_cut,
         !str_sub(geneset, 1, 3) == "TBA") %>%
  mutate(study = "mal067") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

sig_genesets <- unique(mal067$geneset)

sig_translation <- tibble(geneset = sig_genesets,
                          module = str_replace_all(sig_genesets, " ", "\\_") %>%
                            str_replace_all("\\(", "") %>%
                            str_replace_all("\\)", "") %>%
                            str_replace_all("\\,", "") %>%
                            str_replace_all(";", "") %>%
                            str_replace_all("\\_\\&\\_", "\\_") %>%
                            str_replace_all("\\_\\-\\_", "\\_") %>%
                            str_replace_all("\\-", "\\_") %>%
                            str_replace_all("\\.", "\\_") %>%
                            str_replace_all("\\/", "\\_") %>%
                            str_to_lower() %>%
                            str_replace_all("dcs", "d\\_cs") %>%
                            str_replace_all("gtpase", "gt\\_pase"))

vahey <- read_csv(here("data/Vahey_GSEA_disease.csv")) %>%
  filter(geneset %in% sig_genesets,
         time == "T0") %>%
  mutate(study = "vahey") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

mal068_rrr  <- read_csv(here("output/MAL068_RRR_baseline_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal068 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

mal071_rrr  <- read_csv(here("output/MAL071_RRR_baseline_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal071 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

supp_m0_multistudy <- bind_rows(vahey, mal068_rrr, mal071_rrr) %>%
  group_by(study) %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  ungroup() %>%
  bind_rows(mal067) %>%
  left_join(groupings, by = "geneset") %>%
  rename(Annotation = annotation,
         BTM = geneset,
         Study = study) %>%
  mutate(is_significant = FDR <= fdr_cut) %>%
  dplyr::select(Study, Annotation, BTM, Direction, PValue, FDR, is_significant) %>%
  arrange(FDR)

write_csv(supp_m0_multistudy,
          here(out_folder, "supp_tables/table_S9_m0_comparator_multistudy_comparison.csv"))
