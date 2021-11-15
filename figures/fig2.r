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
  library(circlize)
  library(RColorBrewer)
})

## DMSO RTSS vs Comp at M3
dmso_rtss_comp <- read_csv(here("output/vaccine_both_dmso_M3.csv"))
## DMSO RTSS M3 vs M0
dmso_m3_m0 <- read_csv(here("output/M3-M0_rtss_dmso.csv"))
## AMA1 RTSS vs Comp at M3
ama1_rtss_comp <- read_csv(here("output/vaccine_both_ama1_M3.csv"))
## AMA1 RTSS M3 vs M0
ama1_m3_m0 <- read_csv(here("output/M3-M0_rtss_ama1.csv"))
## CSP RTSS vs Comp at M3
csp_rtss_comp <- read_csv(here("output/vaccine_both_csp_M3.csv"))
## CSP RTSS M3 vs M0
csp_m3_m0 <- read_csv(here("output/M3-M0_rtss_csp.csv"))
## HBS RTSS vs Comp at M3
hbs_rtss_comp <- read_csv(here("output/vaccine_both_hbs_M3.csv"))
## HBS RTSS M3 vs M0
hbs_m3_m0 <- read_csv(here("output/M3-M0_rtss_hbs.csv"))

.build_data <- function(data, stim, comp) {
  data %>%
    mutate(stimulation = stim,
           comparison = comp)
}

## put the data together
comp_data <- bind_rows(
  .build_data(dmso_rtss_comp, "DMSO", "RTS,S vs Comp"),
  .build_data(dmso_m3_m0, "DMSO", "M3 vs M0"),
  .build_data(ama1_rtss_comp, "AMA1", "RTS,S vs Comp"),
  .build_data(ama1_m3_m0, "AMA1", "M3 vs M0"),
  .build_data(csp_rtss_comp, "CSP", "RTS,S vs Comp"),
  .build_data(csp_m3_m0, "CSP", "M3 vs M0"),
  .build_data(hbs_rtss_comp, "HBsAg", "RTS,S vs Comp"),
  .build_data(hbs_m3_m0, "HBsAg", "M3 vs M0")
) %>%
  mutate(category = glue("{stimulation}, {comparison}"),
         category2 = glue("{comparison}, {stimulation}"),
         category2 = factor(category2,
                            levels = paste0(c(rep("RTS,S vs Comp", 4), rep("M3 vs M0", 4)),
                                            ", ", c("DMSO", "AMA1", "CSP", "HBsAg"))))

fdr_cut <- 0.2

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

## Get the significant genesets
sig_genesets <- comp_data %>%
  filter(FDR <= fdr_cut,
         !str_detect(geneset, "TBA")) %>%
  pull(geneset) %>%
  unique()

write_csv(comp_data %>%
            filter(comparison == "RTS,S vs Comp",
                   FDR <= fdr_cut,
                   !str_detect(geneset, "TBA")) %>%
            select(geneset) %>%
            distinct(),
          here("output/fig2_sig_genesets_downselection.csv"))

## prepare data for plots
plot_data <- comp_data %>%
  filter(geneset %in% sig_genesets) %>%
  mutate(log10p = pmin(4, pmax(-4, log10(PValue) * if_else(Direction == "Up", -1, 1))),
         log10fdr = log10(FDR) * if_else(Direction == "Up", -1, 1)) %>%
  left_join(groupings %>% select(geneset, annotation), by = "geneset") %>%
  filter(!is.na(annotation)) %>%
  mutate(order = paste(annotation, geneset))

## Set colors for annotations
ann_colors <- groupings %>%
  select(annotation) %>%
  distinct() %>%
  arrange(annotation) %>%
  mutate(col = c(brewer.pal(9, "Set1"),
                 brewer.pal(9, "Pastel1"),
                 brewer.pal(8, "Set2"))[1:n()]) %>%
  deframe()

## heatmap color scheme
col_fun = colorRamp2(c(-2, 0, 2),  rev(RColorBrewer::brewer.pal(3, "RdBu")))

## Set colors for annotations
ann_colors <- groupings %>%
  dplyr::select(annotation) %>%
  distinct() %>%
  arrange(annotation) %>%
  mutate(col = c(brewer.pal(9, "Set1"),
                 brewer.pal(9, "Pastel1"),
                 brewer.pal(8, "Set2"))[1:n()]) %>%
  deframe()

## need the data wide for complex heatmap
wide_data <- plot_data %>%
  select(geneset, annotation, category2, log10fdr) %>%
  pivot_wider(id_cols = c("geneset", "annotation"),
              names_from = category2, values_from = log10fdr) %>%
  arrange(annotation)

## just the heatmap matrix
hm_mat <- as.matrix(wide_data %>%
                      select(`RTS,S vs Comp, DMSO`,
                             `RTS,S vs Comp, CSP`,
                             `RTS,S vs Comp, HBsAg`,
                             `RTS,S vs Comp, AMA1`,
                             `M3 vs M0, DMSO`,
                             `M3 vs M0, CSP`,
                             `M3 vs M0, HBsAg`,
                             `M3 vs M0, AMA1`))
rownames(hm_mat) <- wide_data$geneset

## row annotations
row_ha = rowAnnotation(f = wide_data$annotation,
                       show_legend = F,
                       annotation_label = "",
                       col = list(f = ann_colors[unique(wide_data$annotation)]))

fdr_log_cut <- -log10(fdr_cut)

ht <- Heatmap(hm_mat,

              cluster_columns = F,
              column_labels = rep(c("DMSO", "CSP", "HBS", "AMA1"), 2),
              column_names_rot = 0,
              column_names_centered = TRUE,
              column_split = factor(c(rep("RTS,S vs Comparator, M3", 4),
                                      rep("RTS,S M3 vs M0", 4)),
                                    levels = c("RTS,S vs Comparator, M3",
                                               "RTS,S M3 vs M0")),

              cluster_rows = F,
              show_row_names = F,
              row_split = wide_data$annotation,
              row_title_rot = 0,
              left_annotation = row_ha,

              heatmap_legend_param = list(labels = c("< -2", "-1", "0", "1", "> 2"),
                                          title = "Signed log10 FDR"),
              col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                curfdr <- 10 ^ (-abs(hm_mat[i, j]))
                if (curfdr < fdr_cut) {
                  grid.rect(x = x, y = y, width = width, height = height,
                            gp = gpar(col = "#333333", fill = "transparent"))
                  if (curfdr < 0.01) {
                    grid.text("***", x, y, gp = gpar(fontsize = 10), vjust = "top")
                  } else if (curfdr < 0.05) {
                    grid.text("**", x, y, gp = gpar(fontsize = 10), vjust = "top")
                  } else {
                    grid.text("*", x, y, gp = gpar(fontsize = 10), vjust = "top")
                  }
                }
              })
# draw(ht)
pdf(here("output/figures/Fig2.pdf"), width = 12, height = 9)
ht
dev.off()

png(here("output/figures/Fig2.png"), width = 1000, height = 750)
ht
dev.off()
