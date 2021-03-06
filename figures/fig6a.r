suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(here)
  library(DT)
  library(magrittr)
  library(readr)
  library(glue)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  })

fdr_cut <- 0.2

## significant genesets (down-selection)
sig_genesets <- read_csv(here("output/fig2_sig_genesets_downselection.csv")) %>%
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

mal067 <- read_csv(here("output/disease_rtss_dmso_M3.csv")) %>%
  filter(geneset %in% sig_genesets) %>%
  mutate(study = "mal067") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)
vahey <- read_csv(here("data/Vahey_GSEA_disease.csv")) %>%
  filter(geneset %in% sig_genesets,
         time == "T4") %>%
  mutate(study = "vahey") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)
mal068_rrr  <- read_csv(here("output/MAL068_RRR_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal068 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)
mal071_rrr  <- read_csv(here("output/MAL071_RRR_module_disease.csv")) %>%
  mutate(module = str_to_lower(module),
         study = "mal071 RRR") %>%
  inner_join(sig_translation, by = "module") %>%
  dplyr::select(study, geneset, NGenes, Direction, PValue, FDR)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))



## heatmap color scheme
col_fun = colorRamp2(c(-4, 0, 4),  rev(RColorBrewer::brewer.pal(3, "RdBu")))

## Set colors for annotations
ann_colors <- groupings %>%
  dplyr::select(annotation) %>%
  distinct() %>%
  arrange(annotation) %>%
  mutate(col = c(brewer.pal(9, "Set1"),
                 brewer.pal(9, "Pastel1"),
                 brewer.pal(8, "Set2"))[1:n()]) %>%
  deframe()

mal067_sig <- mal067 %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  filter(FDR <= fdr_cut)

plot_data <- bind_rows(vahey, mal068_rrr,
                       mal071_rrr) %>%
  filter(geneset %in% mal067_sig$geneset) %>%
  group_by(study) %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  bind_rows(mal067_sig) %>%
  left_join(groupings, by = "geneset") %>%
  mutate(log10fdr = log10(FDR) * if_else(Direction == "Up", -1, 1)) %>%
  ungroup()

## need the data wide for complex heatmap
wide_data <- plot_data %>%
  filter(!is.na(annotation)) %>%
  dplyr::select(geneset, annotation, study, log10fdr) %>%
  pivot_wider(id_cols = c("geneset", "annotation"),
              names_from = study, values_from = log10fdr) %>%
  arrange(annotation)

## just the heatmap matrix
hm_mat <- as.matrix(wide_data %>%
                      dplyr::select(mal067, vahey,
                             `mal068 RRR`,
                             `mal071 RRR`))
rownames(hm_mat) <- wide_data$geneset

## row annotations
row_ha = rowAnnotation(f = wide_data$annotation,
                       show_legend = F,
                       annotation_label = "",
                       col = list(f = ann_colors[unique(wide_data$annotation)]))

fdr_log_cut <- -log10(0.2)

fig6a <- Heatmap(hm_mat,

                    cluster_columns = F,
                    column_names_rot = 0,
                    column_names_centered = TRUE,

                    cluster_rows = F,
                    show_row_names = T,
                    row_split = wide_data$annotation,
                    row_title_rot = 0,
                    left_annotation = row_ha,
                    row_names_max_width = max_text_width(
                      rownames(hm_mat),
                      gp = gpar(fontsize = 12)
                    ),

                    heatmap_legend_param = list(labels = c("< -4", "-2", "0", "2", "> 4"),
                                                title = "Signed FDR"),
                    col = col_fun,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      curfdr <- 10 ^ (-abs(hm_mat[i, j]))
                      if (curfdr < fdr_cut) {
                        grid.rect(x = x, y = y, width = width, height = height,
                                  gp = gpar(col = "#333333", fill = "transparent"))
                        if (curfdr < 0.01) {
                          grid.text("***", x, y, gp = gpar(fontsize = 12), vjust = "top")
                        } else if (curfdr < 0.05) {
                          grid.text("**", x, y, gp = gpar(fontsize = 12), vjust = "top")
                        } else {
                          grid.text("*", x, y, gp = gpar(fontsize = 12), vjust = "top")
                        }
                      }
                    })

pdf(here("output/figures/Fig6a.pdf"), width = 12, height = 6)
fig6a
dev.off()
png(here("output/figures/Fig6a.png"), width = 1200, height = 600)
fig6a
dev.off()
