suppressPackageStartupMessages({
  library(Biobase)
  library(tidyverse)
  library(data.table)
  library(here)
  library(DT)
  library(magrittr)
  library(readr)
  library(venn)
  library(heatmap3)
  library(glue)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  })

fdr_cut <- 0.2

rtss <- read_csv(here("output/disease_rtss_dmso_M0.csv")) %>%
  filter(!str_sub(geneset, 1, 3) == "TBA") %>%
  mutate(vacc = "rtss") %>%
  dplyr::select(vacc, geneset, NGenes, Direction, PValue, FDR)

comp <- read_csv(here("output/disease_comp_dmso_M0.csv")) %>%
  filter(!str_sub(geneset, 1, 3) == "TBA") %>%
  mutate(vacc = "comparator") %>%
  dplyr::select(vacc, geneset, NGenes, Direction, PValue, FDR)

both <- read_csv(here("output/disease_combined_dmso_M0.csv")) %>%
  filter(!str_sub(geneset, 1, 3) == "TBA") %>%
  mutate(vacc = "both") %>%
  dplyr::select(vacc, geneset, NGenes, Direction, PValue, FDR)


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

plot_data <- bind_rows(rtss, comp, both) %>%
  mutate(log10fdr = log10(FDR) * if_else(Direction == "Up", -1, 1)) %>%
  left_join(groupings, by = "geneset")

sig_modules <- plot_data %>%
  filter(FDR <= fdr_cut) %>%
  pull(geneset)

plot_data <- plot_data %>%
  filter(geneset %in% sig_modules)

## need the data wide for complex heatmap
wide_data <- plot_data %>%
  filter(!is.na(annotation)) %>%
  dplyr::select(geneset, annotation, vacc, log10fdr) %>%
  pivot_wider(id_cols = c("geneset", "annotation"),
              names_from = vacc, values_from = log10fdr) %>%
  arrange(annotation)

## just the heatmap matrix
hm_mat <- as.matrix(wide_data %>%
                      dplyr::select(rtss, both, comparator))
rownames(hm_mat) <- wide_data$geneset

## row annotations
row_ha = rowAnnotation(f = wide_data$annotation,
                       show_legend = F,
                       annotation_label = "",
                       col = list(f = ann_colors[unique(wide_data$annotation)]))

fdr_log_cut <- -log10(0.2)

fig_m0_dis_comp <- Heatmap(hm_mat,

                    cluster_columns = F,
                    column_names_rot = 0,
                    column_names_centered = TRUE,

                    cluster_rows = F,
                    show_row_names = F,
                    row_split = wide_data$annotation,
                    row_title_rot = 0,
                    left_annotation = row_ha,

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

pdf(here("output/figures/m0_dis_comp.pdf"), width = 9, height = 9)
fig_m0_dis_comp
dev.off()
png(here("output/figures/m0_dis_comp.png"), width = 750, height = 750)
fig_m0_dis_comp
dev.off()
