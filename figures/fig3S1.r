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

## DMSO Cases vs Controls at M3
dmso_case_control <- read_csv(here("output/disease_comp_dmso_M3.csv"))

## significant genesets (down-selection)
sig_genesets <- read_csv(here("output/fig2_sig_genesets_downselection.csv")) %>%
  pull(geneset)

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})"))

.build_data <- function(data, stim, comp) {
  data %>%
    mutate(stimulation = stim,
           comparison = comp)
}

## put the data together
plot_data <- .build_data(dmso_case_control, "DMSO", "Cases vs Controls") %>%
  filter(geneset %in% sig_genesets) %>%
  group_by(stimulation) %>%
  mutate(FDR2 = p.adjust(PValue, method = "fdr"),
    stimulation = factor(stimulation, levels = c("DMSO", "AMA1", "CSP", "HBsAg")),
         log10p = pmin(4, pmax(-4, log10(PValue) * if_else(Direction == "Up", -1, 1))),
         log10fdr = log10(FDR2) * if_else(Direction == "Up", -1, 1)) %>%
  ungroup() %>%
  left_join(groupings %>% dplyr::select(geneset, annotation), by = "geneset") %>%
  mutate(order = paste(annotation, geneset))

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


sig_modules <- plot_data %>%
  filter(!is.na(annotation),
         FDR2 <= 0.2) %>%
  pull(geneset) %>%
  unique()

## need the data wide for complex heatmap
wide_data <- plot_data %>%
  filter(geneset %in% sig_modules) %>%
  dplyr::select(geneset, annotation, stimulation, log10fdr) %>%
  pivot_wider(id_cols = c("geneset", "annotation"),
              names_from = stimulation, values_from = log10fdr) %>%
  arrange(annotation)


## just the heatmap matrix
hm_mat <- as.matrix(wide_data %>%
                      dplyr::select(DMSO))
rownames(hm_mat) <- wide_data$geneset

## row annotations
row_ha = rowAnnotation(f = wide_data$annotation,
                       show_legend = F,
                       annotation_label = "",
                       col = list(f = ann_colors[unique(wide_data$annotation)]))

fdr_cut <- 0.2
fdr_log_cut <- -log10(fdr_cut)


fig3S1 <- Heatmap(hm_mat,

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

              heatmap_legend_param = list(labels = c("< -2", "-1", "0", "1", "> 2"),
                                          title = "Signed log10 FDR"),
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


pdf(here("output/figures/Fig3S1.pdf"), width = 12, height = 9)
draw(fig3S1)
dev.off()
png(here("output/figures/Fig3S1.png"), width = 1000, height = 750)
draw(fig3S1)
dev.off()


