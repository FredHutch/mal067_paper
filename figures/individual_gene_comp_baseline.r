library(here)
library(tidyverse)
library(glue)

mal067_res <- read_csv(here("output/single_gene/mal067_M0_results.csv")) %>%
  mutate(study = "mal067",
         direction = if_else(logFC > 0, "Up", "Down"),
         signed_log_p = -sign(logFC) * log10(P.Value)) %>%
  dplyr::rename(pval = P.Value)
mal068_res <- read_csv(here("output/single_gene/mal068_baseline_results.csv")) %>%
  mutate(study = "mal068",
         direction = if_else(logFC > 0, "Up", "Down"),
         signed_log_p = -sign(logFC) * log10(P.Value)) %>%
  dplyr::rename(pval = P.Value)
mal071_res <- read_csv(here("output/single_gene/mal071_baseline_results.csv")) %>%
  mutate(study = "mal071",
         direction = if_else(logFC > 0, "Up", "Down"),
         signed_log_p = -sign(logFC) * log10(P.Value)) %>%
  dplyr::rename(pval = P.Value)
vahey_res <- read_csv(here("output/single_gene/vahey_baseline_results.csv")) %>%
  mutate(study = "vahey",
         direction = if_else(logFC > 0, "Up", "Down"),
         signed_log_p = -sign(logFC) * log10(P.Value)) %>%
  dplyr::rename(pval = P.Value)


comp <- mal067_res %>%
  dplyr::select(gene_name, direction, pval, logFC, signed_log_p) %>%
  dplyr::rename(direction_067 = direction,
                pval_067 = pval,
                logFC_067 = logFC,
                logp_067 = signed_log_p) %>%
  inner_join(mal068_res %>%
               dplyr::select(gene_name, direction, pval, logFC, signed_log_p) %>%
               dplyr::rename(direction_068 = direction,
                             pval_068 = pval,
                             logFC_068 = logFC,
                             logp_068 = signed_log_p),
             by = "gene_name") %>%
  inner_join(mal071_res %>%
               dplyr::select(gene_name, direction, pval, logFC, signed_log_p) %>%
               dplyr::rename(direction_071 = direction,
                             pval_071 = pval,
                             logFC_071 = logFC,
                             logp_071 = signed_log_p),
             by = "gene_name") %>%
  inner_join(vahey_res %>%
               dplyr::select(gene_name, direction, pval, logFC, signed_log_p) %>%
               dplyr::rename(direction_vahey = direction,
                             pval_vahey = pval,
                             logFC_vahey = logFC,
                             logp_vahey = signed_log_p),
             by = "gene_name")

fdr_cut <- 0.2

## load mal067 data
library(mal067data)
meta_dt <- as_tibble(pData(mal067_eset))

m0_disease <- read_csv(here("output/dmso_M0_rtss_disease.csv")) %>%
  filter(FDR <= fdr_cut,
         !str_sub(geneset, 1, 3) == "TBA") %>%
  dplyr::select(geneset, NGenes, Direction, PValue, FDR)

m0_genesets <- m0_disease %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  filter(FDR <= fdr_cut) %>%
  pull(geneset)



## set up GSEA analysis

library(Biobase)
library(GSEABase)
min_gene_set <- 5
btm_gtm <- getGmt(here("data/BTM_for_GSEA_20131008.gmt"))
genesets_btm <- btm_gtm[m0_genesets]

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


col_fun2 = colorRamp2(c(-2, 0, 2),  rev(RColorBrewer::brewer.pal(3, "RdBu")))

for (i in 1:length(m0_genesets)) {
  btm_name <- names(genesets_btm[i])
  genes <- geneIds(genesets_btm[i])[[1]]
  btm_m <- str_split(btm_name, "\\(")[[1]] %>% last() %>% str_sub(1, -2)

  ngenes <- length(genes)
  fs <- 4
  fsa <- 8
  if (ngenes < 100) {
    fs <- 6
    fsa <- 10
  }
  if (ngenes < 50) {
    fs <- 8
    fsa <- 12
  }

  btm_comp <- comp %>%
    filter(gene_name %in% genes) %>%
    arrange(logp_067)

  hm_mat2 <- btm_comp %>%
    dplyr::select(logp_067, logp_vahey, logp_068, logp_071) %>%
    as.matrix()

  colnames(hm_mat2) <- c("MAL067", "WRAIR 1032", "MAL068", "MAL071")
  rownames(hm_mat2) <- btm_comp$gene_name

  comp_hm2 <- Heatmap(hm_mat2,

                      cluster_columns = F,
                      column_names_rot = 0,
                      column_names_centered = TRUE,

                      cluster_rows = F,
                      show_row_names = T,
                      row_title_rot = 0,
                      row_names_gp = grid::gpar(fontsize = fs),

                      heatmap_legend_param = list(title = "Signed Log P-value"),
                      col = col_fun2,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        curfdr <- 10 ^ (-abs(hm_mat2[i, j]))
                        if (curfdr < fdr_cut) {
                          grid.rect(x = x, y = y, width = width, height = height,
                                    gp = gpar(col = "#333333", fill = "transparent"))
                          if (curfdr < 0.01) {
                            grid.text("***", x, y, gp = gpar(fontsize = fsa), vjust = "top")
                          } else if (curfdr < 0.05) {
                            grid.text("**", x, y, gp = gpar(fontsize = fsa), vjust = "top")
                          } else {
                            grid.text("*", x, y, gp = gpar(fontsize = fsa), vjust = "top")
                          }
                        }
                      })


  pdf(here(glue("output/figures/single_gene/baseline/{btm_m}_baseline_comp_lpv.pdf")), width = 6, height = 6)
  draw(comp_hm2)
  dev.off()
  png(here(glue("output/figures/single_gene/baseline/{btm_m}_baseline_comp_lpv.png")), width = 600, height = 600)
  draw(comp_hm2)
  dev.off()
}



#################################################################
## No asterisks
#################################################################

for (i in 1:length(m0_genesets)) {
  btm_name <- names(genesets_btm[i])
  genes <- geneIds(genesets_btm[i])[[1]]
  btm_m <- str_split(btm_name, "\\(")[[1]] %>% last() %>% str_sub(1, -2)

  ngenes <- length(genes)
  fs <- 4
  if (ngenes < 100) {
    fs <- 6
  }
  if (ngenes < 50) {
    fs <- 8
  }

    btm_comp <- comp %>%
    filter(gene_name %in% genes) %>%
    arrange(logp_067)

  hm_mat2 <- btm_comp %>%
    dplyr::select(logp_067, logp_vahey, logp_068, logp_071) %>%
    as.matrix()

  colnames(hm_mat2) <- c("MAL067", "WRAIR 1032", "MAL068", "MAL071")
  rownames(hm_mat2) <- btm_comp$gene_name

  comp_hm2 <- Heatmap(hm_mat2,

                      cluster_columns = F,
                      column_names_rot = 0,
                      column_names_centered = TRUE,

                      cluster_rows = F,
                      show_row_names = T,
                      row_title_rot = 0,
                      row_names_gp = grid::gpar(fontsize = fs),

                      heatmap_legend_param = list(title = "Signed Log P-value"),
                      col = col_fun2,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        curfdr <- 10 ^ (-abs(hm_mat2[i, j]))
                        if (curfdr < fdr_cut) {
                          grid.rect(x = x, y = y, width = width, height = height,
                                    gp = gpar(col = "#333333", fill = "transparent"))
                        }
                      })


  pdf(here(glue("output/figures/single_gene/baseline/{btm_m}_baseline_comp_lpv_noast.pdf")), width = 6, height = 6)
  draw(comp_hm2)
  dev.off()
  png(here(glue("output/figures/single_gene/baseline/{btm_m}_baseline_comp_lpv_noast.png")), width = 600, height = 600)
  draw(comp_hm2)
  dev.off()
}

## annotations for the genesets
groupings <- read_csv(here("data/btm_annotation_table_LNC.csv")) %>%
  mutate(geneset = glue("{`Module title`} ({ID})")) %>%
  filter(geneset %in% m0_genesets)

for (m in unique(groupings$annotation)) {
  print(m)
  gs <- groupings %>%
    filter(annotation == m) %>%
    pull(geneset)

  genes <- geneIds(genesets_btm[gs]) %>% unlist() %>% unique()
  btm_m <- str_replace_all(m, " ", "_") %>%
    str_replace_all("/", "-")

  ngenes <- length(genes)
  print(ngenes)
  fs <- 2
  if (ngenes < 200) {
    fs <- 3
  }
  if (ngenes < 150) {
    fs <- 4
  }
  if (ngenes < 100) {
    fs <- 6
  }
  if (ngenes < 50) {
    fs <- 8
  }

    btm_comp <- comp %>%
    filter(gene_name %in% genes) %>%
    arrange(logp_067)

  hm_mat2 <- btm_comp %>%
    dplyr::select(logp_067, logp_vahey, logp_068, logp_071) %>%
    as.matrix()

  colnames(hm_mat2) <- c("MAL067", "WRAIR 1032", "MAL068", "MAL071")
  rownames(hm_mat2) <- btm_comp$gene_name

  comp_hm2 <- Heatmap(hm_mat2,

                      cluster_columns = F,
                      column_names_rot = 0,
                      column_names_centered = TRUE,

                      cluster_rows = F,
                      show_row_names = T,
                      row_title_rot = 0,
                      row_names_gp = grid::gpar(fontsize = fs),

                      heatmap_legend_param = list(title = "Signed Log P-value"),
                      col = col_fun2,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        curfdr <- 10 ^ (-abs(hm_mat2[i, j]))
                        if (curfdr < fdr_cut) {
                          grid.rect(x = x, y = y, width = width, height = height,
                                    gp = gpar(col = "#333333", fill = "transparent"))
                        }
                      })


  pdf(here(glue("output/figures/single_gene/baseline/annotations/{btm_m}_baseline_comp_lpv_noast.pdf")),
      width = 6, height = 6)
  draw(comp_hm2)
  dev.off()
  png(here(glue("output/figures/single_gene/baseline/annotations/{btm_m}_baseline_comp_lpv_noast.png")),
      width = 600, height = 600)
  draw(comp_hm2)
  dev.off()
}

m1 <- "Monocytes"
m2 <- "Dendritic Cells"
m3 <- "Cell Cycle"

gs1 <- groupings %>%
  filter(annotation == m1) %>%
  pull(geneset)
genes1 <- geneIds(genesets_btm[gs1]) %>% unlist() %>% unique()

gs2 <- groupings %>%
  filter(annotation == m2) %>%
  pull(geneset)
genes2 <- geneIds(genesets_btm[gs2]) %>% unlist() %>% unique()

gs3 <- groupings %>%
  filter(annotation == m3) %>%
  pull(geneset)
genes3 <- geneIds(genesets_btm[gs3]) %>% unlist() %>% unique()

length(genes1)
## 297
length(genes2)
## 129
length(genes3)
## 335
length(intersect(genes1, genes2))
## 17
length(intersect(genes1, genes3))
## 130
length(intersect(genes3, genes2))
## 8
length(intersect(genes1, genes2) %>% intersect(genes3))
## 8
intersect(genes1, genes2) %>% intersect(genes3)
## "CD86"   "CSF1R"  "IGSF6"  "ITGAX"
## "SIRPA"  "STAB1"  "TLR8"   "TYROBP"

all_genes <- intersect(genes1, genes2) %>% intersect(genes3)

comp %>%
  filter(gene_name %in% all_genes) %>%
  arrange(logp_067) %>%
  dplyr::select(gene_name, pval_067, pval_068, pval_071, pval_vahey)
