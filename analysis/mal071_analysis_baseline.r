suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(janitor)
  library(MAL071)
  library(SummarizedExperiment)
  library(pheatmap)
  library(here)
  library(impute)
  library(qusage)
  library(broom)
})

fdr_cut <- 0.2

# Library to load RNAseq data
data(MAL071_rnaseq)
# Only look at day of challenge timepoint
MAL071_rnaseq <- MAL071_rnaseq[, MAL071_rnaseq$visit_day == -7]
# Convert the visit_day to factor
# This is easier for doing testing at the end
MAL071_rnaseq$visit_day <- as.factor(MAL071_rnaseq$visit_day)
MAL071_rnaseq$infection <- as.factor(MAL071_rnaseq$infection)

# Number of missing genes per sample
nb_missing_genes <- apply(is.na(assay(MAL071_rnaseq)), 2, sum)
missing_all_genes <- names(nb_missing_genes[nb_missing_genes == nrow(MAL071_rnaseq)])

# I remove samples that have too many missing data
MAL071_rnaseq <- MAL071_rnaseq[,
                               !colnames(MAL071_rnaseq) %in% missing_all_genes]


# Impute missing values
imputed_data <- impute.knn(assay(MAL071_rnaseq))
# Replace the old matrix with the imputed one
assay(MAL071_rnaseq) <- imputed_data$data

# Set up btm for camera analysis
btm <- read.gmt(here("data/BTM_for_GSEA_20131008.gmt"))
# cleaning the names
## btm <- clean_names(btm)
names(btm) <- make_clean_names(names(btm))
# remove modules with TBA in name
btm <- btm[!str_detect(names(btm), "tba")]
# Convert gene names to gene indices
gene_ids <- ids2indices(btm, rowData(MAL071_rnaseq)$gene_name)
# Number of genes per module
n_module <- sapply(gene_ids, length)
# Keep modules with at least 5 genes
gene_ids <- gene_ids[n_module > 5]

rrr_data <- MAL071_rnaseq[, MAL071_rnaseq$arm == "RRR"]



design_inf_rrr <- model.matrix(~ age + infection, colData(rrr_data))

cam_inf_rrr <- camera(assay(rrr_data), index = gene_ids,
                  design = design_inf_rrr, contrast = "infection1")
cam_inf_rrr <-as_tibble(cam_inf_rrr, rownames = "module") %>%
  mutate(slqvalue = log10(FDR) * if_else(Direction == "Up", -1, 1),
         is_significant = FDR < fdr_cut)
write_csv(cam_inf_rrr, here("output/MAL071_RRR_baseline_module_disease.csv"))


