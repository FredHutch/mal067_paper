library(data.table)
library(tidyverse)
library(Biobase)
library(GSEABase)
library(magrittr)
library(janitor)
library(MAL068)
library(SummarizedExperiment)
library(pheatmap)
library(here)
library(impute)
library(qusage)
library(broom)


fdr_cut <- 0.2
## significant genesets (down-selection)
sig_genesets <- read_csv(here("output/fig2_sig_genesets_downselection.csv")) %>%
  pull(geneset)

m3_disease <- read_csv(here("output/disease_rtss_dmso_M3.csv")) %>%
  filter(geneset %in% sig_genesets) %>%
  dplyr::select(geneset, NGenes, Direction, PValue, FDR)

m3_genesets <- m3_disease %>%
  mutate(FDR = p.adjust(PValue, "fdr")) %>%
  filter(FDR <= fdr_cut) %>%
  pull(geneset)

## set up GSEA analysis
min_gene_set <- 5
btm_gtm <- getGmt(here("data/BTM_for_GSEA_20131008.gmt"))


genesets_btm <- btm_gtm[m3_genesets]
genesets_ind <- btm_ind[m3_genesets] %>% unlist() %>% unique()

genesets_genes <- geneIds(genesets_btm) %>% unlist() %>% unique()



# Library to load RNAseq data
data(MAL068_rnaseq)
# Remove post challenge timepoints
MAL068_rnaseq <- MAL068_rnaseq[, MAL068_rnaseq$visit_day == 78]
# Convert the visit_day to factor
# This is easier for doing testing at the end
MAL068_rnaseq$visit_day <- as.factor(MAL068_rnaseq$visit_day)
MAL068_rnaseq$infection <- as.factor(MAL068_rnaseq$infection)

# Number of missing genes per sample
nb_missing_genes <- apply(is.na(assay(MAL068_rnaseq)), 2, sum)
missing_all_genes <- names(nb_missing_genes[nb_missing_genes == nrow(MAL068_rnaseq)])

# I remove samples that have too many missing data
MAL068_rnaseq <- MAL068_rnaseq[!colnames(MAL068_rnaseq) %in% missing_all_genes]

rrr_data <- MAL068_rnaseq[rowData(MAL068_rnaseq)$gene_name %in% genesets_genes,
                          MAL068_rnaseq$arm == "RRR"]


design_dis <- model.matrix(~ age + sex + infection,
                           colData(rrr_data))


dupcor_dis <- duplicateCorrelation(assay(rrr_data),
                                   design_dis,
                                   block = rrr_data$ptid)

lm_dis <- lmFit(assay(rrr_data),
                design_dis,
                correlation = dupcor_dis)
eb_dis <- eBayes(lm_dis)


mal068_res <- topTable(eb_dis,
                       coef = "infection1",
                       number = "inf",
                       sort.by = "none") %>%
       arrange(P.Value)

mal068_res <- mal068_res %>%
  mutate(gene_id = rownames(mal068_res)) %>%
  as_tibble() %>%
  left_join(as_tibble(rowData(MAL068_rnaseq)), by = "gene_id")

write_csv(mal068_res, here("output/single_gene/mal068_results.csv"))

