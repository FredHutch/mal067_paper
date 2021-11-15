library(Biobase)
library(DT)
library(GSEABase)
library(limma)
library(tidyverse)
library(data.table)

## load mal067 data
library(mal067data)

## set up GSEA analysis
min_gene_set <- 5
btm_gtm <- getGmt(here("data/BTM_for_GSEA_20131008.gmt"))
btm_gtm <- btm_gtm[!str_detect(names(btm_gtm), "TBA")]
btm_ind <- ids2indices(geneIds(btm_gtm), rownames(mal067_eset))
btm_ind <- btm_ind[sapply(btm_ind, length) > min_gene_set]

fdr_cut <- 0.2

runGSEA <- function(mal_dat, form, coef, btm=btm_ind) {

  mal_dat$pid <- droplevels(as.factor(mal_dat$pid))

  ## voom the data with the linear model
  design_vac <- model.matrix(formula(form), mal_dat)
  mal_voom_sm <- voom(mal_dat, design=design_vac)

  ## run camera
  cam_vac <- as_tibble(camera(mal_voom_sm, btm, design=design_vac, contrast=coef),
                       rownames="geneset")

  return(cam_vac)
} ## runGSEA

runGSEA_block_pid <- function(mal_dat, form, coef, btm=btm_ind, dup_cor_file = NULL) {

  mal_dat$pid <- as.factor(mal_dat$pid)

  design_vac <- model.matrix(formula(form), mal_dat)

  ## calculate duplicate correlation
  dupcor_imm <- NULL
  if (!is.null(dup_cor_file) & file.exists(dup_cor_file)) {
    dupcor_imm <- read_rds(dup_cor_file)
  } else {
    dupcor_imm <- duplicateCorrelation(mal_dat,
                         design_vac,
                         block = mal_dat$pid)
    if (!is.null(dup_cor_file)) {
      write_rds(dupcor_imm, dup_cor_file)
    }
  }

  lm_res <- lmFit(mal_dat,
                  design_vac,
                  correlation = dupcor_imm)
  eb_res <- eBayes(lm_res)

  ## run camera
  cam_res <- as_tibble(cameraPR(eb_res$t[,coef], btm),
                          rownames = "geneset")

  return(cam_res)
} ## runGSEA_block_pid






