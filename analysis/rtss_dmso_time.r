library(here)
library(magrittr)
library(glue)
library(tidyverse)

source(here("code", "runGSEA.r"))

meta_dt <- as_tibble(pData(mal067_eset))

## parameters set by parent: stim, time
stim="dmso"
time="M3"

################### RTS,S M3 vs M0 ########################

form_visit_main <- "~plate + total_reads + age + visit"
coeff <- "visitM3"

## vaccine: select subset of expressionSet
meta_dt %>%
  filter(stimulation == stim,
         vaccine == "rtss") %$%
  col_id ->
  sample_m3m0


## generate subset data and drop extra stimulation levels
mal_m3m0 <- mal067_eset[, sample_m3m0]
mal_m3m0$stimulation <- fct_drop(mal_m3m0$stimulation)
mal_m3m0$case <- fct_drop(mal_m3m0$case)

cam_m3m0 <- runGSEA_block_pid(mal_m3m0,
                   form_visit_main,
                   coef = coeff,
                   dup_cor_file = here(glue("output/dupcor/rtss_dmso_time_dupcor.rds")))

write_csv(cam_m3m0, file.path(here("output/"), "M3-M0_rtss_dmso.csv"))
