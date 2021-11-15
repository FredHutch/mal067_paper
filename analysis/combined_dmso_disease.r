library(here)
library(magrittr)
library(glue)
library(tidyverse)

source(here("code", "runGSEA.r"))

meta_dt <- as_tibble(pData(mal067_eset))

## parameters set by parent: stim, time
stim <- "dmso"
time <- "M0"

form_main_dis_case <- "~ plate + total_reads + age_weeks + case"
coeff_case <- "casecase"

## disease: select subset of expressionSet
meta_dt %>%
  filter(stimulation == stim,
         case == "case" | case == "control",
         visit == time) %$%
  col_id ->
  sample_ids


## generate subset data and drop extra stimulation levels
mal_dis <- mal067_eset[, sample_ids]
mal_dis$stimulation <- fct_drop(mal_dis$stimulation)
mal_dis$case <- fct_drop(mal_dis$case)

cam_dis_rtss <- runGSEA(mal_dis,
                        form_main_dis_case,
                        coef=coeff_case)

suffix <- paste(stim, time, sep="_")
write_csv(cam_dis_rtss, file.path(here("output/"), glue("disease_combined_{suffix}.csv")))
