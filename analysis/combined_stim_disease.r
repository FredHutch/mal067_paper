library(here)
library(magrittr)
library(glue)
library(tidyverse)

library(mal067data)
data(mal067_eset)

source(here("code", "runGSEA.r"))

meta_dt <- as_tibble(pData(mal067_eset))

## parameters
stims <-  c("csp", "hbs", "ama1")
time <- "M0"
coeff_case <- "casecase"

min_gene_set <- 5
form_main_dis_case <- "~ plate + total_reads + age_weeks + stimulation*case"

for (stim in stims) {
  temp <- meta_dt %>%
    filter(stimulation == "dmso" | stimulation == stim,
           visit == time,
           case != "neither") %>%
    dplyr::select(pid, stimulation) %>%
    distinct()

  print(glue("{stim}: total pids - {n_distinct(temp$pid)}"))
  print(glue("stim - {sum(temp$stimulation == stim)}"))
  print(glue("DMSO - {sum(temp$stimulation == 'dmso')}"))
}

for (stim in stims) {
  ## disease: select subset of expressionSet
  meta_dt %>%
    filter(stimulation == "dmso" | stimulation == stim,
           visit == time,
           case != "neither") %$%
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
}
