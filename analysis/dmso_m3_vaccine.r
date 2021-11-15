library(here)
library(magrittr)
library(glue)
library(tidyverse)

source(here("code", "runGSEA.r"))

meta_dt <- as_tibble(pData(mal067_eset))

## parameters set by parent: stim, time
stim="dmso"

################### DMSO M3 RTS,S vs Comparator ########################

form_main_vac_both_m0 <- "~plate + total_reads + age_weeks + vaccine"
form_main_vac_both_m3 <- "~plate + total_reads + age + vaccine"

for (time in c("M0", "M3")) {
  form_vac <- if_else(time == "M0",
                      form_main_vac_both_m0,
                      form_main_vac_both_m3)

  ## vaccine: select subset of expressionSet
  meta_dt %>%
    filter(stimulation == stim,
           visit == time) %>%
    mutate(stimulation = droplevels(stimulation),
           case = droplevels(case)) %$%
    col_id ->
    sample_vac


  ## generate subset data and drop extra stimulation levels
  mal_vac <- mal067_eset[, sample_vac]
  mal_vac$stimulation <- fct_drop(mal_vac$stimulation)
  mal_vac$case <- fct_drop(mal_vac$case)

  cam_vac <- runGSEA(mal_vac,
                     form_vac,
                     coef="vaccinertss")
  suffix <- paste(stim, time, sep="_")
  write_csv(cam_vac, file.path(here("output/"), glue("vaccine_both_{suffix}.csv")))

}
