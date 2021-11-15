library(here)
library(magrittr)
library(glue)
library(tidyverse)

source(here("code", "runGSEA.r"))

meta_dt <- as_tibble(pData(mal067_eset))

stims <- c("csp", "hbs", "ama1")
time <- "M3"


form_visit_main <- "~plate + total_reads + age + visit*stimulation"
################### vaccine and disease  analysis ########################

for (stim in stims) {
  temp <- meta_dt %>%
    filter(stimulation == "dmso" | stimulation == stim,
           vaccine == "comparator") %>%
    dplyr::select(pid, stimulation) %>%
    distinct()

  print(glue("{stim}: total pids - {n_distinct(temp$pid)}"))
  print(glue("stim - {sum(temp$stimulation == stim)}"))
  print(glue("DMSO - {sum(temp$stimulation == 'dmso')}"))
}

for (stim in stims) {
  coeff <- glue("visitM3:stimulation{stim}")

  ## vaccine: select subset of expressionSet
  meta_dt %>%
    filter(stimulation == "dmso" | stimulation == stim,
           vaccine == "comparator") %$%
    col_id ->
    sample_m3m0


  ## generate subset data and drop extra stimulation levels
  mal_m3m0 <- mal067_eset[, sample_m3m0]
  mal_m3m0$stimulation <- fct_drop(mal_m3m0$stimulation)
  mal_m3m0$case <- fct_drop(mal_m3m0$case)

  cam_m3m0 <- runGSEA_block_pid(mal_m3m0,
                                form_visit_main,
                                coef = coeff,
                                dup_cor_file = here(glue("output/dupcor/comp_{stim}_time_dupcor.rds")))

  write_csv(cam_m3m0, file.path(here("output/"), glue("M3-M0_comp_{stim}.csv")))
}
