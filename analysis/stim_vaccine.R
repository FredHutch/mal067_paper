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

min_gene_set <- 5
form_interact_vac_m0 <- "~plate + total_reads + age_weeks + stimulation*vaccine"
form_interact_vac_m3 <- "~plate + total_reads + age + stimulation*vaccine"


for (time in c("M3", "M0")) {
  for (stim in stims) {
  coeff <- glue("stimulation{stim}:vaccinertss")
    temp <- meta_dt %>%
      filter(stimulation == "dmso" | stimulation == stim,
             visit == time) %>%
      dplyr::select(pid, vaccine) %>%
      distinct()

    print(glue("{time} - {stim}: total pids - {nrow(temp)}"))
    print(glue("comp - {sum(temp$vaccine == 'comparator')}"))
    print(glue("rtss - {sum(temp$vaccine == 'rtss')}"))

    form_vac <- if_else(time == "M0",
                        form_interact_vac_m0,
                        form_interact_vac_m3)

    ## vaccine: select subset of expressionSet
    meta_dt %>%
      filter(stimulation == "dmso" | stimulation == stim,
             visit == time) %>%
      mutate(stimulation = droplevels(stimulation),
             case = droplevels(case)) %$%
      col_id ->
      sample_vac


    ## generate subset data and drop extra stimulation levels
    mal_vac <- mal067_eset[, sample_vac]
    mal_vac$stimulation <- fct_drop(mal_vac$stimulation)
    mal_vac$case <- fct_drop(mal_vac$case)

    cam_block <- runGSEA_block_pid(mal_vac,
                       form_vac,
                       coef=coeff,
                       dup_cor_file = here(glue("output/dupcor/{time}_{stim}_vaccine_dupcor.rds")))
    suffix <- paste(stim, time, sep="_")
    write_csv(cam_block, file.path(here("output/"), glue("vaccine_both_{suffix}.csv")))
  }
}
