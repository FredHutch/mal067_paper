suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(glue)
})

out_folder <- "output"

match_data <- read_csv(here("data/RTSSFullPheno_20211001.csv"))

simplified_data <- match_data %>%
  dplyr::select(pid, site, agec, rna_seq_casecon.m12,
         rna_seq_match_id_m12, date.M0,
         m3.5_Start_m3_date_plus_14,
         tomalaria, malaria, sex, age,
         vaccine2c, malaria_before_m3.5, end_m3.5_dec13,
         time_to_followup, malaria_transmision) %>%
  distinct() %>%
  dplyr::rename(casecon = rna_seq_casecon.m12,
         match_id = rna_seq_match_id_m12) %>%
  mutate(m0 = as.Date(date.M0, format = "%d/%m/%Y"),
         end_dec13 = as.Date(end_m3.5_dec13, format = "%d/%m/%Y"),
         m3_plus_14 = as.Date(m3.5_Start_m3_date_plus_14, format = "%d/%m/%Y")) %>%
  arrange(match_id)


## Matched by vaccine, site, and age category (one case per match id)
match_summary <- simplified_data %>%
  filter(!is.na(match_id)) %>%
  group_by(match_id) %>%
  summarize(nSubj = n(),
            vacc = unique(vaccine2c),
            site = unique(site),
            agec = unique(agec),
            n_male = sum(sex == "M"),
            n_female = sum(sex == "F"),
            m0_range = as.numeric(max(m0) - min(m0)),
            end_dec13_range = as.numeric(max(end_dec13) - min(end_dec13)),
            m3_plus_14_range = as.numeric(max(m3_plus_14) - min(m3_plus_14)),
            min_age = min(age),
            max_age = max(age),
            age_range = max_age - min_age,
            min_tomalaria = min(tomalaria),
            max_tomalaria = max(tomalaria),
            tomalaria_range = max_tomalaria - min_tomalaria)


extract_table <- simplified_data %>%
  filter(!is.na(match_id)) %>%
  dplyr::rename(subject = pid,
                ageCategory = agec,
                caseControl = casecon,
                startDate = m0,
                vaccine = vaccine2c) %>%
  dplyr::select(subject, match_id, site, ageCategory, sex,
                vaccine, startDate, time_to_followup, caseControl) %>%
  arrange(match_id, caseControl, startDate)

write_csv(extract_table, here(out_folder, "supp_tables/match_summary.csv"))
