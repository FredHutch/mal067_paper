library(here)
library(glue)
library(tidyverse)
library(Biobase)

## load mal067 data
library(mal067data)
meta_dt <- as_tibble(pData(mal067_eset))

subject_data <- meta_dt %>%
  dplyr::select(pid, visit, site, age, vaccine, case) %>%
  distinct()

subject_data$vaccine[subject_data$pid == 3010] <- "comparator"
subject_data$vaccine[subject_data$pid == 3020] <- "rtss"

subject_data <- distinct(subject_data)

bag_data <- subject_data %>%
  filter(site == "BAGAMOYO")

man_data <- subject_data %>%
  filter(site == "MANHICA")

m0_data <- subject_data %>%
  filter(visit == "M0")

m3_data <- subject_data %>%
  filter(visit == "M3")


bag_data %>%
  filter(visit == "M0") %>%
  dplyr::select(age, vaccine, case) %>%
  table()

man_data %>%
  filter(visit == "M0") %>%
  dplyr::select(age, vaccine, case) %>%
  table()

bag_data %>%
  filter(visit == "M3") %>%
  dplyr::select(age, vaccine, case) %>%
  table()

man_data %>%
  filter(visit == "M3") %>%
  dplyr::select(age, vaccine, case) %>%
  table()
