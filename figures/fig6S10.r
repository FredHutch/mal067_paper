suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(glue)
})

library(mal067data)
data(mal067_meta_flow)

## All varying fields
flow_data <- read_tsv(here("data/170830-RTSS case control phenotyping.txt"), guess_max = 77000) %>%
  dplyr::select(`Batch #`, PTID, VISITNO, TESTDT, SAMP_ORD,
                VIABL1, RECOVR1, VIABL2, RECOVR2, ASSAYID,
                Sample, `$DATE`,
                `WELL ID`, `Sample Order`, `EXPERIMENT NAME`, Label,
                Count, `Include in analysis`, Lineage,
                `Subset of lineage`, Reliable) %>%
  filter(PTID != "FH12 (SAC207295)")  ## remove control samples

## Relevant fields?
flow_data <- read_tsv(here("data/170830-RTSS case control phenotyping.txt"), guess_max = 77000) %>%
  filter(PTID != "FH12 (SAC207295)") %>%
  dplyr::select(PTID, VISITNO, Label,
                Count, `Include in analysis`, Lineage,
                `Subset of lineage`, Reliable) %>%
  rename(ptid = PTID,
         visit = VISITNO,
         label = Label,
         ct = Count,
         include = `Include in analysis`,
         lineage = Lineage,
         lineage_subset = `Subset of lineage`,
         reliable = Reliable)

flow_cts <- flow_data %>%
  group_by(ptid, visit) %>%
  summarize(total_ct = sum(ct, na.rm = T),
            .groups = "drop")

rtss_flow <- flow_data %>%
  filter(lineage %in% c("CD14+CD16+ Intermediate Mono",
                        "CD14+CD16- Classical Mono",
                        "CD14loCD16hi Inflammatory Mono",
                        "Lymphocytes",
                        "Singlets")) %>%
  pivot_wider(c(ptid, visit),
              names_from = lineage,
              values_from = ct) %>%
  dplyr::rename(Intermediate = `CD14+CD16+ Intermediate Mono`,
                Classical = `CD14+CD16- Classical Mono`,
                Inflammatory = `CD14loCD16hi Inflammatory Mono`) %>%
  mutate(mono_ct = Intermediate + Classical + Inflammatory,
         mono_freq = mono_ct / Singlets) %>%
  left_join(mal067_meta_flow %>%
              dplyr::select(sid, case, site, age, vaccine, match) %>%
              distinct(),
            by = c("ptid" = "sid")) %>%
  dplyr::filter(!is.na(visit),
                vaccine == "rtss") %>%
  mutate(case = factor(if_else(case == "case", "Cases", "Controls"),
                       levels = c("Controls", "Cases")),
         visit = if_else(visit == "M0", "Month 0", "Month 3"))

## Monocyte frequency?
mono_box <- ggplot(rtss_flow,
       aes(x = case, y = mono_freq)) +
  geom_boxplot(show.legend = FALSE, outlier.color = NA) +
  geom_jitter(height = 0, width = 0.3) +
  theme_bw() +
  facet_wrap(~visit) +
  labs(x = "", y = "Monocyte Frequency")


## Inflammatory Monocytes
inf_box <- ggplot(rtss_flow,
       aes(x = case, y = Inflammatory / Singlets)) +
  geom_boxplot(show.legend = FALSE, outlier.color = NA) +
  geom_jitter(height = 0, width = 0.3) +
  theme_bw() +
  facet_wrap(~visit) +
  labs(x = "", y = "Inflammatory Monocyte Frequency")

## Inflammatory Monocytes/Lymphocyte Ratio
ratio_box <- ggplot(rtss_flow,
       aes(x = case, y = Inflammatory / Lymphocytes)) +
  geom_boxplot(show.legend = FALSE, outlier.color = NA) +
  geom_jitter(height = 0, width = 0.3) +
  theme_bw() +
  facet_wrap(~visit) +
  labs(x = "", y = "Inflammatory Monocyte/Lymphocyte Ratio")

ggsave(here("output/figures/fig6_supp10A_monocytes.pdf"), mono_box, width = 5, height = 8)
ggsave(here("output/figures/fig6_supp10B_inflammatory.pdf"), inf_box, width = 5, height = 8)
ggsave(here("output/figures/fig6_supp10C_ratio.pdf"), ratio_box, width = 5, height = 8)


library(lmerTest)

m0_data <- rtss_flow %>%
  filter(visit == "Month 0") %>%
  mutate(inf = Inflammatory / Singlets,
         ratio = Inflammatory / Lymphocytes)

m3_data <- rtss_flow %>%
  filter(visit == "Month 3") %>%
  mutate(inf = Inflammatory / Singlets,
         ratio = Inflammatory / Lymphocytes)

m0_mono <- lmer(mono_freq ~ case  + (1 | match), m0_data)
## p-value 0.697
m3_mono <- lmer(mono_freq ~ case  + (1 | match), m3_data)
## p-value 0.981
m0_inf <- lmer(inf ~ case  + (1 | match), m0_data)
## p-value 0.2188
m3_inf <- lmer(inf ~ case  + (1 | match), m3_data)
## p-value 0.114
m0_ratio <- lmer(ratio ~ case  + (1 | match), m0_data)
## p-value 0.1627
m3_ratio <- lmer(ratio ~ case  + (1 | match), m3_data)
## p-value 0.16





