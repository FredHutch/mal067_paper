suppressPackageStartupMessages({
  library(mal067data)
  library(MAL067Package)
  library(tidyverse)
  library(here)
  library(glue)
  library(COMPASS)
  library(cowplot)
})

data("mal067_ics_counts")
mal067_ics_counts = as_tibble(mal067_ics_counts)
mal067_ics_counts %>%
  select(Stim) %>%
  table

CD4_ICS_Primary_Magnitude = mal067_ics_counts %>%
  dplyr::filter(Stim %in% c("negctrl","CSP"),
         Parent == "/S/Exclude/14-/Lv/L/3+/4+",
         Population == "/S/Exclude/14-/Lv/L/3+/4+/IL2\\TNFa\\CD154",
         `Sample Order` != 99)

CD4_ICS_Primary_Magnitude = CD4_ICS_Primary_Magnitude %>%
  mutate(primary_magnitude = (Count+0.1) / (ParentCount))

cd4_mag <- read_rds(here("data/PRIMARY_CD4_MAGNITUDE.rds"))


lowq = CD4_ICS_Primary_Magnitude %>% group_by(pid) %>% summarize(lowtcell = any(ParentCount < 20000))

primary_compass_fit = readRDS(here("data/primary_compass_fit_filtered.rds"))

#merge with low quality scores for T cell counts
primary_scores = scores(primary_compass_fit) %>%
  inner_join(lowq %>%
               select(pid,lowtcell) %>%
               group_by(pid) %>%
               summarize(lowtcell = any(lowtcell)) %>%
               unique())

data("mal067_rx")
primary_scores  = mal067_rx %>% select(PTID, VISITNO,
                                       malaria_transmision,
                                       malaria_before_m3.5) %>%
  inner_join(primary_scores)


fig4a = primary_scores %>%
  dplyr::filter(vaccine2c=="rtss") %>%
  ggplot() +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(height=0) +
  aes(y=PFS,x=VISITNO)+
  labs(x = "Visit", y = "Polyfunctionality score") +
  theme_bw()

ggsave(here("output/figures/Fig4a.png"),
       fig4a,
       height = 6, width = 8)
ggsave(here("output/figures/Fig4a.pdf"),
       fig4a,
       height = 6, width = 8)


plot_data <- inner_join(cd4_mag, lowq)

fig4b = ggplot(plot_data, aes(VISITNO, difference_trunc)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(height = 0) +
  scale_y_sqrt() +
  labs(x = "Visit", y = "% CD4+ expressing IL2+ or TNFa or CD154") +
  theme_bw()

ggsave(here("output/figures/Fig4b.png"),
       fig4b,
       height = 6, width = 8)
ggsave(here("output/figures/Fig4b.pdf"),
       fig4b,
       height = 6, width = 8)




m0_pfs <- primary_scores %>%
  filter(vaccine2c == "rtss", VISITNO == "M0") %>%
  pull(PFS)
m3_pfs <- primary_scores %>%
  filter(vaccine2c == "rtss", VISITNO == "M3") %>%
  pull(PFS)

wilcox.test(m0_pfs, m3_pfs)
## 0.0002287

## library(nlme)
library(lmerTest)

visit_scores <- primary_scores %>%
  filter(vaccine2c == "rtss",
         VISITNO %in% c("M0", "M3")) %>%
  mutate(pid = as.factor(pid))
visit_lmer <- lmer(PFS ~ VISITNO  + (1 | pid), visit_scores)
## p-value 0.00469

inf_scores <- primary_scores %>%
  filter(vaccine2c == "rtss",
         malaria %in% c(0, 1),
         VISITNO == "M3") %>%
  mutate(VISITNO = as.factor(VISITNO))
inf_lmer2 <- lmer(PFS ~ malaria + (1 | match_id_m12), inf_scores)
## p-value 0.126


cd4_scores <- inner_join(cd4_mag, lowq)

visit_cd4 <- cd4_scores %>%
  filter(vaccine2c == "rtss",
         VISITNO %in% c("M0", "M3")) %>%
  mutate(pid = as.factor(pid))
visit_lmer_cd4 <- lmer(difference_trunc ~ VISITNO  + (1 | pid), visit_cd4)
## p-value 0.127

inf_cd4 <- cd4_scores %>%
  filter(vaccine2c == "rtss",
         malaria %in% c(0, 1),
         VISITNO == "M3") %>%
inf_lmer_cd42 <- lmer(difference_trunc ~ malaria + (1 | match_id_m12), inf_cd4)
## p-value 0.083



m0_mag <- cd4_scores %>%
  filter(vaccine2c == "rtss", VISITNO == "M0") %>%
  pull(difference_trunc)
m3_mag <- cd4_scores %>%
  filter(vaccine2c == "rtss", VISITNO == "M3") %>%
  pull(difference_trunc)

wilcox.test(m0_mag, m3_mag)
