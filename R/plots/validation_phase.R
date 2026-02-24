library(tidyverse)
library(paletteer)
library(ggpmisc)
library(ggrepel)

df_r2 <- readRDS("tables/variance_covariates.rds")%>%
  mutate(pfdr = p.adjust(p.value_time_day)) %>%
  filter(pfdr < 0.05)

df_eff <- readRDS("tables/harmonic_models.rds") %>%
  mutate(pfdr = p.adjust(pvalue_h)) %>%
  filter(pfdr < 0.05)

d1 <- inner_join(df_r2, df_eff, by = c("Biomarker", "Label", "Type", "FID")) %>%
  left_join(data.table::fread("data/olink_assay.dat"), by = c("Biomarker" = "Assay"))

set <- d1 %>% filter(amplitude_24hfreq > 0.1 & t_r2_time_day > 0.01)

sum(set$Type == "Proteins")
# [1] 86

### VALIDATION OLINK

ref <- data.table::fread("~/OneDrive - Nexus365/projects/circadian/clara_results_top19.csv") %>%
  inner_join(d1, by = c("Assay" = "Biomarker")) %>%
  mutate(acro_1h_adj = acrophase_24hfreq +
           (((acrophase_hr - acrophase_24hfreq + 12) %% 24) - 12)
         ) %>%
  mutate(dif_amp = amplitude_24hfreq - amplitude,
         dif_phase = acrophase_24hfreq - acrophase_hr)

t <- cor.test(ref$acrophase_24hfreq, ref$acro_1h_adj, method = "pearson")
# Pearson's product-moment correlation
#
# data:  ref$acrophase_24hfreq and ref$acrophase_hr
# t = 11.949, df = 16, p-value = 2.189e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8638553 0.9808869
# sample estimates:
#       cor
# 0.9482773

cor.test(ref$amplitude_24hfreq, ref$amplitude, method = "pearson")

to_rad <- function(x, period = 24) (x %% period) / period * 2*pi

t1 <- to_rad(ref$acrophase_24hfreq)
t2 <- to_rad(ref$acrophase_hr)

x <- cbind(t1, t2)
f <- BAMBI::circ_cor(x, type = "js")


p_olink <- ggplot(ref,
                  aes(x = acrophase_24hfreq, y = acro_1h_adj)) +
  geom_point() +
  geom_label_repel(
    data = ref %>% filter(abs(dif_phase) > 4),
    size = 2.5,
    aes(label = Label),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  ggpmisc::stat_poly_eq(
    mapping    = aes(label = paste("italic(R) ==", round(t$estimate, 2))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("italic(p) ==", formatC(t$p.value, format = "e", digits = 0))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.90,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("n == 18")),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 3.5,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "UKB acrophase (h)", y = "TREASURE acrophase (h)") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray")) +
  ggtitle("C") +
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave("plots/acrophase_validation.png", p_olink, width = 4, height = 4.5)

t_amp <- cor.test(ref$amplitude_24hfreq, ref$amplitude, method = "pearson")

a_olink <- ref %>%
  ggplot(aes(
    x     = amplitude_24hfreq,
    y     = amplitude
  )) +
  geom_point() +
  geom_label_repel(
    data = ref %>% filter(abs(dif_amp) > 0.3),
    size = 2.5,
    aes(label = Label),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  ggpmisc::stat_poly_eq(
    mapping    = aes(label = paste("italic(R) ==", round(t_amp$estimate, 2))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("italic(p) ==", formatC(t_amp$p.value, format = "e", digits = 0))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.90,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("n == 18")),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 3.5,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "UKB amplitude", y = "TREASURE amplitude") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray")) +
  ggtitle("F") +
  theme(
    plot.title = element_text(face = "bold")
  )


a_olink

ggsave("plots/amplitude_validation.png", a_olink, width = 4, height = 4)




# Specht et al

# circ_diff24 <- function(a, b) {
#   # returns difference in hours in [-12,12]
#   delta <- (a - b + 12) %% 24 - 12
#   abs(delta)
# }


specht <- readxl::read_xlsx("data/1-s2.0-S2352721823002401-mmc1.xlsx", skip = 1) %>%
  rename(
    acro_fund   = `acrophase of fundamental harmonic in 2-harmonic fit  (time to DLMO in hours)`,
    acro_1st    = `acrophase of first harmonic in 2-harmonic fit  (time to DLMO in hours)`,
    acro_1h = `acrophase of 1-harmonic fit (time to DLMO in hours)`,
    amp_1h = `amplitude of 1 harmonic fit`) %>%
  inner_join(d1) %>%
  filter(amplitude_24hfreq > 0.1 & t_r2_time_day > 0.01) %>%
  mutate(acro_1h = case_when(acrophase_24hfreq > 18 & acro_1h < 10 ~ acro_1h + 24,
                                  TRUE ~ acro_1h),
         acro_1h_adj = acrophase_24hfreq +
           (((acro_1h - acrophase_24hfreq + 12) %% 24) - 12)) %>%
  group_by(Label) %>% mutate(n = n()) %>%
  mutate(Label2 = case_when(n > 1 ~ paste0(Label,"_", row_number()),
                            TRUE ~ Label)) %>%
  mutate(dif_amp = amplitude_24hfreq - amp_1h,
         dif_phase = acrophase_24hfreq - acro_1h)


t_pha_specht_1h <- cor.test(specht$acrophase_24hfreq, specht$acro_1h_adj, method = "pearson")

t_pha_specht_1h_CIRC <- cor.test(specht$acrophase_24hfreq[specht$`Rhythmic Category` == "circadian"], specht$acro_1h_adj[specht$`Rhythmic Category` == "circadian"], method = "pearson")

t_pha_specht_1h_RHY <- cor.test(specht$acrophase_24hfreq[specht$`Rhythmic Category` == "diurnal"], specht$acro_1h_adj[specht$`Rhythmic Category` == "diurnal"], method = "pearson")

length(unique(specht$Gene))
# [1] 44

p_pha <- specht %>%
  ggplot(aes(x = acrophase_24hfreq, y = acro_1h_adj)) +
  geom_point() +
  geom_label_repel(
    data = specht %>% filter(abs(dif_phase) > 4),
    size = 2.5,
    aes(label = Label2),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  ggpmisc::stat_poly_eq(
    mapping    = aes(label = paste("italic(R) ==", round(t_pha_specht_1h$estimate, 2))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("italic(p) ==", formatC(t_pha_specht_1h$p.value, format = "e", digits = 0))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.90,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("n == 53")),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 3.5,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(x = "Olink/UKB acrophase (h)", y = "SomaScan acrophase (h)") +
  scale_color_manual(
    values = paletteer_dynamic("cartography::green.pal", 2),
    guide = guide_legend(ncol = 2)) +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray"),
        legend.position = "bottom") +
  ggtitle("D") +
  theme(
    plot.title = element_text(face = "bold")
  )

p_pha


t_amp_specht_1h <- cor.test(specht$amplitude_24hfreq, specht$amp_1h, method = "pearson")


p_sp_amp <- specht %>%
  ggplot(aes(x = amplitude_24hfreq, y = amp_1h)) +
  geom_point() +
  geom_label_repel(
    data = specht %>% filter(abs(dif_amp) > 0.3),
    size = 2.5,
    aes(label = Label2),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  ggpmisc::stat_poly_eq(
    mapping    = aes(label = paste("italic(R) ==", round(t_amp_specht_1h$estimate, 2))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("italic(p) ==", formatC(t_amp_specht_1h$p.value, format = "e", digits = 0))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.90,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("n == 53")),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 3.5,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  #scale_y_continuous(limits = c(0, 30)) +
  labs(x = "Olink/UKB amplitude", y = "SomaScan amplitude") +
  scale_color_manual(
    values = paletteer_dynamic("cartography::green.pal", 2),
    guide = guide_legend(ncol = 2)) +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray"),
        legend.position = "bottom") +
  ggtitle("G") +
  theme(
    plot.title = element_text(face = "bold")
  )
p_sp_amp

