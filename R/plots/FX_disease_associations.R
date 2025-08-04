install.packages("forcats")
library(forcats)
library(stringr)
library(dplyr)
library(ggplot2)
library(stringr)

hr_gapabs <- readRDS("results_gap_abs_diseases_1y.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_gap <- readRDS("results_gap_diseases_1y.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_prots <-  readRDS("results_gap_abs_proteins_diseases_1y.rds") %>%
  left_join(list_diseases, by = c("disease" = "p")) %>%
  mutate(model = "Model3")

#
# hr_res <- readRDS("/mnt/project/results_res_diseases.rds") %>%
#   left_join(list_diseases, by = c("disease" = "p"))
#
# hr_resabs <- readRDS("/mnt/project/results_res_abs_diseases.rds") %>%
#   left_join(list_diseases, by = c("disease" = "p"))
#
# hr_gapabs_rint <- readRDS("results_gap_abs_rint_diseases.rds") %>%
#   left_join(list_diseases, by = c("disease" = "p"))

bind_rows(hr_gap, hr_gapabs, hr_prots) %>%
#bind_rows(hr_res, hr_resabs) %>%
#hr_gapabs %>%
  mutate(type = case_when(type == "logistic" ~ "Prevalent",
                          type == "cox" ~ "Incident"),
         type = factor(type, levels = c())
         term = case_when(str_detect(term, "abs") ~ "Dysregulation",
                          TRUE ~ "Acceleration"),
         p_val = p.adjust(p.value, method = "fdr"),
         sig = case_when(p_val < 0.05 ~ "*",
                         TRUE ~ "")) %>%
  filter(cases > 50) %>%
  group_by(names) %>%
  mutate(max_effect = max(effect, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(names = forcats::fct_reorder(names, max_effect, .desc = FALSE),
         model_num = str_extract(model, "Model[123]"),
         model_num = factor(model_num, levels = c("Model1", "Model2", "Model3"), labels = c("Model1: sex + age + 20PCs",
                                                                                  "Model2: Model1 + BMI + smoking + fasting h",
                                                                                  "Model3: Model2 + 20 top proteins"))) %>%
  ggplot(aes(x = effect, y = names, color = model_num, label = sig)) +
  geom_point(position = position_dodge(width = 0.5), alpha = 0.8) +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_text(color = "black", aes(hjust = 10)) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_continuous(limits = c(0.5, 1.7)) +
  facet_grid(class ~  term + type, scales = "free", space = "free") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(1, "lines")
  )



prots <- readRDS("/mnt/project/results_proteins_disease.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))
