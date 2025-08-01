install.packages("forcats")
library(forcats)
library(stringr)

hr_gapabs <- readRDS("results_gap_abs_diseases.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_gap <- readRDS("results_gap_diseases.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_res <- readRDS("results_res_diseases.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_resabs <- readRDS("results_res_abs_diseases.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_gapabs_rint <- readRDS("results_gap_abs_rint_diseases.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_gapabs_rint %>%
  filter(cases > 50) %>%
  group_by(names) %>%
  mutate(max_effect = max(effect, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(names = forcats::fct_reorder(names, max_effect, .desc = FALSE),
         model_num = str_extract(model, "Model[12]"),
         model_num = factor(model_num, levels = c("Model1", "Model2"), labels = c("Model1: sex + age + 20PCs", "Model2: Model1 + BMI + smoking + fasting h"))) %>%
  ggplot(aes(x = effect, y = names, color = model_num)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_continuous(limits = c(0.75, 1.51)) +
  facet_grid(class ~ type, scales = "free", space = "free") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(1, "lines")
  )
