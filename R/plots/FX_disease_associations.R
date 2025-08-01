install.packages("forcats")
library(forcats)
library(stringr)

hr_res <- readRDS("results_gap_abs_diseases.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

hr_res %>%
  filter(cases > 50) %>%
  group_by(names) %>%
  mutate(max_effect = max(effect, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(names = forcats::fct_reorder(names, max_effect, .desc = FALSE),
         model_num = str_extract(model, "Model[12]")) %>%
  ggplot(aes(x = effect, y = names, color = model_num)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_continuous(limits = c(0.7, 1.6)) +
  facet_grid(class ~ type, scales = "free", space = "free") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(1, "lines")
  )
