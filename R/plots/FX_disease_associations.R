install.packages("forcats")
library(forcats)
library(stringr)
library(dplyr)
library(ggplot2)
library(stringr)
install.packages("paletteer")
library(paletteer)


data <- bind_rows(readRDS("results_gap_abs_diseases_1y.rds"),
                  readRDS("results_gap_diseases_1y.rds"),
                  readRDS("results_gap_abs_proteins_diseases_1y.rds"),
                  readRDS("results_gap_proteins_diseases_1y.rds")) %>%
  left_join(list_diseases, by = c("disease" = "p")) %>%
  mutate(type = case_when(type == "logistic" ~ "Prevalent",
                          type == "cox" ~ "Incident"),
         type = factor(type, levels = c("Prevalent", "Incident")),
         term = case_when(str_detect(term, "abs") ~ "Dysregulation",
                          TRUE ~ "Acceleration"),
         p_val = p.adjust(p.value, method = "fdr", n = max(27, n())),
         sig = case_when(p_val < 0.05 ~ "*",
                         TRUE ~ ""))

plot_assoc <- data %>%
  filter(cases > 40) %>%
  filter(model %in% c("Model1", "Model2")) %>%
  mutate(names = forcats::fct_reorder(names, cases, .desc = FALSE),
         model_num = str_extract(model, "Model[12]"),
         model_num = factor(model_num, levels = c("Model2", "Model1"), labels = c( "Model2: Model1 + BMI + smoking + fasting h",
                                                                                     "Model1: sex + age + 20PCs"))) %>%
  ggplot(aes(x = effect, y = names, color = model_num, label = sig)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_text(data = subset(data, model == "Model2"), position = position_dodge(width = 0.5),
            color = "black", hjust = 0, size = 5,
            aes(x = hi95 + 0.05)) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_continuous(limits = c(0.5, 1.7), n.breaks = 4) +
  labs(x = "OR / HR", color = "") +
  guides(
    color = guide_legend(
      nrow       = 2,
      byrow      = TRUE,
      override.aes = list(size = 6),
      reverse = T
    )) +
  scale_color_manual(values = c("#8BA6A9", "#A7CECB")) +
  facet_grid(class ~  term + type, scales = "free", space = "free") +
  theme_pubclean() +
  theme(strip.background    = element_rect(fill   = "white",
                                           colour = "black",
                                           size   = 0.5),
    legend.position = "bottom",
    axis.title.y = element_blank(),
    axis.line.y = element_line(colour = "black"),
    strip.text.y = element_text(angle = 0), strip.text.y.right = element_text(angle = 270),
    panel.spacing.y = unit(1, "lines")
  )

ggsave("plots/associations_diseases.png", plot_assoc, width = 8, height = 8)

prots <- readRDS("/mnt/project/results_proteins_disease.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))



### Plots counts
c_dis <- readRDS("/mnt/project/counts_diseases_circadian.rds") %>%
  left_join(list_diseases, by = c("disease" = "p"))

library(scales)

p_counts <- c_dis %>%
  mutate(
    lab_name = paste0(names, " (n=", total, ")"),
    lab_name = fct_reorder(lab_name, total)
  ) %>%
  pivot_longer(c(prevalent, incident),
               names_to  = "case_type",
               values_to = "count") %>%
  mutate(case_type = factor(case_type,
                            levels = c("prevalent", "incident"),
                            labels = c("Prevalent", "Incident"))) %>%
  ggplot(aes(x = count, y = lab_name, fill = case_type, order = case_type)) +
  geom_col(position = position_fill(reverse = TRUE), width = 1, color = "black") +
  geom_text(aes(label = count),
            position = position_fill(vjust = 0.5, reverse = TRUE),
            color = "black", size = 4) +
  scale_x_continuous() +
  facet_grid(rows = vars(class),
             scales = "free_y",
             space  = "free") +
  labs(
    x    = "Proportion",
    y    = NULL,
    fill = "Case type"
  ) +
  scale_fill_manual(values = c("#FAF3DD", "#C8D5B9")) +
  theme_classic2() +
  theme(
    strip.text.y = element_text(angle = 0),
    strip.text.y.right = element_text(angle = 270),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

ggsave("plots/counts_diseases.png", p_counts, width = 10, height = 8)
