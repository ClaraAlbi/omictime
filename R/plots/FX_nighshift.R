library(tidyverse)

results <- readRDS("data_share/results_associations_phenotypes_CA.rds") %>%
  filter(str_detect(term, "chrono_Nightshift")) %>%
  mutate(FDR = p.adjust(p.value),
         term = str_remove(term, "chrono_Nightshift")) %>%
  separate(term, into = c("chrono", "shift"), sep = "_") %>%
  mutate(chrono = factor(chrono, levels = c("Morning", "Don't know", "Evening")),
         shift = factor(shift, levels = c("Never", "Sometimes", "Always")),
         lower = ifelse(reference, 0, estimate - 1.96 * std.error),
         upper = ifelse(reference, 0, estimate + 1.96 * std.error),
         in_mins = estimate * 60)

point_df <- data.frame(
  chrono   = factor("Morning", levels = levels(results$chrono)),
  estimate = c(0, NA, NA),
  shift    = factor(c("Never", "Sometimes", "Always"),
                    levels = levels(results$shift))
)

ps <- ggplot(
  results,
  aes(x = chrono, y = estimate, fill = shift, alpha = FDR < 0.05)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge(width = 0.8),
    width = 0.15,
    linewidth = 0.6
  ) +
  geom_point(
    data = point_df,
    aes(x = chrono, y = estimate, fill = shift, group = shift),
    position = position_dodge(width = 0.8),
    shape = 21,
    size = 3,
    color = "black",
    stroke = 0.8,
    alpha = 1,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Never"     = "#40BCD8",
      "Sometimes" = "#F39237",
      "Always"    = "#D63230"
    )
  ) +
  labs(
    x = "Chronotype",
    y = "CA (β, 95% CI)",
    fill = "Night shift work"
  ) +
  scale_alpha_manual(
    values = c(`TRUE` = 1, `FALSE` = 0.4),
    breaks = c(TRUE, FALSE),
    labels = c("Yes", "No"),
    na.translate = FALSE,
    guide = guide_legend(
      override.aes = list(
        colour = "black",   # or NA if you want no outline
        fill = "black",     # ensures solid neutral legend key
        shape = 16
      )
    )) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(linewidth = 0.1),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

ggsave(
  "plots/FX_pheno_nightshift_chrono.png",
  ps,
  width = 6,
  height = 4
)
