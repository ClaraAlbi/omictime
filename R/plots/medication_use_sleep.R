


results <- bind_rows(readRDS("data_share/results_associations_sleep_depression_CA.rds") %>% mutate(type = "CA"),
                     readRDS("data_share/results_associations_sleep_depression_CM.rds") %>% mutate(type = "CM"))

plot_data <- results %>%
  add_column(n = c(c(1628, 392, 410, 345, 201, 801, 904, 3229), c(1628, 392, 410, 345, 201, 801, 904, 3229)),
             o = c(c(1, 8, 2, 5, 6, 3, 4, 7), c(1, 8, 2, 5, 6, 3, 4, 7))) %>%
  mutate(lower = ifelse(reference, 0, estimate - 1.96 * std.error),
         upper = ifelse(reference, 0, estimate + 1.96 * std.error),
         term = str_remove(term, "sleep_depression"),
    display_term = ifelse(reference,
                          paste0(term, " (ref)"), term),
    display_term = paste0(display_term, "\n", n),
    display_term = fct_reorder(display_term, -o),
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value)
    )

p_dep <- plot_data %>%
  ggplot(aes(x = estimate, y = display_term, alpha = FDR < 0.05, color = type)) +

  geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.2), size = 1, fatten = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_color_manual(values = c("#2ca02c", "#9467bd")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(title = "Sleep changed during sadness/depression \nor loss of interest",
       x = "Effect size (SE)",
       y = NULL, color = " ", alpha = "FDR < 5%") +
  theme_classic(base_size = 14) +
  theme(
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(), legend.position = "right"
  )

ggsave("plots/F8_depression_sleep.png", p_dep, width = 7, height = 5)
