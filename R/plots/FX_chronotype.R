
install.packages("ggpubr")

all_results <- readRDS("results_chronotype_diseases_1y.rds")

p_chr <- all_results %>%
  filter(cases > 40) %>%
  left_join(list_diseases, by = c("disease" = "p")) %>%
  mutate(q = str_sub(term, -1),
         term = str_sub(term, 1, 3),
         type = case_when(type == "logistic" ~ "Prevalent",
                          type == "cox" ~ "Incident"),
         type = factor(type, levels = c("Prevalent", "Incident")),
         term = case_when(term == "chr" ~ "Self-report",
                          term == "gap" ~ "Predicted"),
         term = factor(term, levels = c("Self-report", "Predicted")),
         p_val = p.adjust(p.value, method = "fdr"),
         sig = case_when(p_val < 0.05 ~ "*",
                         TRUE ~ "")) %>%
  ggplot(aes(x = effect, y = names, color = q, shape = term, label = sig)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(xmin = lo95, xmax = hi95),
                width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_text(position = position_dodge(width = 0.5),
           color = "black", hjust = 0, size = 5,
           aes(x = hi95 + 0.03)) +
  scale_x_continuous(n.breaks = 4) +
  labs(x = "OR / HR", color = "") +
  guides(
    color = guide_legend(title = "x vs 5",
      nrow       = 1,
      byrow      = FALSE,
      override.aes = list(size = 6)
    ),
    shape = "none") +
  facet_grid(class ~  term + type, scales = "free_y", space = "free") +
  ggpubr::theme_pubclean() +
  theme(strip.background    = element_rect(fill   = "white",
                                           colour = "black",
                                           size   = 0.5),
        legend.position = "bottom",
        legend.justification = "left",
        legend.box.just    = "left",
        legend.title = element_text(),
        axis.title.y = element_blank(),
        axis.line.y = element_line(colour = "black"),
        strip.text.y = element_text(angle = 0), strip.text.y.right = element_text(angle = 270),
        panel.spacing.y = unit(1, "lines")
  )


ggsave("plots/associations_chronotype.png", p_chr, width = 8, height = 8)


