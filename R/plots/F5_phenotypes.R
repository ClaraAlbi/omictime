library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
install.packages("ggstatsplot")
library(ggstatsplot)

covs <- readRDS("/mnt/project/biomarkers/covs.rds")
l <- list.files("/mnt/project/biomarkers_3", full.names = T)
preds_olink <- tibble(f = l[str_detect(l, "predictions_olink")]) %>%
  mutate(d = map(f, readRDS)) %>%
  unnest(d) %>%
  mutate(gap = pred_lasso - time_day)
preds_olink$res <- residuals(lm(pred_lasso ~ time_day, data = preds_olink))



data <- preds_olink %>%
  left_join(covs) %>%
  left_join(data.table::fread("/mnt/project/clara/covariates.csv")) %>%
  mutate(Sex = factor(sex, labels = c("Female", "Male"))) %>%
  group_by(p21000_i0) %>%
  mutate(n = n()) %>% filter(n > 200) %>%
  mutate(
    ancestry = factor(
      p21000_i0,
      levels = c("British",  "Any other white background","African",
                 "Caribbean", "Indian", "Irish")
    ), ancestry = forcats::fct_recode(ancestry,"Other white" = "Any other white background")
  ) %>% ungroup()


plot_demo1 <- data %>%
  ggplot(aes(x = age_recruitment, y = res, color = Sex)) + geom_smooth() +
  theme_classic(base_size = 14) +
  labs(x = "Age", y = "Acceleration", color = "Sex") +
  ggtitle("C") +
  paletteer::scale_color_paletteer_d("nbapalettes::cavaliers_retro") +
  theme(legend.position      = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 20, face  = "bold"))


# Models
broom::tidy(glm(res ~ age_recruitment, data = data %>% filter(Sex == "Male")))

gam_mod <- mgcv::gam(
  res ~ s(age_recruitment,         # a smooth function of age
          k = 3       # you can tweak number of knots
  ),
  data = data %>% filter(Sex == "Female"),
  method = "REML"
)
summary(gam_mod)

broom::tidy(lm(res ~ ancestry, data = data))

# ANCESTRY


ggbetweenstats(data %>% slice(1:1000), gap, p21000_i0)


plot_demo2 <- data %>%
  filter(!is.na(ancestry)) %>%
  ggplot(aes(x = ancestry, y = res, fill = ancestry)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = 2) +
  # manual significance bar between x = 1 and x = 2
  annotate("segment",
           x    = 1,    # left group
           xend = 3,    # right group
           y    = y_max + 0.05,
           yend = y_max + 0.05,
           size = 0.4) +
  annotate("text",
           x    = 2,                    # midpoint
           y    = y_max + 0.3,           # a bit above the bar
           label = "p = 7e-16",
           size  = 2) +
  annotate("segment",
           x    = 1,    # left group
           xend = 4,    # right group
           y    = y_max + 0.65,
           yend = y_max + 0.65,
           size = 0.4) +
  annotate("text",
           x    = 2.5,                    # midpoint
           y    = y_max + 0.9,           # a bit above the bar
           label = "p = 2e-12",
           size  = 2) +
  annotate("segment",
           x    = 1,    # left group
           xend = 5,    # right group
           y    = y_max + 1.2,
           yend = y_max + 1.2,
           size = 0.4) +
  annotate("text",
           x    = 3,                    # midpoint
           y    = y_max + 1.5,           # a bit above the bar
           label = "p = 1e-11",
           size  = 2) +
  labs(
    x = "Self-reported ethnicity",
    y = "Acceleration"
  ) +
  ggtitle("D") +
  paletteer::scale_fill_paletteer_d("rcartocolor::Earth", direction = -1) +
  theme_classic(base_size = 14) +
  theme(
    axis.text       = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(size = 20, face = "bold")
  )


part2 <- cowplot::plot_grid(plot_demo1, plot_demo2, ncol = 2)



### Acceleration

acc <- data.table::fread("/mnt/project/accelerometer.csv")

data %>%
  left_join(acc) %>% mutate(a = sum(p40030_i0 != "")) %>% summarise(a)
