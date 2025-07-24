library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
install.packages("ggpubr")
install.packages("ggrepel")
install.packages("paletteer")
install.packages("forcats")
install.packages("cowplot")
install.packages("broom")
install.packages("ggdraw")

covs <- readRDS("/mnt/project/biomarkers/covs.rds")
ancestry <- data.table::fread("/mnt/project/clara/covariates.csv")

l <- list.files("/mnt/project/biomarkers_3", full.names = T)

preds_olink <- tibble(f = l[str_detect(l, "predictions_olink")]) %>%
    mutate(d = map(f, readRDS)) %>%
    unnest(d)
preds_olink$res <- residuals(lm(pred_lasso ~ time_day, data = preds_olink))

mod <- lm(pred_lasso ~ time_day, data = preds_olink)

df2 <- broom::augment(mod, data = preds_olink)

# Identify the two biggest |residuals|
top2 <- df2 %>%
  filter(eid %in% c(5723240, 1218408))
  filter(time_day == 12  | time_day == 18) %>%
  mutate(absres = abs(.resid)) %>%
  slice_max(absres, n = 4, with_ties = FALSE) %>%
  slice_sample(n = 2)

p_ex <- ggplot(df2, aes(x = time_day, y = pred_lasso, color = res)) +
  geom_point() +
  geom_abline(
    intercept = coef(mod)[1],
    slope     = coef(mod)[2],
    color     = "black",
    size      = 1
  ) +
  # highlight the two points
  geom_point(
    data  = top2,
    aes(x = time_day, y = pred_lasso),
    color = "black",
    size  = 3
  ) +
  # draw residual vectors for just those two
  geom_segment(
    data  = top2,
    aes(
      x    = time_day,
      y    = pred_lasso,
      xend = time_day,
      yend = .fitted
    ),
    arrow      = arrow(length = unit(0.15, "cm")),
    color      = "black",
    size       = 1
  ) +
  ggrepel::geom_label_repel(
    data = top2,
    aes(x = time_day, y = pred_lasso,
        label = paste0("res=", round(.resid, 2), "h")),
    nudge_x    = 0.5,       # move right
    nudge_y    = 0,         # no vertical shift
    hjust      = 0,         # left-align text so it extends rightward
    segment.size = 0,
    direction  = "y",       # only repel vertically, to keep all labels at same x-offset
    color      = "black",
    size       = 4,
    box.padding= 0.1
  ) +
  ggtitle("A") +
  labs(
    x     = "Recorded time of day",
    y     = "Predicted proteomic time", color = "Acceleration"
  ) +
  paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, limits = c(-4, 4),
                                     oob    = scales::squish) +
  theme_classic(base_size = 14) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"), plot.title = element_text(size = 20, face  = "bold"))

p_hist <- preds_olink %>%
  ggplot(aes(x = res, fill = res)) +
  geom_histogram(
    aes(fill = ..x..),   # map bin midpoint to fill
    bins = 30,           # or whatever bin count you prefer
    color = "white"      # optional: white borders between bins
  ) +
  labs(x = "Acceleration") +
  ggtitle("B") +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 20, face  = "bold"))


part1 <- cowplot::plot_grid(p_ex, p_hist, ncol = 2, rel_widths = c(1.5, 1))

### COVARIATES

data <- preds_olink %>%
  left_join(covs) %>%
  left_join(ancestry) %>%
  mutate(Sex = factor(sex, labels = c("Female", "Male"))) %>%
  group_by(p21000_i0) %>%
  mutate(n = n()) %>% filter(n > 200) %>%
  mutate(
    ancestry = factor(
      p21000_i0,
      levels = c("British",  "Any other white background","African",
                 "Caribbean", "Indian", "Irish")
    ), ancestry = forcats::fct_recode(ancestry,"Other white" = "Any other white background")
  )



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

broom::tidy(glm(res ~ age_recruitment, data = data %>% filter(Sex == "Male")))

gam_mod <- mgcv::gam(
  res ~ s(age_recruitment,         # a smooth function of age
      k = 10       # you can tweak number of knots
    ),
  data = data %>% filter(Sex == "Female"),
  method = "REML"
)

broom::tidy(lm(res ~ ancestry, data = data))

data %>%
  mutate(bmi = weight / (height/100)^2) %>%
  ggplot(aes(x = bmi, y = res)) + geom_smooth()

y_max <- data %>%
  filter(!is.na(ancestry)) %>%
  pull(res) %>%
  max(na.rm = TRUE)

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


###

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  mutate(gap = pred_lasso - time_day) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  filter(n == 3)

df <- df %>% group_by(i) %>%
  nest() %>%
  mutate(res = map(data, ~residuals(lm(pred_lasso ~ time_day, data = .x)))) %>%
  unnest()

res_wide <- df %>%
  pivot_wider(id_cols = eid,
    names_from  = i,
    values_from = c(res, time_day),
    names_prefix = "i"
  )


make_pair_plot <- function(v1, v2){
  # column names
  x_res  <- paste0("res_i",      v1)
  y_res  <- paste0("res_i",      v2)
  x_time <- paste0("time_day_i", v1)
  y_time <- paste0("time_day_i", v2)

  # 1) compute the two correlations
  r_acc  <- cor(res_wide[[x_res]],  res_wide[[y_res]],  use = "pairwise.complete.obs")
  r_time <- cor(res_wide[[x_time]], res_wide[[y_time]], use = "pairwise.complete.obs")

  # 2) build a little data frame for the two text labels
  label_df <- tibble(
    x     = c(-3.5, -3.5),            # both start at x = –3.5
    y     = c( 3.8,  3.4),            # one at y=3.8, one at y=3.4
    label = c(
      paste0("italic(r)[Acceleration] == ", round(r_acc,  2)),
      paste0("italic(r)[Time~day]       == ", round(r_time, 2))
    ),
    col   = c("#2374AB", "#E85F5C")
  )

  # 3) draw!
  ggplot(res_wide, aes_string(x = x_res, y = y_res)) +
    # residuals
    geom_point(color = "#2374AB", alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "#2374AB") +

    # the two text labels
    geom_text(
      data        = label_df,
      aes(x = x, y = y, label = label, color = col),
      parse       = TRUE,
      hjust       = 0,
      size        = 3.5,
      show.legend = FALSE
    ) +
    scale_x_continuous(limits = c(-4, 4)) +
    scale_y_continuous(limits = c(-4, 4)) +
    scale_color_identity() +

    # zoom to –4…4 without dropping any data/text

    labs(
      title = paste0("i", v1, " vs i", v2),
      x     = paste0("Acceleration ","i", v1),
      y     = paste0("Acceleration ", "i", v2)
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title   = element_text(face = "bold", size = 16),
      axis.title   = element_text(face = "bold")
    )
}

# rebuild your three panels
p1 <- make_pair_plot(0, 2)
p2 <- make_pair_plot(0, 3)
p3 <- make_pair_plot(2, 3)

final_plot <-  cowplot::plot_grid(p1, p2, p3, nrow = 1)


title_grob <- cowplot::ggdraw() +
  cowplot::draw_label(
    "E",
    fontface = "bold",
    x        = 0,        # left-align
    hjust    = -1.5,        # left justification
    size     = 22
  ) +
  theme(
    # strip away any background color
    plot.background  = element_rect(color = "white"),
    panel.background = element_rect(color = "white")
  )

# 2. Stack title + plot
titled_plot <- cowplot::plot_grid(
  title_grob,
  final_plot,
  ncol        = 1,
  rel_heights = c(0.05, 1)  # 10% height for title, 90% for the grid
)


full <- cowplot::plot_grid(part1, part2, titled_plot, nrow = 3)

ggsave("plots/F4_combined.png", full, width = 10, height = 12)


df2 <- df %>%
  filter(n == 3) %>%
  pivot_wider(id_cols = c(eid), names_from = i, values_from = c(res, y), names_prefix = "i") %>%
  rowwise() %>%
  mutate(max_drawdown = max(c_across(starts_with("res"))) - min(c_across(starts_with("res")))) %>%
  ungroup() %>%
  pivot_longer(
    cols       = c(res_i0, res_i2, res_i3, y_i0, y_i2, y_i3),
    names_to   = c(".value", "i"),
    names_pattern = "(.*)_i(.*)"
  ) %>%
  filter(!is.na(y))

hist(df2$max_drawdown)
sum(df2$max_drawdown > 3)

ggplot(df2, aes(x = y, y = res, group = eid)) +
  # first draw all “flat” subjects in light grey
  geom_line(data = filter(df2, max_drawdown < 3),
            colour = "grey80", alpha = 0.5) +
  # then draw the big oscillators on top
  geom_line(data = filter(df2, max_drawdown > 3),
            aes(colour = max_drawdown), size = 0.5) +
  scale_colour_viridis_c(option = "C", name = "max_drawdown") +
  labs(x = "Year sample", y = "Residual") +
  theme_minimal(base_size = 14)



a <- df %>%
  filter(n == 3) %>%
  mutate(y = as.numeric(y)) %>%
  pivot_wider(id_cols = c(eid), names_from = i, values_from = c(res, y), names_prefix = "i") %>%
  rowwise() %>%
  mutate(mean_i03 = mean(c(res_i0, res_i2, res_i3)),
         period = max(c_across(starts_with("y"))) - min(c_across(starts_with("y")))) %>%
  ungroup() %>%
  pivot_longer(
    cols       = c(res_i0, res_i2, res_i3, y_i0, y_i2, y_i3),
    names_to   = c(".value", "i"),
    names_pattern = "(.*)_i(.*)"
  )

ggplot(a, aes(x = mean_i03)) + geom_histogram()

a %>%
  pivot_wider(id_cols = c(eid), names_from = i, values_from = res, names_prefix = "i") %>%
  ggplot(aes(x = i0)) + geom_histogram()

hist(a$mean_i03)


ggplot(df, aes(x = time_day, y = res, color = i)) + geom_point()
