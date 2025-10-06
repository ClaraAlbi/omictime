library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
#install.packages("ggpubr")
install.packages("paletteer")
#install.packages("forcats")
install.packages("cowplot")
install.packages("broom")
#install.packages("ggdraw")
install.packages("ggrepel")
library(ggplot2)

preds_olink <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(i == 0 & !is.na(cv)) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)),
         gap = pred_mean - time_day,
         mod_sd = sd(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  filter(!is.na(time_day))

preds_olink$res <- residuals(lm(pred_mean ~ time_day, data = preds_olink))
mod <- lm(pred_mean ~ time_day, data = preds_olink)

df2 <- broom::augment(mod, data = preds_olink)
df2$gap <- df2$pred_lasso - df2$time_day

# Identify the two biggest |residuals|
top2 <- df2 %>%
  filter(eid %in% c(5723240, 1218408))


p_ex <- ggplot(df2, aes(x = time_day, y = pred_mean, color = res)) +
  geom_pointrange(
    aes(ymin = pred_mean - mod_sd, ymax = pred_mean + mod_sd), size = 0.4) +
  geom_abline(
    intercept = coef(mod)[1],
    slope     = coef(mod)[2],
    color     = "black",
    size      = 1
  ) +
  # highlight the two points
  geom_point(
    data  = top2,
    aes(x = time_day, y = pred_mean),
    color = "black",
    size  = 3
  ) +
  # draw residual vectors for just those two
  geom_segment(
    data  = top2,
    aes(
      x    = time_day,
      y    = pred_mean,
      xend = time_day,
      yend = .fitted
    ),
    arrow      = arrow(length = unit(0.15, "cm")),
    color      = "black",
    size       = 1
  ) +
  ggrepel::geom_label_repel(
    data = top2,
    aes(x = time_day, y = pred_mean,
        label = paste0(round(res, 1), "h")),
    nudge_x    = 0.5,       # move right
    nudge_y    = 0,         # no vertical shift
    hjust      = 0,         # left-align text so it extends rightward
    segment.size = 0,
    direction  = "y",       # only repel vertically, to keep all labels at same x-offset
    color      = "black",
    size       = 5,
    box.padding= 0.1
  ) +
  labs(
    x     = "Recorded time of day",
    y     = "Mean predicted proteomic time", color = "Circadian acceleration"
  ) +
  paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging",
                                     direction = -1,
                                     limits = c(-6, 6)) +
  guides(
    colour = guide_colourbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 14,
      barheight      = 1.2,
      reverse = FALSE
    )
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position   = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12),
        legend.margin = margin(0,0,0,0),
        plot.title = element_text(size = 20, face  = "bold"))


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
  scale_x_continuous(limits = c(-5,5)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 20, face  = "bold"))


#part1 <- cowplot::plot_grid(p_ex, p_hist, ncol = 2, rel_widths = c(1.2, 1))



p_ex %>%
  { . +
      ggtitle("") +
      theme_classic(base_size = 20) +
      theme(
        legend.position   = "top",
        axis.title        = element_text(face = "bold"),
        legend.title      = element_text(face = "bold", size = 20),
        legend.text       = element_text(size = 14),
        legend.margin     = margin(0,0,0,0),
        plot.title        = element_text(size = 20, face = "bold")
      ) } %>%
  { ggsave("slides/acceleration.png", plot = ., width = 10, height = 8) }

p_hist2 <- preds_olink %>%
  ggplot(aes(x = res, fill = res)) +
  geom_histogram(
    aes(fill = ..x..),   # map bin midpoint to fill
    bins = 30,           # or whatever bin count you prefer
    color = "white"      # optional: white borders between bins
  ) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1) +
  scale_x_continuous(limits = c(-5,5), sec.axis = dup_axis(name = "Circadian acceleration", labels = NULL)) +
  scale_y_continuous(
    labels = function(x) paste0(x / 1000, "k"),
    breaks = seq(1000, 5000, 1000), # optional: control tick positions
    expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.title.x.top = element_text(size = 16, face  = "bold"),
        axis.ticks.x.top = element_blank(),
        plot.title = element_text(size = 20, face  = "bold"))

ggsave("slides/hist_acceleration.png", plot = p_hist2, width = 6, height = 5)

p_hist3 <- preds_olink %>%
  ggplot(aes(x = abs(res), fill = abs(res))) +
  geom_histogram(
    aes(fill = ..x..),   # map bin midpoint to fill
    bins = 30,           # or whatever bin count you prefer
    color = "white"      # optional: white borders between bins
  ) +
  paletteer::scale_fill_paletteer_c("ggthemes::Classic Gray", direction = 1) +
  labs(x = "Hours") +
  scale_x_continuous(limits = c(0,5), sec.axis = dup_axis(name = "Circadian misalignment", labels = NULL, expand = c(0,0))) +
  scale_y_continuous(
    labels = function(x) paste0(x / 1000, "k"),
    breaks = seq(1000, 5000, 1000), # optional: control tick positions
    expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x.top = element_text(size = 16, face  = "bold"),
        axis.ticks.x.top = element_blank(),
        plot.title = element_text(size = 20, face  = "bold"))

ggsave("slides/hist_misalignment.png", plot = p_hist3, width = 6, height = 5)

p_hist_comb <- cowplot::plot_grid(p_hist2, p_hist3, nrow = 2, rel_heights = c(0.7, 0.8))

part1 <- cowplot::plot_grid(p_ex, p_hist_comb, ncol = 2, rel_widths = c(1.7, 1), labels = c("A", "B"), label_size = 18)



###

time <- readRDS("/mnt/project/biomarkers/time.rds")

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  left_join(time) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  filter(!is.na(time_day)) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) #%>%
  filter(n == 3)

df$res <- residuals(lm(pred_mean ~ time_day, data = df))

res_wide <- df %>%
  pivot_wider(id_cols = eid,
    names_from  = i,
    values_from = c(res, time_day, y),
    names_prefix = "i"
  )


make_pair_plot <- function(v1, v2){
  # column names
  x_gap  <- paste0("res_i",      v1)
  y_gap  <- paste0("res_i",      v2)
  x_time <- paste0("time_day_i", v1)
  y_time <- paste0("time_day_i", v2)

  # 1) compute correlations
  r_acc  <- cor.test(res_wide[[x_gap]],  res_wide[[y_gap]],  use = "pairwise.complete.obs")
  r_time <- cor.test(res_wide[[x_time]], res_wide[[y_time]], use = "pairwise.complete.obs")
  n <- res_wide %>% filter(!is.na(res_wide[[y_gap]]) & !is.na(res_wide[[x_gap]])) %>% nrow()
  av_y <- mean(as.numeric(res_wide[[paste0("y_i", v2)]]) - as.numeric(res_wide[[paste0("y_i", v1)]]), na.rm = T)
  sd_y <- sd(as.numeric(res_wide[[paste0("y_i", v2)]]) - as.numeric(res_wide[[paste0("y_i", v1)]]), na.rm = T)

  # 2) build a little data frame for the two text labels
  label_df <- tibble(
    x     = c(-5, -5, -5),
    y     = c( 6,  5.3, 4.6),
    label = c(paste0("N == ", n),
      paste0("italic(r)[Acceleration] == ", round(r_acc$estimate,  2), "~(p ==", sprintf("%.0e", r_acc$p.value), ")"),
      paste0("italic(r)[Time~day] == ", round(r_time$estimate, 2), "~(p ==", sprintf("%.0e", r_time$p.value), ")")
    ),
    col   = c("black","#2374AB", "darkgreen" )
  )

  # 3) draw!
  ggplot(res_wide, aes_string(x = x_gap, y = y_gap)) +
    # residuals
    geom_point(color = "#2374AB", alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "#2374AB") +

    # the two text labels
    geom_text(
      data        = label_df,
      aes(x = x, y = y, label = label, color = col),
      parse       = TRUE,
      hjust       = 0,
      size        = 4,
      show.legend = FALSE
    ) +
    scale_x_continuous(limits = c(-5, 6)) +
    scale_y_continuous(limits = c(-5, 6)) +
    scale_color_identity() +

    # zoom to –4…4 without dropping any data/text

    labs(
      title = paste0("i", v1, " vs i", v2),
      subtitle = paste0(round(av_y, 1), " (±", round(sd_y, 1),")", " years"),
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

final_plot <-  cowplot::plot_grid(p3, p1, p2, nrow = 1)


full <- cowplot::plot_grid(part1,  final_plot, nrow = 2, labels = c("", "C"), label_size = 18)

ggsave("plots/F4_combined.png", full, width = 10, height = 8)

