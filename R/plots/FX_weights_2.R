library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
install.packages("patchwork")
library(patchwork)
library(ggplot2)
install.packages("forcats")


cv_sets <- readRDS("/mnt/project/LASSO_exp/cv.olink_14panels_subsets.rds")

data <- readRDS("/mnt/project/LASSO_exp/res_olink_14panels.rds") %>%
  mutate(across(-eid, ~ tidyr::replace_na(.x, 0)))

thresholds <- c(0, 0.01, 0.05, 0.1)

compute_cv_score <- function(k, threshold, data, cv_sets) {

  cv_coef <- readRDS(
    sprintf("/mnt/project/LASSO_exp/coefs_olink_14panels_cv%d.rds", k)
  ) %>%
    filter(model == "LASSO")

  w <- cv_coef$weight[-1]
  names(w) <- cv_coef$feature[-1]

  # apply threshold
  w[abs(w) < threshold] <- 0

  data[cv_sets == k, ] %>%
    select(eid, all_of(names(w))) %>%
    transmute(
      eid = eid,
      score = as.numeric(as.matrix(select(., -eid)) %*% w),
      cv = paste0("cv", k),
      threshold = threshold,
      n = sum(w != 0)
    ) %>%
    as_tibble()
}

res <- crossing(
  k = 1:5,
  threshold = thresholds
) %>%
  pmap_dfr(~ compute_cv_score(..1, ..2, data, cv_sets))


p1 <- res %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, time_day)) %>%
  group_by(cv, threshold) %>%
  summarise(c = cor(time_day, score)^2) %>%
  ggplot(aes(x = as.factor(threshold), y = c)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0.4, 0.8)) +
  labs(x= "LASSO weight threshold", y = "R2") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major.y = element_line(colour = "gray"))

df <- res %>%
  left_join(
    readRDS("/mnt/project/biomarkers/time.rds") %>%
      select(eid, time_day),
    by = "eid"
  ) %>%
  group_by(cv, threshold, n) %>%
  summarise(
    c = cor(time_day, score)^2,
    .groups = "drop"
  )

n_mean_df <- df %>%
  group_by(threshold) %>%
  summarise(
    n_mean = mean(n),
    y = max(n),
    .groups = "drop"
  )

p2 <- df %>%
  ggplot(aes(x = as.factor(threshold), y = n)) +
  geom_boxplot() +
  geom_text(size = 3,
    data = n_mean_df,
    aes(
      x = as.factor(threshold),
      y = y,
      label = paste0(round(n_mean, 1)),
    ),
    vjust = -0.5
  ) +
  coord_cartesian(clip = "off") +
  labs(x= "LASSO weight threshold", y = "Number of features") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major.y = element_line(colour = "gray"))

coef_avg <- map_dfr(1:5, ~
                      readRDS(
                        sprintf("/mnt/project/LASSO_exp/coefs_olink_14panels_cv%d.rds", .x)
                      ) %>%
                      filter(model == "LASSO") %>%
                      filter(feature != "(Intercept)") %>%
                      mutate(cv = .x, feature = toupper(feature))
)

saveRDS(coef_avg, "data_share/LASSO_14pan_weights.rds")

top_features <- LASSO_14pan_weights %>% group_by(feature) %>% summarise(m_w = mean(weight)) %>%
  filter(abs(m_w) > 0.1) %>% pull(feature)


write(top_features, "data_share/top_66_proteins.txt")

p3 <- coef_avg %>%
  group_by(cv) %>% filter(abs(weight) > 0.1) %>%
  ggplot(aes(x = forcats::fct_reorder(feature, desc(abs(weight))), y = abs(weight), fill = weight >0)) +
  geom_boxplot() +
  labs(x= "Protein", y = "abs(weight)") +
  theme_classic(base_size = 10) +
  scale_fill_manual(
    labels = c(`TRUE` = "> 0", `FALSE` = "< 0"),
    values = c(
      `TRUE`  = "#2C7FB8",  # muted blue (positive)
      `FALSE` = "#D7301F"   # muted red (negative)
    ),
    name = "Sign weight"
  ) +
  theme(panel.grid.major.y = element_line(colour = "gray"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = c(0.9, 0.7),
        legend.background = element_rect(
          fill = "white",
          colour = "black"
        ))



pf <- (p1 + p2) / p3 +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 2))

ggsave("plots/F3_LASSO_coef_v2.png", pf, width = 10, height = 7)


