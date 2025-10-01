library(purrr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

l <- c(list.files("/mnt/project/biomarkers_3",
                  pattern = "coefs", full.names = TRUE)[-c(31:35, 1:5, 16:20)],
       list.files("/mnt/project/biomarkers_3/covariate_res/MODELS",
                  pattern = "coefs", full.names = TRUE))


weights <- tibble(f = l[str_detect(l, "coefs_all_")]) %>%
  mutate(mod = map(f, readRDS),
         cv = map_chr(f, ~sub(".*cv([0-9]+)\\.rds$", "\\1", .x))) %>%
  unnest(mod) %>%
  filter(feature != "(Intercept)") %>%
  filter(feature != 0) %>%
  group_by(model, feature) %>%
  summarise(mean_weight = mean(weight),
            sd_weight   = sd(weight), .groups = "drop",
            n = n()) %>%
  filter(n == 5) %>%
  mutate(abs_weight = abs(mean_weight)) %>%
  group_by(model) %>%
  mutate(scaled_importance = abs_weight / max(abs_weight)) %>%
  ungroup()

weights %>%
  group_by(model) %>% count()


imp <- weights %>%
  left_join(fields %>% select(field_id, title) %>% mutate(field_id = as.character(field_id)), by = c("feature" = "field_id")) %>%
  mutate(phen = case_when(is.na(title) ~ feature,
                          TRUE ~ title)) %>% select(-title, -feature) %>%
  filter(phen != "(Intercept)") %>%
  filter(n == 5) %>%
  group_by(model) %>%
  slice_max(order_by = abs(scaled_importance), n = 50) %>%
  ungroup()

importance_avg <- imp %>%
  group_by(phen) %>%
  summarise(avg_importance = mean(scaled_importance, na.rm = TRUE), .groups = "drop")

top20_heat <- imp %>%
  left_join(importance_avg, by = "phen") %>%
  mutate(phen = reorder(phen, avg_importance))

top20_heat <- top20_heat %>%
  mutate(feature_wrapped = str_wrap(phen, width = 30))

plot_import <- ggplot(top20_heat, aes(x = model, y = reorder(feature_wrapped, avg_importance), fill = scaled_importance)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "white", high = "darkgreen", midpoint = 0) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(5, 30, 5, 5)  # extra space on right
  ) +
  labs(
    x = "Model",
    y = "Feature",
    fill = "Importance"
  )

ggsave("plots/FS_weights_rank.png", plot_import, width = 8, height = 12)




####
library()
top_sets <- weights %>%
  group_by(model) %>%
  mutate(rank = rank(-scaled_importance, ties.method = "first")) %>%
  filter(rank <= 500) %>%
  summarise(features = list(feature), .groups = "drop")

# convert to list
feat_list <- setNames(top_sets$features, top_sets$model)

upset_data <- fromList(feat_list)

# drop features appearing in only 1 model
upset_data_filtered <- upset_data[rowSums(upset_data) > 1, ]

# plot upset
png("plots/FS_weights_upset.png", width = 2000, height = 1500, res = 300)
upset(upset_data_filtered,
      nsets = ncol(upset_data_filtered),
      order.by = "freq",
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size")
dev.off()

