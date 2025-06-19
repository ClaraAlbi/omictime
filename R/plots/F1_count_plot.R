
library(tidyverse)


fields <- data.table::fread("data/field.tsv")

df_r2 <- bind_rows(readRDS("data/aov_labs.rds") %>% mutate(type = "Biochemistry") %>%
                  left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                readRDS("data/aov_counts.rds") %>% mutate(type = "Cell_counts") %>%
                  left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                readRDS("data/aov_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                  left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)
  ) %>%
  mutate(color_var = case_when(type == "Proteomics-Olink" ~ "#76B041",
                               type == "Metabolomics-NMR" ~ "#2374AB",
                               type == "Cell_counts" ~ "#8F3985",
                               type == "Biochemistry" ~ "#E85F5C"))

plot_data <- df_r2 %>%
  filter(term == "time_day") %>%
  mutate(nt = case_when(
    #pr2 > 0.10 ~ ">10",
    pr2 >= 0.05 & pr2 <= 0.1 ~ "10-5%",
    pr2 >= 0.01 & pr2 <= 0.05 ~ "5-1%",
    pr2 >= 0.005 & pr2 <= 0.01 ~ "1-0.5%",
    pr2 < 0.005 ~ "<0.5%"
  ),
  nt = fct_relevel(nt, "10-5%", "5-1%", "1-0.5%", "<0.5%")) %>%
  group_by(nt, type, color_var) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type, color_var) %>%
  mutate(prop = n / sum(n)) %>%
  mutate(label_inside = ifelse(prop > 0.01, as.character(n), NA),
         mid_pos = cumsum(prop) - prop / 2)


plot_data2 <- plot_data %>%
  mutate(
    #nt = fct_rev(nt),    # ← reverse the fill‐order
    facet_html = sprintf("<span style='color:%s'>%s</span>", color_var, type),
    facet_html = factor(
      facet_html,
      levels = c(
        "<span style='color:#E85F5C'>Biochemistry</span>",
        "<span style='color:#8F3985'>Cell_counts</span>",
        "<span style='color:#2374AB'>Metabolomics-NMR</span>",
        "<span style='color:#76B041'>Proteomics-Olink</span>"
      )
    )
  )

outside_labels <- plot_data2 %>%
  filter(is.na(label_inside)) %>%
  group_by(type) %>%
  slice_max(prop, n = 1) %>%
  ungroup()

plot_bars_h <- plot_data2 %>%
  ggplot(aes(x = type, y = prop, fill = nt)) +

  # fill the panel (try >1 if you want even thicker)
  geom_bar(
    stat     = "identity",
    width    = 0.8,       # full width
    color    = "black",
    position = position_stack(reverse = TRUE)
  ) +

  # remove the tiny margin on the x‐axis so the bar truly spans edge to edge
  scale_x_discrete(expand = c(0, 0)) +

  geom_text(
    aes(label = label_inside),
    position    = position_stack(vjust = 0.5, reverse = TRUE),
    size        = 4,
    show.legend = FALSE
  ) +

  ggrepel::geom_text_repel(
    data        = outside_labels,
    aes(x      = type, y = mid_pos, label = n),
    inherit.aes = FALSE,
    direction   = "x",
    nudge_x     = 0.7,
    segment.color = "grey20",
    size        = 4,
    show.legend = FALSE
  ) +

  scale_y_reverse(expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~facet_html, scales = "free", ncol = 4) +
  labs(fill = "Time of day R²", y = "Number of biomarkers") +
  scale_fill_paletteer_d("ghibli::LaputaLight", direction = -1) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +

  theme_minimal() +
  theme(
    # **zero** space between panels
    #panel.spacing   = unit(0, "lines"),

    strip.placement  = "inside",
    strip.background = element_blank(),
    strip.text       = element_markdown(size = 14, hjust = 0),

    axis.text        = element_blank(),
    axis.title.x       = element_blank(),
    panel.grid       = element_blank(),

    legend.position  = "right",
    legend.direction = "vertical",
    legend.key.size  = unit(1, "cm"),
    legend.spacing.x = unit(1, "cm"),
    legend.text      = element_text(size = 12),
    legend.title     = element_text(size = 12)
  )

ggsave("plots/plot_vars_v.png", plot_bars_h, height = 6, width = 6)



df_r2 %>%
  filter(term == "time_day") %>%
  summarise(prop = sum(p.value < 0.05)/n(), n = n())

df_r2 %>%
  filter(term == "time_day") %>%
  group_by(type) %>%
  summarise(prop = sum(p.value < 0.05)/n(), n = n()) %>%
  filter(`p.value < 0.05` == TRUE) %>%
  group_by(type) %>% summarise(a = `p.value < 0.05`/n)

df_r2 %>%
  filter(term == "time_day") %>%
  filter(pr2 >= 0.01) %>% count()

max_vs <- df_r2 %>%
  filter(term != "Residuals") %>%
  group_by(phen, title, type) %>%
  slice_max(order_by = pr2, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(phen, title, term, type, pr2) %>%
  filter(pr2 > 0.01) %>%
  filter(term == "time_day")


