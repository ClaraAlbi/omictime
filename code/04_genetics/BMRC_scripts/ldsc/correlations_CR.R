source("config/paths.R")

options(bitmapType='cairo-png')
library(tidyverse)
library(patchwork)
library(ggplot2)

resdir <- bmrc_path("results/v3/")
files <- list.files(resdir, pattern = "\\.log$", full.names = TRUE)

results_CR <- map_dfr(files, extract_ldsc) %>%
  mutate(fdr = p.adjust(p),
         sig2 = ifelse(fdr < 0.05, "•", ifelse(p < 0.05, "·", ""))) 

results_CR %>% filter(p1_name ==  "res_GWAS" | p2_name == "res_GWAS") %>% arrange(desc(abs(rg)))
results_CR %>% filter(p1_name ==  "chrono_v2_GWAS" | p2_name == "chrono_v2_GWAS") %>% arrange(desc(abs(rg)))
results_CR %>% filter(p1_name ==  "RAgroup.linear.munge.sumstats.gz" | p2_name == "RAgroup.linear.munge.sumstats.gz") %>% arrange(desc(abs(rg)))

x_levels <- tribble(~trait,    ~name, ~type,
  "Circadian Acceleration (CA)","res_GWAS",                "Blood",
#  "Chronotype (binary)", "Jones_2018_morning_person", "Self-report",
  "Chronotype", "chrono_v2_GWAS",      "Self-report",
#  "Low relative amplitude", "RAgroup.logistic.munge.sumstats.gz", "Accelerometer",
  "Relative amplitude", "RAgroup.linear.munge.sumstats.gz", "Accelerometer",
  "L5 time","ACC_L5_TIME", "Accelerometer",
  "M10 time","ACC_M10_TIME", "Accelerometer",
  "Sleep midpoint", "ACC_SLEEP_MIDP", "Accelerometer",
  "Diurnal inactivity","ACC_DIURNAL_INACT", "Accelerometer",
  "Sleep duration","ACC_SLEEP_DUR", "Accelerometer",
  "Sleep duration sd","ACC_SLEEP_DUR_SD", "Accelerometer",
  "Sleep efficiency","ACC_SLEEP_EFF", "Accelerometer",
  "Num sleep episodes","ACC_N_SLEEP_EPISODES", "Accelerometer"
)

trait_order  <- x_levels$name
trait_labels <- x_levels$trait

df2 <- results_CR %>%
  filter(p1_name != "RAgroup_full_sample.assoc.logistic_tidy.ma.full.sumstats.gz")  %>%
  filter(p2_name != "RAgroup_full_sample.assoc.logistic_tidy.ma.full.sumstats.gz" ) %>%
  filter(p2_name != "morning_person-5.1" ) %>%
  filter(p2_name != "chrono_b" ) %>%
  mutate(
    i1_raw = match(p1_name, trait_order),
    i2_raw = match(p2_name, trait_order)
  ) %>%
  filter(!is.na(i1_raw), !is.na(i2_raw)) %>%
  mutate(
    i1 = pmin(i1_raw, i2_raw),
    i2 = pmax(i1_raw, i2_raw),
    t1 = trait_order[i1],
    t2 = trait_order[i2]
  ) %>%
  ungroup() %>%              # ← IMPORTANT
  arrange(i1, i2) %>%
  mutate(y = factor(t1,
                    levels = rev(trait_order),
                    labels = rev(trait_labels)),
         x = factor(t2,
                    levels = trait_order,
                    labels = trait_labels))

diag_df <- tibble(
  p1_name = trait_order,
  p2_name = trait_order,
  rg      = NA_real_,
  sig2    = ""
) %>%
  mutate(
    y = factor(p2_name, levels = rev(trait_order), labels = rev(trait_labels)),
    x = factor(p1_name, levels = trait_order, labels = trait_labels)
  )


out2 <- plot_df %>%
  filter(!is.na(t1)) %>%
  select(p1 = x, p2 = y, rg, se, z, p, fdr, h2_p1, h2se_p1, h2_p2, h2se_p2)

plot_df <- bind_rows(df2, diag_df) %>%
  #filter(p < 0.05) %>%
  left_join(x_levels, by = c("p2_name" = "name")) %>%
  mutate(type = case_when(t1 == "res_GWAS"  ~ "Blood",
                          t2 == "RAgroup.linear.munge.sumstats.gz" & t1 == "Jones_2018_morning_person" ~ "Self-report",
                          t2 == "RAgroup.linear.munge.sumstats.gz" & t1 == "chrono_v2_GWAS" ~ "Self-report",
                          TRUE ~ type),
         type = factor(type, levels = c("Blood", "Self-report", "Accelerometer")))
  
saveRDS(plot_df, "OUTPUTS/table_correlations.rds")

p_all <- ggplot(plot_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = rg), color = "black", linewidth = 0.3) +
  geom_text(aes(label = sig2), size = 4) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-1, 1),
    name = "rg"
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
  theme_classic(base_size = 10) +
  theme(
    #panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    axis.title  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
    
    # place legend inside bottom-left of the PANEL
    legend.position = c(0.02, 0.02),   # 0,0 = bottom-left; 1,1 = top-right
    legend.justification  = c("left", "bottom"),
    legend.background     = element_rect(fill = alpha("white", 0.8), color = "black", linewidth = 0.3),
    legend.key.size       = unit(4, "mm"),
    legend.title          = element_text(size = 9),
    legend.text           = element_text(size = 8),
    plot.margin            = margin(t = 20, r = 20, 1, 1),
    strip.text.y = element_text(size = 6)
  )

p_all
ggsave(paste0("OUTPUTS/rg_chrono_lp.png"), p_all, width = 7, height = 7)


t <- tibble::tribble(
  ~trait,                       ~name,                        ~type,            ~n,
  "Diurnal inactivity",         "ACC_DIURNAL_INACT",          "Accelerometer",     2,
  "L5 time",                    "ACC_L5_TIME",                "Accelerometer",     5,
  "M10 time",                   "ACC_M10_TIME",               "Accelerometer",     1,
  "Num sleep episodes",         "ACC_N_SLEEP_EPISODES",       "Accelerometer",    21,
  "Sleep duration",             "ACC_SLEEP_DUR",              "Accelerometer",    10,
  "Sleep efficiency",           "ACC_SLEEP_EFF",              "Accelerometer",     4,
  "Sleep midpoint",              "ACC_SLEEP_MIDP",                   "Accelerometer",    0,
  "Sleep duration sd",           "ACC_SLEEP_DUR_SD",                 "Accelerometer",    0,
  "Relative amplitude",          "RAgroup.linear.munge.sumstats.gz", "Accelerometer",    1,
  #"Chronotype (binary)",        "Jones_2018_morning_person",  "Self-report",      65,
  "Chronotype",  "chrono_v2_GWAS",             "Self-report",     136,
  "Circadian Acceleration (CA)",  "res_GWAS",                        "Blood", 20
)


p_h2 <- bind_rows(map_dfr(files, extract_ldsc),
  map_dfr(files, extract_ldsc) %>% select(-p1_name, -h2_p1, -h2se_p1) %>% rename(p1_name = p2_name, h2_p1 = h2_p2, h2se_p1 =  h2se_p2)) %>%
  filter(p1_name %in% x_levels$name) %>%
  #map_dfr(files, extract_ldsc) %>% select(-p1_name, -h2_p1, -h2se_p1) %>% rename(p1_name = p2_name, h2_p1 = h2_p2, h2se_p1 =  h2se_p2) %>%
  distinct(p1_name, .keep_all = TRUE) %>%
  filter(p1_name != "chrono_b") %>%
  left_join(x_levels, by = c("p1_name" = "name")) %>%
  left_join(t) %>%
  mutate(type = factor(type, levels = c("Blood", "Self-report", "Accelerometer")),
         trait = factor(name, levels = trait_order, labels = trait_labels)) %>%
  ggplot(aes(x = h2_p1, y = fct_rev(trait))) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = h2_p1 - 2*h2se_p1,
                    xmax = h2_p1 + 2*h2se_p1), width = 0.2) +
  facet_grid(rows = vars(type), scales = "free_y", space = "free_y", drop   = FALSE) +
  geom_text(
    aes(x = 0.01, label = n),
    hjust = 0,
    color = "white"
  ) +
  #scale_y_reverse(expand = c(0, 0)) +
  labs(x = "SNP-h2") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(size = 7),
        axis.title.x = element_text(size = 8),
    strip.text.y = element_blank(),   # remove facet labels
    strip.background = element_blank(),
    #axis.line.y = element_blank(),
    axis.title.y = element_blank(),   # remove y-axis title
    axis.text.y  = element_blank(),   # remove y-axis text
    axis.ticks.y = element_blank(),
    plot.margin = margin(1, 15, 1, 1)
  )
p_h2

p_combined <- plot_grid(
  p_all,
  p_h2,
  ncol = 2,
  rel_widths = c(4, 1),
  align = "h",
  axis = "tb"
)

p_combined

ggsave(paste0("OUTPUTS/rg_chrono.png"), p_combined, width = 8, height = 6)



### with manhattan

p_gen <- plot_grid(
  p,
  p_combined,
  nrow = 2, rel_heights = c(0.5, 1), labels = "AUTO"
)

ggsave(paste0("OUTPUTS/Figure_5.png"), p_gen, width = png_w-1, height = png_h + 7, dpi = dpi)





p_both <- p_combined + pcomp +
  plot_layout(widths = c(2, 1))
ggsave(paste0("OUTPUTS/rg_both.png"), p_all, width = 10, height = 7)


p_both <- p_all + p_h2 + pcomp +
  plot_layout(widths = c(2.5, 0.8, 2))+
  plot_annotation(
    tag_levels = "A",
  ) &
  theme(
    plot.tag = element_text(size = 14)
  )
p_both[[2]] <- p_both[[2]] + theme(plot.tag = element_blank())

ggsave(paste0("OUTPUTS/rg_both.png"), p_both, width = 10, height = 5)
