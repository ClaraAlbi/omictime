library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
install.packages("forcats")
#install.packages("ggstatsplot")
#library(ggstatsplot)

covs <- readRDS("/mnt/project/biomarkers/covs.rds")

preds_olink <- readRDS("/mnt/project/olink_int_replication.rds") %>% filter(i == 0 & !is.na(cv)) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)),
         gap = pred_mean - time_day,
         mod_sd = sd(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  left_join(covs) %>%
  mutate(Sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
         time_day_r = round(time_day, 0)) %>%
  filter(age_recruitment > 39)
preds_olink$res <- residuals(lm(pred_mean ~ time_day, data = preds_olink))
mod <- lm(pred_mean ~ time_day, data = preds_olink)

cor(preds_olink$time_day, preds_olink$pred_mean)^2

ggplot(preds_olink, aes(x = time_day, y = pred_mean, color = res)) +
  geom_pointrange(aes(ymin = pred_mean - mod_sd, ymax = pred_mean + mod_sd), position = "dodge")

ggplot(preds_olink, aes(x = time_day, y = pred_mean, color = res)) +
  geom_pointrange(aes(ymin = pred_mean - mod_sd, ymax = pred_mean + mod_sd)) +
  geom_abline(
    intercept = coef(mod)[1],
    slope     = coef(mod)[2],
    color     = "black",
    size      = 1
  ) +
  # draw residual vectors for just those two
  ggtitle("A") +
  labs(
    x     = "Recorded time of day",
    y     = "Mean predicted proteomic time", color = "Acceleration"
  ) +
  guides(
    colour = guide_colourbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 14,
      barheight      = 1.5,
      reverse = TRUE
    )
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position   = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10), legend.margin = margin(0,0,0,0),
        plot.title = element_text(size = 20, face  = "bold"))



# >%
#   group_by(p21000_i0) %>%
#   mutate(n = n()) %>% filter(n > 200) %>%
#   mutate(
#     ancestry = factor(
#       p21000_i0,
#       levels = c("British",  "Any other white background","African",
#                  "Caribbean", "Indian", "Irish")
#     ), ancestry = forcats::fct_recode(ancestry,"Other white" = "Any other white background")
#   ) %>% ungroup()

# Compare res and mod_sd

ggplot(preds_olink, aes(x = age_recruitment, y = res, color = Sex)) + geom_smooth()

ggplot(preds_olink, aes(x = age_recruitment, y = res, color = Sex)) + geom_boxplot() +
  facet_grid(~age_recruitment, scales = "free")


data %>%
  mutate(time_day_r = round(time_day, 0)) %>%
  group_by(Sex, age_recruitment) %>%
  summarise(m_p = mean(res), n = n(), sd_p = sd(res)) %>%
  filter(n > 10) %>%
  ggplot(aes(x = age_recruitment, y = m_p, color = Sex)) +
  geom_pointrange(aes(ymin = m_p - sd_p, ymax = m_p + sd_p), position = "dodge")

preds_olink %>%
  mutate(time_day_r = round(time_day, 0)) %>%
  group_by(Sex, time_day_r) %>%
  summarise(m_p = mean(res), n = n(), sd_p = sd(res)) %>%
  filter(n > 10) %>%
  ggplot(aes(x = time_day_r, y = m_p, color = Sex)) +
  geom_pointrange(aes(ymin = m_p - sd_p, ymax = m_p + sd_p), position = "dodge")


data %>%
  mutate(time_day_r = round(time_day, 0)) %>%
  group_by(Sex, time_day_r) %>%
  summarise(m_p = mean(gap), n = n(), sd_p = sd(gap)) %>%
  filter(n > 10) %>%
  ggplot(aes(x = time_day_r, y = m_p, color = Sex)) +
  geom_pointrange(aes(ymin = m_p - sd_p, ymax = m_p + sd_p), position = "dodge")

data %>%
  ggplot(aes(x = time_day, y = res, color = Sex)) + geom_point()
  geom_pointrange(aes(ymin = m_p - sd_p, ymax = m_p + sd_p), position = "dodge")


t.test(data$res ~ data$Sex)
t.test(data$res, data$age_recruitment)

plot_demo1 <- data %>%
  filter(age_recruitment > 39) %>%
  ggplot(aes(x = age_recruitment, y = time_day, color = Sex)) + geom_smooth() +
  theme_classic(base_size = 14) +
  labs(x = "Age", y = "Recorded time of day", color = "Sex") +
  ggtitle("C") +
  paletteer::scale_color_paletteer_d("nbapalettes::cavaliers_retro") +
  theme(legend.position      = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(size = 20, face  = "bold"))


# Models
broom::tidy(glm(gap ~ Sex, data = data))
broom::tidy(glm(gap ~ age_recruitment, data = data %>% filter(Sex == "Female")))

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


## well being

fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

mh <- data.table::fread("/mnt/project/psychosocial_MH.csv")
phy <- data.table::fread("physical_activity.csv")
sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid, h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`) %>%
  filter(chrono %in% 1:4) %>%
  mutate(super_sleep = case_when(h_sleep > 1 ~ 1,
                                 TRUE ~0))
wb <- data %>%
  left_join(sleep)


summary(lm(gap ~ super_sleep, data = wb))

a <- fields[fields$field_id %in% as.numeric(str_remove(str_remove(colnames(wb), "p"), "_i0")),]


wb %>%
  group_by(p4548_i0) %>%
  summarise(m_gap = mean(gap))
  ggplot(aes(x = p4548_i0, y = gap)) + geom_violin()

health <- wb  %>%
  filter(p4548_i0 != "") %>%
  mutate(score = case_when(
    p4548_i0 == "Extremely happy"   ~  3,
    p4548_i0 == "Very happy" ~ 2,
    p4548_i0 == "Moderately happy"                    ~  1,

    p4548_i0 == "Moderately unhappy"                  ~ -1,
    p4548_i0 == "Very unhappy" ~ -2,
    p4548_i0 == "Extremely unhappy" ~ -3,
    TRUE                                              ~  NA_integer_   # covers "" / “Do not know” / “Prefer not to answer”
)) #%>% group_by(score) %>%
  summarise(m_gap = mean(gap), sd_gap = sd(gap), n = n())

summary(lm(gap ~ score, data = health))

summary(lm(score ~ abs(res) + sex + factor(assessment_centre) + p22189, data = health))

wb %>%
  group_by(p20127_i0) %>%
  summarise(m_gap = mean(gap), sd_gap = sd(gap))


hist(wb$p22189)

wb %>%
  group_by(p4548_i0) %>%
  summarise(m_gap = mean(gap), sd_gap = sd(gap))

ggplot(wb, aes(x = p22189, y = gap)) + geom_smooth()

cor.test(wb$gap, wb$p22189, use = "complete.obs")



