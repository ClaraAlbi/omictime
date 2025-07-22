

install.packages("ggpubr")

df <- readRDS("olink_int_replication.rds") %>%
  mutate(gap = pred_lasso - time_day) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  filter(n ==3)

df <- df %>% group_by(i) %>%
  nest() %>%
  mutate(res = map(data, ~residuals(lm(pred_lasso ~ time_day, data = .x)))) %>%
  unnest()

ggplot(df, aes(x = time_day, y = res, color = i)) + geom_point()

res_wide <- df %>%
  pivot_wider(id_cols = eid,
    names_from  = i,
    values_from = res,
    names_prefix = "i"
  )

make_pair_plot <- function(xvar, yvar, xlab, ylab) {
  ggplot(res_wide, aes_string(x = xvar, y = yvar)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    ggpubr::stat_cor(
      aes(label = ..r.label..),
      method = "pearson",
      size = 3
    ) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 12)
}

p1 <- make_pair_plot(
  "i0",
  "i2",
  "Residuals: Initial assessment",
  "Residuals: Imaging"
)

p2 <- make_pair_plot(
  "i0",
  "i3",
  "Residuals: Initial assessment",
  "Residuals: First repeat imaging"
)
p3 <- make_pair_plot(
  "i2",
  "i3",
  "Residuals: Imaging",
  "Residuals: First repeat imaging"
)

# 6. Arrange with cowplot
final_plot <- plot_grid(
  p1, p2, p3,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_size = 14
)


cor(res_wide$i0, res_wide$i2)


hist(df2$max_drawdown)

df2 <- df %>%
  filter(n == 3) %>%
  pivot_wider(id_cols = c(eid, sex), names_from = i, values_from = res, names_prefix = "i") %>%
  rowwise() %>%
  mutate(max_drawdown = max(c_across(starts_with("i"))) - min(c_across(starts_with("i")))) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("i"), names_to = "time", values_to = "residual")

sum(df2$max_drawdown > 3)

ggplot(df2, aes(x = time, y = residual, group = eid)) +
  # first draw all “flat” subjects in light grey
  geom_line(data = filter(df2, max_drawdown < 3),
            colour = "grey80", alpha = 0.5) +
  # then draw the big oscillators on top
  geom_line(data = filter(df2, max_drawdown > 3),
            aes(colour = max_drawdown), size = 0.5) +
  facet_grid(~sex) +
  scale_colour_viridis_c(option = "C", name = "max_drawdown") +
  labs(x = "Instance", y = "Residual",
       title = "Top Oscillators") +
  theme_minimal(base_size = 14)


cors <- df %>%
  group_by(i) %>%
  nest() %>%
  mutate(cor_matrix = map_dbl(data, ~ cor(.x$time_day, .x$pred_lasso, use = "pairwise.complete.obs")^2))

pp <- df %>%
  left_join(cors %>% select(i, cor_matrix)) %>%
  ggplot(aes(x = time_day, y = pred_lasso)) +
  geom_point() +
  facet_wrap(~paste0("i",i) + paste0("R2 = ",round(cor_matrix, 2))) +
  labs(x = "Time of day", y = "Predicted time") +
  theme_minimal() + theme(axis.title.x = element_blank())

pv <- df %>%
  ggplot(aes(y = as.factor(i), x = time_day)) +
  geom_violin() +
  facet_wrap(~i, scales = "free_y") +  # match facet exactly
  labs(y = "Density", x = "Time of day") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), strip.text.x = element_blank())

pf <- cowplot::plot_grid(pp, pv, nrow = 2, rel_heights = c(1, 0.6), align = "v", axis = "lr")

### times
time_i0 <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  filter(y > 2006)

time_i1 <- data.table::fread("/mnt/project/blood_sampling_instance1.tsv") %>%
  #filter(eid %in% preds_i0_nmr$eid) %>%
  mutate(max_time = pmax(`3166-1.0`,`3166-1.1`,`3166-1.2`,`3166-1.3`,`3166-1.4`, `3166-1.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date) %>%
  filter(y > 2006)

i2_meta <- data.table::fread("/mnt/project/blood_sampling_instance2.tsv") %>%
  filter(!is.na(`3166-2.0`)) %>%
  #filter(eid %in% i2_data$eid) %>%
  mutate(max_time = pmax(`3166-2.0`,`3166-2.1`,`3166-2.2`,`3166-2.3`,`3166-2.4`, `3166-2.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date) %>%
  filter(y > 2006)





### GAP differences

df %>%
  mutate(gap = time_day - pred_lasso) %>%
  group_by(eid) %>%
  summarise(m_gap = mean(abs(gap), na.rm = T),
            m_time = mean(time_day, na.rm = T),
            v_gap = var(abs(gap), na.rm = T),
            v_time = var(time_day, na.rm = T),
            n = n()) %>%
  arrange(desc(v_time)) %>%
  ggplot(aes(x = v_gap, y = v_time)) +
  geom_point()

d <- preds_i0_olink %>%
  mutate(gap = time_day - pred_lasso)

m <- lm(time_day ~ pred_lasso, data = d)
d$res <- residuals(lm(time_day ~ pred_lasso, data = d))

d %>%
  mutate(gap = time_day - pred_lasso) %>%
  ggplot(aes(x = time_day, y = res2)) +
  geom_point() +
  geom_smooth(method = "lm")


## AGExSEX plots
covs <- data.table::fread("/mnt/project/covariates.tsv") %>%
  select(eid, `31-0.0`, `34-0.0`, `21022-0.0`)
colnames(covs) <- c("eid", "Sex", "year_birth", "Age_baseline")

#data_i0 <-
preds_i0 %>%
  left_join(covs) %>%
  mutate(Age = Age_baseline, Sex = as.factor(Sex),
         gap = time_day - pred_lasso,
         absgap = abs(gap)) %>%
  pivot_longer(c(gap, absgap)) %>%
  #group_by(Age, Sex) %>% summarise(m_gap = mean(gap), n = n())
  ggplot(aes(x = Age, y = value, color = Sex)) +
  geom_smooth() +
  facet_wrap(~name, scales = "free")

covs_i2 <- preds_i2 %>%
  left_join(i2_meta) %>%
  left_join(covs) %>%
  mutate(Age = y - year_birth, Sex = as.factor(Sex),
         gap = y_test - pred_lasso,
         absgap = abs(gap))
covs_i2 %>%
  pivot_longer(c(gap, absgap)) %>%
  #group_by(Age, Sex) %>% summarise(m_gap = mean(gap), n = n())
  ggplot(aes(x = Age, y = value, color = Sex)) +
  geom_smooth() +
  facet_wrap(~name, scales = "free")

covs_i3 <- preds_i3 %>%
  left_join(i3_meta) %>%
  left_join(covs) %>%
  mutate(Age = y - year_birth, Sex = as.factor(Sex),
         gap = y_test - pred_lasso,
         absgap = abs(gap))

covs_i3 %>%
  pivot_longer(c(gap, absgap)) %>%
  #group_by(Age, Sex) %>% summarise(m_gap = mean(gap), n = n())
  ggplot(aes(x = Age, y = value, color = Sex)) +
  geom_smooth() +
  facet_wrap(~name, scales = "free")


i2i3 <- inner_join(covs_i3, covs_i2, by = c("eid", "Sex"))

i2i3 %>% ggplot(aes(x = pred_lasso.x, y = pred_lasso.y, color = Sex)) +
  geom_point()

i2i3 %>% ggplot(aes(x = y_test.x, y = y_test.y, color = Sex)) +
  geom_point()

