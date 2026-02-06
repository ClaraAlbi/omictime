

p <- readRDS("data_share/LASSO_14pan_weights.rds") %>%
  distinct(feature, .keep_all = TRUE) %>%
  filter(abs(weight) > 0)

write(p$feature, "data_share/list_prots0.1.txt")
write(p$feature, "data_share/list_prots_all.txt")
