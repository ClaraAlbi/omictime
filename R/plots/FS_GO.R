
library(tidyverse)

pr <- data.table::fread("data/explore_ukb (1).csv")

go <- data.table::fread("data/explore_ukb (1).csv") %>%
  mutate(
    go_bp = str_trim(`Biological Process`))  %>%
  separate_rows(go_bp, sep = ",\\s*") %>%
  filter(go_bp != "")

a <- go %>% filter(Gene %in% pr$Gene) %>% group_by(go_bp) %>% summarise(l = list(Gene))


top <- go %>% group_by(go_bp) %>% count() %>% arrange(desc(n))
top %>% filter(str_detect(go_bp, "circadian")) %>% pull(go_bp) %>% paste0(., collapse = ", ")

go_i <- go %>%
  filter(str_detect(go_bp, "circadian"))
