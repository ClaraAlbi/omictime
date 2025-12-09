library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(scales)

dx download file-J4v3Z48J499kv6FXzyFQg24b

rel <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20", `22006-0.0`, `22021-0.0`, `22000-0.0`)
colnames(rel) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")

table(rel$is_white)



proteins_in_model <-


prot <- readRDS("/mnt/project/olink_instance_0_QC.rds") %>%
  filter(eid %in% rel$eid[rel$rel == 0]) %>%
  #filter(eid %in% rel$eid[rel$is_white == 1]) %>%
  select(-eid)

x <- as.matrix(prot)

x[is.na(x)] <- mean(x, na.rm = TRUE)

pca_res <- prcomp(x)

summary(pca_res)

# Extract:
scores <- pca_res$x                 # principal component scores (observations × PCs)
loadings <- pca_res$rotation        # variable loadings (original variables × PCs)
sdev <- pca_res$sdev

var_explained <- sdev^2 / sum(sdev^2)
cum_var <- cumsum(var_explained)

k <- which(cum_var >= 0.50)[1]

pc_scores_k <- pca_res$x[, 1:k]
pc_loadings_k <- pca_res$rotation[, 1:k]

# Make a small table
pca_var_tbl <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  StdDev = round(sdev, 4),
  Variance = round(var_explained, 4),
  Cumulative = round(cum_var, 4)
)
print(pca_var_tbl)

qplot(seq_along(var_explained)[1:30], var_explained[1:30], geom = "line") +
  geom_point() +
  scale_x_continuous(breaks = seq_along(var_explained)[1:30]) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Principal Component", y = "Proportion of variance explained",
       title = "Scree plot") + theme_minimal()

scores <- cbind(eid = prot$eid, as.data.frame(pc_scores_k)) %>%
  left_join(data.table::fread("ancestry_new.csv"))
ggplot(scores, aes(PC1, PC2, color = p30079)) +
  geom_point() +
  theme_minimal()

saveRDS(scores, "protein_PCA.rds")
