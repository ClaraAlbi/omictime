

olink <- readRDS("/mnt/project/biomarkers_res/res_olink_tech_14panels.rds")

olink_mat <- as.matrix(olink %>% select(-eid))

# Extract C1qtnf5 vector
c1q <- olink_mat[, "c1qtnf1"]

# Correlate C1qtnf5 with all proteins
cors <- apply(
  olink_mat,
  2,
  function(x) cor(c1q, x, use = "pairwise.complete.obs")
)

# Convert to data frame
cor_df <- data.frame(
  protein = names(cors),
  r = as.numeric(cors)
)

# Optional: remove self-correlation
cor_df <- cor_df[cor_df$protein != "c1qtnf1", ]

# View top correlations
cs <- cor_df[order(-cor_df$r), ]
