library(glmnet)

m <- readRDS("data_share/cv.olink_14_lasso_cv1.rds")

cm <- coef(m)

nz <- as.vector(cm) != 0

d <- data.frame(
  variable = rownames(cm)[nz],
  weight   = as.numeric(cm[nz])
)

d %>%
  arrange(desc(abs(weight))) %>%
  slice(1:20)
