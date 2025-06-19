

# Data frame for plotting
observed <- -log10(sort(gwas$P))
expected <- -log10(ppoints(length(gwas$P)))

df <- data.frame(expected = expected, observed = observed)

# QQ plot
p <- ggplot(df, aes(x = expected, y = observed)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(x = expression(Expected~~-log[10](italic(p))),
       y = expression(Observed~~-log[10](italic(p))),
       title = "QQ Plot") +
  theme_minimal()

ggsave("qqplot_absgap.png", p)

chisq <- qchisq(1 - gwas$P, df = 1)
lambda_gc <- median(chisq) / qchisq(0.5, df = 1)
lambda_gc
