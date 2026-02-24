

install.packages("patchwork")
library(patchwork)
pf <- (free(p_v3 + p_v1 + p_v2) / free(p_c)) + (p2 + p3 / p5) +
  plot_layout(heights = c(0.5, 0.4, 1))




ggsave("plots/Figure_3.png", pf, width = 10, height = 15, dpi = 320)
