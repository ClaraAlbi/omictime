

install.packages("patchwork")
library(patchwork)


pc <- (p2 + p3 + p5) + plot_annotation(title = "C") &
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    )
  )

p_rep <- (p_v3 + p_v1 + p_v2 ) +
  plot_annotation(title = "A") &
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    )
  )


p_tre <- p_c + theme(plot.margin = margin(0, 150, 0, 150)) +
  plot_annotation(title = "B") &
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14,
      hjust = 0
    ))



pf <- cowplot::plot_grid(p_rep, p_tre, pc, nrow = 3, rel_heights = c(1, 0.7, 1.8))

pf
ggsave("plots/Figure_4.png", pf, width = 10, height = 12, dpi = 320)
