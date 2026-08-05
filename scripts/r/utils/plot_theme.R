# ============================================================
# scripts/r/utils/plot_theme.R
# Shared plotting style + save helper.
# ============================================================

library(ggplot2)

group_palette <- c("control" = "#0072B2", "IBD" = "#D55E00")

cohort_palette <- c(
  "HMP_2019_ibdmdb" = "#009E73",
  "HallAB_2017"     = "#CC79A7",
  "IjazUZ_2017"     = "#E69F00",
  "NielsenHB_2014"  = "#56B4E9"
)

theme_pipeline <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(color = "grey30", size = base_size - 1),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

save_figure <- function(plot, filename_stem, width = 6, height = 5, dpi = 300) {
  ggsave(here("results", "figures", paste0(filename_stem, ".pdf")), plot, width = width, height = height)
  ggsave(here("results", "figures", paste0(filename_stem, ".png")), plot, width = width, height = height, dpi = dpi)
}
