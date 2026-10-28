# =============================================================================
# analyses/fig_style.R  -  one figure style for the analyses/ scripts
#
# Figures are drawn AT PRINT SIZE (mm) so the point sizes below are the sizes on
# the page. The earlier figures were drawn 9-10 inches wide with 7-8 pt text and
# shrunk to a 174 mm column, which put text near 5 pt.
#
#   FIG_W2 = 174 mm (double column), FIG_W1 = 85 mm (single column)
#   axis tick labels 8 pt; axis titles, legend, strips 9 pt; panel titles 10 pt
#   geom_text: pass pt(8) as size, so labels inside panels are also 8 pt
#
# Usage:  source("analyses/fig_style.R"); ... + theme_fig(); save_fig(p, path, FIG_W2, 150)
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

FIG_W2 <- 174
FIG_W1 <- 85
pt <- function(x) x / ggplot2::.pt          # text size in points, for geom_text

theme_fig <- function(base = 9, grid = FALSE) {
  t <- theme_classic(base_size = base) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      axis.text         = element_text(size = base - 1, colour = "black"),
      axis.title        = element_text(size = base),
      axis.line         = element_line(linewidth = 0.35),
      axis.ticks        = element_line(linewidth = 0.35),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base, face = "bold", hjust = 0,
                                       margin = margin(b = 3)),
      legend.text       = element_text(size = base),
      legend.title      = element_text(size = base),
      legend.key.size   = unit(10, "pt"),
      plot.title        = element_text(size = base + 1, face = "bold", hjust = 0,
                                       margin = margin(b = 4)),
      plot.subtitle     = element_text(size = base, hjust = 0),
      plot.caption      = element_text(size = base - 1, hjust = 0, colour = "grey30"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      panel.spacing     = unit(8, "pt"),
      plot.margin       = margin(5, 8, 5, 5))
  if (grid) t <- t + theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3))
  t
}

# PDF through Quartz on macOS: embeds Helvetica and handles Unicode (beta, >=).
# cairo_pdf needs XQuartz, which is not installed here; base pdf() drops Unicode.
pdf_dev <- function(filename, width, height, ...) {
  if (capabilities("aqua")) {
    grDevices::quartz(file = filename, type = "pdf", width = width, height = height,
                      family = "Helvetica", bg = "white")
  } else {
    grDevices::pdf(filename, width = width, height = height)
  }
}

save_fig <- function(p, path_noext, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(paste0(path_noext, ".pdf"), p, width = w, height = h, device = pdf_dev)
  ggsave(paste0(path_noext, ".png"), p, width = w, height = h, dpi = 300,
         device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png", bg = "white")
  invisible(path_noext)
}
