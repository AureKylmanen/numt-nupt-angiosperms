#!/usr/bin/env Rscript
# 03_plot_pgls_panels.R -- manuscript scatter panels
#
#   Rscript 03_plot_pgls_panels.R
#
# Reads pgls_fits.rds from 02_pgls_models.R and writes 9 cm PDF panels to

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
source(file.path(here, "common.R"))
suppressPackageStartupMessages(library(ggplot2))

COL_NUMT <- "#2c728e"   # <- set to your NUMT colour
COL_NUPT <- "#21918c"   # <- set to your NUPT colour
COL_BOTH <- "#440154"
FIG <- file.path(OUT, "Figures"); dir.create(FIG, showWarnings = FALSE)

fits <- readRDS(file.path(OUT, "pgls_fits.rds"))

fmtP <- function(p) if (p < 0.001) "italic(P) < 0.001" else sprintf("italic(P) == %.3f", p)

panel <- function(key, file, xlab, ylab, col, isometry = FALSE) {
  f <- fits[[key]]
  if (is.null(f)) { cat("skip", key, "(no fit)\n"); return(invisible()) }
  dat <- f$data; xv <- names(dat)[3]; yv <- names(dat)[2]
  xs <- seq(min(dat[[xv]]), max(dat[[xv]]), length.out = 150)
  X  <- cbind(1, xs)
  yh <- as.numeric(X %*% coef(f$fit))
  se <- sqrt(pmax(0, rowSums((X %*% vcov(f$fit)) * X)))
  band <- data.frame(x = xs, y = yh, lo = yh - 1.96 * se, hi = yh + 1.96 * se)
  tb <- f$table[1, ]
  lab <- sprintf("list(italic(b) == %.3f, %s, lambda == %.2f, italic(n) == %d)",
                 tb$slope, fmtP(tb$P), tb$lambda, tb$n)

  g <- ggplot(dat, aes(.data[[xv]], .data[[yv]])) +
    geom_point(colour = col, alpha = 0.45, size = 1, stroke = 0) +
    geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE,
                fill = col, alpha = 0.2) +
    geom_line(data = band, aes(x, y), inherit.aes = FALSE, colour = col, linewidth = 0.8)
  if (isometry)   # slope-1 reference through the bivariate mean
    g <- g + geom_abline(slope = 1, intercept = mean(dat[[yv]]) - mean(dat[[xv]]),
                         linetype = "dashed", colour = "grey40", linewidth = 0.5)
  g <- g +
    annotate("text", x = -Inf, y = Inf, label = lab, parse = TRUE,
             hjust = -0.05, vjust = 1.4, size = 2.6) +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = 9) +
    theme(axis.text = element_text(colour = "black"))
  ggsave(file.path(FIG, file), g, width = 9, height = 8, units = "cm")
  cat("wrote", file, "\n")
}

panel("numt_prop_genome", "numt_merged_prop_vs_genome_size.pdf",
      expression(log[10]~"nuclear genome size (Mbp)"),
      expression(log[10]~"NUMT proportion"), COL_NUMT)
panel("nupt_prop_genome", "nupt_merged_prop_vs_genome_size.pdf",
      expression(log[10]~"nuclear genome size (Mbp)"),
      expression(log[10]~"NUPT proportion"), COL_NUPT)
panel("numt_bp_genome", "numt_bp_vs_genome_size.pdf",
      expression(log[10]~"nuclear genome size (Mbp)"),
      expression(log[10]~"merged NUMT length (bp)"), COL_NUMT, isometry = TRUE)
panel("nupt_bp_genome", "nupt_bp_vs_genome_size.pdf",
      expression(log[10]~"nuclear genome size (Mbp)"),
      expression(log[10]~"merged NUPT length (bp)"), COL_NUPT, isometry = TRUE)
panel("numt_prop_mito", "numt_merged_prop_vs_mito_size.pdf",
      expression(log[10]~"mitochondrial genome size (kbp)"),
      expression(log[10]~"NUMT proportion"), COL_NUMT)
panel("nupt_prop_plastid", "nupt_merged_prop_vs_plastid_size.pdf",
      expression(log[10]~"plastid genome size (kbp)"),
      expression(log[10]~"NUPT proportion"), COL_NUPT)
panel("nupt_prop_numt_prop", "numt_prop_vs_nupt_prop.pdf",
      expression(log[10]~"NUMT proportion"),
      expression(log[10]~"NUPT proportion"), COL_BOTH)
panel("numt_prop_mito_noSconica", "numt_merged_prop_vs_mito_size_noSconica.pdf",
      expression(log[10]~"mitochondrial genome size (kbp)"),
      expression(log[10]~"NUMT proportion"), COL_NUMT)
panel("nupt_bp_numt_bp", "nupt_bp_vs_numt_bp.pdf",
      expression(log[10]~"merged NUMT length (bp)"),
      expression(log[10]~"merged NUPT length (bp)"), COL_BOTH)
