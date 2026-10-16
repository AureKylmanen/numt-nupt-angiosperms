#!/usr/bin/env Rscript
# plot_ks_distribution.R -- Ks (dS) histogram of coding insertions, one class
#
# Usage:
#   Rscript plot_ks_distribution.R numt
#   Rscript plot_ks_distribution.R nupt
#   Rscript plot_ks_distribution.R nupt --outdir Figures --xmax 0.6
#   Rscript plot_ks_distribution.R nupt --combined all_species_nupt_kaks_fixed_combined.tsv
#
# Reads the combined Ka/Ks table of the class (default
#   <base>/all_species_<class>_kaks_fixed_combined.tsv, built by
#   rebuild_combined_kaks.py from the Ka/Ks file chosen per species in
#   final_species_summary_v2.tsv).
# Writes  <outdir>/ks_distribution_<class>.pdf and .png
#
# Figure size: drawn at 6 x 4.5 in with large text, because the manuscript
# places each panel at half the text width (subfigure 0.49\textwidth), which
# shrinks the plot about 1.8-fold; text then prints at roughly 9-10 pt.
# Axis title uses K_S (as in the manuscript), not d_S.
#
# Values above xmax are counted and reported but not drawn; the dashed line
# marks the dcNUxT threshold (Ks = 0.1). Bins are right-closed
# (0, 0.01], (0.01, 0.02], ... with Ks = 0 in the first bin; the axis is zoomed
# with coord_cartesian so no data are removed.

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !args[1] %in% c("numt", "nupt"))
  stop("first argument must be 'numt' or 'nupt'")
cls <- args[1]

opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
base     <- opt("--base",     "/vol/storage/5T_disk/organelles")
combined <- opt("--combined", file.path(base, paste0("all_species_", cls, "_kaks_fixed_combined.tsv")))
outdir   <- opt("--outdir",   base)
xmax     <- as.numeric(opt("--xmax", "0.6"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fill <- c(numt = "#2c728e", nupt = "#21918c")[[cls]]

if (!file.exists(combined)) stop("combined table not found: ", combined)
df <- read.delim(combined, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Ks" %in% names(df)) stop("no 'Ks' column in ", combined)
ks_values <- suppressWarnings(as.numeric(df$Ks))
ks_values <- ks_values[!is.na(ks_values)]
if (!length(ks_values)) stop("no Ks values in ", combined)

n_above <- sum(ks_values > xmax)
cat("Table:", combined, "\n")
cat("Species:", length(unique(df[[1]])), "\n")
cat("Total Ks values:", length(ks_values), "\n")
cat("Ks > 0.1 (dcNUxT):", sum(ks_values > 0.1), "\n")
cat(sprintf("Values > %.1f (not shown): %d\n", xmax, n_above))

p <- ggplot(data.frame(ks = ks_values[ks_values <= xmax]), aes(x = ks)) +
  geom_histogram(binwidth = 0.01, boundary = 0, closed = "right",
                 fill = fill, color = NA) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, xmax, by = 0.1), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, xmax)) +
  labs(x = expression(italic(K)[S]), y = "Count") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.ticks.length = unit(4, "pt"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 17),
    plot.margin = margin(10, 15, 10, 10)
  )

out <- file.path(outdir, paste0("ks_distribution_", cls))
ggsave(paste0(out, ".pdf"), p, width = 6, height = 4.5, bg = "white")
ggsave(paste0(out, ".png"), p, width = 6, height = 4.5, dpi = 300, bg = "white")
cat("Wrote", paste0(out, c(".pdf", ".png")), sep = "\n  ")
