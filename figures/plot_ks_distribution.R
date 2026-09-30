#!/usr/bin/env Rscript
# plot_ks_distribution.R -- Ks (dS) histogram of coding insertions, one class
#
# Replaces numt_ks_distribution.R and nupt_ks_distribution.R, which differed
# only in input suffix, fill colour and output name.
#
# Usage:
#   Rscript plot_ks_distribution.R numt
#   Rscript plot_ks_distribution.R nupt
#   Rscript plot_ks_distribution.R nupt --suffix _kaks_c11_fixed.csv
#   Rscript plot_ks_distribution.R numt --base /vol/storage/5T_disk/organelles --outdir Figures
#
# Reads   <base>/species/<SPECIES>/results/<SPECIES><suffix>
#         (default suffix: _numt_kaks.csv / _nupt_kaks.csv, as written by
#          02_cnuxt_kaks/*_gene_alignment.sh)
# Writes  <outdir>/ks_distribution_<class>.pdf and .png
#
# Values above 0.6 are counted and reported but not drawn; the dashed line
# marks the dcNUxT threshold (Ks = 0.1).
#
# Binning fix relative to the old per-class scripts: those used
# geom_histogram(binwidth = 0.01) with scale_x_continuous(limits = c(0, 0.6)).
# ggplot centres bins on multiples of the binwidth, so the first bin spanned
# -0.005..0.005, fell outside the limits and was dropped without drawing
# (the "Removed 2 rows" warning) - every Ks <= 0.005 was missing and all other
# bars were offset by half a bin. Bins are now (0, 0.01], (0.01, 0.02], ...
# with Ks = 0 included in the first, and the axis is zoomed with
# coord_cartesian so no data are removed. Right-closed bins match the
# dcNUxT definition (Ks > 0.1): the bar ending at 0.1 holds non-dcNUxTs only.

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !args[1] %in% c("numt", "nupt"))
  stop("first argument must be 'numt' or 'nupt'")
cls <- args[1]

opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
base   <- opt("--base",   "/vol/storage/5T_disk/organelles")
suffix <- opt("--suffix", paste0("_", cls, "_kaks.csv"))
outdir <- opt("--outdir", base)
xmax   <- as.numeric(opt("--xmax", "0.6"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fill <- c(numt = "#2c728e", nupt = "#21918c")[[cls]]

ks_values <- numeric(0)
n_files <- 0
for (sp_dir in list.dirs(file.path(base, "species"), recursive = FALSE)) {
  sp <- basename(sp_dir)
  f <- file.path(sp_dir, "results", paste0(sp, suffix))
  if (!file.exists(f)) next
  df <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !"Ks" %in% names(df)) next
  ks <- suppressWarnings(as.numeric(df$Ks))
  ks_values <- c(ks_values, ks[!is.na(ks)])
  n_files <- n_files + 1
}
if (!length(ks_values)) stop("no Ks values found for pattern */results/*", suffix)

n_above <- sum(ks_values > xmax)
cat("Tables read:", n_files, "\n")
cat("Total Ks values:", length(ks_values), "\n")
cat("Ks > 0.1 (dcNUxT):", sum(ks_values > 0.1), "\n")
cat(sprintf("Values > %.1f (not shown): %d\n", xmax, n_above))

p <- ggplot(data.frame(ks = ks_values[ks_values <= xmax]), aes(x = ks)) +
  geom_histogram(binwidth = 0.01, boundary = 0, closed = "right",
                 fill = fill, color = NA) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, xmax, by = 0.1), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, xmax)) +
  labs(x = expression(italic(d)[S] ~ "values"), y = "Counts") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12),
    plot.margin = margin(10, 15, 10, 10)
  )

out <- file.path(outdir, paste0("ks_distribution_", cls))
ggsave(paste0(out, ".pdf"), p, width = 8, height = 5.5, bg = "white")
ggsave(paste0(out, ".png"), p, width = 8, height = 5.5, dpi = 300, bg = "white")
cat("Wrote", paste0(out, c(".pdf", ".png")), sep = "\n  ")
