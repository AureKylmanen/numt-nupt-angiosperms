#!/usr/bin/env Rscript
# 02_pgls_models.R -- PGLS models for the manuscript
#
#   Rscript 02_pgls_models.R
#
# Needs organelle_lengths_final.tsv
# and final_species_summary_v2.tsv
# in $NUXT_OUT (default /vol/storage/5T_disk/organelles/pgls_v2).
# Writes to $NUXT_OUT:
#   pgls_table.tsv        one row per model term (Supplemental Table, tab:pgls)
#   lambda_profiles.tsv   slope and P at every lambda (Supplemental Table S4)
#   pgls_fits.rds         fitted objects, read by 03_plot_pgls_panels.R
#   excluded_species.tsv  who was dropped from which class and why

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
source(file.path(here, "common.R"))

d  <- load_data()
df <- d$df; tr <- d$tree

models <- list(
  # nuclear genome size (proportion and bp versions; bp slope = prop slope + 1)
  numt_prop_genome  = list(log_numt_prop ~ log_genome),
  nupt_prop_genome  = list(log_nupt_prop ~ log_genome),
  numt_bp_genome    = list(log_numt_bp   ~ log_genome),
  nupt_bp_genome    = list(log_nupt_bp   ~ log_genome),
  # donor organelle size
  numt_prop_mito    = list(log_numt_prop ~ log_mito),
  # Silene conica (11.3 Mb, 128 chromosomes) is the extreme point on this axis
  numt_prop_mito_noSconica = list(log_numt_prop ~ log_mito, subset = df$species != "Silene_conica"),
  nupt_prop_plastid = list(log_nupt_prop ~ log_plastid),
  nupt_prop_plastid_noparasites = list(log_nupt_prop ~ log_plastid, subset = !df$parasite),
  # covariation between classes, no shared denominator
  nupt_prop_numt_prop = list(log_nupt_prop ~ log_numt_prop),   # Figure 1c
  nupt_bp_numt_bp   = list(log_nupt_bp ~ log_numt_bp),            # same, no shared denominator
  # are donor size and nuclear size confounded?
  mito_genome       = list(log_mito    ~ log_genome),
  plastid_genome    = list(log_plastid ~ log_genome),
  # joint models: donor effect with nuclear size held constant
  numt_bp_genome_mito    = list(log_numt_bp ~ log_genome + log_mito),
  nupt_bp_genome_plastid = list(log_nupt_bp ~ log_genome + log_plastid)
)

fits <- list()
for (nm in names(models)) {
  m <- models[[nm]]
  cat("fitting", nm, "\n")
  fits[[nm]] <- fit_profile(m[[1]], df, tr, nm, subset = m$subset)
}
fits <- Filter(Negate(is.null), fits)

tab  <- do.call(rbind, lapply(fits, `[[`, "table"))
prof <- do.call(rbind, lapply(fits, `[[`, "profile"))
rownames(tab) <- NULL

write.table(tab,  file.path(OUT, "pgls_table.tsv"),      sep = "\t", quote = FALSE, row.names = FALSE)
write.table(prof, file.path(OUT, "lambda_profiles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(fits, file.path(OUT, "pgls_fits.rds"))

excl <- df[df$exclude_numt | df$exclude_nupt,
           c("species", "exclude_numt", "exclude_nupt", "mito_length", "plastid_length")]
write.table(excl, file.path(OUT, "excluded_species.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fmtP <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))
cat("\n===== PGLS at ML lambda (95% profile CI) | OLS =====\n")
for (i in seq_len(nrow(tab))) with(tab[i, ], cat(sprintf(
  "%-30s %-12s n=%3d  lambda=%.2f (%.2f-%.2f)%s  b=%+.3f SE=%.3f P=%s | P(b=1)=%s | OLS b=%+.3f P=%s | %s\n",
  model, term, n, lambda, lambda_lo, lambda_hi, ifelse(lambda_boundary, "*", ""),
  slope, SE, fmtP(P), fmtP(P_vs_1), OLS_slope, fmtP(OLS_P), robust_to_lambda)))
cat("\n* lambda estimate on the 0/1 boundary\n")
cat("P(b=1) is only meaningful for bp ~ size models (test of isometry).\n")
cat("Excluded species written to excluded_species.tsv\n")
