#!/usr/bin/env Rscript
# mtpt_sensitivity.R
# Does double-counting of sequence hit by both organelle queries (NUMT and
# NUPT calls on the same nuclear bp, mostly via plastid-derived regions of the
# mitogenome = MTPTs) change the comparative results?
#
# Fits the same PGLS models on three versions of the insertion proportions:
#   merged     : merged (10 kb) proportions from final_species_summary_v2.tsv (= main analysis)
#   d0_total   : union of unmerged fragments / genome size
#   d0_excl    : (union of unmerged fragments - shared bp) / genome size
# merged vs d0_total isolates the effect of merging; d0_total vs d0_excl
# isolates the effect of the NUMT/NUPT overlap.
# All models are refitted here with the same lambda profile, so the three
# versions are compared like for like.
#
# Usage (from /vol/storage/5T_disk/organelles):
#   Rscript mtpt_sensitivity.R
# Outputs (in mtpt_sensitivity/):
#   pgls_by_version.tsv     slope, SE, P, ML lambda per model and version
#   ratio_signal.tsv        lambda of log10(NUPT/NUMT) per version
#   species_overlap.tsv     per-species shared fractions and MTPT coverage (supp table)
#   shared_vs_mtpt.pdf      shared fraction of NUPT bp vs plastome fraction in mitogenome

suppressPackageStartupMessages({ library(ape); library(nlme); library(ggplot2) })

BASE    <- Sys.getenv("BASE_DIR", getwd())
SUMMARY <- file.path(BASE, "final_species_summary_v2.tsv")
MTPT    <- file.path(BASE, "species/mtpt_coverage.tsv")
TREE    <- file.path(BASE, "tree/species_names_grafted.nwk")
OUT     <- file.path(BASE, "mtpt_sensitivity")
GRID    <- seq(0, 1, by = 0.02)
dir.create(OUT, showWarnings = FALSE)

pick <- function(df, cands, what) {
  hit <- cands[cands %in% names(df)]
  if (!length(hit)) stop("No column for ", what, " in ", deparse(substitute(df)),
                         ". Columns: ", paste(names(df), collapse = ", "))
  df[[hit[1]]]
}

## ---- data -----------------------------------------------------------------
s  <- read.delim(SUMMARY, check.names = FALSE)
d <- data.frame(
  species    = s$species,
  genome     = pick(s, c("genome_size"), "genome size"),
  mito       = pick(s, c("mito_length", "mito_size"), "mitogenome size"),
  plastid    = pick(s, c("plastid_length", "plastid_size"), "plastome size"),
  numt_merged = pick(s, c("numt_merged_prop"), "NUMT merged proportion"),
  nupt_merged = pick(s, c("nupt_merged_prop"), "NUPT merged proportion"))
# union of fragments per class and sequence detected by both queries
# (numt_nupt_detection.sh -> build_species_summary.py)
d$numt_d0_bp <- s$numt_union_bp
d$nupt_d0_bp <- s$nupt_union_bp
d$shared_bp  <- s$shared_bp
num <- c("genome","mito","plastid","numt_merged","nupt_merged","numt_d0_bp","nupt_d0_bp","shared_bp")
d[num] <- lapply(d[num], function(x) suppressWarnings(as.numeric(x)))

d$numt_d0      <- d$numt_d0_bp / d$genome
d$nupt_d0      <- d$nupt_d0_bp / d$genome
d$numt_d0_excl <- (d$numt_d0_bp - d$shared_bp) / d$genome
d$nupt_d0_excl <- (d$nupt_d0_bp - d$shared_bp) / d$genome
d$shared_frac_nupt <- d$shared_bp / d$nupt_d0_bp
d$shared_frac_numt <- d$shared_bp / d$numt_d0_bp

if (file.exists(MTPT)) {
  m <- read.delim(MTPT)
  m$plastid_len <- suppressWarnings(as.numeric(m$plastid_len))
  m$mtpt_cov_bp <- suppressWarnings(as.numeric(m$mtpt_cov_bp))
  m$mtpt_frac <- m$mtpt_cov_bp / m$plastid_len
  d <- merge(d, m[, c("species", "mtpt_cov_bp", "mtpt_frac")], by = "species", all.x = TRUE)
  cat("Plastome fraction present in the mitogenome: median", round(median(d$mtpt_frac, na.rm = TRUE), 3),
      "range", paste(round(range(d$mtpt_frac, na.rm = TRUE), 3), collapse = "-"), "\n")
  big <- "bp_in_hits_ge10kb_ge99.5"
  if (big %in% names(m)) cat("Species whose mitogenome has >=10 kb segments >=99.5% identical to the plastome:",
                             sum(as.numeric(m[[big]][m$species %in% d$species]) > 0, na.rm = TRUE), "\n")
} else d$mtpt_frac <- NA

## ---- tree -----------------------------------------------------------------
extend_ultrametric <- function(tr) {
  nt <- Ntip(tr); depth <- node.depth.edgelength(tr)[seq_len(nt)]
  term <- match(seq_len(nt), tr$edge[, 2])
  tr$edge.length[term] <- tr$edge.length[term] + (max(depth) - depth); tr
}
tr <- read.tree(TREE)
keep <- intersect(tr$tip.label, d$species)
cat("Species: summary+overlap =", nrow(d), "; also in tree =", length(keep), "\n")
tr <- keep.tip(tr, keep)
if (!is.ultrametric(tr)) tr <- extend_ultrametric(tr)
d <- d[match(keep, d$species), ]

## ---- lambda-profile PGLS ----------------------------------------------------
fit_profile <- function(y, x = NULL, dat) {
  ok <- is.finite(y) & (if (is.null(x)) TRUE else is.finite(x))
  dd <- data.frame(sp = dat$species[ok], y = y[ok])
  if (!is.null(x)) dd$x <- x[ok]
  t2 <- keep.tip(tr, dd$sp)
  f  <- if (is.null(x)) y ~ 1 else y ~ x
  ll <- sapply(GRID, function(l) {
    g <- tryCatch(gls(f, dd, correlation = corPagel(l, t2, form = ~sp, fixed = TRUE),
                      method = "ML"), error = function(e) NULL)
    if (is.null(g)) NA else as.numeric(logLik(g))
  })
  lam <- GRID[which.max(ll)]
  g <- gls(f, dd, correlation = corPagel(lam, t2, form = ~sp, fixed = TRUE), method = "ML")
  ci <- range(GRID[ll >= max(ll, na.rm = TRUE) - 1.92], na.rm = TRUE)
  out <- data.frame(n = nrow(dd), lambda = lam, lambda_lo = ci[1], lambda_hi = ci[2])
  if (!is.null(x)) {
    tt <- summary(g)$tTable
    out$slope <- tt["x", "Value"]; out$se <- tt["x", "Std.Error"]; out$p <- tt["x", "p-value"]
  }
  out
}

L <- function(v) ifelse(v > 0, log10(v), NA)
versions <- list(
  merged   = c(numt = "numt_merged", nupt = "nupt_merged"),
  d0_total = c(numt = "numt_d0",     nupt = "nupt_d0"),
  d0_excl  = c(numt = "numt_d0_excl", nupt = "nupt_d0_excl"))

res <- list(); rat <- list()
for (v in names(versions)) {
  nm <- L(d[[versions[[v]]["numt"]]]); np <- L(d[[versions[[v]]["nupt"]]])
  models <- list(
    "NUMT ~ genome size"  = list(nm, L(d$genome)),
    "NUPT ~ genome size"  = list(np, L(d$genome)),
    "NUMT ~ mito size"    = list(nm, L(d$mito)),
    "NUPT ~ plastid size" = list(np, L(d$plastid)),
    "NUPT ~ NUMT"         = list(np, nm))
  for (mn in names(models)) {
    t0 <- Sys.time()
    r <- fit_profile(models[[mn]][[1]], models[[mn]][[2]], d)
    res[[length(res) + 1]] <- cbind(version = v, model = mn, r)
    cat(sprintf("%-9s %-20s slope=%7.3f p=%.3g lambda=%.2f (%.0fs)\n", v, mn, r$slope, r$p,
                r$lambda, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  rr <- fit_profile(np - nm, NULL, d)
  rat[[v]] <- cbind(version = v, model = "log10(NUPT/NUMT) ~ 1", rr)
  cat(sprintf("%-9s ratio lambda=%.2f (%.2f-%.2f)\n", v, rr$lambda, rr$lambda_lo, rr$lambda_hi))
}
res <- do.call(rbind, res); rat <- do.call(rbind, rat)
write.table(res, file.path(OUT, "pgls_by_version.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rat, file.path(OUT, "ratio_signal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- shared fraction vs MTPT content ----------------------------------------
cat("\nShared fraction of NUPT bp: median", round(median(d$shared_frac_nupt, na.rm = TRUE), 3),
    "IQR", paste(round(quantile(d$shared_frac_nupt, c(.25, .75), na.rm = TRUE), 3), collapse = "-"), "\n")
cat("Shared fraction of NUMT bp: median", round(median(d$shared_frac_numt, na.rm = TRUE), 3),
    "IQR", paste(round(quantile(d$shared_frac_numt, c(.25, .75), na.rm = TRUE), 3), collapse = "-"), "\n")
okm <- is.finite(d$mtpt_frac) & is.finite(d$shared_frac_nupt)
cat("Species with MTPT coverage:", sum(okm), "\n")
if (sum(okm) > 5) {
  print(cor.test(d$mtpt_frac[okm], d$shared_frac_nupt[okm], method = "spearman", exact = FALSE))
  p <- ggplot(d[okm, ], aes(mtpt_frac, shared_frac_nupt)) +
    geom_abline(linetype = 2, colour = "grey50") + geom_point(alpha = .6, size = 1.4) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() +
    labs(x = "Fraction of plastome present in mitogenome",
         y = "Fraction of NUPT sequence also detected as NUMT")
  ggsave(file.path(OUT, "shared_vs_mtpt.pdf"), p, width = 4.5, height = 4.5)
  d$above_line <- d$shared_frac_nupt - d$mtpt_frac
  cat("\nLargest excess over MTPT prediction (check scaffold composition):\n")
  print(head(d[okm, c("species", "mtpt_frac", "shared_frac_nupt", "above_line")][
    order(-d$above_line[okm]), ], 15), row.names = FALSE)
}
write.table(d[, c("species", "numt_d0_bp", "nupt_d0_bp", "shared_bp", "shared_frac_numt",
                  "shared_frac_nupt", "mtpt_cov_bp", "mtpt_frac")[
                    c("species", "numt_d0_bp", "nupt_d0_bp", "shared_bp", "shared_frac_numt",
                      "shared_frac_nupt", "mtpt_cov_bp", "mtpt_frac") %in% names(d)]],
            file.path(OUT, "species_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nDone:", OUT, "\n")
