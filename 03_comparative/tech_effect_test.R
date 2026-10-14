#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Does sequencing technology (long-read / hybrid / short-read) affect the
# NUMT and NUPT content we detect?
#
# Layers, from weakest to strongest control:
#   1. Descriptives + Kruskal-Wallis / Wilcoxon (no covariates)
#   2. OLS: log10(y) ~ tech + log10(genome size)  [+ log10(N50)]
#   3. PGLS (corPagel) on the same models, LRT for the tech term
#   4. Sensitivity: drop assemblies with >50% of insertion bp on
#      organelle-dominated scaffolds (retained organelle contigs inflate
#      content and could mask a short-read deficit)
#   5. Within-genus paired comparison (genera with both short-read and
#      long-read-based assemblies)
#   6. Insertion-level prediction: if short-read assemblers collapse recent
#      insertions, the deficit should be in high-identity and long insertions
#      (needs insertion_metrics.tsv from tech_insertion_metrics.py)
#
# Usage:
#   Rscript tech_effect_test.R [summary] [assembly_stats] [tree] [metrics] [outdir]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ape); library(nlme); library(ggplot2)
})

args  <- commandArgs(trailingOnly = TRUE)
f_sum <- if (length(args) >= 1) args[1] else "final_species_summary_v2.tsv"
f_asm <- if (length(args) >= 2) args[2] else "assembly_technology.tsv"
f_tre <- if (length(args) >= 3) args[3] else "/vol/storage/5T_disk/organelles/tree/species_names_grafted.nwk"
f_met <- if (length(args) >= 4) args[4] else "insertion_metrics.tsv"
outd  <- if (length(args) >= 5) args[5] else "tech_effect"
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) read.delim(f, quote = "", comment.char = "",
                             stringsAsFactors = FALSE, na.strings = c("NA", ""))
sm  <- rd(f_sum)
asm <- rd(f_asm)
d   <- merge(sm, asm[, c("species", "assembly", "n50", "tech_class",
                          "sequencing_tech", "numt_pct_suspect_bp",
                          "nupt_pct_suspect_bp")], by = "species")
if (file.exists(f_met)) d <- merge(d, rd(f_met), by = "species", all.x = TRUE)

cat("Species in summary:", nrow(sm), " merged with assembly stats:", nrow(d), "\n")
cat("\ntech_class counts (all):\n"); print(table(d$tech_class, useNA = "ifany"))

keep <- c("long_read", "hybrid", "short_read")
dropped <- d[!(d$tech_class %in% keep), c("species", "tech_class", "sequencing_tech")]
if (nrow(dropped)) {
  cat("\nExcluded (other/unknown/missing) -- check these by hand:\n")
  print(dropped, row.names = FALSE)
}
d <- d[d$tech_class %in% keep, ]
d$tech  <- factor(d$tech_class, levels = keep)                 # long_read = reference
d$tech2 <- factor(ifelse(d$tech_class == "short_read", "short_read", "long_based"),
                  levels = c("long_based", "short_read"))
d$lgs   <- log10(d$genome_size)
d$ln50  <- log10(d$n50)
d$genus <- sub("_.*", "", d$species)

tree <- read.tree(f_tre)
d <- d[d$species %in% tree$tip.label, , drop = FALSE]
cat("\nAnalysed (in tree, tech known):", nrow(d), "\n")
print(table(d$tech))

# Responses: log10 proportion and log10 merged count; count models always
# carry genome size, so the tech coefficient is a per-bp effect in both.
lg <- function(x) { x[!is.finite(x) | x <= 0] <- NA; log10(x) }
d$y_numt_prop  <- lg(d$numt_merged_prop)
d$y_nupt_prop  <- lg(d$nupt_merged_prop)
d$y_numt_count <- lg(d$numt_merged_count)
d$y_nupt_count <- lg(d$nupt_merged_count)
resp <- c(numt_prop = "y_numt_prop", nupt_prop = "y_nupt_prop",
          numt_count = "y_numt_count", nupt_count = "y_nupt_count")

sink(file.path(outd, "tech_effect_results.txt"), split = TRUE)

# ---- 1. Descriptives and rank tests ---------------------------------------
cat("\n==== 1. Descriptives (medians) ====\n")
desc <- aggregate(cbind(numt_merged_count, nupt_merged_count,
                        numt_merged_prop, nupt_merged_prop,
                        genome_size / 1e6, n50 / 1e6) ~ tech, data = d, FUN = median)
names(desc)[6:7] <- c("genome_Mb", "N50_Mb")
print(desc, digits = 3)
cat("\nGenome size differs by tech? Kruskal-Wallis p =",
    signif(kruskal.test(lgs ~ tech, d)$p.value, 3), "\n")
cat("N50 differs by tech?         Kruskal-Wallis p =",
    signif(kruskal.test(ln50 ~ tech, d)$p.value, 3), "\n")

for (nm in names(resp)) {
  y <- d[[resp[nm]]]
  kw <- kruskal.test(y ~ d$tech)
  wx <- wilcox.test(y ~ d$tech2)
  cat(sprintf("%-11s KW p=%-9.3g  short vs long-based Wilcoxon p=%.3g\n",
              nm, kw$p.value, wx$p.value))
}

# ---- helpers ---------------------------------------------------------------
fit_pgls <- function(form, dd, tr) {
  dd <- dd[complete.cases(dd[, all.vars(form)]), ]
  tr <- keep.tip(tr, dd$species)
  m <- tryCatch(gls(form, data = dd, method = "ML",
                    correlation = corPagel(0.5, tr, form = ~species)),
                error = function(e) NULL)
  lam_how <- "estimated"
  if (!is.null(m)) {                      # lambda outside [0,1] is not interpretable
    l <- as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE))
    if (l < 0 || l > 1) m <- NULL
  }
  if (is.null(m)) {                       # grid-profile lambda when free fit fails
    lam_how <- "grid-profiled"
    grid <- seq(0, 1, by = 0.02); best <- NULL
    for (l in grid) {
      mm <- tryCatch(gls(form, data = dd, method = "ML",
                         correlation = corPagel(l, tr, form = ~species, fixed = TRUE)),
                     error = function(e) NULL)
      if (!is.null(mm) && (is.null(best) || logLik(mm) > logLik(best))) best <- mm
    }
    m <- best
  }
  attr(m, "lam_how") <- lam_how
  attr(m, "n") <- nrow(dd)
  m
}
lam_of <- function(m) as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE))

report_tech <- function(m, m0, label) {
  if (is.null(m) || is.null(m0)) { cat(label, ": fit failed\n"); return(invisible(NULL)) }
  tt  <- summary(m)$tTable
  ci  <- tryCatch(intervals(m, which = "coef")$coef, error = function(e) NULL)
  lr  <- 2 * (as.numeric(logLik(m)) - as.numeric(logLik(m0)))
  dfd <- length(coef(m)) - length(coef(m0))
  cat(sprintf("\n%s  (n=%d, lambda=%.3f %s)\n", label, attr(m, "n"),
              lam_of(m), attr(m, "lam_how")))
  cat(sprintf("  LRT tech term: chi2=%.2f, df=%d, p=%.3g\n",
              lr, dfd, pchisq(lr, dfd, lower.tail = FALSE)))
  for (r in grep("^tech", rownames(tt), value = TRUE)) {
    fc <- 10^tt[r, "Value"]
    lo <- if (!is.null(ci)) 10^ci[r, "lower"] else NA
    hi <- if (!is.null(ci)) 10^ci[r, "upper"] else NA
    cat(sprintf("  %-18s fold vs ref = %.2f [95%% CI %.2f-%.2f], p=%.3g\n",
                r, fc, lo, hi, tt[r, "p-value"]))
  }
  if ("lgs" %in% rownames(tt))
    cat(sprintf("  log10 genome size slope = %.3f (p=%.3g)\n",
                tt["lgs", "Value"], tt["lgs", "p-value"]))
  if ("ln50" %in% rownames(tt))
    cat(sprintf("  log10 N50 slope         = %.3f (p=%.3g)\n",
                tt["ln50", "Value"], tt["ln50", "p-value"]))
  invisible(NULL)
}

run_block <- function(dd, title) {
  cat("\n====", title, "====\n")
  for (nm in names(resp)) {
    y <- resp[nm]
    for (tv in c("tech", "tech2")) {
      for (withN50 in c(FALSE, TRUE)) {
        rhs0 <- if (withN50) "lgs + ln50" else "lgs"
        f1 <- as.formula(paste(y, "~", tv, "+", rhs0))
        f0 <- as.formula(paste(y, "~", rhs0))
        # same rows for both models
        sub <- dd[complete.cases(dd[, c(y, tv, "lgs", "ln50")]), ]
        # OLS
        o1 <- lm(f1, sub); o0 <- lm(f0, sub)
        p_ols <- anova(o0, o1)$`Pr(>F)`[2]
        # PGLS
        m1 <- fit_pgls(f1, sub, tree); m0 <- fit_pgls(f0, sub, tree)
        lab <- sprintf("%s ~ %s + %s", nm, tv, rhs0)
        cat(sprintf("\n[OLS] %s: F-test tech p=%.3g", lab, p_ols))
        report_tech(m1, m0, paste("[PGLS]", lab))
      }
    }
  }
}

# ---- 2-3. OLS and PGLS, full set -------------------------------------------
run_block(d, "2-3. OLS and PGLS, all assemblies")

# ---- 4. Sensitivity: drop assemblies dominated by organelle-like scaffolds -
cat("\n==== 4. Sensitivity: excluding >50% suspect-scaffold bp (per class) ====\n")
resp_all <- resp
for (cls in c("numt", "nupt")) {
  col <- paste0(cls, "_pct_suspect_bp")
  dd <- d[!is.na(d[[col]]) & d[[col]] <= 50, ]
  cat(sprintf("\n%s: %d of %d species retained\n", toupper(cls), nrow(dd), nrow(d)))
  print(table(dd$tech))
  resp <- resp_all[grepl(cls, names(resp_all))]
  run_block(dd, paste("4.", toupper(cls), "low-suspect subset"))
}
resp <- resp_all

# ---- 5. Within-genus paired comparison -------------------------------------
cat("\n==== 5. Within-genus comparison (short_read vs long_based) ====\n")
for (nm in c("numt_prop", "nupt_prop")) {
  y <- resp[nm]
  sub <- d[!is.na(d[[y]]), ]
  sub$res <- resid(lm(as.formula(paste(y, "~ lgs")), sub))   # genome-size-corrected
  g <- aggregate(res ~ genus + tech2, sub, mean)
  w <- reshape(g, idvar = "genus", timevar = "tech2", direction = "wide")
  w <- w[complete.cases(w), ]
  if (nrow(w) < 3) { cat(nm, ": only", nrow(w), "genera with both classes\n"); next }
  diff <- w$res.short_read - w$res.long_based
  wt <- wilcox.test(diff, mu = 0)
  cat(sprintf("%s: %d genera; median short - long = %.3f log10 units (fold %.2f); signed-rank p=%.3g\n",
              nm, nrow(w), median(diff), 10^median(diff), wt$p.value))
  cat("  genera:", paste(w$genus, collapse = ", "), "\n")
}

# ---- 6. Insertion-level prediction -----------------------------------------
cat("\n==== 6. Insertion-level: high-identity and long insertions ====\n")
lgt <- function(p) { p <- pmin(pmax(p, 1e-4), 1 - 1e-4); qlogis(p) }
ins <- c(numt_frac_bp_ge99 = "numt_frac_bp_ge99", nupt_frac_bp_ge99 = "nupt_frac_bp_ge99",
         numt_frac_bp_ge5kb = "numt_frac_bp_ge5kb", nupt_frac_bp_ge5kb = "nupt_frac_bp_ge5kb")
ins <- ins[ins %in% names(d)]
if (!length(ins)) {
  cat("insertion_metrics.tsv not found / no metric columns; skipped\n")
} else {
  for (nm in names(ins)) {
    d[[paste0("z_", nm)]] <- lgt(as.numeric(d[[ins[nm]]]))
    y <- paste0("z_", nm)
    med <- tapply(as.numeric(d[[ins[nm]]]), d$tech, median, na.rm = TRUE)
    cat(sprintf("\n%s medians: %s\n", nm,
                paste(sprintf("%s=%.3f", names(med), med), collapse = "  ")))
    f1 <- as.formula(paste(y, "~ tech2 + lgs")); f0 <- as.formula(paste(y, "~ lgs"))
    sub <- d[complete.cases(d[, c(y, "tech2", "lgs")]), ]
    m1 <- fit_pgls(f1, sub, tree); m0 <- fit_pgls(f0, sub, tree)
    if (!is.null(m1) && !is.null(m0)) {
      tt <- summary(m1)$tTable
      lr <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
      cat(sprintf("  PGLS logit(%s) ~ tech2 + lgs: short_read coef=%.3f (odds ratio %.2f), LRT p=%.3g, lambda=%.3f\n",
                  nm, tt["tech2short_read", "Value"], exp(tt["tech2short_read", "Value"]),
                  pchisq(lr, 1, lower.tail = FALSE), lam_of(m1)))
    }
  }
}
sink()

# ---- Figure ----------------------------------------------------------------
long <- rbind(
  data.frame(species = d$species, tech = d$tech, class = "NUMT", y = d$y_numt_prop),
  data.frame(species = d$species, tech = d$tech, class = "NUPT", y = d$y_nupt_prop))
long$res <- NA
for (cl in c("NUMT", "NUPT")) {
  i <- long$class == cl & !is.na(long$y)
  long$res[i] <- resid(lm(long$y[i] ~ d$lgs[match(long$species[i], d$species)]))
}
lab_n <- sapply(split(d$tech, d$tech), length)
levels(long$tech) <- sprintf("%s\n(n=%d)", levels(long$tech), lab_n[levels(long$tech)])

p1 <- ggplot(long[!is.na(long$y), ], aes(tech, y)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.4) +
  facet_wrap(~class) +
  labs(x = NULL, y = expression(log[10]~"proportion of assembly (merged)"),
       title = "Raw") + theme_bw(base_size = 10)
p2 <- ggplot(long[!is.na(long$res), ], aes(tech, res)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.4) +
  facet_wrap(~class) +
  labs(x = NULL, y = "residual after log10 genome size",
       title = "Genome-size corrected") + theme_bw(base_size = 10)
pdf(file.path(outd, "tech_effect_boxplots.pdf"), width = 7, height = 7)
print(p1); print(p2)
invisible(dev.off())
cat("\nWrote", file.path(outd, "tech_effect_results.txt"), "and",
    file.path(outd, "tech_effect_boxplots.pdf"), "\n")
