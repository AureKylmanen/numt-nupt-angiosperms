#!/usr/bin/env Rscript
# assembly_quality.R
# Does assembly quality explain the NUMT/NUPT insertion proportions, and does
# it explain the difference in phylogenetic signal between them?
#
#   Rscript assembly_quality.R
#
# Inputs ($NUXT_BASE, default /vol/storage/5T_disk/organelles):
#   final_species_summary_v2.tsv   merged NUMT/NUPT proportions
#   tree/assembly_stats.tsv        00_data/make_assembly_stats.py (+ technology column)
#   tree/species_names_grafted.nwk TimeTree backbone with grafted species; pruned to
#                                  the analysed species and made ultrametric by
#                                  extending terminal branches, as in common.R
#
# Frames assembly quality as a measured covariate rather than asserting
# contamination: every species is retained and the effect is quantified.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(nlme)
})

BASE <- Sys.getenv("NUXT_BASE", Sys.getenv("BASE_DIR", getwd()))
DIR  <- file.path(BASE, "tree")
QC   <- file.path(DIR, "assembly_stats.tsv")

## ---- column mapping: EDIT to match your file ------------------------------
COL <- list(
  species   = "species",
  n50       = "scaffold_n50",     # scaffold N50 in bp
  n_scaff   = "n_scaffolds",      # optional, set NA to skip
  tech      = "sequencing_tech",  # long / hybrid / short
  org_frac  = "organelle_scaffold_frac"  # optional; fraction of insertion bp on
)                                        # organelle-dominated scaffolds

tr   <- read.tree(file.path(DIR, "species_names_grafted.nwk"))
extend_ultrametric <- function(tr) {   # = phytools::force.ultrametric(method = "extend")
  nt <- Ntip(tr); depth <- node.depth.edgelength(tr)[seq_len(nt)]
  term <- match(seq_len(nt), tr$edge[, 2])
  tr$edge.length[term] <- tr$edge.length[term] + (max(depth) - depth); tr
}
summ <- read.delim(file.path(BASE, "final_species_summary_v2.tsv"))
numt <- data.frame(species = summ$species, NUMTs = summ$numt_merged_prop)
nupt <- data.frame(species = summ$species, NUPTs = summ$nupt_merged_prop)

if (!file.exists(QC)) stop("Need ", QC, " with at least: ",
                           paste(unlist(COL[c("species","n50","tech")]), collapse=", "))
qc <- read.delim(QC)
cat("Columns in", basename(QC), ":\n  ", paste(names(qc), collapse=", "), "\n\n")
if (!COL$n50 %in% names(qc))
  stop("Need at least a scaffold-N50 column; got: ", paste(names(qc), collapse=", "))
HAS_TECH <- !is.na(COL$tech) && COL$tech %in% names(qc)
if (!HAS_TECH) cat("NOTE: no sequencing-technology column; using contiguity only.\n",
                   "      To add it, rerun make_assembly_stats.py with --tech-col,\n",
                   "      or merge it in from combined_dataset.csv.\n\n")

d <- Reduce(function(a,b) merge(a,b,by="species"),
            list(numt[,c("species","NUMTs")], nupt[,c("species","NUPTs")], qc))
d <- d[d$species %in% tr$tip.label, ]
tr <- extend_ultrametric(keep.tip(tr, d$species)); d <- d[match(tr$tip.label, d$species), ]
rownames(d) <- d$species
cat("Species with tree + traits + assembly stats:", nrow(d), "\n")

lg <- function(x) { x <- as.numeric(x); x[x <= 0] <- NA; log10(x) }
d$lNUMT <- lg(d$NUMTs); d$lNUPT <- lg(d$NUPTs)
d$lN50 <- lg(d[[COL$n50]])

## quality predictors, from whatever the stats file actually has
QCOV <- "lN50"
if (!is.na(COL$n_scaff) && COL$n_scaff %in% names(d)) {
  d$lNSCAF <- lg(d[[COL$n_scaff]]); QCOV <- c(QCOV, "lNSCAF")
}
if ("frac_bp_in_top10" %in% names(d)) {
  d$ftop10 <- as.numeric(d$frac_bp_in_top10); QCOV <- c(QCOV, "ftop10")
}
if (HAS_TECH) {
  d$tech <- factor(d[[COL$tech]]); QCOV <- c(QCOV, "tech")
  cat("Technology:", paste(names(table(d$tech)), table(d$tech),
                           sep="=", collapse="  "), "\n")
}
QF <- paste(QCOV, collapse = " + ")
cat("Quality predictors:", QF, "\n")
cat("Correlations among them:\n")
print(round(cor(d[, intersect(QCOV, c("lN50","lNSCAF","ftop10"))],
                use = "pairwise.complete.obs"), 3))
cat("\n")

lam_of <- function(m) {
  l <- tryCatch(coef(m$modelStruct$corStruct, unconstrained = FALSE), error = function(e) NA)
  as.numeric(l)[1]
}

## corPagel's optimiser is sensitive to the starting lambda; retry across a
## few starts before declaring failure.
fit_pagel <- function(f, sub, t2, starts = c(1, 0.8, 0.5, 0.3, 0.1, 0.95)) {
  for (s0 in starts) {
    m <- try(suppressWarnings(
           gls(f, data = sub, method = "ML",
               correlation = corPagel(s0, t2, form = ~species))),
         silent = TRUE)
    if (!inherits(m, "try-error")) {
      l <- tryCatch(as.numeric(coef(m$modelStruct$corStruct,
                                    unconstrained = FALSE))[1],
                    error = function(e) NA)
      if (is.finite(l)) return(m)
    }
  }
  NULL
}

pgls <- function(f, dat, label) {
  sub <- dat[complete.cases(dat[, all.vars(f)]), ]
  t2  <- keep.tip(tr, sub$species)
  m <- fit_pagel(f, sub, t2)
  if (is.null(m)) { cat(sprintf("%-46s FAILED (no start converged)\n", label)); return(NULL) }
  s <- summary(m)$tTable
  cat(sprintf("%-46s n=%3d  lambda=%.3f\n", label, nrow(sub), lam_of(m)))
  if (nrow(s) > 1) for (i in seq_len(nrow(s))[-1])
    cat(sprintf("       %-28s beta=%7.3f  P=%.4g\n", rownames(s)[i], s[i,1], s[i,4]))
  invisible(m)
}

cat("=== T0: is assembly quality itself phylogenetically structured? =========\n")
x <- setNames(d$lN50, d$species); x <- x[is.finite(x)]
k <- phylosig(keep.tip(tr, names(x)), x, method="K", test=TRUE, nsim=1000)
l <- phylosig(keep.tip(tr, names(x)), x, method="lambda", test=TRUE)
cat(sprintf("log10 N50:  K=%.3f (P=%.3g)   lambda=%.3f (P=%.3g)\n", k$K,k$P,l$lambda,l$P))
cat("  lambda near 0  -> quality noise is diagonal, biases trait lambda DOWN\n")
cat("  lambda high    -> quality is clade-structured, can bias trait lambda UP\n\n")

cat("=== T1: does assembly quality predict insertion content? ================\n")
pgls(lNUMT ~ lN50,        d, "NUMT ~ N50")
pgls(lNUPT ~ lN50,        d, "NUPT ~ N50")
cat("\n  -- one predictor at a time (joint model is collinear) --\n")
for (v in intersect(QCOV, c("lN50","lNSCAF","ftop10"))) {
  pgls(as.formula(paste("lNUMT ~", v)), d, paste("NUMT ~", v))
  pgls(as.formula(paste("lNUPT ~", v)), d, paste("NUPT ~", v))
}
cat("\n  -- all together --\n")
pgls(as.formula(paste("lNUMT ~", QF)), d, "NUMT ~ quality")
pgls(as.formula(paste("lNUPT ~", QF)), d, "NUPT ~ quality")
cat("\n  Contamination predicts NEGATIVE slope on N50 and POSITIVE on scaffold\n")
cat("  count (stray organellar contigs fragment the assembly and add bp).\n")
cat("  Opposite signs argue against contamination as the driver.\n")
cat("\n")

cat("=== T2: does the lambda gap survive conditioning on quality? ============\n")
m0n <- pgls(lNUMT ~ 1,           d, "NUMT ~ 1                 (baseline)")
m0p <- pgls(lNUPT ~ 1,           d, "NUPT ~ 1                 (baseline)")
m1n <- pgls(as.formula(paste("lNUMT ~", QF)), d, "NUMT ~ quality   (adjusted)")
m1p <- pgls(as.formula(paste("lNUPT ~", QF)), d, "NUPT ~ quality   (adjusted)")
cat(sprintf("\n  NUMT lambda: %.3f -> %.3f   NUPT lambda: %.3f -> %.3f\n",
            lam_of(m0n), lam_of(m1n), lam_of(m0p), lam_of(m1p)))
cat(sprintf("  gap: %.3f -> %.3f\n", lam_of(m0n)-lam_of(m0p), lam_of(m1n)-lam_of(m1p)))
cat("  gap unchanged -> the difference is not an assembly artefact\n")
cat("  gap shrinks   -> quantify and report as a limitation\n\n")

if (!is.na(COL$org_frac) && COL$org_frac %in% names(d)) {
  cat("=== T2b: organelle-dominated scaffold fraction as covariate =============\n")
  d$orgf <- as.numeric(d[[COL$org_frac]])
  m2n <- pgls(as.formula(paste("lNUMT ~", QF, "+ orgf")), d, "NUMT ~ quality + org.frac")
  m2p <- pgls(as.formula(paste("lNUPT ~", QF, "+ orgf")), d, "NUPT ~ quality + org.frac")
  cat(sprintf("\n  lambda: NUMT %.3f   NUPT %.3f   gap %.3f\n\n",
              lam_of(m2n), lam_of(m2p), lam_of(m2n)-lam_of(m2p)))
}

cat("=== T4: log-ratio test (cancels anything acting equally on both) ========\n")
cat("  Response = log10(NUPT prop / NUMT prop). Genome size, contiguity,\n")
cat("  coverage and the shared denominator all cancel. A non-zero slope means\n")
cat("  fragmentation affects the two classes UNEQUALLY.\n")
d$lRATIO <- d$lNUPT - d$lNUMT
pgls(lRATIO ~ lN50,        d, "NUPT:NUMT ratio ~ N50")
pgls(as.formula(paste("lRATIO ~", QF)), d, "NUPT:NUMT ratio ~ quality")
xr <- setNames(d$lRATIO, d$species); xr <- xr[is.finite(xr)]
lr <- phylosig(keep.tip(tr, names(xr)), xr, method = "lambda", test = TRUE)
cat(sprintf("  ratio lambda = %.3f (P = %.3g)\n\n", lr$lambda, lr$P))

cat("=== T5: phylogeny vs fragmentation - variance partition ================\n")
part <- function(trait, label) {
  f0 <- as.formula(paste(trait, "~ 1"))
  f1 <- as.formula(paste(trait, "~", QF))
  sub <- d[complete.cases(d[, c(trait, QCOV)]), ]
  t2  <- keep.tip(tr, sub$species)
  # OLS (no phylogeny) vs PGLS (phylogeny), with and without fragmentation
  o0 <- gls(f0, data = sub, method = "ML")
  o1 <- gls(f1, data = sub, method = "ML")
  p0 <- fit_pagel(f0, sub, t2)
  p1 <- fit_pagel(f1, sub, t2)
  if (is.null(p0) || is.null(p1)) {
    cat(sprintf("\n  %s: PGLS did not converge, skipped\n", label)); return(invisible(NULL))
  }
  ll <- sapply(list(o0, o1, p0, p1), logLik)
  ai <- sapply(list(o0, o1, p0, p1), AIC)
  cat(sprintf("\n  %s (n=%d)\n", label, nrow(sub)))
  cat(sprintf("    %-26s logLik=%9.2f  AIC=%9.2f\n",
              c("null (neither)", "fragmentation only",
                "phylogeny only", "both"), ll, ai))
  cat(sprintf("    delta-AIC vs null:  frag %+.1f   phylo %+.1f   both %+.1f\n",
              ai[2]-ai[1], ai[3]-ai[1], ai[4]-ai[1]))
  cat(sprintf("    fragmentation gain GIVEN phylogeny: %+.1f AIC\n", ai[4]-ai[3]))
  cat(sprintf("    phylogeny gain GIVEN fragmentation: %+.1f AIC\n", ai[4]-ai[2]))
  cat(sprintf("    lambda: %.3f (alone) -> %.3f (adjusted)\n",
              lam_of(p0), lam_of(p1)))
}
part("lNUMT", "NUMT proportion")
part("lNUPT", "NUPT proportion")
cat("\n  Negative delta-AIC favours the model. If 'phylogeny gain GIVEN\n")
cat("  fragmentation' stays strongly negative, the signal is phylogenetic and\n")
cat("  not an assembly artefact.\n\n")

cat("=== T3: lambda on high-quality assemblies only ==========================\n")
for (thr in c(1e6, 1e7)) {
  s <- d[is.finite(d$lN50) & d[[COL$n50]] >= thr, ]
  if (nrow(s) < 50) { cat(sprintf("N50 >= %.0e: only %d species, skipped\n", thr, nrow(s))); next }
  t2 <- keep.tip(tr, s$species)
  ln <- phylosig(t2, setNames(s$lNUMT, s$species), method="lambda", test=TRUE)
  lp <- phylosig(t2, setNames(s$lNUPT, s$species), method="lambda", test=TRUE)
  cat(sprintf("N50 >= %.0e  n=%3d   lambda NUMT=%.3f  NUPT=%.3f  gap=%.3f\n",
              thr, nrow(s), ln$lambda, lp$lambda, ln$lambda - lp$lambda))
}
cat("\nFull-dataset reference: NUMT 0.682 (0.539-0.798), NUPT 0.381 (0.202-0.564)\n")
