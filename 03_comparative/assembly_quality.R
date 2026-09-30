#!/usr/bin/env Rscript
# assembly_quality.R
# Does assembly quality explain the NUMT/NUPT insertion proportions, and does
# it explain the difference in phylogenetic signal between them?
#
#   Rscript assembly_quality.R
#
# Frames assembly quality as a measured covariate rather than asserting
# contamination: every species is retained and the effect is quantified.

suppressPackageStartupMessages({
  library(ape); library(phytools); library(nlme)
})

DIR <- "/vol/storage/5T_disk/organelles/tree"
QC  <- file.path(DIR, "assembly_stats.tsv")

## ---- column mapping: EDIT to match your file ------------------------------
COL <- list(
  species   = "species",
  n50       = "scaffold_n50",     # scaffold N50 in bp
  n_scaff   = "n_scaffolds",      # optional, set NA to skip
  tech      = "sequencing_tech",  # long / hybrid / short
  org_frac  = "organelle_scaffold_frac"  # optional; fraction of insertion bp on
)                                        # organelle-dominated scaffolds

tr   <- read.tree(file.path(DIR, "tree_analysis.nwk"))
numt <- read.delim(file.path(DIR, "numt_merged_proportions.tsv"))
nupt <- read.delim(file.path(DIR, "nupt_merged_proportions.tsv"))

if (!file.exists(QC)) stop("Need ", QC, " with at least: ",
                           paste(unlist(COL[c("species","n50","tech")]), collapse=", "))
qc <- read.delim(QC)
cat("Columns in", basename(QC), ":\n  ", paste(names(qc), collapse=", "), "\n\n")
miss <- setdiff(unlist(COL[c("species","n50","tech")]), names(qc))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse=", "),
                       "\n  Edit COL at the top of this script.")

d <- Reduce(function(a,b) merge(a,b,by="species"),
            list(numt[,c("species","NUMTs")], nupt[,c("species","NUPTs")], qc))
d <- d[d$species %in% tr$tip.label, ]
tr <- keep.tip(tr, d$species); d <- d[match(tr$tip.label, d$species), ]
rownames(d) <- d$species
cat("Species with tree + traits + assembly stats:", nrow(d), "\n")

lg <- function(x) { x <- as.numeric(x); x[x <= 0] <- NA; log10(x) }
d$lNUMT <- lg(d$NUMTs); d$lNUPT <- lg(d$NUPTs)
d$lN50  <- lg(d[[COL$n50]])
d$tech  <- factor(d[[COL$tech]])
cat("Technology:", paste(names(table(d$tech)), table(d$tech), sep="=", collapse="  "), "\n\n")

lam_of <- function(m) {
  l <- tryCatch(coef(m$modelStruct$corStruct, unconstrained = FALSE), error = function(e) NA)
  as.numeric(l)[1]
}

pgls <- function(f, dat, label) {
  sub <- dat[complete.cases(dat[, all.vars(f)]), ]
  t2  <- keep.tip(tr, sub$species)
  m <- try(gls(f, data = sub, correlation = corPagel(1, t2, form = ~species),
               method = "ML"), silent = TRUE)
  if (inherits(m, "try-error")) { cat(sprintf("%-46s FAILED\n", label)); return(NULL) }
  s <- summary(m)$tTable
  cat(sprintf("%-46s n=%3d  lambda=%.3f\n", label, nrow(sub), lam_of(m)))
  for (i in 2:nrow(s))
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
pgls(lNUMT ~ lN50 + tech, d, "NUMT ~ N50 + technology")
pgls(lNUPT ~ lN50 + tech, d, "NUPT ~ N50 + technology")
cat("\n")

cat("=== T2: does the lambda gap survive conditioning on quality? ============\n")
m0n <- pgls(lNUMT ~ 1,           d, "NUMT ~ 1                 (baseline)")
m0p <- pgls(lNUPT ~ 1,           d, "NUPT ~ 1                 (baseline)")
m1n <- pgls(lNUMT ~ lN50 + tech, d, "NUMT ~ N50 + technology  (adjusted)")
m1p <- pgls(lNUPT ~ lN50 + tech, d, "NUPT ~ N50 + technology  (adjusted)")
cat(sprintf("\n  NUMT lambda: %.3f -> %.3f   NUPT lambda: %.3f -> %.3f\n",
            lam_of(m0n), lam_of(m1n), lam_of(m0p), lam_of(m1p)))
cat(sprintf("  gap: %.3f -> %.3f\n", lam_of(m0n)-lam_of(m0p), lam_of(m1n)-lam_of(m1p)))
cat("  gap unchanged -> the difference is not an assembly artefact\n")
cat("  gap shrinks   -> quantify and report as a limitation\n\n")

if (!is.na(COL$org_frac) && COL$org_frac %in% names(d)) {
  cat("=== T2b: organelle-dominated scaffold fraction as covariate =============\n")
  d$orgf <- as.numeric(d[[COL$org_frac]])
  m2n <- pgls(lNUMT ~ lN50 + tech + orgf, d, "NUMT ~ N50 + tech + org.frac")
  m2p <- pgls(lNUPT ~ lN50 + tech + orgf, d, "NUPT ~ N50 + tech + org.frac")
  cat(sprintf("\n  lambda: NUMT %.3f   NUPT %.3f   gap %.3f\n\n",
              lam_of(m2n), lam_of(m2p), lam_of(m2n)-lam_of(m2p)))
}

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
