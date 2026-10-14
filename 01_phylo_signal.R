#!/usr/bin/env Rscript
# 01_phylo_signal.R -- phylogenetic signal in insertion content
#
#   Rscript 01_phylo_signal.R            (PHYLOSIG_NSIM=10000 for finer K P-values)
#
# Traits (log10, on the same tree and tip set as the PGLS models):
#   NUMT proportion, NUPT proportion, NUPT:NUMT ratio, scaffold N50 (if the
#   assembly-stats table has an N50 column)
#
# For each trait:
#   Blomberg's K with a tip-randomisation test (phytools::phylosig)
#   Pagel's lambda with a likelihood-ratio test against lambda = 0 (phylosig)
#   and against lambda = 1 (Brownian motion; geiger::fitContinuous, lambda vs BM)
#
# Writes to $NUXT_OUT: phylo_signal_table.tsv
# Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 exported in the shell:
# fitContinuous forks, and forked R with multi-threaded BLAS can deadlock.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
source(file.path(here, "common.R"))

need <- c("phytools", "geiger")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({ library(phytools); library(geiger) })
# single-threaded linear algebra: forked workers + threaded BLAS can hang
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

NSIM <- as.integer(Sys.getenv("PHYLOSIG_NSIM", "1000"))
cat("Blomberg K randomisations per trait:", NSIM, "\n")

d  <- load_data()
df <- d$df; tr <- d$tree

traits <- list(
  "NUMT proportion"  = df$log_numt_prop,
  "NUPT proportion"  = df$log_nupt_prop,
  "NUPT:NUMT ratio"  = df$log_nupt_prop - df$log_numt_prop
)

# scaffold N50, if the assembly statistics table carries it
if (file.exists(STATS_F)) {
  st <- read.delim(STATS_F, stringsAsFactors = FALSE)
  n50_col <- grep("n50", names(st), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(n50_col)) {
    v <- as.numeric(st[[n50_col]][match(df$species, st$species)])
    traits[[paste0("log10 ", n50_col)]] <- ifelse(v > 0, log10(v), NA)
  } else cat("No N50 column in", STATS_F, "- N50 row skipped\n")
}

rows <- list()
for (nm in names(traits)) {
  x <- setNames(traits[[nm]], df$species)
  x <- x[is.finite(x)]
  t2 <- keep.tip(tr, names(x))
  x <- x[t2$tip.label]
  t0 <- Sys.time()
  k   <- phylosig(t2, x, method = "K", test = TRUE, nsim = NSIM)
  lam <- phylosig(t2, x, method = "lambda", test = TRUE)
  # ncores = 1: fitContinuous otherwise forks one process per optimisation
  # start (parallel::mclapply), which can deadlock with a multithreaded BLAS
  fl  <- fitContinuous(t2, x, model = "lambda", ncores = 1)
  fb  <- fitContinuous(t2, x, model = "BM", ncores = 1)
  lr1 <- 2 * (fl$opt$lnL - fb$opt$lnL)
  rows[[nm]] <- data.frame(
    trait = nm, n = length(x),
    K = k$K, K_P = k$P,
    lambda = lam$lambda, lambda_LR_vs0 = 2 * (lam$logL - lam$logL0), lambda_P_vs0 = lam$P,
    lambda_LR_vs1 = lr1, lambda_P_vs1 = pchisq(max(lr1, 0), df = 1, lower.tail = FALSE))
  cat(sprintf("%-22s n=%d  K=%.3f (P=%.3g)  lambda=%.3f (P vs 0=%.3g, vs 1=%.3g)  [%.0f s]\n",
              nm, length(x), k$K, k$P, lam$lambda, lam$P, rows[[nm]]$lambda_P_vs1,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  # write after every trait, so finished rows survive an interrupted run
  write.table(do.call(rbind, rows), file.path(OUT, "phylo_signal_table.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
tab <- do.call(rbind, rows)
write.table(tab, file.path(OUT, "phylo_signal_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nwrote", file.path(OUT, "phylo_signal_table.tsv"), "\n")
