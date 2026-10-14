# common.R -- shared by 01_phylo_signal.R, 02_pgls_models.R, 03_plot_pgls_panels.R
#
# Data: insertion content and organelle sizes per species from
# final_species_summary_v2.tsv (01_detection/build_species_summary.py); nuclear
# genome size from the genome FASTA (tree/assembly_stats.tsv,
# 00_data/make_assembly_stats.py) where available.
# Parasitic species with reduced plastomes are flagged (genus list below) for
# the plastid-size model without parasites.
#
# Fitting: lambda profiled on a grid (corPagel(fixed = FALSE) fails to
# converge for several models); report fit at ML lambda with a 95% profile
# CI from the 1.92 log-likelihood drop, and OLS (lambda = 0) alongside.

suppressPackageStartupMessages({
  library(ape); library(nlme); library(dplyr)
})

BASE <- Sys.getenv("NUXT_BASE", "/vol/storage/5T_disk/organelles")
OUT  <- Sys.getenv("NUXT_OUT",  file.path(BASE, "pgls_v2"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TREE_F  <- file.path(BASE, "tree/species_names_grafted.nwk")
SUMM_F  <- file.path(BASE, "final_species_summary_v2.tsv")
STATS_F <- file.path(BASE, "tree/assembly_stats.tsv")
PARASITE_GENERA <- c("Cuscuta", "Orobanche", "Phelipanche", "Aeginetia",
                     "Epifagus", "Striga", "Rafflesia", "Sapria")

GRID <- seq(0, 1, by = 0.02)

# Same as phytools::force.ultrametric(method = "extend"): lengthen each
# terminal branch so every tip reaches the height of the deepest tip.
# Written out here so the PGLS scripts need ape only, not phytools.
extend_ultrametric <- function(tr) {
  nt <- Ntip(tr)
  depth <- node.depth.edgelength(tr)[seq_len(nt)]
  term <- match(seq_len(nt), tr$edge[, 2])
  tr$edge.length[term] <- tr$edge.length[term] + (max(depth) - depth)
  tr
}

load_data <- function() {
  tree <- read.tree(TREE_F)
  df <- read.delim(SUMM_F, stringsAsFactors = FALSE)[
    , c("species", "genome_size", "numt_merged_prop", "nupt_merged_prop",
        "mito_length", "plastid_length")] %>%
    filter(!is.na(genome_size), !is.na(numt_merged_prop), !is.na(nupt_merged_prop)) %>%
    mutate(parasite = sub("_.*", "", species) %in% PARASITE_GENERA,
           # insertion bp are independent of which genome size is used
           numt_bp = numt_merged_prop * genome_size,
           nupt_bp = nupt_merged_prop * genome_size)

  # FASTA-derived genome size where available (fixes truncated Hordeum etc.)
  if (file.exists(STATS_F)) {
    st <- read.delim(STATS_F)[, c("species", "total_bp")]
    df <- left_join(df, st, by = "species") %>%
      mutate(genome_size = ifelse(!is.na(total_bp), total_bp, genome_size))
  }

  df <- df %>% mutate(
    numt_merged_prop = numt_bp / genome_size,
    nupt_merged_prop = nupt_bp / genome_size,
    log_genome  = log10(genome_size / 1e6),        # Mbp
    log_mito    = log10(mito_length / 1e3),        # kbp
    log_plastid = log10(plastid_length / 1e3),     # kbp
    log_numt_prop = log10(numt_merged_prop),
    log_nupt_prop = log10(nupt_merged_prop),
    log_numt_bp   = log10(numt_bp),
    log_nupt_bp   = log10(nupt_bp))

  overlap <- intersect(tree$tip.label, df$species)
  tr <- keep.tip(tree, overlap)
  tr$edge.length[tr$edge.length < 0] <- 0
  if (!is.ultrametric(tr)) tr <- extend_ultrametric(tr)
  df <- df[match(tr$tip.label, df$species), ]
  rownames(df) <- df$species
  cat(sprintf("species in tree and data: %d  (parasites %d)\n", nrow(df), sum(df$parasite)))
  list(df = df, tree = tr)
}

## fit one model across the lambda grid ----------------------------------------
fit_profile <- function(formula, df, tree, label, subset = NULL) {
  vars <- all.vars(formula)
  ok <- complete.cases(df[, vars]) & apply(is.finite(as.matrix(df[, vars])), 1, all)
  if (!is.null(subset)) ok <- ok & subset
  dat <- df[ok, c("species", vars)]
  t2  <- keep.tip(tree, dat$species)
  dat <- dat[match(t2$tip.label, dat$species), ]

  fits <- lapply(GRID, function(l) {
    ll <- min(max(l, 1e-6), 1 - 1e-6)
    m <- try(suppressWarnings(gls(formula, data = dat, method = "ML",
             correlation = corPagel(ll, phy = t2, form = ~species, fixed = TRUE))),
             silent = TRUE)
    if (inherits(m, "try-error")) NULL else m
  })
  keep <- !vapply(fits, is.null, logical(1))
  if (!any(keep)) { cat(label, ": all fits failed\n"); return(NULL) }
  lls  <- vapply(fits[keep], function(m) as.numeric(logLik(m)), numeric(1))
  lams <- GRID[keep]
  best <- fits[keep][[which.max(lls)]]
  lhat <- lams[which.max(lls)]
  inCI <- lls >= max(lls) - 1.92
  ols  <- summary(lm(formula, data = dat))$coefficients

  s <- summary(best)$tTable
  terms <- rownames(s)[-1]
  res <- do.call(rbind, lapply(terms, function(tm) data.frame(
    model = label, term = tm, n = nrow(dat),
    lambda = lhat, lambda_lo = min(lams[inCI]), lambda_hi = max(lams[inCI]),
    lambda_boundary = lhat %in% c(0, 1),
    slope = s[tm, "Value"], SE = s[tm, "Std.Error"], P = s[tm, "p-value"],
    # t-test against isometry (slope = 1); meaningful for bp ~ size models
    P_vs_1 = 2 * pt(-abs((s[tm, "Value"] - 1) / s[tm, "Std.Error"]),
                    df = nrow(dat) - length(terms) - 1),
    OLS_slope = ols[tm, "Estimate"], OLS_P = ols[tm, "Pr(>|t|)"],
    stringsAsFactors = FALSE)))

  # slope and P of every term at every lambda (Supplemental Table S4)
  sens <- do.call(rbind, lapply(seq_along(lams), function(i) {
    tt <- summary(fits[keep][[i]])$tTable
    data.frame(model = label, term = terms, lambda = lams[i], logLik = lls[i],
               slope = tt[terms, "Value"], P = tt[terms, "p-value"])
  }))
  res$robust_to_lambda <- vapply(terms, function(tm) {
    st <- sens[sens$term == tm, ]
    if (all(st$P < 0.05)) "significant at all lambda"
    else if (all(st$P >= 0.05)) "never significant"
    else sprintf("significant for lambda in [%.2f, %.2f]",
                 min(st$lambda[st$P < 0.05]), max(st$lambda[st$P < 0.05]))
  }, character(1))

  list(table = res, profile = sens, fit = best, data = dat, lambda = lhat)
}
