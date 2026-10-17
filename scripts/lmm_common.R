# Shared set-up of the Gaussian LMMs (Fig 4c-d, ED6): libraries, interaction matrix, per-strain targets, inverse relatedness matrix, prior.
# Sourced by fig4_cd.R and fig_ed6.R, which define ROOT first.
.libPaths(c(file.path(ROOT, ".Rlib"), .libPaths()))
suppressMessages({library(MCMCglmm); library(ape); library(parallel)})
set.seed(1)
M <- read.table(file.path(ROOT, "data/interactions/interaction_matrix.csv"), sep=";", header=TRUE, check.names=FALSE, row.names=1)
strain_targets <- function() data.frame(bacteria = rownames(M), n_lytic = rowSums(M > 0, na.rm=TRUE),
  avg_mlc = apply(M, 1, function(x) { x <- x[!is.na(x) & x > 0]; if (length(x)) mean(x) else NA }))
# d: data frame with a `bacteria` column. Returns the rows that are in the ultrametric tree, plus the inverse relatedness matrix.
prepare_lmm <- function(d) {
  d$animal <- make.names(d$bacteria)
  tree <- read.tree(file.path(ROOT, "data/genomics/bacteria/panacota/tree/370+host_ultrametric_tree_root=B992.nwk"))
  tree$tip.label <- make.names(tree$tip.label)
  tree <- drop.tip(tree, setdiff(tree$tip.label, d$animal))
  list(d = d[d$animal %in% tree$tip.label, ], Ainv = inverseA(tree, scale=TRUE)$Ainv)
}
prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))
