# SI Text 7: number of lytic interactions of a phage against the number of RBPs it encodes, Gaussian LMM (MCMCglmm), phage phylogeny as random effect
# (the VIRIDIC UPGMA tree of data/genomics/phages/tree), 1,000,000 iterations, burn-in 50,000, thin 250, prior V = 1, nu = 0.002 (settings of the Fig 4c-d model).
#   Rscript scripts/si_rbp_lmm.R  (after si_analyses.py)  ->  results/si/rbp_lmm.csv
ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value=TRUE))), ".."))
.libPaths(c(file.path(ROOT, ".Rlib"), .libPaths()))
suppressMessages({library(MCMCglmm); library(ape); library(parallel)})
set.seed(1)
d <- read.csv(file.path(ROOT, "results/si/rbp_count_vs_lytic.csv"), stringsAsFactors = FALSE)
tree <- read.tree(file.path(ROOT, "data/genomics/phages/tree/96_viridic_distance_phylogenetic_tree_algo=upgma.nwk"))
d$animal <- make.names(d$phage); tree$tip.label <- make.names(tree$tip.label)
stopifnot(all(d$animal %in% tree$tip.label))
tree <- drop.tip(tree, setdiff(tree$tip.label, d$animal)); tree$edge.length[tree$edge.length <= 0] <- 1e-8
if (!is.ultrametric(tree)) tree <- chronos(tree, quiet = TRUE) |> multi2di()
Ainv <- inverseA(tree, scale = TRUE)$Ainv
prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))
fit1 <- function(dd, label) {
  m <- MCMCglmm(n_lytic ~ n_rbp, random = ~ animal, data = dd, ginverse = list(animal = Ainv), family = "gaussian", prior = prior, nitt = 1000000, burnin = 50000, thin = 250, verbose = FALSE)
  s <- summary(m)$solutions["n_rbp", ]
  data.frame(set = label, n_phages = nrow(dd), beta = s["post.mean"], ci_low = s["l-95% CI"], ci_high = s["u-95% CI"], pMCMC = s["pMCMC"], row.names = NULL)
}
res <- do.call(rbind, mclapply(list(list(d, "all 96 phages (the 7 listed without RBP count as 0)"), list(d[d$n_rbp > 0, ], "89 phages with at least one RBP")), function(a) fit1(a[[1]], a[[2]]), mc.cores = 2))
res$paper_beta <- 28.11; res$paper_ci <- "-33.09 to 83.58"; res$paper_pMCMC <- 0.33
write.csv(res, file.path(ROOT, "results/si/rbp_lmm.csv"), row.names = FALSE)
print(res)
