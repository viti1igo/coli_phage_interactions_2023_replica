# ED6 A and C: association of bacterial genome size (per Mb, from our assemblies) with the number of lytic interactions (A) and the average MLC (C).
# Same Gaussian LMM as fig4_cd.R. AMR-gene panels (B, D) need an AMR annotation that is not shipped and was not run.
#   Rscript scripts/fig_ed6.R  ->  results/fig4_summary/ed6.csv
ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value=TRUE))), ".."))
source(file.path(ROOT, "scripts/lmm_common.R"))
gs <- read.csv(file.path(ROOT, "data/derived/genome_size.csv"), row.names=1)
d <- strain_targets(); d$mb <- gs[rownames(M), "genome_size"] / 1e6
lmm <- prepare_lmm(d); d <- lmm$d[!is.na(lmm$d$mb), ]; Ainv <- lmm$Ainv
fit1 <- function(y) {
  dd <- d[!is.na(d[[y]]), ]; dd$y <- dd[[y]]
  m <- MCMCglmm(y ~ mb, random = ~ animal, data = dd, ginverse = list(animal = Ainv), family = "gaussian", prior = prior,
                nitt = 1000000, burnin = 50000, thin = 250, verbose = FALSE)
  s <- summary(m)$solutions["mb", ]
  data.frame(panel = if (y == "n_lytic") "A" else "C", target = y, n_strains = nrow(dd), beta_per_Mb = s["post.mean"], ci_low = s["l-95% CI"], ci_high = s["u-95% CI"], pMCMC = s["pMCMC"], ESS = s["eff.samp"], row.names = NULL)
}
res <- do.call(rbind, mclapply(c("n_lytic", "avg_mlc"), fit1, mc.cores = 2))
write.csv(res, file.path(ROOT, "results/fig4_summary/ed6.csv"), row.names=FALSE)
write.csv(d[, c("bacteria", "n_lytic", "avg_mlc", "mb")], file.path(ROOT, "results/fig4_summary/ed6_data.csv"), row.names=FALSE)
print(res)
