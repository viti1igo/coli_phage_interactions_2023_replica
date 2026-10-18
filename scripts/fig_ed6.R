# ED6: association of bacterial genome size (A, C) and number of AMR genes (B, D) with the number of lytic interactions (A, B) and the average MLC (C, D).
# Gaussian LMM as in fig4_cd.R. Genome size is from our assemblies; the AMR gene counts are the paper's (Source Data ED6, written by compare_source_data.py
# to data/derived/paper_amr_counts.csv) because no AMR annotation is shipped. A and C are repeated with the paper's genome sizes as a check.
#   Rscript scripts/fig_ed6.R  ->  results/fig4_summary/ed6.csv, ed6_data.csv
ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value=TRUE))), ".."))
source(file.path(ROOT, "scripts/lmm_common.R"))
gs <- read.csv(file.path(ROOT, "data/derived/genome_size.csv"), row.names=1)
pa <- read.csv(file.path(ROOT, "data/derived/paper_amr_counts.csv"), row.names=1, colClasses=c(bacteria="character"))
d <- strain_targets()
d$mb <- gs[rownames(M), "genome_size"] / 1e6
d$mb_paper <- pa[rownames(M), "genome_size_paper"] / 1e6
d$amr <- pa[rownames(M), "n_amr_genes"]
lmm <- prepare_lmm(d); d <- lmm$d; Ainv <- lmm$Ainv
jobs <- data.frame(panel = c("A", "B", "C", "D", "A_paper_size", "C_paper_size"),
                   y = c("n_lytic", "n_lytic", "avg_mlc", "avg_mlc", "n_lytic", "avg_mlc"),
                   x = c("mb", "amr", "mb", "amr", "mb_paper", "mb_paper"), stringsAsFactors = FALSE)
fit1 <- function(i) {
  j <- jobs[i, ]; dd <- d[!is.na(d[[j$y]]) & !is.na(d[[j$x]]), ]; dd$y <- dd[[j$y]]; dd$x <- dd[[j$x]]
  m <- MCMCglmm(y ~ x, random = ~ animal, data = dd, ginverse = list(animal = Ainv), family = "gaussian", prior = prior,
                nitt = 1000000, burnin = 50000, thin = 250, verbose = FALSE)
  s <- summary(m)$solutions["x", ]
  data.frame(panel = j$panel, target = j$y, covariate = j$x, n_strains = nrow(dd), beta = s["post.mean"], ci_low = s["l-95% CI"], ci_high = s["u-95% CI"], pMCMC = s["pMCMC"], ESS = s["eff.samp"], row.names = NULL)
}
res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), fit1, mc.cores = 3))
write.csv(res, file.path(ROOT, paste0("results/fig4_summary/ed6", SUFFIX, ".csv")), row.names=FALSE)
write.csv(d[, c("bacteria", "n_lytic", "avg_mlc", "mb", "mb_paper", "amr")], file.path(ROOT, paste0("results/fig4_summary/ed6_data", SUFFIX, ".csv")), row.names=FALSE)
print(res)
