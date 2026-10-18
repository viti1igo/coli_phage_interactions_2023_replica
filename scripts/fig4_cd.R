# Fig 4C/D: association of the number of defence systems with (C) the number of lytic interactions of a strain and (D) its average MLC upon infection.
# Gaussian LMM (MCMCglmm), phylogeny as a random effect through the inverse relatedness matrix, 1,000,000 iterations, burn-in 50,000, thin 250,
# inverse-Gamma prior (V=1, nu=0.002), as in the paper's Methods. Run with the shipped defence counts and with ours (T22).
#   Rscript scripts/fig4_cd.R      ->  results/fig4_summary/fig4cd.csv
ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value=TRUE))), ".."))
source(file.path(ROOT, "scripts/lmm_common.R"))
pic <- read.table(file.path(ROOT, "data/genomics/bacteria/picard_collection.csv"), sep=";", header=TRUE, row.names=1, na.strings="", check.names=FALSE)
reg <- read.table(file.path(ROOT, "replication/T25_features/picard_collection_regenerated.csv"), sep=";", header=TRUE, row.names=1, na.strings="", check.names=FALSE)
d <- strain_targets()
d$ds_shipped <- pic[rownames(M), "n_defense_systems"]; d$ds_ours <- reg[rownames(M), "n_defense_systems"]
# the paper's Source Data for panel d uses a smaller defense count than panel c (written by compare_source_data.py when the supplement is downloaded)
f4d <- file.path(ROOT, "data/derived/paper_fig4d_defense_counts.csv")
d$ds_paper4d <- if (file.exists(f4d)) read.csv(f4d, row.names=1, colClasses=c(bacteria="character"))[rownames(M), "n_defense_4d"] else NA
lmm <- prepare_lmm(d); d <- lmm$d; Ainv <- lmm$Ainv
jobs <- expand.grid(model = c("C: n lytic interactions", "D: average MLC upon infection"), ds = c("ds_shipped", "ds_ours", "ds_paper4d"), stringsAsFactors = FALSE)
jobs <- jobs[!(jobs$ds == "ds_paper4d" & (startsWith(jobs$model, "C") | all(is.na(d$ds_paper4d)))), ]
fit1 <- function(i) {
  j <- jobs[i, ]; y <- if (startsWith(j$model, "C")) "n_lytic" else "avg_mlc"
  dd <- d[!is.na(d[[y]]) & !is.na(d[[j$ds]]), ]; dd$x <- dd[[j$ds]]; dd$y <- dd[[y]]
  m <- MCMCglmm(y ~ x, random = ~ animal, data = dd, ginverse = list(animal = Ainv), family = "gaussian", prior = prior,
                nitt = 1000000, burnin = 50000, thin = 250, verbose = FALSE)
  s <- summary(m)$solutions["x", ]
  data.frame(model = j$model, defence_count = j$ds, n_strains = nrow(dd), beta = s["post.mean"], ci_low = s["l-95% CI"], ci_high = s["u-95% CI"], pMCMC = s["pMCMC"], ESS = s["eff.samp"], row.names = NULL)
}
res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), fit1, mc.cores = 4))
res$paper_beta <- ifelse(startsWith(res$model, "C"), -0.22, -0.04)
res$paper_ci <- ifelse(startsWith(res$model, "C"), "-0.74 to 0.34", "-0.05 to -0.01")
res$paper_pMCMC <- ifelse(startsWith(res$model, "C"), 0.41, 0.004)
dir.create(file.path(ROOT, "results/fig4_summary"), showWarnings=FALSE, recursive=TRUE)
write.csv(res, file.path(ROOT, "results/fig4_summary/fig4cd.csv"), row.names=FALSE)
print(res)
