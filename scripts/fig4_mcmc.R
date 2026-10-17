# Fig 4 (per-phage binomial/threshold MCMCglmm with phylogenetic random effect), adapted from the authors' Rmd.
# Data prep is verbatim; only paths, iteration counts, and parallel/outputs differ.
#   Rscript scripts/fig4_mcmc.R <out_dir> <nitt> <thin> <burnin> <cores> [phage1,phage2,...]
ROOT <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value=TRUE))), ".."))
.libPaths(c(file.path(ROOT, ".Rlib"), .libPaths()))
args <- commandArgs(trailingOnly=TRUE)
out_dir <- args[1]; NITT <- as.integer(args[2]); THIN <- as.integer(args[3]); BURN <- as.integer(args[4]); CORES <- as.integer(args[5])
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
library(nlme)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(MCMCglmm)
library(tibble)
library(tidyverse)
library(r2glmm)
library(performance)

options(max.print=9999999)

bact_features = read.table(file.path(ROOT, "data/genomics/bacteria/picard_collection.csv"), sep=";", header=TRUE, na.strings = "")
phage_features = read.table(file.path(ROOT, "data/genomics/phages/guelin_collection.csv"), sep=";", header=TRUE)
interaction_matrix = read.table(file.path(ROOT, "data/interactions/interaction_matrix.csv"), sep=";", header=TRUE)
coli_receptors = read.table(file.path(ROOT, "data/genomics/bacteria/outer_membrane_proteins/blast_results_cured_clusters=99_wide.tsv"), sep="\t", header=TRUE)
coli_ds = read.table(file.path(ROOT, "data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv"), sep=";", header=TRUE)

interactions_long = interaction_matrix %>% pivot_longer(cols = -c("bacteria"), names_to="phage", values_to="interaction")

bact_features$bacteria = make.names(bact_features$bacteria)
interactions_long$bacteria = make.names(interactions_long$bacteria)
interactions_long$interaction = as.integer(interactions_long$interaction > 0)
rownames(bact_features) = bact_features$bacteria
rownames(phage_features) = phage_features$phage
rownames(interaction_matrix) = interaction_matrix$bacteria

# Prepare categorical data
bact_features[is.na(bact_features[, "Klebs_capsule_type"]),"Klebs_capsule_type"] = "NA"
bact_features[is.na(bact_features[, "ABC_serotype"]), "ABC_serotype"] = "NA"
bact_features$Clermont_Phylo = as.factor(bact_features$Clermont_Phylo)
bact_features$ST_Warwick = as.factor(bact_features$ST_Warwick)
bact_features$O.type = as.factor(bact_features$O.type)
bact_features$H.type = as.factor(bact_features$H.type)
bact_features$LPS_type = as.factor(bact_features$LPS_type)
bact_features$Klebs_capsule_type = as.factor(bact_features$Klebs_capsule_type)
bact_features$ABC_serotype = as.factor(bact_features$ABC_serotype)

# Prepare receptors and DS

coli_receptors$bacteria = make.names(coli_receptors$bacteria)  # FIX vs authors' Rmd: names were not make.names()'d, so merge() silently dropped 136/403 strains
coli_receptors[is.na(coli_receptors),] = "NA"
coli_receptors = lapply(coli_receptors, as.factor)
coli_ds[,-1] = coli_ds[,-1] > 0 # to boolean
coli_ds$bacteria = make.names(coli_ds$bacteria)
n_obs_sys = colSums(coli_ds[,-1])

filter_out_non_frequent_sys = TRUE
if (filter_out_non_frequent_sys) {
  print(dim(coli_ds))
  print("Filtering out non frequent defense systems.")
  thresh_sys = 15
  print(paste("Threshold :", thresh_sys))
  coli_ds = data.frame(coli_ds$bacteria, coli_ds[,-1][, n_obs_sys > thresh_sys])
  colnames(coli_ds)[1] = "bacteria"
  print(dim(coli_ds))
}
bact_features = merge(bact_features, coli_receptors, by="bacteria")
bact_features = merge(bact_features, coli_ds, by="bacteria")
bact_features$bacteria = bact_features$bacteria
phage_features$phage = make.names(phage_features$phage)
rownames(bact_features) = bact_features$bacteria

bact_features = bact_features[,c("bacteria", "O.type", "H.type", "LPS_type", "Klebs_capsule_type", "ABC_serotype", "BTUB", "FADL", "FHUA", "LAMB", "LPTD", "NFRA", "OMPA", "OMPC", "OMPF", "TOLC", "TSX", "YNCD", colnames(coli_ds[,-1]))]

# Load phylogenetic trees
use_ultrametric_tree = TRUE
if (use_ultrametric_tree) {
  bact_phylo_tree = ape::read.tree(file.path(ROOT, "data/genomics/bacteria/panacota/tree/370+host_ultrametric_tree_root=B992.nwk"))
  bact_phylo_tree$tip.label = make.names(bact_phylo_tree$tip.label)
  # Root tree and scale to ultrametric tree
  bact_phylo_tree = ape::root(bact_phylo_tree, "B992", resolve.root=TRUE)
  # Get inverse of phylogenetic covariance matrix
  bact_Ainv = inverseA(bact_phylo_tree, scale=TRUE)$Ainv
} else {
  bact_phylo_tree = ape::read.tree(file.path(ROOT, "data/genomics/bacteria/panacota/tree/370+host_tree_root=albertii.nwk"))
  bact_phylo_tree = ape::root(bact_phylo_tree, "B992", resolve.root=TRUE)
  bact_phylo_tree$edge.length = abs(bact_phylo_tree$edge.length)
  bact_phylo_tree$edge.length = bact_phylo_tree$edge.length + 5e-9
  bact_Ainv = inverseA(bact_phylo_tree, scale=FALSE)$Ainv
}
phage_phylo_tree = ape::read.tree(file.path(ROOT, "data/genomics/phages/tree/96_viridic_distance_phylogenetic_tree_algo=upgma.nwk"))
bact_phylo_tree$edge.length = bact_phylo_tree$edge.length + 5e-9
phage_phylo_tree$edge.length = phage_phylo_tree$edge.length + 5e-9
phage_phylo_tree$tip.label = make.names(phage_phylo_tree$tip.label)
phage_Ainv = inverseA(phage_phylo_tree, scale=FALSE)$Ainv

# Check that all tree nodes are in the bact_features dataframe
i = 0
for (name in bact_phylo_tree$tip.label) {
  if (!(name %in% rownames(bact_features))) {
    # tibble::add_row(bact_features, bacteria=name)
    print(name)
    bact_phylo_tree = ape::drop.tip(bact_phylo_tree, name)
    i = i + 1
  }
}
print("Removed tips:")
print(i)

i = 0
for (name in unique(interactions_long$bacteria)) {
  if (!(name %in% bact_phylo_tree$tip.label)) {
    interactions_long = subset(interactions_long, bacteria != name)
    i = i + 1
  }
}
print("Removed bacteria:")
print(i)


avoid_low_cardinal_levels = TRUE
if (avoid_low_cardinal_levels) {
  print("Lumping low cardinality levels for O.type and H.type")
  thresh = 15
  
  bact_features = bact_features %>% mutate(H.type=if_else(H.type == "H-", "Other", H.type), ABC_serotype=if_else(ABC_serotype %in% c("1", "3", "6", "Unknown", "NA"), "Other", ABC_serotype), Klebs_capsule_type=if_else(Klebs_capsule_type %in% c("K10", "K127", "K16", "K39", "K54", "K25", "K63", "K9", "NA"), "Other", Klebs_capsule_type))
  bact_features$H.type = as.factor(bact_features$H.type)
  bact_features$ABC_serotype = as.factor(bact_features$ABC_serotype)
  
  bact_features = bact_features %>% mutate(O.type=fct_lump_min(O.type, thresh, other_level="Other"), H.type=fct_lump_min(H.type, thresh, other_level="Other"), BTUB=fct_lump_min(BTUB, thresh, other_level="Other"), FADL=fct_lump_min(FADL, thresh, other_level="Other"), FHUA=fct_lump_min(FHUA, thresh, other_level="Other"), LAMB=fct_lump_min(LAMB, thresh, other_level="Other"), LPTD=fct_lump_min(LPTD, thresh, other_level="Other"), NFRA=fct_lump_min(NFRA, thresh, other_level="Other"), OMPA=fct_lump_min(OMPA, thresh, other_level="Other"), OMPC=fct_lump_min(OMPC, thresh, other_level="Other"), OMPF=fct_lump_min(OMPF, thresh, other_level="Other"), TOLC=fct_lump_min(TOLC, thresh, other_level="Other"), TSX=fct_lump_min(TSX, thresh, other_level="Other"), YNCD=fct_lump_min(YNCD, thresh, other_level="Other")) #ABC_serotype=fct_lump_min(ABC_serotype, thresh, other_level="Other"))  # Avoid overparametrization by only taking into account O.type and H.type with more than n observations
  print("Min size threshold: ")
  print(thresh)
}


library(parallel)
phages <- phage_features$phage
if (length(args) >= 6) phages <- make.names(strsplit(args[6], ",")[[1]])
fit_one <- function(phage_of_interest) {
  out <- file.path(out_dir, paste0(phage_of_interest, ".tsv"))
  if (file.exists(out)) return(invisible(NULL))
  d <- merge(subset(interactions_long, phage == phage_of_interest), bact_features, by="bacteria")
  Ainv <- inverseA(bact_phylo_tree, scale=FALSE)$Ainv
  prior <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(R1 = list(V = 1, nu = 0.002)))
  formula <- as.formula(paste("interaction ~ ", paste(colnames(d[,-(1:3)]), collapse=" + ")))
  set.seed(1)
  fit <- function(prior) MCMCglmm(formula, random = ~ bacteria, data = d, ginverse = list(bacteria=Ainv), family="threshold",
                                  scale=scale_tree, nitt=NITT, thin=THIN, burnin=BURN, prior=prior, verbose=FALSE)
  m0 <- tryCatch(fit(prior), error=function(e) e)
  if (inherits(m0, "error") && grepl("singular", conditionMessage(m0))) {
    # DEVIATION from the paper, only for phages that fail with "Mixed model equations singular": weakly informative N(0, 10^2) prior on the fixed effects
    p <- ncol(model.matrix(formula, d)); prior$B <- list(mu=rep(0, p), V=diag(p) * 100)
    m0 <- fit(prior)
    cat(phage_of_interest, "\n", file=file.path(out_dir, "STRONG_PRIOR_PHAGES.txt"), append=TRUE)
  }
  if (inherits(m0, "error")) stop(conditionMessage(m0))
  sols <- summary(m0)$solutions
  write.table(data.frame(feature=row.names(sols), sols, row.names=NULL), out, sep="\t", quote=FALSE, row.names=FALSE)
}
res <- mclapply(phages, function(p) tryCatch(fit_one(p), error=function(e) message(p, ": ", conditionMessage(e))), mc.cores=CORES, mc.preschedule=FALSE)
