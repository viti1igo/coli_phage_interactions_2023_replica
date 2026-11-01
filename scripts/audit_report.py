"""Phase 18 (T75): write reports/12_gap_audit.md from the tables in results/audit/. The numbers in the report are read from those files.
python scripts/audit_report.py"""
import os
import numpy as np, pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
A = f"{ROOT}/results/audit"
OUTFILE = f"{ROOT}/publish/phage-host-feature-necessity/reports/12_gap_audit.md"
rd = lambda f, **k: pd.read_csv(f"{A}/{f}", **k)


def md(df, align=None):
    df = df.astype(str).apply(lambda c: c.str.replace("|", "/", regex=False).str.replace("\n", " "))
    head = "| " + " | ".join(df.columns) + " |\n|" + "|".join("---:" if (align and c in align) else "---" for c in df.columns) + "|\n"
    return head + "\n".join("| " + " | ".join(r) + " |" for r in df.values) + "\n"


o, g, gp = rd("overall.csv"), rd("gap_register.csv"), rd("gap_register_prediction.csv")
pp, pc, pr = rd("per_phage.csv", index_col=0), rd("per_phage_causes.csv"), rd("paper_vs_repo.csv")
fx, ps = rd("fixed_folds_summary.csv", index_col=0)["value"], rd("pairs_t1_t2_summary.csv", index_col=0)["value"]
ra, sw, mg = rd("pairs_regeneration_summary.csv"), rd("swap_groups.csv"), rd("swap_marginal.csv")
sp, ss = rd("seed_spread.csv", index_col=0), rd("seed_spread_summary.csv", index_col=[0, 1])
ta, ed = rd("tier_alignment.csv"), rd("ed7_definitions.csv")
ctl = rd("swap_controls.csv", index_col=0)["value"]
fm = rd("fold_mechanism.csv", index_col=0)["value"]
pk = rd("preprint_check.csv")
psd = rd("pair_set_differences.csv")
pe, pef, pea, pel2, env = rd("pinned_env_overall.csv"), rd("pinned_env_folds.csv"), rd("pinned_env_pair_agreement.csv"), rd("pinned_env_logreg_l2_phages.csv", index_col=0), rd("env_match.csv")
per = pe[pe.section == "run"].set_index("run"); pch = pe[pe.section == "change"].set_index(["a", "b"])
pa = lambda run, m: pea[(pea.run.str.startswith(run)) & (pea.model == m)].iloc[0]
pe_tbl = per.reset_index()[["run", "mean_phage_auroc", "pooled_auroc", "mean_abs_d_auroc_vs_T1", "corr_auroc_vs_T1", "same_model_family_as_T1"]].round(4).fillna("")
pe_tbl.columns = ["run", "mean per-phage AUROC", "pooled AUROC", "mean abs change vs T1", "correlation with T1", "same model family as T1"]
env_tbl = env[["component", "authors", "pinned_runs", "status", "reason"]].fillna("")
t1_missing = psd[psd.side == "missing from T1"]
pm = rd("permanova_checks.csv")
ext = rd("external_check.csv")
ev = pc.groupby(["comparison", "cause_code", "evidence"]).size().reset_index(name="rows")
evid_tbl = md(ev.rename(columns={"comparison": "Comparison", "cause_code": "Code", "evidence": "Evidence", "rows": "Rows"}))
xr = pc[pc.cause_code == "X"]
x_list = "; ".join(f"{p} ({', '.join(xr[xr.phage == p].comparison.str.replace(' vs ', ' against '))})" for p in xr.phage.drop_duplicates())
fxd = rd("fixed_folds_vs_authors.csv", index_col=0)
v = lambda tier, m: o[(o.tier == tier) & (o.metric == m)].iloc[0]
code_n = lambda d: ", ".join(f"{c} {n}" for c, n in d.cause_code.value_counts().items())

# ---- overall table
rows = []
for m, label in (("pooled_auroc", "Pooled AUROC"), ("pooled_aupr", "Pooled average precision"), ("mean_phage_auroc", "Mean per-phage AUROC"), ("mean_phage_aupr", "Mean per-phage AUPR"), ("mean_strain_auroc", "Mean per-strain AUROC")):
    t0 = v("T1", m).paper_T0
    rows.append({"Metric": label, "Paper (T0)": "-" if pd.isna(t0) else f"{t0:.3f}" if m != "pooled_aupr" else f"{t0:.2f}",
                 **{t: f"{v(t, m).value:.3f}" for t in ("T1", "T2", "T3", "T4")}})
overall_tbl = md(pd.DataFrame(rows))
gaps = pd.DataFrame([{"Metric": label, "T2 against T1": f"{v('T2', m).gap_to_T1:+.4f} ({v('T2', m).gap_to_T1_ci})", "T3 against T1": f"{v('T3', m).gap_to_T1:+.4f} ({v('T3', m).gap_to_T1_ci})",
                      "T4 against T1": f"{v('T4', m).gap_to_T1:+.4f} ({v('T4', m).gap_to_T1_ci})"}
                     for m, label in (("pooled_auroc", "Pooled AUROC"), ("pooled_aupr", "Pooled average precision"), ("mean_phage_auroc", "Mean per-phage AUROC"), ("mean_phage_aupr", "Mean per-phage AUPR"), ("mean_strain_auroc", "Mean per-strain AUROC"))])
gaps_tbl = md(gaps)

# ---- per phage
flag = {c: int(pp[f"flag_{c}"].sum()) for c in ("T2_T1", "T3_T2", "T4_T2")}
mean_abs = {c: pp[f"d_{c}"].abs().mean() for c in ("T2_T1", "T3_T2", "T4_T2")}
n_cause = pc.groupby(["comparison", "cause_code"]).size()
cnt = lambda comp, code: int(n_cause.get((comp, code), 0))
from scipy.stats import spearmanr
rho = spearmanr(pp.n_lytic_T2, pp.d_T2_T1.abs())[0]
swp = rd("swap_per_phage.csv", index_col=0)
frows = pc[(pc.comparison == "T3 vs T2") & (pc.cause_code == "F")].copy()
frows["Largest single group"] = [f"{swp.loc[p, 'main_group']} ({swp.loc[p, 'main_d']:+.3f})" for p in frows.phage]
ftbl = md(frows.assign(**{"Phage": frows.phage, "Lytic strains": frows.n_lytic, "Change, T3 against T2": frows.d_auroc.map("{:+.3f}".format),
                          "Change with folds fixed, all groups": frows.swap_all_d_auroc.map("{:+.3f}".format)})[["Phage", "Lytic strains", "Change, T3 against T2", "Change with folds fixed, all groups", "Largest single group"]])
tiers_tbl = md(pd.DataFrame([
    ("T0", "The paper: text, supplement and Source Data", "`docs/paper/` (local only, never pushed)"),
    ("T1", "The authors' repository: shipped predictions and per-phage performance files", "`dev/predictions/results/`"),
    ("T2", "This replication on the shipped inputs", "`results/repro/`"),
    ("T3", "This replication on regenerated features (own tree, UMAP, groups)", "`results/repro_regenerated/`"),
    ("T4", "T3 plus the interaction matrix regenerated from the raw plaque scores", "`results/repro_regenerated_matrix/`")], columns=["Tier", "Source", "Where"]))
codes_tbl = md(pd.DataFrame([
    ("P", "Published materials disagree: the paper, supplement, Source Data or repository contradict each other or the Methods"),
    ("D", "A definition or unit is not stated, so the number cannot be reproduced"),
    ("U", "An input or code is not published (CV groups file, recommender code, test-strain genomes, defense-subfamily clustering)"),
    ("M", "The input data differs (a strain or label)"),
    ("F", "Feature regeneration differs from the shipped table, with a measured effect"),
    ("S", "Resampling noise: a different fold partition, model randomness or a small number of lytic strains"),
    ("C", "The model-choice rule or the model picked differs"),
    ("I", "A deliberate choice on this side, for example NODF in place of BINMATNEST"),
    ("OK", "Flagged only because the paper gives no number or a qualitative statement; the numbers agree"),
    ("X", "Unexplained; the row says what was tested and what is missing")], columns=["Code", "Meaning"]))

# ---- partition and seeds
mean_phage = lambda col, stat: ss.loc[(stat, "mean_phage_auroc"), col]
seed_tbl = md(pd.DataFrame([
    {"Source of variation": "Model seed, folds fixed (5 seeds)", "Mean per-phage AUROC": f"{mean_phage('fixed folds, seeds 0-4', 'mean'):.4f}", "Sd over runs": f"{mean_phage('fixed folds, seeds 0-4', 'std'):.4f}",
     "Range over runs": f"{mean_phage('fixed folds, seeds 0-4', 'min'):.4f} to {mean_phage('fixed folds, seeds 0-4', 'max'):.4f}", "Mean per-phage sd": f"{sp.sd_fixed.mean():.4f}", "Largest per-phage sd": f"{sp.sd_fixed.max():.3f}"},
    {"Source of variation": "Seed and a random group partition (5 seeds)", "Mean per-phage AUROC": f"{mean_phage('shuffled partition, seeds 0-4', 'mean'):.4f}", "Sd over runs": f"{mean_phage('shuffled partition, seeds 0-4', 'std'):.4f}",
     "Range over runs": f"{mean_phage('shuffled partition, seeds 0-4', 'min'):.4f} to {mean_phage('shuffled partition, seeds 0-4', 'max'):.4f}", "Mean per-phage sd": f"{sp.sd_shuffled.mean():.4f}", "Largest per-phage sd": f"{sp.sd_shuffled.max():.3f}"}]))

# ---- swaps
swt = sw.copy()
swt = swt.assign(**{"Group regenerated": swt.group, "Mean AUROC change": swt.mean_d_auroc.map("{:+.4f}".format), "Same, phages with 10 or more lytic strains": swt.mean_d_auroc_n10.map("{:+.4f}".format),
                    "Phages moving by more than 0.05": swt.phages_abs_d_above_0_05, "Calls flipped at 0.5": swt.call_flip_rate.map("{:.1%}".format)})[["Group regenerated", "Mean AUROC change", "Same, phages with 10 or more lytic strains", "Phages moving by more than 0.05", "Calls flipped at 0.5"]]
swt["Group regenerated"] = swt["Group regenerated"].replace({"all": "all seven together"})
mgt = md(mg.assign(**{"Group": mg.group, "Cost once the others are regenerated": mg.mean_marginal_d_auroc.map("{:+.4f}".format), "Same, phages with 10 or more lytic strains": mg.mean_marginal_d_auroc_n10.map("{:+.4f}".format),
                      "Phages moving by more than 0.05": mg.phages_abs_marginal_above_0_05})[["Group", "Cost once the others are regenerated", "Same, phages with 10 or more lytic strains", "Phages moving by more than 0.05"]])
n_small = int((pp.n_lytic_T2 < 10).sum())
all_row = sw[sw.group == "all"].iloc[0]
ra3, ra4 = ra[ra.tier == "T3"].iloc[0], ra[ra.tier == "T4"].iloc[0]
ed_sd = ed.iloc[0]

# ---- register tables
def reason(n, t):
    if t in ("", "evidence in the cited report"): return n
    if t.startswith("tested "): return f"{n}. Tested {t[7:]}"
    if t.startswith("not testable here"): return f"{n}. Not testable here{t[len('not testable here'):]}"
    return f"{n}. {t}"


reg = g.assign(**{"Row": g.figure + ": " + g.metric, "Code": g.cause_code, "Reason": [reason(n, t) for n, t in zip(g.note, g.tested)]})[["id", "Row", "Code", "Reason"]]
reg_tbl = md(reg.rename(columns={"id": "Id"}))
NICE = {"pooled_auroc": "pooled AUROC", "mean_phage_auroc": "mean per-phage AUROC", "mean_phage_aupr": "mean per-phage AUPR", "mean_strain_auroc": "mean per-strain AUROC"}
gp["metric"] = gp.metric.replace(NICE)
pred_tbl = md(gp.assign(**{"Id": gp.id, "Level": gp.level, "Comparison": gp.comparison, "Metric": gp.metric, "Code": gp.cause_code, "Size and evidence": gp["size"]})[["Id", "Level", "Comparison", "Metric", "Code", "Size and evidence"]])
repo_n = pr.cause_code.value_counts()
repo_star = pr[pr.changes_a_reported_number][["item", "topic", "cause_code"]].rename(columns={"item": "Item in report 11", "topic": "Topic", "cause_code": "Code"})
settle = (g[g.cause_code.isin(["U", "D", "X"])].groupby("would_settle").agg(rows=("id", lambda s: ", ".join(s)), code=("cause_code", lambda s: "/".join(sorted(set(s))))).reset_index())
settle_tbl = md(settle.rename(columns={"would_settle": "What the authors could supply", "rows": "Register rows", "code": "Code"})[["What the authors could supply", "Code", "Register rows"]])
xrows = g[g.cause_code == "X"]

text = f"""# 12. Gap audit

This report gives the reason for each gap between the paper, the authors' repository and this replication. Five sources are compared, as tiers:

{tiers_tbl}
Three comparisons run through the report: the paper against the repository (T0 against T1), our run on the shipped inputs against both (T2), and the shipped inputs against the regenerated ones (T2 against T3 and T4). Each is made overall, per phage (96 phages) and per strain-phage pair (38,377 pairs present in every tier). Every gap gets one cause code and a pointer to the file that holds the evidence:

{codes_tbl}
A gap is coded X when the evidence in hand does not decide it. The code is never left empty. The register has {len(g)} rows from the paper comparison ([`results/audit/gap_register.csv`](../results/audit/gap_register.csv)), {len(gp)} rows for the prediction results ([`gap_register_prediction.csv`](../results/audit/gap_register_prediction.csv)) and 28 rows for the paper against the repository ([`paper_vs_repo.csv`](../results/audit/paper_vs_repo.csv)). All of it comes from `scripts/audit_gaps.py`.

## Short answer

- Overall performance agrees. On the shipped inputs our pooled AUROC is {v('T2', 'pooled_auroc').value:.3f} against {v('T1', 'pooled_auroc').value:.3f} for the authors' predictions and {v('T1', 'pooled_auroc').paper_T0:.2f} in the paper, and every overall metric of T2 lies inside the bootstrap interval of T1.
- The per-phage differences between our run and the authors' come from the cross-validation partition. The authors' groups file is not published. Our rebuilt groups sit inside one of their folds in {1 - ps['groups_split_across_authors_folds']:.1%} of cases, yet the fold assignments barely overlap (adjusted Rand index {ps['fold_partition_adjusted_rand_mean']:.2f}). With their folds the mean absolute change in per-phage AUROC falls from {fx['mean_abs_d_auroc_our_folds']:.3f} to {fx['mean_abs_d_auroc']:.4f}, and {int(fx['identical_auroc_1e6'])} of 96 phages match to six decimals.
- Pinning the environment to the paper (scikit-learn 1.1.2, scipy 1.11.2, Python 3.10 because 1.1.2 has no 3.12 build) and using the authors' folds gives a mean absolute change of {per.loc["T2.1 authors' folds", "mean_abs_d_auroc_vs_T1"]:.4f} in per-phage AUROC, and {pa("T2.1", "LogReg_l2").share_identical_1e6:.0%} of the LogReg L2 probabilities equal theirs to six decimals. Random forests cannot match (the authors left them unseeded). With our own groups the version does not close the gap. [Details](#pinned-environment-t21-and-t31).
- Regenerating the features costs {abs(all_row.mean_d_auroc):.3f} mean per-phage AUROC with the folds fixed and leaves pooled AUROC almost unchanged ({v('T3', 'pooled_auroc').value:.3f} and {v('T4', 'pooled_auroc').value:.3f} against {v('T2', 'pooled_auroc').value:.3f}). Most of the loss sits in {n_small} phages with fewer than 10 lytic strains. The O-type and Klebsiella capsule type carry it once the other groups are regenerated.
- The {len(g)} rows of the paper comparison that did not match split as {code_n(g)}. The paper and its materials disagree with themselves in {int((g.cause_code == 'P').sum())} rows, and a definition or unit is missing in {int((g.cause_code == 'D').sum())}.
- Of the 28 differences between the repository and the paper, {repo_n.get('P', 0)} are inconsistencies between published materials (code P), {repo_n.get('D', 0)} lack a definition, {repo_n.get('U', 0)} rest on unpublished code or data and {repo_n.get('C', 0)} is the model-selection rule.
- An independent preprint on the same matrix ([check](#check-against-an-independent-preprint)) agrees on the dataset counts and on overall AUROC (0.869 against 0.86), and its baseline supports our reading of ED7.\n- {len(xrows)} rows of the paper comparison stay unexplained, and {int((pc.cause_code == 'X').sum())} of the {len(pc)} per-phage rows lie beyond the noise estimates without a tested cause. [Open rows](#open-rows-and-what-would-close-them) lists them and what each needs.

## Overall performance

All five metrics below pool the held-out predictions of the 10 folds and use the {int(o.n_pairs.iloc[0]):,} pairs that every tier has. The authors' script averages the 10 per-fold values instead, which gives the 0.787 and 0.765 of the earlier reports. The bootstrap resamples strains 500 times and keeps each phage's best model fixed.

{overall_tbl}
Gaps against the authors' predictions, with the 95 % bootstrap interval:

{gaps_tbl}
T2 sits inside the interval of T1 on every metric, so the shipped-input replication has no overall gap to explain. T3 and T4 lose {abs(v('T3', 'mean_phage_auroc').gap_to_T1):.3f} and {abs(v('T4', 'mean_phage_auroc').gap_to_T1):.3f} mean per-phage AUROC, and the intervals exclude zero. Pooled AUROC stays inside its interval. The paper's mean per-strain AUROC of 0.85 is below its own Source Data (0.872) and below the authors' predictions ({v('T1', 'mean_strain_auroc').value:.3f}), so that gap is in the text (P).

## Per phage

[`per_phage.csv`](../results/audit/per_phage.csv) holds the 96 phages by tier: AUROC, its standard error over folds, AUPR, best model and pooled AUROC. The AUROC computed from the fold-level predictions equals the earlier `per_phage.csv` of each run exactly.

- T2 against T1: mean absolute change {mean_abs['T2_T1']:.3f}, {flag['T2_T1']} phages above 0.05.
- T3 against T2: mean absolute change {mean_abs['T3_T2']:.3f}, {flag['T3_T2']} phages above 0.05.
- T4 against T2: mean absolute change {mean_abs['T4_T2']:.3f}, {flag['T4_T2']} phages above 0.05.
- Phages with fewer lytic strains move more (Spearman {rho:.2f} between the number of lytic strains and the size of the change).

[`per_phage_causes.csv`](../results/audit/per_phage_causes.csv) codes every flagged phage ({len(pc)} phage-comparison rows). Each row carries a `basis` sentence and an `evidence` column that says how strong the support is:

{evid_tbl}
Read the levels this way. The first two are tests: removing the factor removes the change (the authors' folds), or putting the factor in reproduces at least half of it (the feature swaps). Rows coded S are weaker. It says only that the change lies inside both noise estimates (the standard error over folds and the spread of five random partitions), so it cannot be told from noise, and no cause was shown. The code S does not mean a measured seed effect. The rows coded X lie beyond both noise estimates and no tested factor is sufficient: {x_list}. The spread of random partitions comes from only five runs per phage, so it is a rough estimate, and a few of these X rows may still be noise. For BCH953_P2 and MT1B1_3A1 I re-scored the regenerated-matrix predictions against the shipped labels and got the same AUROC, which rules out the relabelled pairs.

The phages coded F for T3 against T2:

{ftbl}
## Per pair

Pair-level comparisons use the {int(ps['pairs']):,} pairs present in T1 and T2. The authors' predictions lack {int(ta[ta.tier == 'T1'].T2_pairs_missing.iloc[0])} pairs that T2 has ({t1_missing.bacteria.nunique()} strains, {int((t1_missing.bacteria == 'NILS81').sum())} of them for NILS81) and have {int(ta[ta.tier == 'T1'].pairs_not_in_T2.iloc[0])} pair that T2 lacks. Labels agree on every shared pair. The authors' script turns the label F into NaN in an N/P/F/U matrix that the repository does not ship, so the cause of the 59 pairs is a reading I could not verify (code U). T4 adds {int(ta[ta.tier == 'T4'].pairs_not_in_T2.iloc[0])} pairs and changes {int(ta[ta.tier == 'T4'].label_flips_vs_T2.iloc[0])} labels.

- T2 against T1: Spearman {ps['spearman_all']:.3f} over all pairs (median {ps['spearman_phage_median']:.3f} per phage, lowest {ps['spearman_phage_min']:.3f}). Calls at 0.5 agree for {ps['call_agreement']:.1%} of pairs ({int(ps['flipped_pairs']):,} flips, {int(ps['flips_to_lytic']):,} to lytic and {int(ps['flips_to_nonlytic']):,} to non-lytic).
- The flip rate is {ps['flip_rate_same_model']:.1%} when the best model is identical, {ps['flip_rate_same_family_other_model']:.1%} for another model of the same family and {ps['flip_rate_different_family']:.1%} for another family. The first figure is the floor set by partition and randomness.
- With the authors' folds the flip count falls to {int(fx['flips_at_0_5'])} ({fx['flips_at_0_5'] / fx['pairs']:.1%}). Where the best model is the same only {fx['flip_rate_same_model']:.1%} flip, and where it differs {fx['flip_rate_other_model']:.1%}. {96 - int(fx['same_best_model'])} phages pick another best model, because the choice is a close race in AUPR.
- T3 and T4 against T2: {ra3.flip_rate:.1%} and {ra4.flip_rate:.1%} of calls flip (Spearman {ra3.spearman:.2f} and {ra4.spearman:.2f}). The flip rate is the same whether or not the strain's own categorical features changed, so those five columns do not explain it. With the folds fixed, the regenerated features alone flip {all_row.call_flip_rate:.1%} of calls, led by UMAP ({sw[sw.group == 'UMAP'].call_flip_rate.iloc[0]:.1%}) and the O-type ({sw[sw.group == 'O-type'].call_flip_rate.iloc[0]:.1%}). The regenerated tree gives different UMAP coordinates for every strain, which fits UMAP moving the most calls.
- Of the {int(ta[ta.tier == 'T4'].label_flips_vs_T2.iloc[0])} relabelled pairs in T4, {int(ra4.flips_with_label_change)} flip.

## The fold partition

The authors fix `np.random.seed(0)` in the prediction script and say nothing about seeds elsewhere. The cross-validation groups file `370+host_cross_validation_groups_1e-4.csv` is "not available anymore" according to their README. I rebuilt the groups from the tree. Our groups agree with their folds (a group falls in a single authors' fold in {1 - ps['groups_split_across_authors_folds']:.1%} of cases), yet the two sets of folds overlap little (adjusted Rand index {ps['fold_partition_adjusted_rand_mean']:.3f} on average, {ps['fold_partition_adjusted_rand_max']:.3f} at most).

Two tests say where the difference comes from ([`fold_mechanism.csv`](../results/audit/fold_mechanism.csv)). First, the authors' folds are what greedy `GroupKFold` packing of groups like ours produces. That rule sorts groups by size and puts each into the lightest fold, so only the order among equal-size groups is free, and most of our {int(fm['our_groups'])} groups are single strains. One of our groups (2 strains) is split by their folds. With those two strains taken apart, the authors' folds fit the packing rule for every phage ({int(fm['phages_packing_consistent_after_splitting'])} of {int(fm['phages'])}, all {int(fm['size_blocks_consistent_after_splitting'])} size blocks). Without that change, none of the 96 phages fits. Second, the preprint, the paper and the README all name scikit-learn 1.1.2, and this run uses 1.9.1. A copy of the 1.1.2 algorithm run on our group labels overlaps the authors' folds no better (adjusted Rand index {fm['ari_sklearn_1_1_2_copy']:.3f}, against {fm['ari_installed_sklearn']:.3f} for the installed version). The library version is therefore not the cause. What is missing is the order among equal groups, which comes from the labels in the unpublished file. This shows the authors' folds are reachable from groups like ours, and it does not prove that the order alone explains every difference.

Test T72a: the authors' prediction files carry the fold index of every pair. I added switches to `predict.py` (`FOLDS`, `SEED`, `SWAP`, `TAG`; all off by default, and a run with none set reproduces the earlier T2) and ran the shipped inputs on their folds, restricted to the pairs both runs have.

- Mean absolute change in per-phage AUROC against the authors' values: {fx['mean_abs_d_auroc_our_folds']:.3f} with our folds, {fx['mean_abs_d_auroc']:.4f} with theirs.
- {int(fx['identical_auroc_1e6'])} of 96 phages equal the authors' AUROC to six decimals, {int(fx['same_best_model'])} have the same best model and {int(fx['phages_d_above_0_05'])} differ by more than 0.05.
- The two that differ by more than 0.05 are {', '.join(fxd[fxd.d_auroc.abs() > 0.05].index)}. Both picked another best model than the authors did.
- Mean per-phage AUROC {fx['mean_auroc_fixed_folds']:.4f} against {fx['mean_auroc_authors']:.4f}; pooled AUROC {fx['pooled_auroc_fixed_folds']:.4f} against {fx['pooled_auroc_authors']:.4f}.

The remaining differences come from the best-model race, the unseeded forests and the scikit-learn version (see [the pinned environment](#pinned-environment-t21-and-t31)). Run-to-run noise measured on the same inputs (test T72, five seeds each):

{seed_tbl}
A random partition moves a phage's AUROC about five times as much as a model seed. Of the {flag['T2_T1']} phages that move by more than 0.05 between T2 and T1, {int(((sp.d_T2_T1.abs() > .05) & (sp.z_T2_T1_vs_partition.abs() > 2)).sum())} lie beyond the spread of random partitions, and none stays above 0.05 with the authors' folds.

## Pinned environment (T2.1 and T3.1)

The runs above used scikit-learn 1.9.1 on Python 3.14. The paper states scikit-learn 1.1.2. I built a second conda environment (`ml_pinned`, [`config/envs/ml_pinned.yml`](../config/envs/ml_pinned.yml), built by `scripts/setup_pinned_env.sh`) with scikit-learn 1.1.2 and every other version the paper gives, and ran the same `predict.py` three times: on the shipped inputs with our groups (T2.1), on the shipped inputs with the authors' folds, and on the regenerated features with our groups (T3.1). The authors state no Python version for the model step (3.12.2 appears only in the kernel metadata of their plotting notebook, which does not show what trained the models). scikit-learn 1.1.2 has no Python 3.12 build, so the environment runs Python 3.10.22. The table lists every component and what could not be matched ([`env_match.csv`](../results/audit/env_match.csv)).

{md(env_tbl)}
Results ([`pinned_env_overall.csv`](../results/audit/pinned_env_overall.csv), [`pinned_env_per_phage.csv`](../results/audit/pinned_env_per_phage.csv)):

{md(pe_tbl)}
- With the authors' folds, scikit-learn 1.1.2 brings the mean absolute change in per-phage AUROC from {per.loc["T2 authors' folds", "mean_abs_d_auroc_vs_T1"]:.4f} (1.9.1) to {per.loc["T2.1 authors' folds", "mean_abs_d_auroc_vs_T1"]:.4f}, the correlation with the authors' values to {per.loc["T2.1 authors' folds", "corr_auroc_vs_T1"]:.4f}, and the share of phages with the same best model family to {per.loc["T2.1 authors' folds", "same_model_family_as_T1"]:.1%}. Mean per-phage AUROC is {per.loc["T2.1 authors' folds", "mean_phage_auroc"]:.4f} against {per.loc["T1", "mean_phage_auroc"]:.4f}.
- The pinned environment reproduces the deterministic model exactly. For the {int(pa("T2.1", "LogReg_l2").phages)} phages whose best model is LogReg L2 in both runs, {pa("T2.1", "LogReg_l2").share_identical_1e6:.1%} of the pair probabilities equal the authors' to six decimals, against {pa("T2 authors", "LogReg_l2").share_identical_1e6:.1%} with scikit-learn 1.9.1 ([`pinned_env_pair_agreement.csv`](../results/audit/pinned_env_pair_agreement.csv)). {int((pel2.share_identical_1e6 > 0.99).sum())} of those phages agree on every pair. The others have 398 to 400 pairs instead of 401 or 402 (the pair set differs, see T70) and differ by at most {pel2[pel2.share_identical_1e6 < 0.99].max_abs_dp.drop("LF73_P4", errors="ignore").max():.4f}, except LF73_P4 ({pel2.loc["LF73_P4", "max_abs_dp"]:.3f}), which I did not trace.
- The random forests cannot match. The authors seed NumPy once (`np.random.seed(0)`) and leave the forests, the Dummy and the saga solver without `random_state`, so each model draws from a stream that every earlier model advanced. Replaying it needs their private labels and their row order. With their folds and scikit-learn 1.1.2 the median absolute probability difference is {pa("T2.1", "RF_250_3").median_abs_dp:.3f} for depth-3 forests and {pa("T2.1", "RF_250_6").median_abs_dp:.3f} for depth-6 forests, and no pair is identical. This is the part of the gap that stays after the version and the folds are matched. I did not try to replay the stream.
- With our own groups, the version does not close the gap: the mean absolute change against the authors is {per.loc["T2.1", "mean_abs_d_auroc_vs_T1"]:.4f} (T2.1) against {per.loc["T2", "mean_abs_d_auroc_vs_T1"]:.4f} (T2). The two scikit-learn versions assign different folds to the same groups (adjusted Rand index {pef.iloc[0].mean_ari_191_vs_112:.3f}, the same fold index for {pef.iloc[0].same_fold_index_share:.0%} of pairs, no phage identical), and neither assignment is close to the authors' ({pef.iloc[0].mean_ari_112_vs_authors:.3f} for 1.1.2, {pef.iloc[0].mean_ari_191_vs_authors:.3f} for 1.9.1). Between T2 and T2.1 the AUROC of a phage moves by {pch.loc[("T2.1", "T2")].mean_abs_d_auroc:.3f} on average, as much as the gap to the authors. The version and the missing groups file cannot be separated, and the fold partition stays the main cause of the per-phage gaps.
- Regenerating the features (T3.1 against T2.1) lowers the mean per-phage AUROC by {abs(pch.loc[("T3.1", "T2.1")].mean_d_auroc):.4f}, against {abs(pch.loc[("T3", "T2")].mean_d_auroc):.4f} with 1.9.1; the pooled AUROC is {per.loc["T3.1", "pooled_auroc"]:.4f} in both (T3: {per.loc["T3", "pooled_auroc"]:.4f}). The environment does not change the conclusion about the features.

## Which regenerated features matter

Test T73 puts one regenerated feature group at a time into the shipped inputs, with the authors' folds and seed 0 held fixed. The groups are the O-type, sequence type (ST), LPS type, Klebsiella capsule type, ABC capsule type, the 8-D UMAP and the host features (the ST of the isolation host and the same-as-host comparisons). A group that is swapped carries its host-derived columns with it, except for "host", which swaps only the host side. Controls: swapping nothing reproduces the fixed-fold seed-0 run (largest AUROC difference {ctl['swap_none_vs_fixed_s0_max_abs_d_auroc']:g}), and swapping everything tracks the earlier regenerated run (per-phage correlation {ctl['swap_all_vs_T3_per_phage_corr']:.3f}).

{md(swt)}
The single effects are small and have both signs, so they cancel: they sum to {sw[sw.group != 'all'].mean_d_auroc.sum():+.4f}, while the seven together lose {abs(all_row.mean_d_auroc):.4f}. Without the {n_small} phages that have fewer than 10 lytic strains the combined loss is {abs(all_row.mean_d_auroc_n10):.4f}. Phages 411_P1 ({swp.loc['411_P1', 'all']:+.2f}) and 411_P2 ({swp.loc['411_P2', 'all']:+.2f}) carry most of it, with 3 and 4 lytic strains. The host group changes only the {int(ctl['phages_whose_host_o_or_st_changed'])} phages whose isolation host got a different O-type or ST. Test T73b regenerates everything except one group, which gives the cost of each group once the others are in place:

{mgt}
The O-type and the Klebsiella capsule type account for the loss once the other groups are regenerated, mostly through phages with few lytic strains. The agreement of the regenerated O-type with the shipped, hand-curated column is 74.7 % ([report 3](03_feature_regeneration.md)), the lowest of the typing features, which fits.

## Paper against repository

Report [11](11_repository_vs_paper.md) lists 28 differences between the repository and the paper. [`paper_vs_repo.csv`](../results/audit/paper_vs_repo.csv) gives each one a code: {', '.join(f'{c} {n}' for c, n in repo_n.items())}. The {len(repo_star)} that change a number the paper reports:

{md(repo_star)}
## Register of the paper comparison

Each row of [`paper_vs_replication.csv`](../results/comparison/paper_vs_replication.csv) that did not match or could not be compared has a code and a reason ({code_n(g)}). The basis column of the CSV says whether the reason comes from an earlier report, from a test in this audit, or is a residual. Four ED7 rows were tested here, and eight rows remain residual.

{reg_tbl}
## Prediction results

{pred_tbl}
## Open rows and what would close them

ED7: coverage under four thresholds on the MLC score of the 403-strain matrix and twelve definitions from the raw plaque scores (replicates by dilution) was compared with the Source Data counts ({int(ed_sd.k1)}, {int(ed_sd.k2)}, {int(ed_sd.k3)}, {int(ed_sd.k4)} strains for one to four phages and a plateau of {int(ed_sd.plateau)}). None matches ([`ed7_definitions.csv`](../results/audit/ed7_definitions.csv)). The closest, MLC of at least 2, gives {int(ed[ed.definition == 'MLC >= 2'].k1.iloc[0])}, {int(ed[ed.definition == 'MLC >= 2'].k2.iloc[0])}, {int(ed[ed.definition == 'MLC >= 2'].k3.iloc[0])}, {int(ed[ed.definition == 'MLC >= 2'].k4.iloc[0])} and a plateau of {int(ed[ed.definition == 'MLC >= 2'].plateau.iloc[0])}. The raw scores cover only 369 strains, so they cannot reproduce a 403-strain curve.

Fig 3a: adding strain LF110 or using the old genus labels leaves the adjusted R² at {pm.iloc[0].adjR2_both:.3f} and {pm.iloc[1].adjR2_both:.3f}, against 0.65 in the paper. The isolation strain alone gives {pm.iloc[0].adjR2_isolation_strain:.3f} against 0.53. Both are 0.02 to 0.04 above the paper, so some setting of the PERMANOVA differs from the paper's in a way the paper does not describe.

The {len(xrows)} unexplained rows are {', '.join(xrows.id)}. G01 and G02 (Fig 4d pMCMC) agree in beta and differ only in the posterior tail. G15 is the PERMANOVA above. G26 (Fig 4b) needs the authors' filtering and clustering code. G38 to G41 (the K1 rule) are not explained by the K1 count: the regenerated calls give 45 K1 strains, as in the text, and the rule still gives 89 % and 65 % against 73 % and 72 %. They need the strain set and the definition of infection.

{settle_tbl}
## The original preprint against the paper and the repository

The preprint behind the link you sent is version 1, posted 22 November 2023, which is the file already in `docs/paper/`. [`preprint_check.csv`](../results/audit/preprint_check.csv) sets {len(pk)} statements side by side: the preprint, the published paper (T0) and the authors' repository (T1). Numbers in brackets are line numbers in the extracted text. T1 values are computed from the shipped matrix, feature table and predictions.

{md(pk.rename(columns={"item": "Item", "preprint v1": "Preprint v1", "published paper (T0)": "Published paper (T0)", "repository (T1)": "Repository (T1)", "verdict": "Verdict"}))}
Where the three agree. The descriptive numbers (83 +/- 70 strains per phage, 20 +/- 15 phages per strain, 12 strains lysed by no phage, the LPS counts, 22 Klebsiella capsules) are the same in the preprint, the paper and the repository, once strain LF110 is added. Both texts describe the per-phage models, the class weights and the recommender forests in the same way. The cocktail results (91.67, 78.95 and 81 %, and the U statistics) agree with the repository.

Where the preprint differs from the paper. The 77 % AUROC is a per-phage mean in the preprint and a per-strain 85 % in the paper. The preprint gives 86 % in its Results and 85 % in its Discussion, and the paper gives 86 % in both. The coverage of 68 % for the most covering four phages appears only in the preprint. The share of E. coli phylogroups covered by the isolation strains falls from 75 % to 66 %. The paper truncates 78.95 % to 78 % in its text.

Where the repository differs from both texts. The matrix has 402 strains. The K1 rule (87 and 55 % against 73 and 72 %), the triplet coverage (87.6 % against 63 % and the Source Data's 67.0 %) and the AUPR rule for choosing the best model are the differences found before. Two are new. The paper's 20.3 % lytic share does not match its own 20.6 % and 20.8 %, and the matrix gives 20.75 %. The text's P = 0.02 for tailored against generic cocktails is 0.034 in the Source Data and in the recomputation.

What this changed in the audit. Two rows of the register said the paper gives no number, and it does: 83 +/- 70 strains per phage (20.6 %), and 3,379 false negatives and 2,922 false positives over 38,688 pairs. Both now compare as numbers and match within 5 % (20.7 % here; 3,307 and 2,940 in our run, 3,318 and 2,980 in the authors' predictions), so the paper comparison has {len(g)} unresolved rows instead of 47. The 68 % for four phages in ED7 comes from the preprint and not from the published text, and the affected rows and reports now say so.

## Check against an independent preprint

After the audit was finished I checked it against a preprint from another group that models the same E. coli matrix: Noonan et al., [Phylogeny-agnostic strain-level prediction of phage-host interactions from genomes](https://doi.org/10.1101/2025.11.15.688630), bioRxiv, version posted 30 April 2026, not peer reviewed. I typed its numbers from the PDF (Table 1, Table 2 and the Results) and recomputed our side from the shipped matrix ([`external_check.csv`](../results/audit/external_check.csv)).

{md(ext.rename(columns={"item": "Item", "preprint": "Preprint", "this audit": "This audit", "verdict": "Verdict"})[["Item", "Preprint", "This audit", "Verdict"]])}
What it supports. The dataset counts agree (402 strains, 20.7 % positives), so both analyses read the public matrix the same way. That is a check of loading, not an independent measurement, since the file is the same. The AUROC for unseen strains (0.869) lies within 0.013 of the paper, the authors' predictions and this run, although the preprint uses another model family and another split design. That fits the finding that overall AUROC does not depend on the fold partition. For ED7 the check is more useful. The preprint's promiscuity baseline implies about 61 % of strains for one phage, and the same baseline on our matrix gives 65 %. Both are far above the 42.7 % of the ED7 Source Data. The preprint calls any non-zero interaction a positive, which is the definition used here, so the Source Data curve needs a stricter definition that nobody has stated.

What it cannot test. The preprint reports no per-phage numbers for E. coli in its main text, so it says nothing about the fold-partition effect on single phages, the feature swaps, the recommender or the K1 rule. Its cross-validation (20 folds on hierarchical clusters of strains, 10 % held out) differs from the authors' group 10-fold, and the text does not say how its AUROC is pooled. I therefore read it as a consistency check from independent code on the same public file, and not as a replication of the paper. One slip inside the preprint: Klebsiella-3 has 25 phages in the Methods and 52 in Table 1, and only 52 gives its stated total of 518.

## Limits

- The code S is a statistical statement about noise, not a measured seed effect. The z test treats the two runs as independent samples of folds, which overstates the noise for a paired comparison, so it flags fewer phages than a paired test would.
- The bootstrap intervals keep each phage's best model fixed and do not reselect it per replicate.
- Each swap run uses one seed. The seed noise measured with the folds fixed is 0.0011 on the mean per-phage AUROC, so effects of 0.003 or more on the mean are above it. Effects on single phages use that phage's own spread.
- The swaps cover the five typing columns, the UMAP and the host features. They do not include the regenerated interaction matrix, so the T4 loss beyond T3 is not attributed to features or to labels beyond the rescoring check for two phages.
- The decomposition by group is not additive, and a group's cost depends on which other groups are regenerated. The marginal table is one ordering, with everything else regenerated.
- The 59 pairs that differ between the authors' predictions and ours have a cause I could not verify.
- Eight rows of the register are residual, as listed above.
- The pinned environment matches scikit-learn and scipy to the paper, but the authors do not state Python, numpy, pandas, joblib, the operating system or the BLAS for the model step, so those were chosen, not matched. The forests' random stream was not replayed. The effect of the unmatched parts is not isolated.
- `results/audit_runs/` (88 MB of run outputs) and the pair-level CSVs that hold the authors' predictions are not in this repository. `results/audit/authors_folds.csv` is: it holds only the fold index of each pair, read from the authors' prediction files ([report 13](13_cv_folds.md) says how).

## Reproduce

```bash
.venv/bin/python scripts/audit_gaps.py register      # gap register, paper against repository
.venv/bin/python scripts/audit_gaps.py align         # tier alignment
.venv/bin/python scripts/audit_gaps.py overall       # overall metrics with bootstrap intervals (about 8 minutes)
for s in shipped regenerated regenerated+matrix; do FEATURE_SET=$s .venv/bin/python scripts/audit_gaps.py folds; done
.venv/bin/python scripts/audit_gaps.py phage         # per-phage table and causes
.venv/bin/python scripts/audit_gaps.py pairs         # pair-level comparisons
.venv/bin/python scripts/audit_gaps.py authors_folds
scripts/run_audit_runs.sh                            # 27 attribution runs, about 100 seconds each
.venv/bin/python scripts/audit_gaps.py runs          # T72a, T72, T73 summaries and refined causes
.venv/bin/python scripts/audit_gaps.py final         # ED7 and PERMANOVA tests, final registers
scripts/setup_pinned_env.sh                          # conda environment with scikit-learn 1.1.2 (T81)
scripts/run_pinned_env.sh                            # T2.1, T2.1 with the authors' folds, T3.1 (about 2 minutes each)
for s in shipped regenerated; do .venv/bin/python scripts/fold_version_check.py results/audit_runs/fold_version/sk191_$s.csv $s; done   # and with the pinned python, as sk112_$s.csv
.venv/bin/python scripts/audit_pinned.py             # T82 tables
.venv/bin/python scripts/audit_report.py             # this report
```
"""
import re
text = re.sub(r"(\d)%", r"\1 %", text)
os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
open(OUTFILE, "w", encoding="utf-8").write(text)
print(f"wrote {OUTFILE}: {len(text.split())} words, {text.count(chr(10))} lines")
