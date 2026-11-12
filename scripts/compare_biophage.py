"""Phase 21 (T96): biophage-reproduce (B1 as written, B2 paper-like env) against T2.1 and T1, all on the authors' shipped inputs.
Needs replication/T93_biophage/<env>/training_results/{per_phage_perf.csv, predictions_named.csv, cv_groups_named.csv}.
Output: results/audit/biophage_compare.csv (one row per run) and biophage_compare_per_phage.csv.   .venv/bin/python scripts/compare_biophage.py"""
from pathlib import Path
import pandas as pd
from sklearn.metrics import adjusted_rand_score as ari, roc_auc_score

R = Path(__file__).resolve().parents[1]
A = R / "results/audit"
pp = pd.read_csv(A / "pinned_env_per_phage.csv").set_index("phage")
overall = pd.read_csv(A / "pinned_env_overall.csv").set_index("run")
ours = pd.read_csv(A / "our_folds.csv")
auth = pd.read_csv(A / "authors_folds.csv")
groups_ours = ours.drop_duplicates("bacteria").set_index("bacteria").group
t21 = pd.read_csv(R / "results/audit_runs/T2_1_shipped/per_phage.csv").set_index("phage")


def fold_ari(f, ref, col):
    m = f.merge(ref, on=["phage", "bacteria"], suffixes=("", "_ref"))
    return m.groupby("phage").apply(lambda d: ari(d.fold, d[col])).mean()


rows, per = [], {"T1": pp.T1, "T2.1": pp["T2.1"]}
for o in ["T1", "T2.1"]:
    rows.append(dict(run=o, what=overall.loc[o, "what"], mean_phage_auroc=overall.loc[o, "mean_phage_auroc"],
                     mean_phage_aupr=overall.loc[o, "mean_phage_aupr"], pooled_auroc=overall.loc[o, "pooled_auroc"]))
for tag, env, what, ours_col in [("B1", "bp_asis", "biophage as written: its training.yml (scikit-learn 1.9.1, pandas 3.1.1)", "fold_sklearn_1_9_1"),
                                 ("B2", "bp_pinned", "biophage code in our ml_pinned + biopython 1.79 (scikit-learn 1.1.2, pandas 1.5.3)", "fold_sklearn_1_1_2")]:
    d = R / "replication/T93_biophage" / env / "training_results"
    s = pd.read_csv(d / "per_phage_perf.csv").set_index("phage")
    pr = pd.read_csv(d / "predictions_named.csv")
    best = pr.merge(s.model.rename("best").reset_index(), on="phage").query("model == best")
    known = best[~best.label_was_nan]
    g = pd.read_csv(d / "cv_groups_named.csv", index_col=0).group.reindex(groups_ours.index)
    folds = best[["phage", "bacteria", "fold"]]
    per[tag] = s.auc
    rows.append(dict(run=tag, what=what, mean_phage_auroc=s.auc.mean(), mean_phage_aupr=s.avg_precision.mean(),
                     pooled_auroc=roc_auc_score(best.y_true, best.y_pred_proba),
                     pooled_auroc_known_labels=roc_auc_score(known.y_true, known.y_pred_proba),
                     nan_cells_as_negative=int(best.label_was_nan.sum()),
                     mean_abs_d_auroc_vs_T1=(s.auc - pp.T1).abs().mean(), corr_auroc_vs_T1=s.auc.corr(pp.T1),
                     mean_abs_d_auroc_vs_T2_1=(s.auc - pp["T2.1"]).abs().mean(),
                     best_family_same_as_T2_1=(s.model.str[:2] == t21.best_model.reindex(s.index).str[:2]).mean(),
                     groups_ari_vs_ours=ari(g, groups_ours), n_groups=g.nunique(),
                     fold_ari_vs_authors=fold_ari(folds, auth.rename(columns={"fold": "f_ref"}), "f_ref"),
                     fold_ari_vs_ours_same_sklearn=fold_ari(folds, ours[["phage", "bacteria", ours_col]], ours_col)))
out = pd.DataFrame(rows)
out.to_csv(A / "biophage_compare.csv", index=False)
pd.DataFrame(per).assign(n_lytic=pp.n_lytic).to_csv(A / "biophage_compare_per_phage.csv")
pd.set_option("display.width", 250, "display.max_columns", 30)
print(out.drop(columns="what").round(4).T)
