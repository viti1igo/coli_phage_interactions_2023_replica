"""T25: compare the regenerated bacterial features with the shipped picard_collection.csv, step by step.
Writes results/replication/feature_agreement.csv (one row per feature) and feature_mismatches.csv (every disagreement)
and replication/T25_features/picard_collection_regenerated.csv.   python scripts/compare_features.py"""
import os, re
import numpy as np, pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..")
R = f"{ROOT}/replication"
SHIP = pd.read_csv(f"{ROOT}/data/genomics/bacteria/picard_collection.csv", sep=";", index_col=0)
stem = lambda s: re.sub(r"\.fasta$", "", os.path.basename(str(s)))
norm = lambda x: np.nan if pd.isna(x) or str(x).strip() in ("", "-", "nan", "None") else str(x).strip()
exists = lambda p: os.path.exists(f"{R}/{p}")
ours = {}  # feature -> Series indexed by strain


EXTRA = {}  # columns saved in the regenerated table but not compared directly
RAW = []  # comparisons with the authors' RAW tool outputs, where they ship them


def load_ectyper():
    d = pd.read_csv(f"{R}/T18_ectyper/output.tsv", sep="\t")
    d = d[~d.Name.str.startswith("._")].set_index("Name")
    raw = pd.read_csv(f"{ROOT}/data/genomics/bacteria/o_type/output.tsv", sep="\t").set_index("Name")
    for col in ("O-type", "H-type"):  # raw ECTyper output (incl. "-" and multi-calls) vs the authors' raw output
        j = d[col].to_frame("ours").join(raw[col].rename("ship"), how="inner")
        RAW.append(dict(feature=f"{col} (raw ECTyper output)", n_strains=len(j), n_agree=int((j.ours == j.ship).sum()),
                        agreement_pct=round(100 * (j.ours == j.ship).mean(), 2), shipped_missing=int((j.ship == "-").sum()), ours_missing=int((j.ours == "-").sum())))
    ours["O-type"], ours["H-type"] = d["O-type"].map(norm), d["H-type"].map(norm)


def load_mlst():
    d = pd.read_csv(f"{R}/T16_mlst/mlst.tsv", sep="\t", header=None)
    ours["ST_Warwick"] = pd.Series(d[2].map(norm).values, index=d[0].map(stem))


def load_clermont():
    d = pd.read_csv(f"{R}/T17_clermont/all/all_phylogroups.txt", sep="\t", header=None)
    ours["Clermont_Phylo"] = pd.Series(d[4].map(norm).values, index=d[0].map(stem))


def load_kaptive():
    d = pd.read_csv(f"{R}/T21_kaptive/kaptive_results_table.txt", sep="\t").set_index("Assembly")
    d.index = d.index.map(stem)
    ok = d["Match confidence"].isin(["Perfect", "Very high", "High", "Good"])  # the shipped K types correspond to confident hits only
    ours["Klebs_capsule_type"] = d["Best match locus"].where(ok).map(lambda x: np.nan if pd.isna(x) else "K" + re.sub(r"\D", "", x))


CAPS = {}  # strain -> (model, wholeness) of the best ABC-capsule system; also saved as a cross-table with the shipped serotype


def load_capsule():
    names = pd.read_csv(f"{R}/T14_panacota/annot/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain
    best = {}
    for gem, strain in names.items():
        f = f"{R}/T20_capsule/runs/{gem}/best_solution.tsv"
        if not os.path.exists(f):
            continue
        try:
            d = pd.read_csv(f, sep="\t", comment="#")
        except pd.errors.EmptyDataError:  # no capsule system found: the file holds only comment lines
            best[strain] = np.nan
            continue
        if len(d):
            top = d.sort_values("sys_wholeness", ascending=False).iloc[0]
            best[strain] = top.model_fqn.split("/")[-1]
        else:
            best[strain] = np.nan
    ours["ABC_serotype (model name)"] = pd.Series(best)


def clermont_norm(x):  # ClermonTyping writes cladeI, albertii; the shipped table writes 'Clade I', 'E. albertii'
    return np.nan if pd.isna(x) else re.sub(r"^e\.|[\s.]", "", str(x).lower())


rows, mism = [], []
for step, fn in [("T18", load_ectyper), ("T16", load_mlst), ("T17", load_clermont), ("T21", load_kaptive), ("T20", load_capsule)]:
    try:
        fn()
    except (FileNotFoundError, KeyError, pd.errors.EmptyDataError) as e:
        print(f"{step}: skipped ({type(e).__name__}: {e})")

if "ABC_serotype (model name)" in ours:
    c = ours.pop("ABC_serotype (model name)").reindex(SHIP.index)
    EXTRA["ABC_model (best system)"] = c  # saved in the table, not compared one-to-one (the shipped labels differ in naming)
    sh = SHIP["ABC_serotype"].map(norm)
    pres = (sh.notna()) == (c.notna())
    named = sh.notna() & c.notna() & sh.str.match(r"^K\d+$") & c.str.match(r"^K\d+")
    RAW.append(dict(feature="ABC capsule present (shipped has a serotype vs we found a system)", n_strains=len(sh), n_agree=int(pres.sum()), agreement_pct=round(100 * pres.mean(), 2), shipped_missing=int(sh.isna().sum()), ours_missing=int(c.isna().sum())))
    ok = (sh == c.str.replace(r"_like$", "", regex=True))[named]
    RAW.append(dict(feature="ABC serotype, strains with a K-type in both", n_strains=int(named.sum()), n_agree=int(ok.sum()), agreement_pct=round(100 * ok.mean(), 2), shipped_missing=0, ours_missing=0))
    os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
    pd.crosstab(sh.fillna("none"), c.fillna("none"), rownames=["shipped ABC_serotype"], colnames=["our best model"]).to_csv(f"{ROOT}/results/replication/capsule_crosstab.csv")

for feat, s in ours.items():
    a = SHIP[feat].map(norm)
    if feat == "ST_Warwick":
        a = a.where(a != "new")  # the shipped table writes "new" for a novel ST; mlst reports no ST ("-")
    b = s.reindex(a.index)
    if feat == "Clermont_Phylo":
        a, b = a.map(clermont_norm), b.map(clermont_norm)
    if feat == "ST_Warwick":
        b = b.map(lambda x: x if pd.isna(x) else x.split("*")[0].lstrip("~"))
    both_na = a.isna() & b.isna()
    agree = (a == b) | both_na
    rows.append(dict(feature=feat, n_strains=len(a), n_agree=int(agree.sum()), agreement_pct=round(100 * agree.mean(), 2),
                     shipped_missing=int(a.isna().sum()), ours_missing=int(b.isna().sum())))
    for k in a.index[~agree]:
        mism.append(dict(feature=feat, strain=k, shipped=a[k], ours=b[k]))
rows = RAW + rows
os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
pd.DataFrame(rows).to_csv(f"{ROOT}/results/replication/feature_agreement.csv", index=False)
pd.DataFrame(mism).to_csv(f"{ROOT}/results/replication/feature_mismatches.csv", index=False)
os.makedirs(f"{R}/T25_features", exist_ok=True)
pd.DataFrame({**ours, **EXTRA}).rename_axis("bacteria").to_csv(f"{R}/T25_features/picard_collection_regenerated.csv", sep=";")
print(pd.DataFrame(rows).to_string(index=False))
