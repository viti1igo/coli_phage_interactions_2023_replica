"""Re-implementation of the paper's four-step cocktail recommender from its Methods (the authors' code is not public) and comparison with the
cocktails they recommended (dev/cocktails/data/cocktails_step.csv).  python scripts/recommender.py   -> results/recommender/

Steps (README sections 2.3-2.6 of the authors' repository; the Methods p.2859 differ slightly, see SPEC below): (1) random forest (depth 4, 80 trees,
positive class weight 1.5) on phage identity + booleans 'query has the same O / ST / ABC / Klebsiella / LPS as the isolation strain of the phage'; (2) strains of the Picard collection with the same haplotype (ST, O, H): phages that
infect the majority; (3) the same random forest on phage identity + all bacterial features; (4) the most generalist phages.  Phages with p > 0.5
are recommended in order of probability, steps in sequence, until three phages; two phages of the same genus and isolation host never co-occur.

Limits: the test strains' genomes are not shipped.  Supplementary Table 6 gives ST, O, H, LPS and phylogroup only, so ABC serotype and Klebsiella
capsule are unknown for them, and their UMAP coordinates are recomputed from the shipped distance matrix (together with the Picard strains)."""
import os, sys
import numpy as np, pandas as pd
import umap
from sklearn.ensemble import RandomForestClassifier
import predict as P
import supp

ROOT = P.ROOT
OUT = f"{ROOT}/results/recommender" + os.environ.get("VARIANT", "")
os.makedirs(OUT, exist_ok=True)
C = f"{ROOT}/dev/cocktails/data"
# SPEC: "readme" = the authors' README (step 1 compares O, ST, ABC, Klebsiella, LPS; step 2 keeps phages infecting at least half of the group);
#       "paper" = the Methods text (step 1 compares ST, O, H, LPS, ABC; step 2 keeps more than half)
SPEC = os.environ.get("SPEC", "readme")
STEP1_KEYS = ("ST_Warwick", "O-type", "Klebs_capsule_type", "LPS_type", "ABC_serotype") if SPEC == "readme" else ("ST_Warwick", "O-type", "H-type", "LPS_type", "ABC_serotype")
ENC = os.environ.get("PHAGE_ENCODING", "onehot")            # onehot | ordinal (phage identity as one integer column)
MAXF = "sqrt" if os.environ.get("MAX_FEATURES", "all") == "sqrt" else None   # all features at each split (default); "sqrt" never gives p > 0.5
rf = lambda: RandomForestClassifier(max_depth=4, n_estimators=80, class_weight={0: 1, 1: 1.5}, max_features=MAXF, random_state=0, n_jobs=-1)

Y = P.Y
BACT = P.BACT.copy()
BACT.index = BACT.index.astype(str)
guelin = pd.read_csv(f"{ROOT}/data/genomics/phages/guelin_collection.csv", sep=";", index_col=0)
phages = list(Y.columns)
train = [s for s in Y.index if s in BACT.index]
t6 = supp.sheet(4, "Supplementary Table 6").set_index("bacteria")
t6["ST_Warwick"] = t6.ST_Warwick.astype(str)
BACT["ST_Warwick"] = BACT.ST_Warwick.astype(str)

# ---- UMAP of Picard + test strains from the shipped distance matrix
D = pd.read_csv(f"{C}/picard+test_collection_phylogenetic_distances.tsv", sep="\t", index_col=0)
D.index = D.index.astype(str); D.columns = D.columns.astype(str)
emb = pd.DataFrame(umap.UMAP(n_components=8, metric="precomputed", n_neighbors=15, min_dist=0.1, random_state=0).fit_transform(D.values),
                   index=D.index, columns=[f"UMAP{i}" for i in range(8)])
test_names = list(t6.index)
test = ["T:" + s for s in test_names]  # one test strain (AN03) has the name of a Picard strain
emb_t = emb.reindex([s for s in test_names if s not in train]).rename(index=lambda s: "T:" + s)
emb_all = pd.concat([emb.loc[[s for s in emb.index if s in set(train)]], emb_t])
no_umap = [s for s in test if s not in emb_all.index]
print(f"test strains {len(test)}; without UMAP coordinates (not in the distance matrix): {no_umap}")

# ---- strain feature table for Picard (from the shipped table) and test strains (from Supplementary Table 6)
F = pd.concat([BACT.loc[train, ["O-type", "H-type", "ST_Warwick", "LPS_type", "Klebs_capsule_type", "ABC_serotype"]],
               t6.reindex(test_names).set_axis(test)[["O-type", "H-type", "ST_Warwick", "LPS_type"]].assign(Klebs_capsule_type=np.nan, ABC_serotype=np.nan)])
F = F.join(emb_all, how="left")
for c in ("O-type", "ST_Warwick"):  # lump values seen fewer than 3 times in the Picard collection
    keep = F.loc[train, c].value_counts(); keep = set(keep[keep >= 3].index)
    F[c + "_l"] = F[c].where(F[c].isin(keep) | F[c].isna(), "Other")
host = {p: guelin.loc[p, "Phage_host"] for p in phages}
HF = {p: (BACT.loc[h] if h in BACT.index else None) for p, h in host.items()}


def same(a, b):
    return float(pd.notna(a) and pd.notna(b) and a == b)


def phage_cols(n):
    ids = pd.Series(np.tile(range(len(phages)), n))
    if ENC == "ordinal":
        return ids.to_frame("phage").astype(float)
    return pd.get_dummies(ids.astype(pd.CategoricalDtype(range(len(phages)))), prefix="ph").astype(float)


FD = F[list(STEP1_KEYS)].to_dict("index")
HD = {p: (h[list(STEP1_KEYS)].to_dict() if h is not None else None) for p, h in HF.items()}


def pairs_step1(strains):
    rows = [[same(FD[s][k], HD[p][k]) if HD[p] is not None else 0. for k in STEP1_KEYS] for s in strains for p in phages]
    X = pd.DataFrame(rows, columns=[f"same_{k}" for k in STEP1_KEYS])
    return pd.concat([X, phage_cols(len(strains))], axis=1)


num = [f"UMAP{i}" for i in range(8)]
mu, sd = F.loc[train, num].mean(), F.loc[train, num].std()


def pairs_step3(strains):
    cat = pd.get_dummies(F[["O-type_l", "ST_Warwick_l", "LPS_type", "Klebs_capsule_type", "ABC_serotype"]].astype(object), dtype=float)
    z = (F[num] - mu) / sd
    S = pd.concat([z, cat], axis=1).loc[strains].fillna(0.)
    X = S.loc[np.repeat(strains, len(phages))].reset_index(drop=True)
    return pd.concat([X, phage_cols(len(strains))], axis=1), cat.columns


lab = Y.loc[train].stack().rename("y").reset_index()
lab.columns = ["bacteria", "phage", "y"]
lab = lab.dropna(subset=["y"])
idx_tr = [(s, p) for s, p in zip(lab.bacteria, lab.phage)]
y_tr = lab.y.values


def build(kind):
    if kind == 1:
        Xa = pairs_step1(train)
    else:
        Xa, _ = pairs_step3(train)
    Xa.index = pd.MultiIndex.from_product([train, phages])
    return Xa


Xtr1, Xtr3 = build(1), build(3)
sel = pd.MultiIndex.from_tuples(idx_tr)
m1 = rf().fit(Xtr1.loc[sel].values, y_tr)
m3 = rf().fit(Xtr3.loc[sel].values, y_tr)
ok = [s for s in test if s not in no_umap]
P1 = pd.DataFrame(m1.predict_proba(pairs_step1(test).values)[:, 1].reshape(len(test), -1), index=test, columns=phages)
X3, _ = pairs_step3(ok)
P3 = pd.DataFrame(m3.predict_proba(X3.values)[:, 1].reshape(len(ok), -1), index=ok, columns=phages).reindex(test)
generalists = list(Y.loc[train].gt(0).sum().sort_values(ascending=False, kind="stable").index)
txt = lambda x: "nan" if pd.isna(x) else str(x)
hap = F.loc[train, ["ST_Warwick", "O-type", "H-type"]].apply(lambda r: "|".join(txt(x) for x in r), axis=1).where(F.loc[train, ["ST_Warwick", "O-type", "H-type"]].notna().all(axis=1))


def step2(s):
    if F.loc[s, ["ST_Warwick", "O-type", "H-type"]].isna().any():  # a haplotype needs ST, O and H
        return pd.Series(dtype=float), 0
    h = "|".join(txt(x) for x in F.loc[s, ["ST_Warwick", "O-type", "H-type"]])
    g = hap[hap == h].index
    if len(g) == 0:
        return pd.Series(dtype=float), 0
    fr = Y.loc[g].mean()
    return fr[(fr >= 0.5) if SPEC == "readme" else (fr > 0.5)].sort_values(ascending=False, kind="stable"), len(g)


def allowed(p, chosen):
    return not any(guelin.loc[p, "Genus"] == guelin.loc[c, "Genus"] and host[p] == host[c] for c in chosen)


rec = []
for s in test:
    chosen, steps = [], []
    s2, ng = step2(s)
    cands = [("0_same_as_host", P1.loc[s][P1.loc[s] > 0.5].sort_values(ascending=False, kind="stable")), ("1_neighbours", s2),
             ("3_misc", P3.loc[s][P3.loc[s] > 0.5].sort_values(ascending=False, kind="stable") if s in ok else pd.Series(dtype=float)),
             ("4_rescue", pd.Series(1.0, index=generalists))]
    for name, c in cands:
        for p in c.index:
            if len(chosen) == 3: break
            if p not in chosen and allowed(p, chosen):
                chosen.append(p); steps.append(name)
    rec.append(dict(bacteria=s, p1=chosen[0], p2=chosen[1], p3=chosen[2], step1=steps[0], step2=steps[1], step3=steps[2], n_neighbours=ng))
ours = pd.DataFrame(rec).set_index("bacteria")
ours.rename(index=lambda x: x[2:]).to_csv(f"{OUT}/our_cocktails.csv")

# ---- compare with the authors' recommendations
sh = pd.read_csv(f"{C}/cocktails_step.csv", sep=";")
shc = sh.groupby("bacteria").agg(phages=("phage", lambda x: sorted(x)), steps=("chosen_at_step", lambda x: list(x)))
rows = []
for s in test:
    n = s[2:]
    if n not in shc.index: continue
    a, b = sorted(ours.loc[s, ["p1", "p2", "p3"]]), shc.loc[n, "phages"]
    rows.append(dict(bacteria=n, overlap=len(set(a) & set(b)), exact=a == b, ours=" + ".join(a), authors=" + ".join(b),
                     authors_steps=" / ".join(sorted(set(shc.loc[n, "steps"]))), our_steps=" / ".join(sorted(set(ours.loc[s, ["step1", "step2", "step3"]])))))
cmp = pd.DataFrame(rows).set_index("bacteria")
cmp.to_csv(f"{OUT}/comparison_with_authors.csv")
print(f"strains compared: {len(cmp)}; identical cocktails: {int(cmp.exact.sum())}; mean overlap {cmp.overlap.mean() / 3:.2f} of 3 phages")
print("overlap distribution:", cmp.overlap.value_counts().sort_index().to_dict())
sh["rank_step"] = sh.chosen_at_step
agree = {}
for st in sorted(sh.chosen_at_step.unique()):
    sub = sh[sh.chosen_at_step == st]
    hit = sum(p in set(ours.loc["T:" + b, ["p1", "p2", "p3"]]) for b, p in zip(sub.bacteria, sub.phage) if "T:" + b in ours.index)
    agree[st] = f"{hit} of {len(sub)}"
print("authors' recommended phages also in our cocktail, by the authors' step:", agree)
for k, v in {"ours": ours[["step1", "step2", "step3"]].stack().value_counts(), "authors": sh.chosen_at_step.value_counts()}.items():
    print(k, (v / v.sum() * 100).round(1).to_dict())
pd.Series(agree, name="authors_phages_in_our_cocktail").rename_axis("authors_step").to_csv(f"{OUT}/agreement_by_step.csv")
