"""Fig 4A/B from MCMC outputs: count significant (pMCMC < 1e-3, the paper's threshold) adsorption-factor vs defense-system terms,
and compare with the AF-only results shipped in dev/inference_bacteria.   python scripts/fig4_summarize.py"""
import glob, os
import pandas as pd

ROOT = os.path.join(os.path.dirname(__file__), "..")
AF_PREFIX = ("O.type", "H.type", "LPS_type", "Klebs_capsule_type", "ABC_serotype",
             "BTUB", "FADL", "FHUA", "LAMB", "LPTD", "NFRA", "OMPA", "OMPC", "OMPF", "TOLC", "TSX", "YNCD")  # all non-DS covariates in the model


def load(pattern, name):
    rows = []
    for f in glob.glob(pattern):
        d = pd.read_csv(f, sep="\t")
        d.columns = ["feature", "post_mean", "lo", "hi", "eff_samp", "pMCMC"]
        d = d[d.feature != "(Intercept)"].assign(phage=name(f))
        rows.append(d)
    d = pd.concat(rows)
    d["group"] = ["Adsorption factor" if f.startswith(AF_PREFIX) else "Defense system" for f in d.feature]
    d["term"] = d.feature.str.extract(r"^([A-Za-z]+(?:\.type|_type|_serotype|_capsule_type)?)")[0]
    d["sig"] = d.pMCMC < 1e-3
    return d


mine = load(f"{ROOT}/results/fig4/*.tsv", lambda f: os.path.basename(f)[:-4].removeprefix("X"))
ship = load(f"{ROOT}/dev/inference_bacteria/phage_per_phage/results/AF_niter=500000_minfreq=8/*.res",
            lambda f: os.path.basename(f).split("__")[0].removeprefix("X"))
os.makedirs(f"{ROOT}/results/fig4_summary", exist_ok=True)

print(f"phages fitted here: {mine.phage.nunique()}   shipped (AF-only): {ship.phage.nunique()}")
print("\nSignificant terms (pMCMC<1e-3), AF+DS model re-run (Fig 4A):")
a = mine[mine.sig].groupby("group").size()
print(a.to_string())
print("\nPer covariate family, re-run vs shipped AF-only (Fig 4B):")
cmp = pd.DataFrame({"rerun_AF+DS": mine[mine.sig & (mine.group == "Adsorption factor")].groupby("term").size(),
                    "shipped_AF_only": ship[ship.sig].groupby("term").size()}).fillna(0).astype(int)
print(cmp.sort_values("rerun_AF+DS", ascending=False).to_string())
cmp.to_csv(f"{ROOT}/results/fig4_summary/af_hits_by_family.csv")
a.to_csv(f"{ROOT}/results/fig4_summary/af_vs_ds_hits.csv")

# phage-level agreement of the sign/significance of shared terms
m = mine.merge(ship, on=["phage", "feature"], suffixes=("_new", "_ship"))
print(f"\nshared (phage, term) pairs: {len(m)}; both significant: {(m.sig_new & m.sig_ship).sum()}; "
      f"sig in shipped only: {(~m.sig_new & m.sig_ship).sum()}; sig in re-run only: {(m.sig_new & ~m.sig_ship).sum()}")
print(f"Spearman(pMCMC) over shared pairs: {m.pMCMC_new.corr(m.pMCMC_ship, method='spearman'):.3f}")
