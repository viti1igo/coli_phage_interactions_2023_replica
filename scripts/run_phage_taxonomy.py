"""T27: phage genus / species from whole-genome similarity. The paper used VIRIDIC 1r3.6 (genus >= 70%, species >= 95% intergenomic similarity);
VIRIDIC itself is not installed, so its measure is approximated: BLASTn all-vs-all, non-overlapping HSPs chosen greedily by bit score,
similarity = 2 x identical aligned nucleotides / (len A + len B). Single-linkage clusters are compared with the shipped Genus and Species columns.
python scripts/run_phage_taxonomy.py   ->   results/replication/phage_taxonomy_agreement.csv"""
import glob, os, re, subprocess
from collections import defaultdict
import numpy as np, pandas as pd
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = f"{ROOT}/replication/T27_phage_taxonomy"
os.makedirs(OUT, exist_ok=True)
env = dict(os.environ, MAMBA_ROOT_PREFIX=os.path.expanduser("~/micromamba-phage"))
run = lambda *c: subprocess.run(["micromamba", "run", "-n", "blast", *c], env=env, check=True, capture_output=True, text=True)

files = sorted(f for f in glob.glob(f"{ROOT}/data/genomics/phages/FNA/*.fna") if not os.path.basename(f).startswith("._"))
lens, rec = {}, []
with open(f"{OUT}/all.fna", "w") as o:
    for f in files:
        p = os.path.basename(f)[:-4]; seq = "".join(l.strip() for l in open(f) if not l.startswith(">"))
        lens[p] = len(seq); o.write(f">{p}\n{seq}\n")
if not os.path.exists(f"{OUT}/hits.tsv"):
    run("makeblastdb", "-in", f"{OUT}/all.fna", "-dbtype", "nucl", "-out", f"{OUT}/db")
    run("blastn", "-query", f"{OUT}/all.fna", "-db", f"{OUT}/db", "-task", "blastn", "-evalue", "1e-3", "-dust", "no", "-num_threads", "4",
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore", "-out", f"{OUT}/hits.tsv")
h = pd.read_csv(f"{OUT}/hits.tsv", sep="\t", names=["q", "s", "pid", "len", "qs", "qe", "ss", "se", "bits"])
h = h[h.q != h.s]
names = list(lens); ix = {n: i for i, n in enumerate(names)}
ident = np.zeros((len(names), len(names)))
for (q, s), g in h.groupby(["q", "s"]):
    taken_q, taken_s, tot = [], [], 0.0
    for r in g.sort_values("bits", ascending=False).itertuples():
        qa, qb = sorted((r.qs, r.qe)); sa, sb = sorted((r.ss, r.se))
        if any(qa <= b and a <= qb for a, b in taken_q) or any(sa <= b and a <= sb for a, b in taken_s): continue
        taken_q.append((qa, qb)); taken_s.append((sa, sb)); tot += r.len * r.pid / 100
    ident[ix[q], ix[s]] = tot
ident = (ident + ident.T) / 2
L = np.array([lens[n] for n in names]); sim = 100 * 2 * ident / (L[:, None] + L[None, :]); np.fill_diagonal(sim, 100)
ship = pd.read_csv(f"{ROOT}/data/genomics/phages/guelin_collection.csv", sep=";", index_col=0).reindex(names)
print("phages without a shipped Genus:", int(ship.Genus.isna().sum()), "| without Species:", int(ship.Species.isna().sum()), "(labelled NA)")
ship = ship.fillna("NA")
rows = []
for rank, thr in (("Genus", 70), ("Species", 95)):
    lab = connected_components(sim >= thr, directed=False)[1]
    rows.append(dict(rank=rank, threshold_pct=thr, phages=len(names), our_clusters=len(set(lab)), shipped_groups=ship[rank].nunique(), adjusted_Rand_index=round(adjusted_rand_score(ship[rank], lab), 3)))
    pd.DataFrame({"phage": names, "ours": lab, "shipped": ship[rank].values}).to_csv(f"{OUT}/{rank.lower()}_clusters.csv", index=False)
os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
pd.DataFrame(rows).to_csv(f"{ROOT}/results/replication/phage_taxonomy_agreement.csv", index=False)
pd.DataFrame(sim, index=names, columns=names).round(2).to_csv(f"{OUT}/similarity_matrix.csv")
print(pd.DataFrame(rows).to_string(index=False))
