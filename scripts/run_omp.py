"""T23: variants of the outer-membrane proteins (Fig 1h, Fig 4). Method is not specified in the paper beyond "gene clusters with >99% identity".
Here: BLASTP of Prokka's bundled UniProt entries for each gene against all proteomes, best hit per genome and gene, then mmseqs
clustering of the gene nucleotide sequences at 99% identity; clusters are compared with the shipped partition (adjusted Rand index).
YncD is not in Prokka's database; its K-12 reference was downloaded separately (config/references.yaml).   python scripts/run_omp.py"""
import os, re, subprocess
from collections import defaultdict
import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = os.path.join(os.path.dirname(__file__), "..")
OUT = f"{ROOT}/replication/T23_omp"
ANN = f"{ROOT}/replication/T14_panacota/annot"
SPROT = os.path.expanduser("~/micromamba-phage/envs/pan/db/kingdom/Bacteria/sprot")
GENES = ["btuB", "fadL", "fhuA", "lamB", "lptD", "nfrA", "ompA", "ompC", "ompF", "tolC", "tsx", "yncD"]
SHIP = pd.read_csv(f"{ROOT}/data/genomics/bacteria/outer_membrane_proteins/blast_results_cured_clusters=99_wide.tsv", sep="\t", index_col=0)
env = dict(os.environ, MAMBA_ROOT_PREFIX=os.path.expanduser("~/micromamba-phage"))
run = lambda e, *cmd: subprocess.run(["micromamba", "run", "-n", e, *cmd], env=env, check=True, capture_output=True, text=True)
os.makedirs(OUT, exist_ok=True)

# 1. queries
q = []
for h, seq in re.findall(r"^>(.*)\n([^>]*)", open(SPROT).read(), flags=re.M):
    m = re.search(r"~~~([A-Za-z0-9]+)~~~", h)
    if m and m.group(1) in GENES:
        q.append(f">{m.group(1)}|{h.split()[0]}\n{seq.replace(chr(10), '')}\n")
# YncD is not in Prokka's database: use the K-12 entry downloaded into replication/T19_lps/refs (UniProt P76115, gene pqqU / former yncD)
q.append(">yncD|P76115\n" + "".join(open(f"{ROOT}/replication/T19_lps/refs/yncD.fasta").read().split("\n", 1)[1].split()) + "\n")
open(f"{OUT}/queries.faa", "w").write("".join(q))
print("queries:", len(q))

# 2. BLASTP against all proteins (concatenated by the pangenome step)
if not os.path.exists(f"{OUT}/hits.tsv"):
    db = f"{ANN}/Proteins/ESCO.All.prt"
    run("blast", "makeblastdb", "-in", db, "-dbtype", "prot", "-out", f"{OUT}/db")
    run("blast", "blastp", "-query", f"{OUT}/queries.faa", "-db", f"{OUT}/db", "-evalue", "1e-30", "-num_threads", "2",
        "-outfmt", "6 qseqid sseqid pident qcovs bitscore slen", "-out", f"{OUT}/hits.tsv")
h = pd.read_csv(f"{OUT}/hits.tsv", sep="\t", names=["q", "s", "pident", "qcov", "bits", "slen"])
h = h[(h.qcov >= 80) & (h.pident >= 50)]
h["gene"] = h.q.str.split("|").str[0]
h["genome"] = h.s.str.extract(r"^(ESCO\.\d+\.\d+)\.")[0]
best = h.sort_values("bits", ascending=False).drop_duplicates("s")           # each protein -> its best query (separates OmpC/OmpF/PhoE)
best = best.sort_values("bits", ascending=False).drop_duplicates(["genome", "gene"])  # one protein per genome and gene
best.to_csv(f"{OUT}/best_hits.tsv", sep="\t", index=False)

# 3. nucleotide sequences of the chosen genes, cluster at 99% identity
names = pd.read_csv(f"{ANN}/LSTINFO-list.lst", sep="\t").assign(strain=lambda d: d.orig_name.str.replace(r"\.fasta$", "", regex=True)).set_index("gembase_name").strain
res = {}
for g in GENES:
    sub = best[best.gene == g]
    want = set(sub.s)
    with open(f"{OUT}/{g}.fna", "w") as f:
        for gem in sub.genome.unique():
            seqs = dict(re.findall(r"^>(\S+).*\n([^>]*)", open(f"{ANN}/Genes/{gem}.gen").read(), flags=re.M))
            for pid in want & set(seqs):
                f.write(f">{pid}\n{seqs[pid].replace(chr(10), '')}\n")
    if not os.path.exists(f"{OUT}/{g}_cluster.tsv"):
        run("pan", "mmseqs", "easy-cluster", f"{OUT}/{g}.fna", f"{OUT}/{g}", f"{OUT}/tmp_{g}", "--min-seq-id", "0.99", "-c", "0.8", "--cov-mode", "1", "-v", "1")
    cl = pd.read_csv(f"{OUT}/{g}_cluster.tsv", sep="\t", names=["rep", "member"])
    cl["genome"] = cl.member.str.extract(r"^(ESCO\.\d+\.\d+)\.")[0]
    res[g] = cl.set_index("genome")["rep"].rename(g)
    subprocess.run(["rm", "-rf", f"{OUT}/tmp_{g}"])

ours = pd.DataFrame(res)
ours.index = ours.index.map(names)
ours.rename_axis("bacteria").to_csv(f"{OUT}/omp_clusters_regenerated.tsv", sep="\t")
rows = []
for g in GENES:
    col = g.upper()
    s = SHIP[col].astype("string") if col in SHIP else None
    if s is None:
        continue
    j = pd.concat([ours[g].rename("o"), s.rename("s")], axis=1)
    j["o"] = j.o.fillna("none"); j["s"] = j.s.fillna("none")
    rows.append(dict(gene=g, n_strains=len(j), ours_clusters=j.o.nunique(), shipped_clusters=j.s.nunique(), ARI=round(adjusted_rand_score(j.s, j.o), 3),
                     genome_with_gene_ours=int((j.o != "none").sum()), genome_with_gene_shipped=int((j.s != "none").sum())))
os.makedirs(f"{ROOT}/results/replication", exist_ok=True)
pd.DataFrame(rows).to_csv(f"{ROOT}/results/replication/omp_agreement.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
