"""Reader for the authors' phage genomes (data/genomics/phages/FNA). Some files hold an empty first record or comment lines; the longest record is the genome."""
import os

FNA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "genomics", "phages", "FNA")


def read_genome(phage):
    recs, name = {}, None
    for line in open(f"{FNA}/{phage}.fna"):
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            name = line[1:].split()[0]; recs.setdefault(name, []).append("")
        else:
            recs[name][-1] += line.upper()
    seq = max((s for v in recs.values() for s in v), key=len)
    return seq


def quirks():
    out = []
    for f in sorted(os.listdir(FNA)):
        if f.startswith("._") or not f.endswith(".fna"):
            continue
        txt = open(f"{FNA}/{f}").read()
        n = txt.count(">")
        if n != 1 or any(l.startswith(";") for l in txt.split("\n")):
            out.append((f, n, any(l.startswith(";") for l in txt.split("\n"))))
    return out
