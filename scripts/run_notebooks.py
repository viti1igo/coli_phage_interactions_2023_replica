"""Execute the four notebooks of the authors' repository on a scratch copy of it (their relative paths stay valid, the clone is not touched) and
record, per notebook, every cell that fails and why.
  python scripts/run_notebooks.py shipped    # as shipped
  python scripts/run_notebooks.py patched    # with the substitutions in PATCHES (public equivalents of private/Windows paths)
Outputs: notebooks/<mode>/<name>.ipynb (executed), notebooks/<mode>_errors.csv, files the notebooks wrote -> notebooks/<mode>_files/"""
import json, os, re, shutil, subprocess, sys, tempfile
import nbformat
from nbclient import NotebookClient
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
MODE = sys.argv[1] if len(sys.argv) > 1 else "shipped"
NBS = ["data/interactions/bokeh_draw_interaction_matrix.ipynb", "dev/predictions/visualize_predictions.ipynb",
       "dev/cocktails/analyze_cocktail_interactions_tailored_vs_generic_vs_baseline.ipynb",
       "data/genomics/bacteria/defense_finder/build_bacteria_arsenale_pairwise_distance_table.ipynb"]
# in 'patched' mode the (regex, replacement) pairs of notebook_patches.py are applied to the code cells
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from notebook_patches import PATCHES as _P
PATCHES = _P if MODE == "patched" else {}

work = tempfile.mkdtemp(prefix="nbrun_", dir=os.environ.get("NB_TMP", "/private/tmp"))
for d in ("data", "dev"):  # only the authors' directories (the notebooks read nothing else)
    subprocess.run(["rsync", "-a", "--exclude=assemblies", "--exclude=._*", "--exclude=*.rar", f"{ROOT}/{d}/", f"{work}/repo/{d}/"], check=True)
rows = []
out_nb = f"{ROOT}/notebooks/{MODE}"; os.makedirs(out_nb, exist_ok=True)
for rel in NBS:
    name = os.path.basename(rel)
    nb = nbformat.read(f"{work}/repo/{rel}", as_version=4)
    for c in nb.cells:
        if c.cell_type == "code":
            c.outputs, c.execution_count = [], None
            for pat, rep in PATCHES.get(name, []):
                c.source = re.sub(pat, rep, c.source)
    before = set(os.path.join(dp, f) for dp, _, fs in os.walk(f"{work}/repo") for f in fs)
    client = NotebookClient(nb, timeout=900, kernel_name="python3", allow_errors=True, resources={"metadata": {"path": os.path.dirname(f"{work}/repo/{rel}")}})
    client.execute()
    nbformat.write(nb, f"{out_nb}/{name}")
    n_code = sum(c.cell_type == "code" for c in nb.cells)
    for i, c in enumerate([c for c in nb.cells if c.cell_type == "code"]):
        for o in c.outputs:
            if o.output_type == "error":
                rows.append(dict(notebook=name, code_cell=i + 1, of=n_code, error=o.ename, message=re.sub(r"\s+", " ", o.evalue)[:220], first_line=c.source.strip().split("\n")[0][:90]))
    new = [f for dp, _, fs in os.walk(f"{work}/repo") for f in [os.path.join(dp, x) for x in fs] if f not in before and "/._" not in f]
    for f in new:
        d = f"{ROOT}/notebooks/{MODE}_files/{name[:-6]}"; os.makedirs(d, exist_ok=True); shutil.copy(f, d)
    print(f"{name}: {n_code} code cells, {sum(r['notebook'] == name for r in rows)} errors, {len(new)} new files")
pd.DataFrame(rows).to_csv(f"{ROOT}/notebooks/{MODE}_errors.csv", index=False)
shutil.rmtree(work, ignore_errors=True)
