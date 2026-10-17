"""Fig 5b: pooled ROC and precision-recall curves (best model per phage, all 96 phages), ours vs the authors' shipped predictions."""
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve
import predict as P
from figlib import SURF, INK, MUTED, BLUE, ORANGE, style

ours = pd.read_csv(f"{P.ROOT}/results/repro/predictions.csv")
pub = P.published()[1]
SERIES = [("Re-run", ours, BLUE), ("Published predictions", pub, ORANGE)]  # palette slots 1-2

fig, ax = plt.subplots(1, 2, figsize=(10, 4.6), facecolor=SURF)
for name, d, c in SERIES:
    fpr, tpr, _ = roc_curve(d.y, d.p)
    pr, rc, _ = precision_recall_curve(d.y, d.p)
    ax[0].plot(fpr, tpr, color=c, lw=2, label=f"{name} (AUROC {roc_auc_score(d.y, d.p):.3f})")
    ax[1].plot(rc, pr, color=c, lw=2, label=f"{name} (AP {average_precision_score(d.y, d.p):.3f})")
ax[0].plot([0, 1], [0, 1], color=MUTED, lw=1, ls=(0, (4, 3)), label="Chance")
prev = ours.y.mean()
ax[1].axhline(prev, color=MUTED, lw=1, ls=(0, (4, 3)), label=f"Chance (prevalence {prev:.2f})")
for a, title, xl, yl in [(ax[0], "ROC", "False positive rate", "True positive rate"), (ax[1], "Precision-recall", "Recall", "Precision")]:
    a.set_title(title, color=INK, loc="left", fontsize=12)
    a.set_xlabel(xl, color=MUTED); a.set_ylabel(yl, color=MUTED)
    a.set_xlim(0, 1); a.set_ylim(0, 1.02)
    style(a, "both")
    a.legend(frameon=False, labelcolor=INK, fontsize=9, loc="lower right" if title == "ROC" else "upper right")
fig.suptitle("Pooled predictions, 96 phages, group 10-fold CV (best model per phage)", color=INK, x=0.01, ha="left", fontsize=11)
fig.tight_layout()
out = f"{P.ROOT}/results/repro/fig5_curves.png"
fig.savefig(out, dpi=160, facecolor=SURF)
print(out)
