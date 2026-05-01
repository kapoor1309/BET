"""Density plot of side-chain angles in tripeptide (XRX) helices,
grouped by size class of the left neighbour.

Mirrors plot_angles.R but uses matplotlib so the pipeline works
without an R installation.
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

INPUT = sys.argv[1] if len(sys.argv) > 1 else "final/angles.tsv"
OUTPUT = sys.argv[2] if len(sys.argv) > 2 else "final/ss_profile_HHH_for_arg_with_valid_runs.png"

LEVELS = ["Tiny", "Small", "Intermediate", "Large", "Bulky"]
COLORS = {
    "Tiny": "#e6e2d3",
    "Small": "#f2c879",
    "Intermediate": "#f28e5c",
    "Large": "#e4572e",
    "Bulky": "#b11226",
}

df = pd.read_csv(INPUT, sep="\t")

df["angle"] = ((df["angle"] + 180) % 360) - 180
df = df[df["size_class"].isin(LEVELS)]

n = len(df)

fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor("#bdbdbd")
ax.set_facecolor("#bdbdbd")

x_grid = np.linspace(-180, 180, 1000)

for level in LEVELS:
    sub = df[df["size_class"] == level]["angle"].values
    if len(sub) < 2:
        continue
    kde = gaussian_kde(sub, bw_method="silverman")
    y = kde(x_grid)
    ax.plot(x_grid, y, color=COLORS[level], linewidth=1.6, label=level)

ax.set_xlim(-180, 180)
ax.set_xticks(np.arange(-180, 181, 50))

for x in np.arange(-180, 181, 50):
    ax.axvline(x, color="blue", linestyle=":", linewidth=0.5, alpha=0.7)

ax.set_xlabel(r"Angle between adjacent C-$\alpha$ $\to$ Centroid vectors [°]")
ax.set_ylabel("Norm. Freq. [A.U.]")
ax.set_title(f"Tripeptide (XRX) in Helix (n = {n})")

leg = ax.legend(loc="upper left", facecolor="white", edgecolor="black")
leg.get_frame().set_linewidth(0.8)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

os.makedirs(os.path.dirname(OUTPUT) or ".", exist_ok=True)
plt.tight_layout()
plt.savefig(OUTPUT, dpi=300, facecolor=fig.get_facecolor())
print(f"Saved: {OUTPUT}")
print(f"n = {n}")
