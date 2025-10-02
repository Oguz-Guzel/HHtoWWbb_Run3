# This script generates a heatmap of the branching ratios of the
# Higgs boson decaying into various final states.
# The values are taken from
# https://e-publishing.cern.ch/index.php/CYRM/issue/view/32
# Author: Oguz Guzel oguz.guzel@uclouvain.be
# Date: 2025-05-16

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# branching ratios
bb = 5.809e-1
tautau = 6.256e-2
mumu = 2.171e-4
cc = 2.884e-2
gg = 8.180e-2
gammagamma = 2.270e-3
Zgamma = 1.541e-3
WW = 2.152e-1
ZZ = 2.641e-2

# List of branching ratios in the same order as labels
brs = [bb, WW, gg, tautau, cc, ZZ, gammagamma, Zgamma, mumu]

# Create the symmetric matrix with off-diagonal terms multiplied by 2
size = len(brs)
full_matrix = np.zeros((size, size))
for i in range(size):
    for j in range(i + 1):
        if i == j:
            full_matrix[i, j] = brs[i] * brs[j]
        else:
            full_matrix[i, j] = 2 * brs[i] * brs[j]

full_matrix = full_matrix
plt.rcParams["text.usetex"] = True

# Labels
labels = [
    "bb",
    "WW",
    "gg",
    r"$\tau\tau$",
    "cc",
    "ZZ",
    r"$\gamma\gamma$",
    r"Z$\gamma$",
    r"$\mu\mu$",
]
fig, ax = plt.subplots(figsize=(8, 6))

# Color scale (logarithmic)
norm = mcolors.LogNorm(vmin=1e-8, vmax=full_matrix.max())
cmap = plt.cm.Blues

# Mask upper triangle
mask = np.triu(np.ones_like(full_matrix, dtype=bool), k=1)
masked_data = np.ma.masked_array(full_matrix, mask)

# Plot
c = ax.imshow(masked_data, norm=norm, cmap=cmap)

# Add text annotations
for i in range(len(labels)):
    for j in range(i + 1):  # lower triangle + diagonal
        val = full_matrix[i, j]
        if val < 1e-2:
            mantissa, exp = f"{val:.2e}".split("e")
            exp = int(exp)
            text = rf"{mantissa}" + "\n" + rf"$\cdot 10^{{{exp}}}$"
        else:
            text = f"{val:.4f}"
        ax.text(
            j,
            i,
            text,
            ha="center",
            va="center",
            color=(
                "lime"
                if val >= 1e-2
                else ("white" if val < 1e-2 and val >= 1e-4 else "black")
            ),
            fontsize=12,
            font="serif",
            weight="bold",
        )

# Axes formatting
ax.set_xlabel(r"H$\rightarrow$ YY", fontsize=14, fontfamily="serif")
ax.set_ylabel(r"H$\rightarrow$ XX", fontsize=14, fontfamily="serif")
ax.set_xticks(np.arange(len(labels)))
ax.set_yticks(np.arange(len(labels)))
ax.set_xticklabels(labels, ha="right", fontsize=14, fontfamily="serif")
ax.set_yticklabels(labels, fontsize=14, fontfamily="serif")
ax.tick_params(top=False, bottom=True, left=True, right=False)
# Create a custom legend for the color scale
import matplotlib.patches as mpatches

legend_elements = [
    mpatches.Patch(facecolor="lime", label=r"$\mathcal{BR}\geq 10^{-2}$"),
    mpatches.Patch(
        facecolor="white",
        edgecolor="black",
        label=r"$10^{-2}>\mathcal{BR}>10^{-4}$"
    ),
    mpatches.Patch(facecolor="black", label=r"$\mathcal{BR}\leq 10^{-4}$"),
]

ax.legend(
    handles=legend_elements,
    loc="upper left",
    bbox_to_anchor=(0.5, 0.9),
    fontsize=12,
    title_fontsize=14,
    frameon=True,
)

# Color bar
cbar = fig.colorbar(c, ax=ax, fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=14)

# Add a custom text annotation on the plot
ax.text(
    0.95,
    0.9,
    r"$\mathcal{BR}$(HH$\rightarrow$ XXYY)",
    transform=ax.transAxes,
    fontsize=24,
    color="black",
    ha="right",
    va="bottom",
    fontfamily="serif",
)

plt.tight_layout()
plt.savefig("./br_hh_matrix.svg", format="svg")
plt.show()
