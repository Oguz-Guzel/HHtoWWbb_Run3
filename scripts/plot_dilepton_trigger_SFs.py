import json
import numpy as np
import mplhep as hep
import matplotlib.pyplot as plt

hep.style.use("CMS")
plt.rcParams["figure.dpi"] = 400
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['TeX Gyre Heros']

for channel in ['ee', 'mm', 'mixed']:
    sf = json.load(open(f"data/sf_{channel}_trg_lepton0_pt-trg_lepton1_pt-trig_idsV4.json"))
    x_bins = sf["corrections"][0]["data"]["edges"][0]
    y_bins = sf["corrections"][0]["data"]["edges"][1]
    values = sf["corrections"][0]["data"]["content"]

    # reshape flat list of values into 2D array with shape (ny, nx)
    nx = len(x_bins) - 1
    ny = len(y_bins) - 1
    sf_ee_array = np.array(values).reshape((ny, nx))

    fig, ax = plt.subplots(figsize=(12,10))

    mesh = ax.pcolormesh(x_bins, y_bins, sf_ee_array, cmap='viridis', shading='flat')
    fig.colorbar(mesh, ax=ax, label='Scale Factor')
    ax.set_xlabel('Leading Lepton pT [GeV]', fontsize=20)
    ax.set_ylabel('Subleading Lepton pT [GeV]', fontsize=20)
    ch = "El-El" if channel == 'ee' else "Mu-Mu" if channel == 'mm' else "Mixed (El-Mu)" if channel == 'mixed' else None
    ax.text(
        0, 1.05,
        "Private work (CMS simulation)",
        fontsize=22,
        verticalalignment='top',
        fontproperties="Tex Gyre Heros:italic",
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.text(
        0.83, 1.05,
        "(13.6 TeV)",
        fontsize=22,
        verticalalignment='top',
        fontproperties="Tex Gyre Heros:regular",
        transform=ax.transAxes,
        clip_on=False,
    )

    plt.grid(False)

    # Add text annotations for each bin
    for i in range(ny):
        for j in range(nx):
            x_center = (x_bins[j] + x_bins[j+1]) / 2
            y_center = (y_bins[i] + y_bins[i+1]) / 2
            # plt.text(x_center, y_center, f'{sf_ee_array[i, j]:.3f}', ha='center', va='center', fontsize=6, color='white')

    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.95, top=0.95)
    plt.savefig(f"dilepton_trg_sf_{channel}.png", dpi=300)
