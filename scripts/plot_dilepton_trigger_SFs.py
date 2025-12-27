import json
import numpy as np
import matplotlib.pyplot as plt

for channel in ['ee', 'mm', 'mixed']:
    sf = json.load(open(f"data/sf_{channel}_trg_lepton0_pt-trg_lepton1_pt-trig_idsV4.json"))
    x_bins = sf["corrections"][0]["data"]["edges"][0]
    y_bins = sf["corrections"][0]["data"]["edges"][1]
    values = sf["corrections"][0]["data"]["content"]

    # reshape flat list of values into 2D array with shape (ny, nx)
    nx = len(x_bins) - 1
    ny = len(y_bins) - 1
    sf_ee_array = np.array(values).reshape((ny, nx))

    fig = plt.figure(figsize=(8,6))

    plt.pcolormesh(x_bins, y_bins, sf_ee_array, cmap='viridis', shading='flat')
    plt.colorbar(label='Scale Factor')
    plt.xlabel('Leading Lepton pT [GeV]')
    plt.ylabel('Subleading Lepton pT [GeV]')
    ch = "El-El" if channel == 'ee' else "Mu-Mu" if channel == 'mm' else "Mixed (El-Mu)" if channel == 'mixed' else None
    plt.title(f'Trigger Scale Factors for {ch} Channel')
    plt.grid(False)

    # Add text annotations for each bin
    for i in range(ny):
        for j in range(nx):
            x_center = (x_bins[j] + x_bins[j+1]) / 2
            y_center = (y_bins[i] + y_bins[i+1]) / 2
            # plt.text(x_center, y_center, f'{sf_ee_array[i, j]:.3f}', ha='center', va='center', fontsize=6, color='white')
    plt.savefig(f"dilepton_trg_sf_{channel}.png", dpi=300)
