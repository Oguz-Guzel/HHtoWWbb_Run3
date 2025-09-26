import os
import ROOT
import json
import shutil
import argparse

# Run ROOT in batch mode (no GUI) to avoid display issues
ROOT.gROOT.SetBatch(True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate 2D trigger scale factors and export to JSON"
    )
    parser.add_argument(
        "--bamboo_output",
        required=True,
        help="Path to bamboo directory (don't point to results dir, it will be added automatically.)",
    )
    return parser.parse_args()


def mean2D(h, nonempty_only=True):
    """Compute mean of 2D histogram bin contents.
    If nonempty_only=True, use bins with content>0 or error>0."""
    s = 0.0
    n = 0
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            c = h.GetBinContent(i, j)
            if nonempty_only:
                if c > 0.0 or h.GetBinError(i, j) > 0.0:
                    s += c
                    n += 1
            else:
                s += c
                n += 1
    return (s / n if n else 0.0, n)


def calculate_2d_scale_factors_and_export(
    bamboo_output,
    data_files,
    mc_files,
    output_json_name="trigger_scale_factors.json",
):
    """
    Calculate 2D trigger scale factors from data and MC ROOT histograms

    Parameters:
    - data_files: list of paths to data ROOT files
    - mc_files: list of paths to MC ROOT files
    - output_json_name: output JSON filename
    """
    print(f"\nProcessing {len(data_files)} data files and {len(mc_files)} MC files...")

    # Helper to sum histograms across files for a given pair of names
    def sum_histograms(files, hnum_name, hden_name):
        h_num = None
        h_den = None
        for fpath in files:
            tf = ROOT.TFile.Open(fpath, "READ")
            if not tf:
                print(f"   Warning: Could not open {fpath}")
                continue
            curr_num = tf.Get(hnum_name)
            curr_den = tf.Get(hden_name)
            if not curr_num or not curr_den:
                print(
                    f"   Warning: Missing '{hnum_name}' or '{hden_name}' in {os.path.basename(fpath)}"
                )
                tf.Close()
                continue
            if h_num is None:
                h_num = curr_num.Clone(f"{hnum_name}_combined")
                h_den = curr_den.Clone(f"{hden_name}_combined")
                h_num.SetDirectory(0)
                h_den.SetDirectory(0)
            else:
                h_num.Add(curr_num)
                h_den.Add(curr_den)
            tf.Close()
        return h_num, h_den

    # Process each channel
    results = {}  # ch -> dict with histograms and arrays
    x_edges = None
    y_edges = None
    first_channel_for_return = None

    min_denominator = 1e-6  # avoid division by zero

    # get numerator and denominator histogram names from the first data file
    first_file = ROOT.TFile.Open(data_files[0])
    num_names = [
        hist.GetName()
        for hist in first_file.GetListOfKeys()
        if hist.GetName().startswith("num")
    ]
    den_names = [den.replace("num_", "den_") for den in num_names]

    # Derive channel labels from numerator names: "num_ee" -> "ee"
    channels = [n.replace("num_", "") for n in num_names]

    for idx, ch in enumerate(channels):
        nname = num_names[idx]
        dname = den_names[idx]
        # print(f"\n-- Channel '{ch}' using histograms: {nname} / {dname}")

        # Sum data and MC histograms over files
        h_data_num, h_data_den = sum_histograms(data_files, nname, dname)
        h_mc_num, h_mc_den = sum_histograms(mc_files, nname, dname)

        # Stabilize denominators and sanitize MC negative bins
        neg_mc_bins = 0
        for i in range(1, h_data_den.GetNbinsX() + 1):
            for j in range(1, h_data_den.GetNbinsY() + 1):
                # Clamp denominators to avoid 0 division
                if h_data_den.GetBinContent(i, j) < min_denominator:
                    h_data_den.SetBinContent(i, j, min_denominator)
                if h_mc_den.GetBinContent(i, j) < min_denominator:
                    h_mc_den.SetBinContent(i, j, min_denominator)
                # Zero-out negative MC numerator bins (from negative weights)
                if h_mc_num.GetBinContent(i, j) < 0:
                    h_mc_num.SetBinContent(i, j, 0.0)
                    neg_mc_bins += 1
        # if neg_mc_bins:
        #     print(f"  [{ch}] MC negative numerator bins set to 0: {neg_mc_bins}")

        # Efficiencies
        h_eff_data = h_data_num.Clone(f"h_efficiency_data_{ch}")
        # Binomial stats OK for data (integer counts)
        h_eff_data.Divide(h_data_num, h_data_den, 1, 1, "B")

        h_eff_mc = h_mc_num.Clone(f"h_efficiency_mc_{ch}")
        # Using standard ratio for MC (weights/negative weights incompatible with 'B')
        # h_mc_num.Sumw2(True)
        # h_mc_den.Sumw2(True)
        # h_eff_mc.Sumw2(True)
        h_eff_mc.Divide(h_mc_num, h_mc_den, 1.0, 1.0, "")  # no 'B'

        # Clip MC efficiency to physical range [0,1]
        clipped = 0
        for i in range(1, h_eff_mc.GetNbinsX() + 1):
            for j in range(1, h_eff_mc.GetNbinsY() + 1):
                v = h_eff_mc.GetBinContent(i, j)
                if v < 0:
                    h_eff_mc.SetBinContent(i, j, 0.0)
                    clipped += 1
                elif v > 1:
                    h_eff_mc.SetBinContent(i, j, 1.0)
                    clipped += 1
        # if clipped:
        # print(f"  [{ch}] MC efficiency bins clipped to [0,1]: {clipped}")

        # Scale factors
        h_sf = h_eff_data.Clone(f"h_scale_factors_{ch}")
        h_sf.Divide(h_eff_data, h_eff_mc, 1, 1)

        # Print means for this channel
        mean_data, n_data = mean2D(h_eff_data)
        mean_mc, n_mc = mean2D(h_eff_mc)
        mean_sf, n_sf = mean2D(h_sf)
        print(
            f"  [{ch}] Means: \n"
            f"         Data eff = {mean_data:.4f} (N={n_data}), \n"
            f"         MC eff = {mean_mc:.4f} (N={n_mc}), \n"
            f"         SF = {mean_sf:.4f} (N={n_sf})"
        )

        nbx = h_sf.GetNbinsX()
        nby = h_sf.GetNbinsY()

        # Set common bin edges (assume consistent binning across channels)
        if x_edges is None or y_edges is None:
            x_edges = [h_sf.GetXaxis().GetBinLowEdge(i) for i in range(1, nbx + 1)]
            x_edges.append(h_sf.GetXaxis().GetBinUpEdge(nbx))
            y_edges = [h_sf.GetYaxis().GetBinLowEdge(j) for j in range(1, nby + 1)]
            y_edges.append(h_sf.GetYaxis().GetBinUpEdge(nby))

        # Extract values and errors
        vals = []
        errs = []
        suspicious = 0
        for i in range(1, nbx + 1):
            row_v = []
            row_e = []
            for j in range(1, nby + 1):
                v = h_sf.GetBinContent(i, j)
                e = h_sf.GetBinError(i, j)
                row_v.append(float(v))
                row_e.append(float(e))
            vals.append(row_v)
            errs.append(row_e)
        if suspicious > 0:
            print(f"  [{ch}] Suspicious bins replaced: {suspicious}")

        results[ch] = dict(
            h_data_num=h_data_num,
            h_data_den=h_data_den,
            h_mc_num=h_mc_num,
            h_mc_den=h_mc_den,
            h_eff_data=h_eff_data,
            h_eff_mc=h_eff_mc,
            h_sf=h_sf,
            values=vals,
            errors=errs,
        )
        if first_channel_for_return is None:
            first_channel_for_return = ch

    # Build CorrectionLib JSON with a 'channel' category for nominal and uncertainties
    def make_y_node(values_row):
        return {
            "nodetype": "binning",
            "input": "pt_subleading",
            "edges": y_edges,
            "flow": "clamp",
            "content": [float(v) for v in values_row],
        }

    def make_xy_node(matrix_rows):
        return {
            "nodetype": "binning",
            "input": "pt_leading",
            "edges": x_edges,
            "flow": "clamp",
            "content": [make_y_node(row) for row in matrix_rows],
        }

    correction_json = {
        "schema_version": 2,
        "corrections": [
            {
                "name": "trigger_scale_factors_2d",
                "version": 1,
                "description": "2D trigger scale factors by channel with systematic variations",
                "inputs": [
                    {
                        "name": "channel",
                        "type": "string",
                        "description": "final state",
                    },
                    {
                        "name": "systematic",
                        "type": "string",
                        "description": "Systematic variation (nominal, up, down)",
                    },
                    {
                        "name": "pt_leading",
                        "type": "real",
                        "description": "Leading lepton pT [GeV]",
                    },
                    {
                        "name": "pt_subleading",
                        "type": "real",
                        "description": "Subleading lepton pT [GeV]",
                    },
                ],
                "output": {
                    "name": "scale_factor",
                    "type": "real",
                    "description": "Trigger efficiency SF (data/MC)",
                },
                "data": {
                    "nodetype": "category",
                    "input": "channel",
                    "content": [
                        {
                            "key": ch,
                            "value": {
                                "nodetype": "category",
                                "input": "systematic",
                                "content": [
                                    {
                                        "key": "nominal",
                                        "value": make_xy_node(results[ch]["values"]),
                                    },
                                    {
                                        "key": "up",
                                        "value": make_xy_node(
                                            [
                                                [
                                                    min(v + e, 1.5)
                                                    for v, e in zip(row_v, row_e)
                                                ]
                                                for row_v, row_e in zip(
                                                    results[ch]["values"],
                                                    results[ch]["errors"],
                                                )
                                            ]
                                        ),
                                    },
                                    {
                                        "key": "down",
                                        "value": make_xy_node(
                                            [
                                                [
                                                    max(v - e, 0.0)
                                                    for v, e in zip(row_v, row_e)
                                                ]
                                                for row_v, row_e in zip(
                                                    results[ch]["values"],
                                                    results[ch]["errors"],
                                                )
                                            ]
                                        ),
                                    },
                                ],
                            },
                        }
                        for ch in channels
                    ],
                },
            }
        ],
    }

    # Write JSON
    with open(output_json_name, "w") as f:
        json.dump(correction_json, f, indent=2)

    shutil.move(output_json_name, os.path.join(bamboo_output, output_json_name))
    print(
        f"\n  Scale factors JSON saved as: {os.path.join(bamboo_output, output_json_name)}"
    )

    # Return all histos for any downstream plotting
    h_sf_plt = []
    h_eff_data_plt = []
    h_eff_mc_plt = []
    for ch in channels:
        h_sf_plt.append(results[ch]["h_sf"])
        h_eff_data_plt.append(results[ch]["h_eff_data"])
        h_eff_mc_plt.append(results[ch]["h_eff_mc"])
    return h_sf_plt, h_eff_data_plt, h_eff_mc_plt


def create_comparison_plots(
    bamboo_output, h_data_list, h_mc_list, h_sf_list, output_name="comparison.png"
):
    """Create comparison plots for data efficiency, MC efficiency, and scale factors"""
    if not all([h_data_list, h_mc_list, h_sf_list]):
        print("Warning: Empty histogram list provided to create_comparison_plots.")
        return

    # Clone the first histogram to have a base for summation
    h_data = h_data_list[0].Clone("h_data_combined")
    h_mc = h_mc_list[0].Clone("h_mc_combined")
    h_sf = h_sf_list[0].Clone("h_sf_combined")

    # Add the rest of the histograms in the lists
    for i in range(1, len(h_data_list)):
        h_data.Add(h_data_list[i])
    for i in range(1, len(h_mc_list)):
        h_mc.Add(h_mc_list[i])
    for i in range(1, len(h_sf_list)):
        h_sf.Add(h_sf_list[i])

    canvas = ROOT.TCanvas("c_comparison", "Trigger Efficiency Comparison", 1200, 800)
    canvas.Divide(3, 1)

    canvas.cd(1)
    h_data.SetTitle(
        "Data Efficiency;Leading Lepton p_{T} [GeV];Subleading Lepton p_{T} [GeV]"
    )
    h_data.SetMinimum(0.0)
    h_data.SetMaximum(1.0)
    h_data.SetStats(0)
    h_data.Draw("COLZ")

    canvas.cd(2)
    h_mc.SetTitle(
        "MC Efficiency;Leading Lepton p_{T} [GeV];Subleading Lepton p_{T} [GeV]"
    )
    h_mc.SetMinimum(0.0)
    h_mc.SetMaximum(1.0)
    h_mc.SetStats(0)
    h_mc.Draw("COLZ")

    canvas.cd(3)
    h_sf.SetTitle(
        "Scale Factors (Data/MC);Leading Lepton p_{T} [GeV];Subleading Lepton p_{T} [GeV]"
    )
    h_sf.SetMinimum(0.5)
    h_sf.SetMaximum(1.5)
    h_sf.SetStats(0)
    h_sf.Draw("COLZ")

    canvas.SaveAs(output_name)

    shutil.move(output_name, os.path.join(bamboo_output, output_name))
    print(f"\n  Comparison plots saved as: {os.path.join(bamboo_output, output_name)}")


# Example usage
if __name__ == "__main__":

    args = parse_args()
    bamboo_output_dir = args.bamboo_output
    bamboo_results_dir = os.path.join(bamboo_output_dir, "results")

    # List ROOT files in the bamboo results directory
    data_files = [
        f
        for f in os.listdir(bamboo_results_dir)
        if (f.startswith("Muon") or f.startswith("EGamma")) and f.endswith(".root")
    ]
    mc_files = [
        f
        for f in os.listdir(bamboo_results_dir)
        if f not in data_files
        and f[:2] != "__"
        and not f.startswith("ggH")
        and f.endswith(".root")
    ]
    # add bamboo_results_dir to file paths
    data_files = [os.path.join(bamboo_results_dir, f) for f in data_files]
    mc_files = [os.path.join(bamboo_results_dir, f) for f in mc_files]

    h_scale_factors, h_eff_data, h_eff_mc = calculate_2d_scale_factors_and_export(
        bamboo_output_dir,
        data_files,
        mc_files,
        "di_lepton_trigger_scale_factors.json",
    )

    create_comparison_plots(
        bamboo_output_dir,
        h_eff_data,
        h_eff_mc,
        h_scale_factors,
        "trigger_comparison.png",
    )
