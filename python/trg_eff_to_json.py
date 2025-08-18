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
        "--bamboo_results", required=True, help="Path to bamboo results directory"
    )
    return parser.parse_args()


def calculate_2d_scale_factors_and_export(
    bamboo_results,
    data_files,
    mc_files,
    num_hist_name,
    den_hist_name,
    output_json_name="trigger_scale_factors.json",
):
    """
    Calculate 2D trigger scale factors from data and MC ROOT histograms

    Parameters:
    - data_files: list of paths to data ROOT files
    - mc_files: list of paths to MC ROOT files
    - num_hist_name: name(s) of numerator histogram(s) (str or list[str], e.g. ["num_ee","num_mumu","num_emu"])
    - den_hist_name: name(s) of denominator histogram(s) (str or list[str]), same as num_hist_name
    - output_json_name: output JSON filename
    """
    print(f"\nProcessing {len(data_files)} data files and {len(mc_files)} MC files...")

    # Normalize inputs to lists of equal length
    if isinstance(num_hist_name, str):
        num_names = [num_hist_name]
    else:
        num_names = list(num_hist_name)

    # Derive channel labels from numerator names: "num_ee" -> "ee"
    channels = [n.replace("num_", "") for n in num_names]

    if isinstance(den_hist_name, list):
        den_names = list(den_hist_name)
    else:
        den_names = list(den_hist_name)

    if len(den_names) != len(num_names):
        raise ValueError("num_hist_name and den_hist_name must have the same length")

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

    for idx, ch in enumerate(channels):
        nname = num_names[idx]
        dname = den_names[idx]
        print(f"\n-- Channel '{ch}' using histograms: {nname} / {dname}")

        # Sum data and MC histograms over files
        h_data_num, h_data_den = sum_histograms(data_files, nname, dname)
        h_mc_num, h_mc_den = sum_histograms(mc_files, nname, dname)

        if not h_data_num or not h_data_den:
            raise ValueError(
                f"[{ch}] Could not find required histograms in any DATA files"
            )
        if not h_mc_num or not h_mc_den:
            raise ValueError(
                f"[{ch}] Could not find required histograms in any MC files"
            )

        # Stabilize denominators
        for i in range(1, h_data_den.GetNbinsX() + 1):
            for j in range(1, h_data_den.GetNbinsY() + 1):
                if h_data_den.GetBinContent(i, j) < min_denominator:
                    h_data_den.SetBinContent(i, j, min_denominator)
                if h_mc_den.GetBinContent(i, j) < min_denominator:
                    h_mc_den.SetBinContent(i, j, min_denominator)

        # Efficiencies
        h_eff_data = h_data_num.Clone(f"h_efficiency_data_{ch}")
        h_eff_data.Divide(h_data_num, h_data_den, 1, 1, "B")
        h_eff_mc = h_mc_num.Clone(f"h_efficiency_mc_{ch}")
        h_eff_mc.Divide(h_mc_num, h_mc_den, 1, 1, "B")

        # Scale factors
        h_sf = h_eff_data.Clone(f"h_scale_factors_{ch}")
        h_sf.Divide(h_eff_data, h_eff_mc, 1, 1)

        print(f"  [{ch}] Data eff mean: {h_eff_data.GetMean():.4f}")
        print(f"  [{ch}] MC eff mean:   {h_eff_mc.GetMean():.4f}")
        print(f"  [{ch}] SF mean:       {h_sf.GetMean():.4f}")

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
                if v <= 0 or v > 1.5 or e < 0 or e > 0.5:
                    suspicious += 1
                    v = 1.0
                    e = 0.0
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
                "description": "2D trigger scale factors by channel (leading/subleading lepton pT)",
                "inputs": [
                    {
                        "name": "channel",
                        "type": "string",
                        "description": "ee, mumu, emu",
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
                        {"key": ch, "value": make_xy_node(results[ch]["values"])}
                        for ch in channels
                    ],
                },
            },
            {
                "name": "trigger_scale_factors_2d_unc",
                "version": 1,
                "description": "Absolute uncertainty for 2D trigger scale factors by channel",
                "inputs": [
                    {
                        "name": "channel",
                        "type": "string",
                        "description": "ee, mumu, emu",
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
                    "name": "uncertainty",
                    "type": "real",
                    "description": "Absolute uncertainty on SF",
                },
                "data": {
                    "nodetype": "category",
                    "input": "channel",
                    "content": [
                        {"key": ch, "value": make_xy_node(results[ch]["errors"])}
                        for ch in channels
                    ],
                },
            },
        ],
    }

    # Write JSON
    with open(output_json_name, "w") as f:
        json.dump(correction_json, f, indent=2)

    shutil.move(output_json_name, os.path.join(bamboo_results, "..", output_json_name))
    print(f"\n  Scale factors JSON saved as: {output_json_name} in {bamboo_results}")

    # Return first channel histos for any downstream plotting that expects single hists
    ch0 = first_channel_for_return or channels[0]
    return results[ch0]["h_sf"], results[ch0]["h_eff_data"], results[ch0]["h_eff_mc"]


def create_comparison_plots(h_data, h_mc, h_sf, output_name="comparison.png"):
    """Create comparison plots for data efficiency, MC efficiency, and scale factors"""
    canvas = ROOT.TCanvas("c_comparison", "Trigger Efficiency Comparison", 1200, 800)
    canvas.Divide(3, 1)

    canvas.cd(1)
    h_data.SetTitle(
        "Data Efficiency;Leading Lepton p_{T} [GeV];Subleading Lepton p_{T} [GeV]"
    )
    h_data.SetMinimum(0.0)
    h_data.SetMaximum(1.0)
    h_data.SetStats(0)
    h_data.Draw("COLZ")  # no TEXT for speed

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

    shutil.move(output_name, os.path.join(bamboo_results, "..", output_name))
    print(f"\n  Comparison plots saved as: {output_name} in {bamboo_results}")


# Example usage
if __name__ == "__main__":

    args = parse_args()
    bamboo_results = args.bamboo_results

    # List ROOT files in the bamboo results directory
    data_files = [
        f
        for f in os.listdir(bamboo_results)
        if f.startswith("Muon") or f.startswith("EGamma")
    ]
    mc_files = [
        f for f in os.listdir(bamboo_results) if f not in data_files and f[:2] != "__"
    ]
    # add bamboo_results to file paths
    data_files = [os.path.join(bamboo_results, f) for f in data_files]
    mc_files = [os.path.join(bamboo_results, f) for f in mc_files]

    # Histogram names (should be same in both files)
    num_hist = ["num_ee", "num_mumu", "num_emu"]  # Numerator histogram names
    den_hist = ["den_ee", "den_mumu", "den_emu"]  # Denominator histogram names

    h_scale_factors, h_eff_data, h_eff_mc = calculate_2d_scale_factors_and_export(
        bamboo_results,
        data_files,
        mc_files,
        num_hist,
        den_hist,
        "trigger_scale_factors_run3.json",
    )

    create_comparison_plots(
        h_eff_data, h_eff_mc, h_scale_factors, "trigger_comparison.png"
    )
