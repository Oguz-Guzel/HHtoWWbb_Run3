from bamboo import treefunctions as op


def fillSampleTemplate(template, selEras=None):
    from copy import deepcopy

    outTemplate = {}

    # Expand eras
    for name, sample in template.items():
        if "dbs" in sample:
            for era, das in sample["dbs"].items():
                era = str(era)
                if selEras is not None and era not in selEras:
                    continue
                thisSample = deepcopy(sample)
                if "syst" in thisSample:
                    syst, nom = thisSample["syst"]
                    newName = f"{nom}__{era}__{syst}"
                    thisSample["syst"][1] = f"{name}__{era}"
                else:
                    newName = f"{name}__{era}"
                thisSample["db"] = das
                thisSample["era"] = era
                thisSample.pop("dbs")
                outTemplate[newName] = thisSample
        else:
            outTemplate[name] = sample

    return outTemplate


def labeler(label):
    return {"labels": [{"text": label, "position": [0.16, 0.91], "size": 20}]}


def custom_Plotit(
    cfgName, workdir, inDir, outDir, counterReader, config, plotIt, verbose
):
    import collections
    import shutil
    import os
    import subprocess
    import logging

    logger = logging.getLogger(__name__)

    from bamboo.root import gbl as ROOT

    def openFileAndGet(path, mode="read"):
        """Open ROOT file in a mode, check if open properly, and return TFile handle."""
        tf = ROOT.TFile.Open(path, mode)
        if not tf or not tf.IsOpen():
            raise Exception("Could not open file {}".format(path))
        return tf

    to_hadd = collections.defaultdict(list)
    hadd_cfg = collections.defaultdict(dict)
    keep_cfg = collections.defaultdict(dict)

    _gp = list(config["plotIt"]["groups"].keys())
    if "signal" in _gp:
        _gp.remove("signal")

    lumiCFG = {}
    for smp, smpCfg in config["samples"].items():

        era = smpCfg["era"]
        lumi = config["eras"][smpCfg["era"]]["luminosity"]
        lumiCFG[era] = lumi

        smpNm = smp.split(era)[0]
        mergedHists = {}

        file_to_copy = os.path.join(outDir, f"{smp}.root")
        logger.info(f"Copying samples to {outDir}")
        if smpCfg.get("group") in _gp:
            keep_cfg[smp] = smpCfg
            if os.path.isfile(file_to_copy):
                continue
            shutil.copyfile(
                os.path.join(inDir, f"{smp}.root"), os.path.join(outDir, f"{smp}.root")
            )

        else:
            resultsFile = openFileAndGet(
                os.path.join(inDir, f"{smp}.root"), mode="READ"
            )
            logger.info("Stacking signal sample: {0}".format(smp))
            xsc = smpCfg["cross-section"]
            gevt = counterReader(resultsFile)[smpCfg["generated-events"]]
            br = smpCfg["branching-ratio"]
            smpScale = (lumi * xsc * br) / gevt

            for hk in resultsFile.GetListOfKeys():
                hist = hk.ReadObj()
                if not hist.InheritsFrom("TH1"):
                    continue
                hist.Scale(smpScale)
                name = hist.GetName()
                if name not in mergedHists.keys():
                    mergedHists[name] = hist.Clone()
                    mergedHists[name].SetDirectory(0)
                else:
                    mergedHists[name].Add(hist)
            resultsFile.Close()

            normalizedFile = openFileAndGet(
                os.path.join(outDir, f"{smp}.root"), "recreate"
            )
            for hist in mergedHists.values():
                hist.Write()
            normalizedFile.Close()

            to_hadd[smpNm].append(os.path.join(outDir, f"{smp}.root"))
            hadd_cfg[smpNm].update(smpCfg)

    for smp, val in to_hadd.items():
        sum_f = f"{smp}full.root"
        haddCmd = ["hadd", "-f", os.path.join(outDir, sum_f)] + val
        try:
            if verbose:
                logger.info("running {}".format(" ".join(haddCmd)))
                subprocess.check_call(haddCmd)
            else:
                subprocess.check_call(haddCmd, stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            logger.error("Failed to run {0}".format(" ".join(haddCmd)))

    with open(cfgName, "r") as inf:
        with open(os.path.join(workdir, "plots_full.yml"), "w+") as outf:
            outf.write("configuration:\n")
            outf.write("  blinded-range-fill-color: '#29556270'\n")
            outf.write("  blinded-range-fill-style: 1001\n")
            outf.write("  eras:\n")
            for era in lumiCFG.keys():
                outf.write(f"  - {era}\n")
            outf.write("  error-fill-color: '#ee556270'\n")
            outf.write("  error-fill-style: 3154\n")
            outf.write("  experiment: ' '\n")
            outf.write("  extra-label: Private work (CMS simulation)\n")
            outf.write("  height: 800\n")
            outf.write("  luminosity:\n")
            for era, lumi in lumiCFG.items():
                outf.write(f"    {era}: {lumi}\n")
            outf.write("  luminosity-label: '%1$.2f fb^{-1} (13.6 TeV)'\n")
            outf.write("  margin-bottom: 0.1\n")
            outf.write("  margin-left: 0.125\n")
            outf.write("  margin-right: 0.03\n")
            outf.write("  margin-top: 0.05\n")
            outf.write("  ratio-fit-error-fill-color: '#aa556270'\n")
            outf.write("  ratio-fit-error-fill-style: 1001\n")
            outf.write("  ratio-fit-line-color: '#0B486B'\n")
            outf.write(f"  root: {workdir}/results/normalizedSummedSignal\n")
            outf.write("  width: 800\n")
            outf.write("  yields-table-align: v\n")
            outf.write("  yields-table-text-align: l\n")
            outf.write("files:\n")
            line_found = False

            for smp, cfg in hadd_cfg.items():
                outf.write(f"  {smp}full.root:\n")
                # smpFile = openFileAndGet(
                # os.path.join(outDir, f"{smp}full.root"), mode="READ")
                # gevt = counterReader(smpFile)[smpCfg["generated-events"]]
                # outf.write(f"    generated-events: {gevt}\n")
                for k, v in cfg.items():
                    if k in [
                        "type",
                        "group",
                        "line-width",
                        "line-type",
                        "legend",
                        "line-color",
                        "branching-ratio",
                        "cross-section",
                        "era",
                    ]:
                        outf.write(f"    {k}: {v}\n")

            for smp, cfg in keep_cfg.items():
                outf.write(f"  {smp}.root:\n")
                if cfg["type"] != "data":
                    smpFile = openFileAndGet(
                        os.path.join(inDir, f"{smp}.root"), mode="READ"
                    )
                    gevt = counterReader(smpFile)[smpCfg["generated-events"]]
                    outf.write(f"    generated-events: {gevt}\n")
                else:
                    None
                for k, v in cfg.items():
                    if k in [
                        "type",
                        "group",
                        "era",
                        "cross-section",
                        "branching-ratio",
                    ]:
                        outf.write(f"    {k}: {v}\n")
                if cfg["type"] == "data":
                    outf.write("    cross-section: 1.0\n")

            for line in inf:
                if "groups:" in line:
                    line_found = True
                if line_found:
                    outf.write(line)

    if not os.path.isdir(os.path.join(workdir, "plots_full")):
        os.makedirs(os.path.join(workdir, "plots_full"))

    plotitCmd = [
        plotIt,
        "-o",
        f"{workdir}/plots_full",
        "--",
        f"{workdir}/plots_full.yml",
    ]
    try:
        if verbose:
            logger.info("running {}".format(" ".join(plotitCmd)))
            subprocess.check_call(plotitCmd)
        else:
            subprocess.check_call(plotitCmd, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        logger.error("Failed to run {0}".format(" ".join(plotitCmd)))


def runPDF(workdir, channel="DL", era=None, plotsDir="plots"):
    import os

    plots_dir = [os.path.join(workdir, plotsDir)]
    if era:
        plots_dir.append(os.path.join(workdir, f"plots_{era}"))
    else:
        era = "all"
    for plots_dir in plots_dir:
        return f"""
        cp scripts/empty.pdf {plots_dir}
        cp scripts/controlPlotter_{channel}.tex {plots_dir}
        cd {plots_dir}
        pdflatex -interaction=nonstopmode controlPlotter_{channel}.tex > /dev/null 2>&1
        mv controlPlotter_{channel}.pdf ../controlPlotter_{channel}_{era}.pdf
        cd - > /dev/null
        # pdflatex yields.tex
        # cd ../..
        """


def ml_input_var_binning(var_name):
    "Function to return binning, min and max values for the ML input feature plots."
    from bamboo.plots import EquidistantBinning as EqBin

    if "_Px" in var_name or "_Py" in var_name:
        N, mn, mx = 100, -1000, 1000
    elif "_Pz" in var_name:
        N, mn, mx = 200, -4000, 4000
    elif "_E" in var_name or "di_lepton_dijet_met_mass" in var_name:
        N, mn, mx = 100, 0, 2500
    elif "_charge" in var_name:
        N, mn, mx = 5, -2.5, 2.5
    elif "_btag" in var_name:
        N, mn, mx = 50, 0, 1
    elif "_pdgId" in var_name:
        N, mn, mx = 30, -15, 15
    elif "_pT" in var_name or "di_bjet_mass" in var_name or "di_lepton_met_mass" in var_name:
        N, mn, mx = 100, 0, 1000
    elif "_eta" in var_name:
        N, mn, mx = 30, -3, 3
    elif "_N" in var_name:
        N, mn, mx = 10, 0, 10
    elif "dR_" in var_name:
        N, mn, mx = 100, 0, 10
    elif "InvM_" in var_name:
        N, mn, mx = 100, 0, 1000
    elif "_LD"  in var_name or "di_lepton_mass" in var_name:
        N, mn, mx = 100, 0, 600
    elif "abs_dphi_" in var_name:
        N, mn, mx = 64, 0, 3.2
    elif "_tau" in var_name or "_msoftdrop" in var_name:
        N, mn, mx = 100, 0, 1
    elif "_tag" in var_name:
        N, mn, mx = 2, 0, 1
    elif "HT" in var_name:
        N, mn, mx = 100, 0, 2000
    elif "run_year" in var_name:
        N, mn, mx = 2, 2021.5, 2023.5
    return EqBin(N, mn, mx)


def electron_sc_eta(el_eta, el_phi, PV_x, PV_y, PV_z):
    """Derivation of the super cluster eta, taken from
    https://twiki.cern.ch/twiki/bin/view/CMS/EgammaNanoAOD#How_to_get_photon_supercluster_e.
    Starting from nanoAOD v15 the variable is made available in the tree
    and for v12 it's calculated as
    electron_sc_eta = electron.eta + electron.deltaEtaSC"""
    electron_isScEtaEB = op.switch(
        el_eta < 1.479, op.c_bool(1), op.c_bool(0)
    )  # double check this since isScEtaEB branch is only available for photons
    electron_isScEtaEE = op.switch(
        op.in_range(1.479, el_eta, 3.0), op.c_bool(1), op.c_bool(0)
    )
    # double check this since isScEtaEE branch is only available for photons
    tg_theta_over_2 = op.exp(-el_eta)
    tg_theta = 2 * tg_theta_over_2 / (1 - op.pow(tg_theta_over_2, 2))
    pi = 3.14159265359
    R = 130
    angle_x0_y0 = op.multiSwitch(
        (PV_x > 0, op.atan(PV_y / PV_x)),
        (PV_x < 0, pi + op.atan(PV_y / PV_x)),
        (PV_y > 0, pi / 2),
        -pi / 2,
    )
    alpha = angle_x0_y0 + (pi - el_phi)
    sin_beta = op.sqrt(op.pow(PV_x, 2) + op.pow(PV_y, 2)) / R * op.sin(alpha)
    beta = op.abs(op.asin(sin_beta))
    gamma = pi / 2 - alpha - beta
    l = op.sqrt(
        op.pow(R, 2)
        + op.pow(PV_x, 2)
        + op.pow(PV_y, 2)
        - 2 * R * op.sqrt(op.pow(PV_x, 2) + op.pow(PV_y, 2)) * op.cos(gamma)
    )
    z0_zSC = l / tg_theta
    intersection_z = op.switch(el_eta > 0, 310, -310)
    base = intersection_z - PV_z
    r = base * tg_theta
    crystalX = PV_x + r * op.cos(el_phi)
    crystalY = PV_y + r * op.sin(el_phi)
    tg_sctheta = op.multiSwitch(
        (electron_isScEtaEB, R / (PV_z + z0_zSC)),
        (
            electron_isScEtaEE,
            op.sqrt(op.pow(crystalX, 2) + op.pow(crystalY, 2)) / intersection_z,
        ),
        op.c_float(1.0),
    )
    sctheta = op.atan(tg_sctheta)
    sctheta = op.switch(sctheta < 0, sctheta + pi, sctheta)
    tg_sctheta_over_2 = op.tan(sctheta / 2)
    SCEta = -op.log(tg_sctheta_over_2)
    return op.switch(op.OR(electron_isScEtaEB, electron_isScEtaEE), SCEta, el_eta)
