from copy import deepcopy


def fillSampleTemplate(template, selEras=None):
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
    return {'labels': [{'text': label, 'position': [0.235, 0.9], 'size': 24}]}


def custom_Plotit(cfgName, workdir, inDir, outDir, counterReader, config, plotIt, verbose):
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
    if 'signal' in _gp:
        _gp.remove('signal')

    lumiCFG = {}
    for smp, smpCfg in config["samples"].items():

        era = smpCfg["era"]
        lumi = config["eras"][smpCfg["era"]]["luminosity"]
        lumiCFG[era] = lumi

        smpNm = smp.split(era)[0]
        mergedHists = {}

        if smpCfg.get("group") in _gp:
            shutil.copyfile(os.path.join(
                inDir, f"{smp}.root"), os.path.join(outDir, f"{smp}.root"))
            keep_cfg[smp] = smpCfg
        else:
            resultsFile = openFileAndGet(
                os.path.join(inDir, f"{smp}.root"), mode="READ")
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
                os.path.join(outDir, f"{smp}.root"), "recreate")
            for hist in mergedHists.values():
                hist.Write()
            normalizedFile.Close()

            to_hadd[smpNm].append(os.path.join(outDir, f"{smp}.root"))
            hadd_cfg[smpNm].update(smpCfg)

    for smp, val in to_hadd.items():
        sum_f = f"{smp}full.root"
        haddCmd = ["hadd", "-f", os.path.join(outDir, sum_f)]+val
        try:
            if verbose:
                logger.info("running {}".format(" ".join(haddCmd)))
                subprocess.check_call(haddCmd)
            else:
                subprocess.check_call(haddCmd , stdout=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            logger.error("Failed to run {0}".format(" ".join(haddCmd)))

    with open(cfgName, 'r') as inf:
        with open(os.path.join(workdir, 'plots_full.yml'), 'w+') as outf:
            outf.write('configuration:\n')
            outf.write("  blinded-range-fill-color: '#29556270'\n")
            outf.write('  blinded-range-fill-style: 1001\n')
            outf.write('  eras:\n')
            for era in lumiCFG.keys():
                outf.write(f'  - {era}\n')
            outf.write("  error-fill-color: '#ee556270'\n")
            outf.write('  error-fill-style: 3154\n')
            outf.write('  experiment: CMS\n')
            outf.write('  extra-label: Run 3 (2022) - Work in progress\n')
            outf.write('  height: 600\n')
            outf.write('  luminosity:\n')
            for era, lumi in lumiCFG.items():
                outf.write(f'    {era}: {lumi}\n')
            outf.write("  luminosity-label: '%1$.2f fb^{-1} (13.6 TeV)'\n")
            outf.write('  margin-bottom: 0.15\n')
            outf.write('  margin-left: 0.2\n')
            outf.write('  margin-right: 0.03\n')
            outf.write('  margin-top: 0.05\n')
            outf.write("  ratio-fit-error-fill-color: '#aa556270'\n")
            outf.write('  ratio-fit-error-fill-style: 1001\n')
            outf.write("  ratio-fit-line-color: '#0B486B'\n")
            outf.write(f'  root: {workdir}/results/normalizedSummedSignal\n')
            outf.write('  width: 800\n')
            outf.write('  yields-table-align: v\n')
            outf.write("  yields-table-text-align: l\n")
            outf.write('files:\n')
            line_found = False

            for smp, cfg in hadd_cfg.items():
                outf.write(f'  {smp}full.root:\n')
                # smpFile = openFileAndGet(
                # os.path.join(outDir, f"{smp}full.root"), mode="READ")
                # gevt = counterReader(smpFile)[smpCfg["generated-events"]]
                # outf.write(f"    generated-events: {gevt}\n")
                for k, v in cfg.items():
                    if k in ['type', 'group', 'line-width', 'line-type', 'legend', 'line-color', 'branching-ratio', 'cross-section', 'era']:
                        outf.write(f"    {k}: {v}\n")

            for smp, cfg in keep_cfg.items():
                outf.write(f'  {smp}.root:\n')
                if cfg['type'] != 'data':
                    smpFile = openFileAndGet(
                    os.path.join(inDir, f"{smp}.root"), mode="READ")
                    gevt = counterReader(smpFile)[smpCfg["generated-events"]]
                    outf.write(f"    generated-events: {gevt}\n")
                else: None
                for k, v in cfg.items():
                    if k in ['type', 'group', 'era', 'cross-section', 'branching-ratio']:
                        outf.write(f"    {k}: {v}\n")
                if cfg['type'] == 'data':
                    outf.write("    cross-section: 1.0\n")

            for line in inf:
                if 'groups:' in line:
                    line_found = True
                if line_found:
                    outf.write(line)

    if not os.path.isdir(os.path.join(workdir, "plots_full")):
        os.makedirs(os.path.join(workdir, "plots_full"))

    plotitCmd = [plotIt,
                 "-o", f'{workdir}/plots_full', "--", f"{workdir}/plots_full.yml"]
    try:
        if verbose:
            logger.info("running {}".format(" ".join(plotitCmd)))
            subprocess.check_call(plotitCmd)
        else:
            subprocess.check_call(plotitCmd, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        logger.error("Failed to run {0}".format(" ".join(plotitCmd)))

