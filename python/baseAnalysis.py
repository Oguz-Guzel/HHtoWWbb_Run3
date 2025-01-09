import os
import re
import logging
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.analysismodules import NanoAODModule, HistogramsModule

import utils

logger = logging.getLogger(__name__)

jsonPathBase = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"

PU_JSONFiles = {
    "2022": (jsonPathBase + "LUM/2022_Summer22/puWeights.json.gz", "Collisions2022_355100_357900_eraBCD_GoldenJson"),
    "2022EE": (jsonPathBase + "LUM/2022_Summer22EE/puWeights.json.gz", "Collisions2022_359022_362760_eraEFG_GoldenJson"),
    "2023": (jsonPathBase + "LUM/2023_Summer23/puWeights.json.gz", "Collisions2023_366403_369802_eraBC_GoldenJson"),
    "2023BPix": (jsonPathBase + "LUM/2023_Summer23BPix/puWeights.json.gz", "Collisions2023_369803_370790_eraD_GoldenJson"),
}

JECTags = {
    "2022": {
        "MC": "Summer22_22Sep2023_V2_MC",
        "C": "Summer22_22Sep2023_RunCD_V2_DATA",
        "D": "Summer22_22Sep2023_RunCD_V2_DATA"
    },
    "2022EE": {
        "MC": "Summer22EE_22Sep2023_V2_MC",
        "E": "Summer22EE_22Sep2023_RunE_V2_DATA",
        "F": "Summer22EE_22Sep2023_RunF_V2_DATA",
        "G": "Summer22EE_22Sep2023_RunG_V2_DATA"
    },
    "2023": {
        "MC": "Summer23Prompt23_V1_MC",
        "Cv123": "Summer23Prompt23_RunCv123_V1_DATA",
        "Cv4": "Summer23Prompt23_RunCv4_V1_DATA"
    },
    "2023BPix": {
        "MC": "Summer23BPixPrompt23_V1_MC",
        "D": "Summer23BPixPrompt23_RunD_V1_DATA"
    }
}

JERTags = {
    "2022": "Summer22_22Sep2023_JRV1_MC",
    "2022EE": "Summer22EE_22Sep2023_JRV1_MC",
    "2023": "Summer23Prompt23_RunCv1234_JRV1_MC",
    "2023BPix": "Summer23BPixPrompt23_RunD_JRV1_MC"
}

JEC_JSONFiles = {
    "2022": {
        "AK4": jsonPathBase + "JME/2022_Summer22/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22/fatJet_jerc.json.gz"},
    "2022EE": {
        "AK4": jsonPathBase + "JME/2022_Summer22EE/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22EE/fatJet_jerc.json.gz"},
    "2023": {
        "AK4": jsonPathBase + "JME/2023_Summer23/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2023_Summer23/fatJet_jerc.json.gz"},
    "2023BPix": {
        "AK4": jsonPathBase + "JME/2023_Summer23BPix/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2023_Summer23BPix/fatJet_jerc.json.gz"},
}


def getDataRunEra(sample):
    """Return run era (A/B/...) and the following digits for data sample"""
    result = re.search(r'Run20\d{2}([A-Z]\w*)', sample)
    return result.group(1) if result else None


class NanoBaseHHWWbb(NanoAODModule, HistogramsModule):
    """ Base module for HH->WWbb analysis """

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("-c", "--channel",
                            dest="channel",
                            type=str,
                            default="DL",
                            help='Channel to be selected between SL and DL')
        parser.add_argument("--mvaModels",
                            dest="mvaModels",
                            type=str,
                            default="./DNN/",
                            help="Path to MVA models and Evaluate DNN")
        parser.add_argument("--backend", type=str, default="dataframe",
                            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")
        parser.add_argument("--sync", action="store_true", default=False,
                            help="Run synchronisation")

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):

        # Define the git project's directory
        self.git_project_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..'))
        self.era = sampleCfg["era"] if sampleCfg else None
        self.is_MC = self.isMC(sample)

        from bamboo.plots import CutFlowReport
        self.yields = CutFlowReport(
            "yields", recursive=False, printInLog=True)

        # Decorate the tree
        from bamboo.treedecorators import NanoAODDescription, nanoFatJetCalc, CalcCollectionsGroups
        metName = "PuppiMET"
        nanoJetMETCalc_both = CalcCollectionsGroups(
            Jet=("pt", "mass"), changes={metName: (f"{metName}T1", f"{metName}T1Smear")},
            **{metName: ("pt", "phi")})
        nanoJetMETCalc_data = CalcCollectionsGroups(
            Jet=("pt", "mass"), changes={metName: (f"{metName}T1",)},
            **{metName: ("pt", "phi")})
        systVariations = (([nanoFatJetCalc])
                          + [nanoJetMETCalc_both if self.is_MC else nanoJetMETCalc_data])
        tree, noSel, be, lumiArgs = super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v12", year=self.era[:4], isMC=self.is_MC, systVariations=systVariations),
            backend=backend)

        # JEC/JER
        jecTag = JECTags[self.era]["MC" if self.is_MC else getDataRunEra(
            sample)]
        logger.info(f"JEC tag for sample {sample} is {jecTag}")
        # smearing is for MC only
        smearTag = JERTags[self.era] if self.is_MC else None
        if smearTag:
            logger.info(f"JER tag for sample {sample} is {smearTag}")

        jecArgs = {
            "jsonFile": JEC_JSONFiles[self.era]["AK4"],
            "jec": jecTag,
            # "smear": smearTag,
            "splitJER": True,
            "jsonFileSmearingTool": jsonPathBase+'JME/jer_smear.json.gz',
            "jesUncertaintySources": (["Total"] if self.is_MC else None),
            "isMC": self.is_MC,
            "backend": be,
        }

        from bamboo.analysisutils import configureJets, configureType1MET
        configureJets(tree._Jet, jetType="AK4PFPuppi", **jecArgs)
        metName = "PuppiMET"
        configureType1MET(
            getattr(tree, f"_{metName}T1"),
            enableSystematics=(
                (lambda v: not v.startswith("jer")) if self.is_MC else None),
            **jecArgs)
        jecArgs.update({"jsonFile": JEC_JSONFiles[self.era]["AK8"], })
        jecArgs.update({"jetAlgoSubjet": "AK4PFPuppi", })
        jecArgs.update({"jecSubjet": jecTag, })
        jecArgs.update({"jsonFileSubjet": JEC_JSONFiles[self.era]["AK4"], })
        configureJets(tree._FatJet, jetType="AK8PFPuppi", **jecArgs)
        logger.info("Applying Jet energy and resolution corrections")

        # Number of events before any processing
        self.yields.add(noSel, "noSel")

        # MC weight
        if self.is_MC:
            logger.info("Applying genWeight")
            noSel = noSel.refine('genWeight', weight=tree.genWeight)
        else:
            noSel = noSel.refine('genWeight', weight=op.c_float(1.))
        self.yields.add(noSel, "genWeight")

        # PU weight
        if self.is_MC:
            from bamboo.analysisutils import makePileupWeight
            pileupWeight = makePileupWeight(
                PU_JSONFiles[self.era], tree.Pileup_nTrueInt, systName="pileup", sel=noSel)
            logger.info("Applying PU weight")
            noSel = noSel.refine('puWeight', weight=pileupWeight)
        else:
            noSel = noSel.refine('puWeight', weight=op.c_float(1.))
        self.yields.add(noSel, "puWeight")

        # Triggers
        self.triggers_per_PD = {}

        def addHLTPath(PD, HLT):
            if PD not in self.triggers_per_PD.keys():
                self.triggers_per_PD[PD] = []
            try:
                self.triggers_per_PD[PD].append(
                    getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, cross check!" % HLT)

        if self.era == '2022':
            if sample.startswith("SingleMuon_") or sample.startswith("DoubleMuon_"):
                addHLTPath("DoubleMuon_",
                           "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
                addHLTPath("SingleMuon_", "IsoMu24")
            else:
                addHLTPath("Muon_",
                           "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
                addHLTPath("Muon_", "IsoMu24")

        else:
            addHLTPath("Muon_", "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            addHLTPath("Muon_", "IsoMu24")

        addHLTPath("EGamma_", "Ele30_WPTight_Gsf")
        addHLTPath("EGamma_", "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
        addHLTPath("MuonEG_", "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
        addHLTPath("MuonEG_", "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
        addHLTPath("MuonEG_", "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")

        if self.is_MC:
            noSel = noSel.refine('triggers',  cut=(
                op.OR(*chain.from_iterable(self.triggers_per_PD.values()))))
        else:
            noSel = noSel.refine('triggers', cut=makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggers_per_PD))

        self.yields.add(noSel, "triggers")

        return tree, noSel, be, lumiArgs

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Postprocess: run plotIt

        The list of plots is created if needed (from a representative file,
        this enables rerunning the postprocessing step on the results files),
        and then plotIt is executed
        """
        import os
        if not self.plotList:
            self.plotList = self.getPlotList(
                resultsdir=resultsdir, config=config)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [
            ap for ap in self.plotList if isinstance(ap, CutFlowReport)]
        plotList_plotIt = [ap for ap in self.plotList
                           if (isinstance(ap, Plot) or isinstance(ap, DerivedPlot))
                           and len(ap.binnings) == 1]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            from bamboo.analysisutils import printCutFlowReports
            printCutFlowReports(
                config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt and not self.args.sync:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(
                config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, plotDefaults=self.plotDefaults,
                vetoFileAttributes=self.__class__.CustomSampleAttributes)
            runPlotIt(
                cfgName, workdir=workdir, plotIt=self.args.plotIt, eras=(
                    eraMode, eras),
                verbose=self.args.verbose)
        plotsDir = 'plots'
        # # to be automatised soon
        # # activate the following lines to combine signal samples
        # # hadd signal files and create another plots.yml called plots_full.yml
        # plotsDir = 'plots_full'
        # import os
        # import shutil
        # outDir = os.path.join(resultsdir, "normalizedSummedSignal")
        # if os.path.isdir(outDir):
        #     shutil.rmtree(outDir)
        # os.makedirs(outDir)
        # utils.custom_Plotit(cfgName, workdir, resultsdir, outDir, self.readCounters,
        #                     config, plotIt=self.args.plotIt, verbose=self.args.verbose)
        # # end of merging signal samples

        # create pdf presentation
        if not self.mvaModels and not self.args.sync:
            try:
                for era in eras:
                    os.system(utils.runPDF(workdir=workdir,
                              channel=self.args.channel, era=era))
                    logger.info(f"PDF presentation created for era {era}.\n")
                os.system(utils.runPDF(workdir=workdir,
                          channel=self.args.channel, plotsDir=plotsDir))
                logger.info(
                    f"PDF presentation created for all eras combined.\n")
            except Exception as e:
                logger.info(e)

        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]

        from bamboo.analysisutils import loadPlotIt
        _, samples, _, _, _ = loadPlotIt(
            config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)

        # create sync skims if asked
        if self.args.sync:
            if skims:
                from bamboo.root import gbl
                import pandas as pd
                sync_dfs = []
                for skim in skims:
                    frames = []
                    for smp in samples:
                        for cb in (smp.files if hasattr(smp, "files") else [smp]):
                            tree = cb.tFile.Get(skim.treeName)
                            if not tree:
                                logger.info("WARNING: skim tree %s not found in file %s" % (
                                    skim.treeName, cb.tFile.GetName()))
                                logger.info("         skipping...")
                            else:
                                N = tree.GetEntries()
                                cols = gbl.ROOT.RDataFrame(tree).AsNumpy()
                                if "sync" not in skim.name:
                                    cols["weight"] *= cb.scale
                                    cols["process"] = [smp.name] * \
                                        len(cols["weight"])
                                frames.append(pd.DataFrame(cols))
                    df = pd.concat(frames)
                    df = df[self.order]
                    sync_dfs.append(df)
                df = pd.concat(sync_dfs)
                df = df.sort_values(by='event_no')
                syncFileName = f"{self.channel}_sync.csv"
                df.to_csv(os.path.join(resultsdir, syncFileName))
                logger.info(f"Saved dataframe for sync to {syncFileName}")
            else:
                logger.warning(
                    "No skims are found, hence sync file is not produced.")
