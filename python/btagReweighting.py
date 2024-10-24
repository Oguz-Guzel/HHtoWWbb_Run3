
from bamboo.plots import Plot, Skim
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection

import definitions as defs

import os
from ROOT import TFile
from itertools import chain
import logging
logger = logging.getLogger(__name__)

jsonPathBase = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"

PU_JSONFiles = {
    "2022": (jsonPathBase + "LUM/2022_Summer22/puWeights.json.gz", "Collisions2022_355100_357900_eraBCD_GoldenJson"),
    "2022EE": (jsonPathBase + "LUM/2022_Summer22EE/puWeights.json.gz", "Collisions2022_359022_362760_eraEFG_GoldenJson"),
    "2023": (jsonPathBase + "LUM/2023_Summer23/puWeights.json.gz", "Collisions2023_366403_369802_eraBC_GoldenJson"),
    "2023BPix": (jsonPathBase + "LUM/2023_Summer23BPix/puWeights.json.gz", "Collisions2023_369803_370790_eraD_GoldenJson"),
}

BTV_SF_JSONFiles = {
    "2022": jsonPathBase + "BTV/2022_Summer22/btagging.json.gz",
    "2022EE": jsonPathBase + "BTV/2022_Summer22EE/btagging.json.gz",
    "2023": jsonPathBase + "BTV/2023_Summer23/btagging.json.gz",
    "2023BPix": jsonPathBase + "BTV/2023_Summer23BPix/btagging.json.gz",
}


class _base(NanoAODModule, HistogramsModule):

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("-c", "--channel",
                            dest="channel",
                            type=str,
                            default="DL",
                            help='Channel to be selected between SL and DL')
        parser.add_argument("--backend", type=str, default="dataframe",
                            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'")

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):

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
        systVars = (([nanoFatJetCalc])
                    + [nanoJetMETCalc_both if self.is_MC else nanoJetMETCalc_data])
        tree, noSel, be, lumiArgs = super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v12", year=self.era[:4], isMC=self.is_MC, systVariations=systVars),
            backend=self.args.backend or backend)

        # Number of events before any processing
        self.yields.add(noSel, "noSel")

        # MC weight
        if self.is_MC:
            logger.info("Applying genWeight")
            noSel = noSel.refine('genWeight', weight=tree.genWeight)
        else:
            noSel = noSel.refine('genWeight', weight=op.c_float(1.))
        self.yields.add(noSel, "genWeight")

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
                addHLTPath("SingleMuon_", "IsoMu27")
            else:
                addHLTPath("Muon_",
                           "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
                addHLTPath("Muon_", "IsoMu24")
                addHLTPath("Muon_", "IsoMu27")

        else:
            addHLTPath("Muon_", "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            addHLTPath("Muon_", "IsoMu24")
            addHLTPath("Muon_", "IsoMu27")

        addHLTPath("EGamma_", "Ele32_WPTight_Gsf")
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


class btagReweighting(_base):
    """ Class to apply reweighting that should be done before applying btag SF."""

    def __init__(self, args):
        super(btagReweighting, self).__init__(args)
        self.channel = self.args.channel

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # cutflow report
        plots.append(self.yields)

        # define objects
        defs.defineObjects(self, tree)

        preprend = '_before'

        plots.append(
            Skim("WeightsBeforeBtagSF",
                 {
                     "jet_pt" + preprend: self.ak4Jets[0].p4.Pt(),
                     "jet_eta" + preprend: self.ak4Jets[0].p4.Eta(),
                     "jetMultiplicity" + preprend: op.rng_len(self.ak4Jets),
                     "Sel_weight" + preprend: noSel.weight,
                 },
                 noSel
                 ))

        preprend = '_after'

        # btagging SF
        if self.is_MC:
            from bamboo.scalefactors import get_bTagSF_itFit, makeBtagWeightItFit
            logger.info("Applying btagging SF")
            def btvSF(flav): return get_bTagSF_itFit(
                BTV_SF_JSONFiles[self.era], "particleNet", "btagPNetB", flav, sel=noSel, decorr_eras=True, era=self.era)
            btvWeight = makeBtagWeightItFit(self.ak4Jets, btvSF)
            btagSF = noSel.refine("btagSF", weight=btvWeight)
        else:
            btagSF = noSel.refine("btagSF", weight=op.c_float(1.))

        plots.append(
            Skim("WeightsAfterBtagSF",
                 {
                     "jet_pt" + preprend: self.ak4Jets[0].p4.Pt(),
                     "jet_eta" + preprend: self.ak4Jets[0].p4.Eta(),
                     "jetMultiplicity" + preprend: op.rng_len(self.ak4Jets),
                     "Sel_weight" + preprend: btagSF.weight,
                 },
                 noSel
                 ))

        plots.extend([
            Plot.make1D("noSel_n_ak8", op.rng_len(self.ak8Jets), noSel, EqBin(
                10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet"),
        ])

        return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """PostProcessing step."""
        super().postProcess(taskList, config=config, workdir=workdir,
                            resultsdir=resultsdir)

        if not self.plotList:
            self.plotList = self.getPlotList(
                resultsdir=resultsdir, config=config)

        from bamboo.plots import CutFlowReport, DerivedPlot, Plot, Skim

        plotList_plotIt = [ap for ap in self.plotList if (
            isinstance(ap, Plot) or isinstance(ap, DerivedPlot))]
        for plots in plotList_plotIt:
            logger.info("plots in plotList_plotIt: {0}".format(
                plots.name))
        plotList_cutFlow = [
            ap for ap in self.plotList if isinstance(ap, CutFlowReport)]

        for cutflow in plotList_cutFlow:
            logger.info("cutflow in plotList_cutFlow: {0}".format(cutflow))

        skim_list = [ap for ap in self.plotList if isinstance(ap, Skim)]

        eraMode, eras = self.args.eras
        logger.info("EraMode: {0} and eras: {1}".format(eraMode, eras))
        if eras is None:
            eras = list(config["eras"].keys())

        if plotList_cutFlow:
            from bamboo.analysisutils import printCutFlowReports
            printCutFlowReports(config, plotList_cutFlow, workdir=workdir, resultsdir=resultsdir,
                                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)

        dic_weights = {}
        if skim_list:

            mc_sumWeigth_before = 0
            mc_sumWeigth_after = 0

            mc_sumWeigth_2022_before = 0
            mc_sumWeigth_2022_after = 0

            mc_sumWeigth_2022EE_before = 0
            mc_sumWeigth_2022EE_after = 0

            sumWeight_jet_multiplicity_0_2022_before = 0
            sumWeight_jet_multiplicity_1_2022_before = 0
            sumWeight_jet_multiplicity_2_2022_before = 0
            sumWeight_jet_multiplicity_3_2022_before = 0
            sumWeight_jet_multiplicity_4_2022_before = 0
            sumWeight_jet_multiplicity_5_2022_before = 0
            sumWeight_jet_multiplicity_6_2022_before = 0
            sumWeight_jet_multiplicity_7_2022_before = 0
            sumWeight_jet_multiplicity_8_2022_before = 0
            sumWeight_jet_multiplicity_9_2022_before = 0
            sumWeight_jet_multiplicity_10_2022_before = 0

            sumWeight_jet_multiplicity_0_2022_after = 0
            sumWeight_jet_multiplicity_1_2022_after = 0
            sumWeight_jet_multiplicity_2_2022_after = 0
            sumWeight_jet_multiplicity_3_2022_after = 0
            sumWeight_jet_multiplicity_4_2022_after = 0
            sumWeight_jet_multiplicity_5_2022_after = 0
            sumWeight_jet_multiplicity_6_2022_after = 0
            sumWeight_jet_multiplicity_7_2022_after = 0
            sumWeight_jet_multiplicity_8_2022_after = 0
            sumWeight_jet_multiplicity_9_2022_after = 0
            sumWeight_jet_multiplicity_10_2022_after = 0

            sumWeight_jet_multiplicity_0_2022EE_before = 0
            sumWeight_jet_multiplicity_1_2022EE_before = 0
            sumWeight_jet_multiplicity_2_2022EE_before = 0
            sumWeight_jet_multiplicity_3_2022EE_before = 0
            sumWeight_jet_multiplicity_4_2022EE_before = 0
            sumWeight_jet_multiplicity_5_2022EE_before = 0
            sumWeight_jet_multiplicity_6_2022EE_before = 0
            sumWeight_jet_multiplicity_7_2022EE_before = 0
            sumWeight_jet_multiplicity_8_2022EE_before = 0
            sumWeight_jet_multiplicity_9_2022EE_before = 0
            sumWeight_jet_multiplicity_10_2022EE_before = 0

            sumWeight_jet_multiplicity_0_2022EE_after = 0
            sumWeight_jet_multiplicity_1_2022EE_after = 0
            sumWeight_jet_multiplicity_2_2022EE_after = 0
            sumWeight_jet_multiplicity_3_2022EE_after = 0
            sumWeight_jet_multiplicity_4_2022EE_after = 0
            sumWeight_jet_multiplicity_5_2022EE_after = 0
            sumWeight_jet_multiplicity_6_2022EE_after = 0
            sumWeight_jet_multiplicity_7_2022EE_after = 0
            sumWeight_jet_multiplicity_8_2022EE_after = 0
            sumWeight_jet_multiplicity_9_2022EE_after = 0
            sumWeight_jet_multiplicity_10_2022EE_after = 0

            for proc, smpCfg in config["samples"].items():

                if smpCfg.get("group") == "data":

                    mc_sumWeigth_before += 0.0
                    mc_sumWeigth_after += 0.0
                else:
                    logger.info("Sample name to be process: {}".format(proc))

                    def _openFileAndGet(path, mode="read"):
                        """Open ROOT file in a mode, check if open properly, and return TFile handle"""
                        tf = TFile.Open(path, mode)
                        if not tf or not tf.IsOpen():
                            raise Exception(
                                "Could not open file {}".format(path))
                        return tf

                    sample_rootfile = _openFileAndGet(os.path.join(
                        resultsdir, proc + ".root"), "read")  # already the TFile

                    genEvents = self.readCounters(sample_rootfile)[
                        smpCfg["generated-events"]]
                    lumi = config["eras"][smpCfg["era"]]["luminosity"]
                    Xsection = smpCfg["cross-section"]
                    smpScale = lumi*Xsection/genEvents
                    era = smpCfg["era"]

                    jet_multiplicity_0_before = 0
                    jet_multiplicity_1_before = 0
                    jet_multiplicity_2_before = 0
                    jet_multiplicity_3_before = 0
                    jet_multiplicity_4_before = 0
                    jet_multiplicity_5_before = 0
                    jet_multiplicity_6_before = 0
                    jet_multiplicity_7_before = 0
                    jet_multiplicity_8_before = 0
                    jet_multiplicity_9_before = 0
                    jet_multiplicity_10_before = 0

                    jet_multiplicity_0_after = 0
                    jet_multiplicity_1_after = 0
                    jet_multiplicity_2_after = 0
                    jet_multiplicity_3_after = 0
                    jet_multiplicity_4_after = 0
                    jet_multiplicity_5_after = 0
                    jet_multiplicity_6_after = 0
                    jet_multiplicity_7_after = 0
                    jet_multiplicity_8_after = 0
                    jet_multiplicity_9_after = 0
                    jet_multiplicity_10_after = 0

                    for skim in skim_list:
                        tree = sample_rootfile.Get(skim.treeName)

                        if not tree:
                            logger.info("Warning: skim tree %s not found in file %s" % (
                                skim.treeName, sample_rootfile.GetName()))
                            continue

                        else:
                            branches = tree.GetListOfBranches()

                            if "Sel_weight_after" in [branch.GetName() for branch in branches]:
                                afterweight_sum = 0

                                for entry in tree:
                                    afterweight_sum += entry.Sel_weight_after

                                    if entry.jetMultiplicity_after == 0:
                                        jet_multiplicity_0_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 1:
                                        jet_multiplicity_1_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 2:
                                        jet_multiplicity_2_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 3:
                                        jet_multiplicity_3_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 4:
                                        jet_multiplicity_4_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 5:
                                        jet_multiplicity_5_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 6:
                                        jet_multiplicity_6_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 7:
                                        jet_multiplicity_7_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 8:
                                        jet_multiplicity_8_after += entry.Sel_weight_after

                                    elif entry.jetMultiplicity_after == 9:
                                        jet_multiplicity_9_after += entry.Sel_weight_after

                                    else:
                                        jet_multiplicity_10_after += entry.Sel_weight_after

                            elif "Sel_weight_before" in [branch.GetName() for branch in branches]:
                                beforeweight_sum = 0

                                for entry in tree:

                                    beforeweight_sum += entry.Sel_weight_before

                                    if entry.jetMultiplicity_before == 0:
                                        jet_multiplicity_0_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 1:
                                        jet_multiplicity_1_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 2:
                                        jet_multiplicity_2_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 3:
                                        jet_multiplicity_3_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 4:
                                        jet_multiplicity_4_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 5:
                                        jet_multiplicity_5_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 6:
                                        jet_multiplicity_6_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 7:
                                        jet_multiplicity_7_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 8:
                                        jet_multiplicity_8_before += entry.Sel_weight_before

                                    elif entry.jetMultiplicity_before == 9:
                                        jet_multiplicity_9_before += entry.Sel_weight_before

                                    else:
                                        jet_multiplicity_10_before += entry.Sel_weight_before

                    logger.info("Sum of the Weights, Before {0} and After {1}".format(
                        beforeweight_sum, afterweight_sum))

                    if era == "2022":
                        # sum of the weights before and after applying the b-tag SF: Total weight
                        mc_sumWeigth_2022_before += beforeweight_sum * smpScale
                        mc_sumWeigth_2022_after += afterweight_sum * smpScale
                        # sum of the weights before and after applying the b-tag SF: split by multiplicity
                        sumWeight_jet_multiplicity_0_2022_before += smpScale * jet_multiplicity_0_before
                        sumWeight_jet_multiplicity_1_2022_before += smpScale * jet_multiplicity_1_before
                        sumWeight_jet_multiplicity_2_2022_before += smpScale * jet_multiplicity_2_before
                        sumWeight_jet_multiplicity_3_2022_before += smpScale * jet_multiplicity_3_before
                        sumWeight_jet_multiplicity_4_2022_before += smpScale * jet_multiplicity_4_before
                        sumWeight_jet_multiplicity_5_2022_before += smpScale * jet_multiplicity_5_before
                        sumWeight_jet_multiplicity_6_2022_before += smpScale * jet_multiplicity_6_before
                        sumWeight_jet_multiplicity_7_2022_before += smpScale * jet_multiplicity_7_before
                        sumWeight_jet_multiplicity_8_2022_before += smpScale * jet_multiplicity_8_before
                        sumWeight_jet_multiplicity_9_2022_before += smpScale * jet_multiplicity_9_before
                        sumWeight_jet_multiplicity_10_2022_before += smpScale * jet_multiplicity_10_before

                        sumWeight_jet_multiplicity_0_2022_after += smpScale * jet_multiplicity_0_after
                        sumWeight_jet_multiplicity_1_2022_after += smpScale * jet_multiplicity_1_after
                        sumWeight_jet_multiplicity_2_2022_after += smpScale * jet_multiplicity_2_after
                        sumWeight_jet_multiplicity_3_2022_after += smpScale * jet_multiplicity_3_after
                        sumWeight_jet_multiplicity_4_2022_after += smpScale * jet_multiplicity_4_after
                        sumWeight_jet_multiplicity_5_2022_after += smpScale * jet_multiplicity_5_after
                        sumWeight_jet_multiplicity_6_2022_after += smpScale * jet_multiplicity_6_after
                        sumWeight_jet_multiplicity_7_2022_after += smpScale * jet_multiplicity_7_after
                        sumWeight_jet_multiplicity_8_2022_after += smpScale * jet_multiplicity_8_after
                        sumWeight_jet_multiplicity_9_2022_after += smpScale * jet_multiplicity_9_after
                        sumWeight_jet_multiplicity_10_2022_after += smpScale * jet_multiplicity_10_after

                    elif era == "2022EE":
                        mc_sumWeigth_2022EE_before += beforeweight_sum * smpScale
                        mc_sumWeigth_2022EE_after += afterweight_sum * smpScale

                        # sum of the weights before and after applying the b-tag SF: split by multiplicity 2022 EE
                        sumWeight_jet_multiplicity_0_2022EE_before += smpScale * jet_multiplicity_0_before
                        sumWeight_jet_multiplicity_1_2022EE_before += smpScale * jet_multiplicity_1_before
                        sumWeight_jet_multiplicity_2_2022EE_before += smpScale * jet_multiplicity_2_before
                        sumWeight_jet_multiplicity_3_2022EE_before += smpScale * jet_multiplicity_3_before
                        sumWeight_jet_multiplicity_4_2022EE_before += smpScale * jet_multiplicity_4_before
                        sumWeight_jet_multiplicity_5_2022EE_before += smpScale * jet_multiplicity_5_before
                        sumWeight_jet_multiplicity_6_2022EE_before += smpScale * jet_multiplicity_6_before
                        sumWeight_jet_multiplicity_7_2022EE_before += smpScale * jet_multiplicity_7_before
                        sumWeight_jet_multiplicity_8_2022EE_before += smpScale * jet_multiplicity_8_before
                        sumWeight_jet_multiplicity_9_2022EE_before += smpScale * jet_multiplicity_9_before
                        sumWeight_jet_multiplicity_10_2022EE_before += smpScale * jet_multiplicity_10_before

                        sumWeight_jet_multiplicity_0_2022EE_after += smpScale * jet_multiplicity_0_after
                        sumWeight_jet_multiplicity_1_2022EE_after += smpScale * jet_multiplicity_1_after
                        sumWeight_jet_multiplicity_2_2022EE_after += smpScale * jet_multiplicity_2_after
                        sumWeight_jet_multiplicity_3_2022EE_after += smpScale * jet_multiplicity_3_after
                        sumWeight_jet_multiplicity_4_2022EE_after += smpScale * jet_multiplicity_4_after
                        sumWeight_jet_multiplicity_5_2022EE_after += smpScale * jet_multiplicity_5_after
                        sumWeight_jet_multiplicity_6_2022EE_after += smpScale * jet_multiplicity_6_after
                        sumWeight_jet_multiplicity_7_2022EE_after += smpScale * jet_multiplicity_7_after
                        sumWeight_jet_multiplicity_8_2022EE_after += smpScale * jet_multiplicity_8_after
                        sumWeight_jet_multiplicity_9_2022EE_after += smpScale * jet_multiplicity_9_after
                        sumWeight_jet_multiplicity_10_2022EE_after += smpScale * jet_multiplicity_10_after

                    mc_sumWeigth_before += smpScale * beforeweight_sum
                    mc_sumWeigth_after += smpScale * afterweight_sum

                    logger.info("---------------------")

                    sample_rootfile.Close()
                    # output_rootfile.Close()

            mc_weights_dic = {}
            weights_jet_multiplicity = {}
            weights_jet_multiplicity_2022 = {}
            weights_jet_multiplicity_2022EE = {}

            mc_weights_dic["2022"] = mc_sumWeigth_2022_before / \
                mc_sumWeigth_2022_after
            mc_weights_dic["2022EE"] = mc_sumWeigth_2022EE_before / \
                mc_sumWeigth_2022EE_after
            mc_weights_dic["CombEras"] = mc_sumWeigth_before / \
                mc_sumWeigth_after

            weights_jet_multiplicity_2022["0"] = 1 if sumWeight_jet_multiplicity_0_2022_after == 0 else sumWeight_jet_multiplicity_0_2022_before / \
                sumWeight_jet_multiplicity_0_2022_after
            weights_jet_multiplicity_2022["1"] = 1 if sumWeight_jet_multiplicity_1_2022_after == 0 else sumWeight_jet_multiplicity_1_2022_before / \
                sumWeight_jet_multiplicity_1_2022_after
            weights_jet_multiplicity_2022["2"] = 1 if sumWeight_jet_multiplicity_2_2022_after == 0 else sumWeight_jet_multiplicity_2_2022_before / \
                sumWeight_jet_multiplicity_2_2022_after
            weights_jet_multiplicity_2022["3"] = 1 if sumWeight_jet_multiplicity_3_2022_after == 0 else sumWeight_jet_multiplicity_3_2022_before / \
                sumWeight_jet_multiplicity_3_2022_after
            weights_jet_multiplicity_2022["4"] = 1 if sumWeight_jet_multiplicity_4_2022_after == 0 else sumWeight_jet_multiplicity_4_2022_before / \
                sumWeight_jet_multiplicity_4_2022_after
            weights_jet_multiplicity_2022["5"] = 1 if sumWeight_jet_multiplicity_5_2022_after == 0 else sumWeight_jet_multiplicity_5_2022_before / \
                sumWeight_jet_multiplicity_5_2022_after
            weights_jet_multiplicity_2022["6"] = 1 if sumWeight_jet_multiplicity_6_2022_after == 0 else sumWeight_jet_multiplicity_6_2022_before / \
                sumWeight_jet_multiplicity_6_2022_after
            weights_jet_multiplicity_2022["7"] = 1 if sumWeight_jet_multiplicity_7_2022_after == 0 else sumWeight_jet_multiplicity_7_2022_before / \
                sumWeight_jet_multiplicity_7_2022_after
            weights_jet_multiplicity_2022["8"] = 1 if sumWeight_jet_multiplicity_8_2022_after == 0 else sumWeight_jet_multiplicity_8_2022_before / \
                sumWeight_jet_multiplicity_8_2022_after
            weights_jet_multiplicity_2022["9"] = 1 if sumWeight_jet_multiplicity_9_2022_after == 0 else sumWeight_jet_multiplicity_9_2022_before / \
                sumWeight_jet_multiplicity_9_2022_after
            weights_jet_multiplicity_2022["10"] = 1 if sumWeight_jet_multiplicity_10_2022_after == 0 else sumWeight_jet_multiplicity_10_2022_before / \
                sumWeight_jet_multiplicity_10_2022_after

            weights_jet_multiplicity_2022EE["0"] = 1 if sumWeight_jet_multiplicity_0_2022EE_after == 0 else sumWeight_jet_multiplicity_0_2022EE_before / \
                sumWeight_jet_multiplicity_0_2022EE_after
            weights_jet_multiplicity_2022EE["1"] = 1 if sumWeight_jet_multiplicity_1_2022EE_after == 0 else sumWeight_jet_multiplicity_1_2022EE_before / \
                sumWeight_jet_multiplicity_1_2022EE_after
            weights_jet_multiplicity_2022EE["2"] = 1 if sumWeight_jet_multiplicity_2_2022EE_after == 0 else sumWeight_jet_multiplicity_2_2022EE_before / \
                sumWeight_jet_multiplicity_2_2022EE_after
            weights_jet_multiplicity_2022EE["3"] = 1 if sumWeight_jet_multiplicity_3_2022EE_after == 0 else sumWeight_jet_multiplicity_3_2022EE_before / \
                sumWeight_jet_multiplicity_3_2022EE_after
            weights_jet_multiplicity_2022EE["4"] = 1 if sumWeight_jet_multiplicity_4_2022EE_after == 0 else sumWeight_jet_multiplicity_4_2022EE_before / \
                sumWeight_jet_multiplicity_4_2022EE_after
            weights_jet_multiplicity_2022EE["5"] = 1 if sumWeight_jet_multiplicity_5_2022EE_after == 0 else sumWeight_jet_multiplicity_5_2022EE_before / \
                sumWeight_jet_multiplicity_5_2022EE_after
            weights_jet_multiplicity_2022EE["6"] = 1 if sumWeight_jet_multiplicity_6_2022EE_after == 0 else sumWeight_jet_multiplicity_6_2022EE_before / \
                sumWeight_jet_multiplicity_6_2022EE_after
            weights_jet_multiplicity_2022EE["7"] = 1 if sumWeight_jet_multiplicity_7_2022EE_after == 0 else sumWeight_jet_multiplicity_7_2022EE_before / \
                sumWeight_jet_multiplicity_7_2022EE_after
            weights_jet_multiplicity_2022EE["8"] = 1 if sumWeight_jet_multiplicity_8_2022EE_after == 0 else sumWeight_jet_multiplicity_8_2022EE_before / \
                sumWeight_jet_multiplicity_8_2022EE_after
            weights_jet_multiplicity_2022EE["9"] = 1 if sumWeight_jet_multiplicity_9_2022EE_after == 0 else sumWeight_jet_multiplicity_9_2022EE_before / \
                sumWeight_jet_multiplicity_9_2022EE_after
            weights_jet_multiplicity_2022EE["10"] = 1 if sumWeight_jet_multiplicity_10_2022EE_after == 0 else sumWeight_jet_multiplicity_10_2022EE_before / \
                sumWeight_jet_multiplicity_10_2022EE_after

            weights_jet_multiplicity["2022"] = weights_jet_multiplicity_2022
            weights_jet_multiplicity["2022EE"] = weights_jet_multiplicity_2022EE

            dic_weights["MC"] = mc_weights_dic
            dic_weights["JetMulti"] = weights_jet_multiplicity

        dic_MC = weights_jet_multiplicity.copy()

        def _convert_DictToJSON(self, leptontau_btag_rescale_weights_shape):
            import os.path
            from correctionlib.schemav2 import VERSION, Correction, Variable, Category, CorrectionSet
            from correctionlib.JSONEncoder import write

            if not os.path.isdir(os.path.join(workdir, 'data')):
                os.makedirs(os.path.join(workdir, 'data'))

            ratio_correction_btagg_path = os.path.join(workdir, 'data')

            inputs = [
                Variable(name="year", type="string",
                         description="Year: 2022(preEE),2022EE(posEE)"),
                Variable(name="jet_multiplicity", type="string",
                         description="Jet Mulitplicity")
            ]

            output = Variable(name="ratio", type="real",
                              description="Ratio to correct the b-tag SF shape")

            def _get_DataContent(leptontau_btag_rescale_weights_shape):

                data_content = Category.parse_obj({
                    "nodetype": "category",
                    "input": "year",
                    "content": [
                        {
                            "key": year,
                            "value": Category.parse_obj({
                                "nodetype": "category",
                                "input": "jet_multiplicity",
                                "content": [
                                    {
                                        "key": multiplicity,
                                        "value": ratio,
                                    } for multiplicity, ratio, in multiplicity_ratio_dic.items()
                                ]
                            }),

                        } for year, multiplicity_ratio_dic in leptontau_btag_rescale_weights_shape.items()
                    ]
                })

                return data_content

            corr = Correction.parse_obj({
                "version": 1,
                "name": "Ratio_btagSF_shape",
                "description": "Ratio correction for the b-tag SF shape",
                "inputs": inputs,
                "output": output,
                "data": _get_DataContent(leptontau_btag_rescale_weights_shape)
            })

            correction_set = CorrectionSet(
                schema_version=VERSION, corrections=[corr])
            write(correction_set, f"{ratio_correction_btagg_path}/RatioCorr_btagShapeSF_{self.channel}ch.json.gz",
                  sort_keys=True, indent=2, maxlistlen=25, maxdictlen=3, breakbrackets=False)

        _convert_DictToJSON(self, dic_MC)
