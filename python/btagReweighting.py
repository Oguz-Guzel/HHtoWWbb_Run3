
from bamboo.plots import Skim, CutFlowReport
from bamboo import treefunctions as op

from bamboo.analysismodules import NanoAODModule, HistogramsModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection

import definitions as defs
from selections import makeDLSelection
from scalefactors import ScaleFactors as sf

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

        # Define the git project's directory
        self.git_project_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..'))

        self.era = sampleCfg["era"] if sampleCfg else None
        self.is_MC = self.isMC(sample)

        # Decorate the tree
        from bamboo.treedecorators import NanoAODDescription
        tree, noSel, be, lumiArgs = super().prepareTree(
            tree, sample=sample, sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v12", year=self.era[:4], isMC=self.is_MC),
            backend=self.args.backend or backend)

        # MC weight
        genWeight = tree.genWeight if self.is_MC else op.c_float(1.)
        noSel = noSel.refine('genWeight', weight=genWeight)

        # initialise CFR
        self.yields = CutFlowReport(
            "yields", recursive=False, printInLog=True)

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

        return tree, noSel, be, lumiArgs


class btagReweighting(_base):
    """ Class to create the ratio that sould be applied to as a weight 
        before making any btagged jet selection, see
        https://btv-wiki.docs.cern.ch/PerformanceCalibration/shapeCorrectionSFRecommendations/#effect-on-event-yields"""

    def __init__(self, args):
        super(btagReweighting, self).__init__(args)
        self.channel = self.args.channel

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # cutflow report
        plots.append(self.yields)

        # define objects
        defs.defineObjects(self, tree)

        # common scalefactors
        noSel = sf.top_pT_reweight(self, tree, noSel, sample)

        _, pre_sels = makeDLSelection(
            self, noSel)

        # pre_sels is [DL_boosted_pre_ee, DL_boosted_pre_mumu, DL_boosted_pre_emu, DL_resolved_pre_ee, DL_resolved_pre_mumu, DL_resolved_pre_emu]

        for ps in pre_sels:
            # add weights before btagSF
            plots.append(
                Skim("WeightsBeforeBtagSF_"+ps.name,
                     {"jetMultiplicity_before": op.rng_len(self.ak4Jets),
                         "Sel_weight_before": ps.weight},
                     ps
                     )
            )

            # apply btagSF
            ps = sf.btagSF(self, ps)

            # add weights after btagSF
            plots.append(
                Skim("WeightsAfterBtagSF_"+ps.name,
                     {"jetMultiplicity_after": op.rng_len(self.ak4Jets),
                      "Sel_weight_after": ps.weight},
                     ps
                     )
            )
        
        return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super().postProcess(taskList, config=config, workdir=workdir, resultsdir=resultsdir)

        if not self.plotList:
            self.plotList = self.getPlotList(
                resultsdir=resultsdir, config=config)

        skim_list = [ap for ap in self.plotList if isinstance(ap, Skim)]

        # Fallback if no eras are specified
        _, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())

        # Initialize dictionaries for each era
        sumW_perEra = {era: {"before": [0]*11,
                             "after": [0]*11} for era in eras}

        def accumulate_weights(tree, weight_branch, multiplicity_branch, accumulator):
            for entry in tree:
                jm = min(getattr(entry, multiplicity_branch), 10)
                accumulator[jm] += getattr(entry, weight_branch)

        def _openFileAndGet(path, mode="read"):
            tf = TFile.Open(path, mode)
            if not tf or not tf.IsOpen():
                raise Exception(f"Could not open file {path}")
            return tf

        for proc, smpCfg in config["samples"].items():
            if smpCfg.get("group") == "data":
                continue

            sample_rootfile = _openFileAndGet(
                os.path.join(resultsdir, f"{proc}.root"), "read")
            genEvents = self.readCounters(sample_rootfile)[
                smpCfg["generated-events"]]
            lumi = config["eras"][smpCfg["era"]]["luminosity"]
            Xsection = smpCfg["cross-section"]
            smpScale = lumi * Xsection / genEvents
            era = smpCfg["era"]

            nJet_before = [0]*11
            nJet_after = [0]*11

            for skim in skim_list:
                tree = sample_rootfile.Get(skim.treeName)
                if not tree:
                    logger.info(
                        f"Warning: skim tree {skim.treeName} not found in file {sample_rootfile.GetName()}")
                    continue

                branch_names = [branch.GetName()
                                for branch in tree.GetListOfBranches()]

                if "Sel_weight_after" in branch_names:
                    accumulate_weights(
                        tree, "Sel_weight_after", "jetMultiplicity_after", nJet_after)
                elif "Sel_weight_before" in branch_names:
                    accumulate_weights(
                        tree, "Sel_weight_before", "jetMultiplicity_before", nJet_before)

            for i in range(11):
                sumW_perEra[era]["before"][i] += nJet_before[i] * smpScale
                sumW_perEra[era]["after"][i] += nJet_after[i] * smpScale

            sample_rootfile.Close()

        weights_per_nJet = {
            era: {
                i: 1 if sumW_perEra[era]["after"][i] == 0 else sumW_perEra[era]["before"][i] /
                sumW_perEra[era]["after"][i]
                for i in range(11)
            }
            for era in eras
        }

        dict_MC = weights_per_nJet.copy()
        self._convert_DictToJSON(dict_MC, workdir)

    def _convert_DictToJSON(self, ratio_dict, workdir):
        from correctionlib.schemav2 import VERSION, Correction, Variable, Category, CorrectionSet, Binning
        from correctionlib.JSONEncoder import write

        data_dir = os.path.join(workdir, 'data')
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir)

        inputs = [
            Variable(name="year", type="string", description="Year-based era"),
            Variable(name="jet_multiplicity", type="real",
                     description="Jet Multiplicity")
        ]
        output = Variable(name="ratio", type="real",
                          description="Ratio to correct the b-tag SF shape")

        def _get_DataContent(ratio_dict):
            def dict_to_list(dict_):
                lm = []
                lr = []
                for m, r in dict_.items():
                    lm.append(m-0.5)
                    lr.append(r)
                lm.append(list(dict_)[-1]+0.5)
                return [lm, lr]

            return Category.parse_obj({
                "nodetype": "category",
                "input": "year",
                "content": [
                    {
                        "key": year,
                        "value": Binning.parse_obj({
                            "nodetype": "binning",
                            "input": "jet_multiplicity",
                            "edges": dict_to_list(multiplicity_ratio_dict)[0],
                            "content": dict_to_list(multiplicity_ratio_dict)[1],
                            "flow": "clamp"
                        })
                    } for year, multiplicity_ratio_dict in ratio_dict.items()
                ]
            })

        corr = Correction.parse_obj({
            "version": 1,
            "name": "Ratio_btagSF_shape",
            "description": "Ratio correction for the b-tag SF shape",
            "inputs": inputs,
            "output": output,
            "data": _get_DataContent(ratio_dict=ratio_dict)
        })

        correction_set = CorrectionSet(
            schema_version=VERSION, corrections=[corr])
        output_file = os.path.join(data_dir, f"btagSF_rescaling.json.gz")
        write(correction_set, output_file, sort_keys=True, indent=2,
              maxlistlen=25, maxdictlen=3, breakbrackets=False)
        logger.info(f'written to: {output_file}')
