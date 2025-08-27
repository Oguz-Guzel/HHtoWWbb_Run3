import os
import re
import logging
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.analysismodules import NanoAODHistoModule

from bamboo.plots import EquidistantBinning as EqBin

import utils

logger = logging.getLogger(__name__)

jsonPathBase = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"

PU_JSONFiles = {
    "2022": (
        jsonPathBase + "LUM/2022_Summer22/puWeights.json.gz",
        "Collisions2022_355100_357900_eraBCD_GoldenJson",
    ),
    "2022EE": (
        jsonPathBase + "LUM/2022_Summer22EE/puWeights.json.gz",
        "Collisions2022_359022_362760_eraEFG_GoldenJson",
    ),
    "2023": (
        jsonPathBase + "LUM/2023_Summer23/puWeights.json.gz",
        "Collisions2023_366403_369802_eraBC_GoldenJson",
    ),
    "2023BPix": (
        jsonPathBase + "LUM/2023_Summer23BPix/puWeights.json.gz",
        "Collisions2023_369803_370790_eraD_GoldenJson",
    ),
}

JECTags = {
    "2022": {
        "MC": "Summer22_22Sep2023_V2_MC",
        "C": "Summer22_22Sep2023_RunCD_V2_DATA",
        "D": "Summer22_22Sep2023_RunCD_V2_DATA",
    },
    "2022EE": {
        "MC": "Summer22EE_22Sep2023_V2_MC",
        "E": "Summer22EE_22Sep2023_RunE_V2_DATA",
        "F": "Summer22EE_22Sep2023_RunF_V2_DATA",
        "G": "Summer22EE_22Sep2023_RunG_V2_DATA",
    },
    "2023": {
        "MC": "Summer23Prompt23_V2_MC",
        "C": "Summer23Prompt23_V2_DATA",
    },
    "2023BPix": {
        "MC": "Summer23BPixPrompt23_V3_MC",
        "D": "Summer23BPixPrompt23_V3_DATA",
    },
}

JERTags = {
    "2022": "Summer22_22Sep2023_JRV1_MC",
    "2022EE": "Summer22EE_22Sep2023_JRV1_MC",
    "2023": "Summer23Prompt23_RunCv1234_JRV1_MC",
    "2023BPix": "Summer23BPixPrompt23_RunD_JRV1_MC",
}

JEC_JSONFiles = {
    "2022": {
        "AK4": jsonPathBase + "JME/2022_Summer22/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22/fatJet_jerc.json.gz",
    },
    "2022EE": {
        "AK4": jsonPathBase + "JME/2022_Summer22EE/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22EE/fatJet_jerc.json.gz",
    },
    "2023": {
        "AK4": jsonPathBase + "JME/2023_Summer23/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2023_Summer23/fatJet_jerc.json.gz",
    },
    "2023BPix": {
        "AK4": jsonPathBase + "JME/2023_Summer23BPix/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2023_Summer23BPix/fatJet_jerc.json.gz",
    },
}

EGamma_SS_SF_JSONFiles = {
    "2022": (
        jsonPathBase + "EGM/2022_Summer22/electronSS_EtDependent.json.gz",
        "EGMScale_Compound_Ele_2022preEE",
        "EGMSmearAndSyst_ElePTsplit_2022preEE",
    ),
    "2022EE": (
        jsonPathBase + "EGM/2022_Summer22EE/electronSS_EtDependent.json.gz",
        "EGMScale_Compound_Ele_2022postEE",
        "EGMSmearAndSyst_ElePTsplit_2022postEE",
    ),
    "2023": (
        jsonPathBase + "EGM/2023_Summer23/electronSS_EtDependent.json.gz",
        "EGMScale_Compound_Ele_2023preBPIX",
        "EGMSmearAndSyst_ElePTsplit_2023preBPIX",
    ),
    "2023BPix": (
        jsonPathBase + "EGM/2023_Summer23BPix/electronSS_EtDependent.json.gz",
        "EGMScale_Compound_Ele_2023postBPIX",
        "EGMSmearAndSyst_ElePTsplit_2023postBPIX",
    ),
}


def getDataRunEra(sample):
    """Return run era (A/B/...) and the following digits for data sample"""
    result = re.search(r"Run20\d{2}([A-Z]\w*)", sample)
    return result.group(1) if result else None


class _base(NanoAODHistoModule):
    """Base module for HH->WWbb analysis"""

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument(
            "--backend",
            type=str,
            default="dataframe",
            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'",
        )

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):

        # Define the git project's directory
        self.git_project_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..")
        )
        self.era = sampleCfg["era"] if sampleCfg else None
        self.is_MC = self.isMC(sample)

        from bamboo.plots import CutFlowReport

        self.yields = CutFlowReport("yields", recursive=True, printInLog=False)

        # Decorate the tree
        from bamboo.treedecorators import (
            NanoAODDescription,
            nanoFatJetCalc,
            CalcCollectionsGroups,
        )

        metName = "PuppiMET"
        nanoJetMETCalc_both = CalcCollectionsGroups(
            Jet=("pt", "mass"),
            changes={metName: (f"{metName}T1", f"{metName}T1Smear")},
            **{metName: ("pt", "phi")},
        )
        nanoJetMETCalc_data = CalcCollectionsGroups(
            Jet=("pt", "mass"),
            changes={metName: (f"{metName}T1",)},
            **{metName: ("pt", "phi")},
        )
        nanoElectronCalc = CalcCollectionsGroups(Electron=("nElectron"))
        systVariations = ([nanoFatJetCalc, nanoElectronCalc]) + [
            nanoJetMETCalc_both if self.is_MC else nanoJetMETCalc_data
        ]
        tree, noSel, be, lumiArgs = super().prepareTree(
            tree,
            sample=sample,
            sampleCfg=sampleCfg,
            description=NanoAODDescription.get(
                "v12", year=self.era[:4], isMC=self.is_MC, systVariations=systVariations
            ),
            backend=backend,
        )

        # Number of events before any processing
        self.yields.add(noSel, "noSel")

        # JEC/JER
        jecTag = JECTags[self.era]["MC" if self.is_MC else getDataRunEra(sample)]
        logger.info(f"JEC tag for sample {sample} is {jecTag}")
        # smearing is for MC only
        smearTag = JERTags[self.era] if self.is_MC else None
        if smearTag:
            logger.info(f"JER tag for sample {sample} is {smearTag}")

        jecArgs = {
            "jsonFile": JEC_JSONFiles[self.era]["AK4"],
            "jec": jecTag,
            "smear": smearTag,
            "splitJER": True,
            "jsonFileSmearingTool": jsonPathBase + "JME/jer_smear.json.gz",
            "jesUncertaintySources": (["Total"] if self.is_MC else None),
            "isMC": self.is_MC,
            "backend": be,
        }

        from bamboo.analysisutils import (
            configureJets,
            configureType1MET,
            configureElectrons,
        )

        configureJets(tree._Jet, jetType="AK4PFPuppi", **jecArgs)

        metName = "PuppiMET"
        configureType1MET(
            getattr(tree, f"_{metName}T1"),
            enableSystematics=(
                (lambda v: not v.startswith("jer")) if self.is_MC else None
            ),
            **jecArgs,
        )

        jecArgs.update(
            {
                "jsonFile": JEC_JSONFiles[self.era]["AK8"],
                "jetAlgoSubjet": "AK4PFPuppi",
                "jecSubjet": jecTag,
                "jsonFileSubjet": JEC_JSONFiles[self.era]["AK4"],
            }
        )
        configureJets(tree._FatJet, jetType="AK8PFPuppi", **jecArgs)
        logger.info("Applying Jet energy and resolution corrections")

        jsonFileRandomGenerator = os.path.join(
            self.git_project_dir, "../bamboo/tests/data/randomNumbers.json.gz"
        )
        configureElectrons(
            tree._Electron,
            paramsFile=EGamma_SS_SF_JSONFiles[self.era][0],
            scale=EGamma_SS_SF_JSONFiles[self.era][1],
            smearing=EGamma_SS_SF_JSONFiles[self.era][2],
            jsonFileRandomGenerator=jsonFileRandomGenerator,
            addSystematics=True if self.is_MC else False,
            isMC=self.is_MC,
            backend=be,
        )
        logger.info("Applying Electron scale and smear (SS) corrections")

        # MC weight
        if self.is_MC:
            logger.info("Applying genWeight")
            noSel = noSel.refine("genWeight", weight=tree.genWeight)
        else:
            noSel = noSel.refine("genWeight", weight=op.c_float(1.0))
        self.yields.add(noSel, "genWeight")

        # PU weight
        if self.is_MC:
            from bamboo.analysisutils import makePileupWeight

            pileupWeight = makePileupWeight(
                PU_JSONFiles[self.era],
                tree.Pileup_nTrueInt,
                systName="pileup",
                sel=noSel,
            )
            logger.info("Applying PU weight")
            noSel = noSel.refine("puWeight", weight=pileupWeight)
        else:
            noSel = noSel.refine("puWeight", weight=op.c_float(1.0))
        self.yields.add(noSel, "puWeight")

        return tree, noSel, be, lumiArgs


class TriggerEff(_base):
    """
    Trigger efficiency module for di-lepton triggers in Run 3.
    This module calculates the trigger efficiency for di-lepton events
    and produces scale factors that can be used to correct event weights.
    It supports multiple lepton channels (ee, mumu, emu).
    """

    def __init__(self, args):
        super().__init__(args)

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot
        import definitions as defs
        from selections import makeDLSelection

        # call defined objects
        defs.defineObjects(self, tree)

        # get DL selections

        final_state_selections = makeDLSelection(self, noSel, tree, sample, trigger_study=True)

        # Triggers
        self.triggers_per_PD = {}

        def addHLTPath(PD, HLT):
            if PD not in self.triggers_per_PD.keys():
                self.triggers_per_PD[PD] = []
            try:
                self.triggers_per_PD[PD].append(getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, cross check!" % HLT)

        addHLTPath("Muon_", "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
        addHLTPath("EGamma_", "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
        addHLTPath("EGamma_", "DoubleEle33_CaloIdL_MW")
        addHLTPath("MuonEG_", "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
        addHLTPath("MuonEG_", "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
        addHLTPath("MuonEG_", "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")

        from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection

        if self.is_MC:
            trigger_cut = op.OR(*chain.from_iterable(self.triggers_per_PD.values()))
        else:
            trigger_cut = makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggers_per_PD
            )

        plots = []

        # --- add leading/subleading pT 2D efficiency histograms ---
        pt_bins_l1 = EqBin(20, 0, 300)  # 20 GeV bins up to 300 GeV
        pt_bins_l2 = EqBin(20, 0, 300)
        binning = (pt_bins_l1, pt_bins_l2)

        el1, el2 = self.tightElectrons[0], self.tightElectrons[1]
        mu1, mu2 = self.tightMuons[0], self.tightMuons[1]

        emu_leading_pt = op.switch(el1.pt > mu1.pt, el1.pt, mu1.pt)
        emu_subleading_pt = op.switch(el1.pt > mu1.pt, mu1.pt, el1.pt)

        di_leptons = {
            "_ee": (el1.pt, el2.pt),
            "_mumu": (mu1.pt, mu2.pt),
            "_emu": (emu_leading_pt, emu_subleading_pt),
        }
        fs = ""
        for v in ["den", "num"]:
            for selection in final_state_selections:
                for fs in di_leptons.keys():
                    if fs in selection.name:
                        fs = fs
                        break
                plots.append(
                    Plot.make2D(
                        f"{v}_{selection.name}",
                        di_leptons[fs],
                        (
                            selection
                            if v == "den"
                            else selection.refine(
                                f"{selection.name}_triggers{fs}", cut=trigger_cut
                            )
                        ),
                        binning,
                        title=f"{fs} " + "Denominator" if v == "den" else "Numerator",
                    )
                )

        return plots
