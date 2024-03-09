import os
import re
import logging
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection

logger = logging.getLogger(__name__)

JECTagDatabase = {
    "2022": {
        "MC": "Summer22_22Sep2023_V2_MC",
        "C": "Summer22_22Sep2023_RunCD_V2_DATA",
        "D": "Summer22_22Sep2023_RunCD_V2_DATA"},
    "2022EE": {
        "MC": "Summer22EE_22Sep2023_V2_MC",
        "E": "Summer22EE_22Sep2023_RunE_V2_DATA",
        "F": "Summer22EE_22Sep2023_RunF_V2_DATA",
        "G": "Summer22EE_22Sep2023_RunG_V2_DATA"},
}

JERTagDatabase = {
    "2022": "Summer22EEPrompt22_JRV1_MC",
    "2022EE": "Summer22EEPrompt22_JRV1_MC",
}

jsonPathBase = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"

puWeightsTuple = {
    "2022": (jsonPathBase + "LUM/2022_Summer22/puWeights.json.gz", "Collisions2022_355100_357900_eraBCD_GoldenJson"),
    "2022EE": (jsonPathBase + "LUM/2022_Summer22EE/puWeights.json.gz", "Collisions2022_359022_362760_eraEFG_GoldenJson"),
}

JEC_JSONFiles = {
    "2022": {
        "AK4": jsonPathBase + "JME/2022_Summer22/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22/fatJet_jerc.json.gz"},
    "2022EE": {
        "AK4": jsonPathBase + "JME/2022_Summer22EE/jet_jerc.json.gz",
        "AK8": jsonPathBase + "JME/2022_Summer22EE/fatJet_jerc.json.gz"},
}

BTV_SF_JSONFiles = {
    "2022": jsonPathBase + "BTV/2022_Summer22/btagging.json.gz",
    "2022EE": jsonPathBase + "BTV/2022_Summer22EE/btagging.json.gz",
}

MUO_SF_JSONFiles = {
    "2022": jsonPathBase + "MUO/2022_27Jun2023/muon_Z.json.gz",
    "2022EE": jsonPathBase + "MUO/2022EE_27Jun2023/muon_Z.json.gz",
}

EL_SF_JSONFileDirs = {
    "2022": "2022Re-recoBCD",
    "2022EE": "2022Re-recoE+PromptFG",
}


def getRunEra(sample):
    """Return run era (A/B/...) for data sample"""
    result = re.search(r'Run20..([A-Z]?)', sample)
    if result is None:
        return "MC"
    else:
        return result.group(1)


class ScaleFactors():
    """Class to handle scale factors"""

    def commonSF(self, tree, sel, sample):
        # MC weight
        if self.is_MC:
            logger.info("Applying genWeight")
            sel = sel.refine('genWeight', weight=tree.genWeight)
        self.yields.add(sel, "genWeight")

        # PU weight
        if self.is_MC:
            from bamboo.analysisutils import makePileupWeight
            pileupWeight = makePileupWeight(
                puWeightsTuple[self.era], tree.Pileup_nTrueInt, systName="pileup", sel=sel)
            logger.info("Applying PU weight")
            sel = sel.refine('puWeight', weight=pileupWeight)
        self.yields.add(sel, "puWeight")
        # Triggers
        self.triggersPerPrimaryDataset = {}

        def addHLTPath(PD, HLT):
            if PD not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[PD] = []
            try:
                self.triggersPerPrimaryDataset[PD].append(
                    getattr(tree.HLT, HLT))
            except AttributeError:
                print("Couldn't find branch tree.HLT.%s, cross check!" % HLT)

        if self.era == '2022':
            if sample.startswith("SingleMuon_") or sample.startswith("DoubleMuon_"):
                addHLTPath("DoubleMuon_", "TripleMu_12_10_5")
                addHLTPath("DoubleMuon_",
                           "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
                addHLTPath("SingleMuon_", "IsoMu24")
                addHLTPath("SingleMuon_", "IsoMu27")
            else:
                addHLTPath("Muon_", "TripleMu_12_10_5")
                addHLTPath("Muon_",
                           "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
                addHLTPath("Muon_", "IsoMu24")
                addHLTPath("Muon_", "IsoMu27")
            addHLTPath("EGamma_", "Ele32_WPTight_Gsf")
            addHLTPath("EGamma_", "Ele16_Ele12_Ele8_CaloIdL_TrackIdL")
            addHLTPath("EGamma_", "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG_",
                       "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG_", "Mu8_DiEle12_CaloIdL_TrackIdL")
            addHLTPath("MuonEG_", "DiMu9_Ele9_CaloIdL_TrackIdL_DZ")

        if self.era == '2022EE':
            if sample.startswith("SingleMuon_") or sample.startswith("DoubleMuon_"):
                addHLTPath("DoubleMuon_", "TripleMu_12_10_5")
                addHLTPath("DoubleMuon_",
                           "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
                addHLTPath("SingleMuon_", "IsoMu24")
                addHLTPath("SingleMuon_", "IsoMu27")
            else:
                addHLTPath("Muon_", "TripleMu_12_10_5")
                addHLTPath("Muon_",
                           "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
                addHLTPath("Muon_", "IsoMu24")
                addHLTPath("Muon_", "IsoMu27")
            addHLTPath("EGamma_", "Ele16_Ele12_Ele8_CaloIdL_TrackIdL")
            addHLTPath("MuonEG_", "DiMu9_Ele9_CaloIdL_TrackIdL_DZ")
            addHLTPath("MuonEG_", "Mu8_DiEle12_CaloIdL_TrackIdL")
            addHLTPath("EGamma_", "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG_",
                       "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("EGamma_", "Ele32_WPTight_Gsf")

        if self.is_MC:
            sel = sel.refine('trigger',  cut=(
                op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values()))))
        else:
            sel = sel.refine('trigger', cut=makeMultiPrimaryDatasetTriggerSelection(
                sample, self.triggersPerPrimaryDataset))

        self.yields.add(sel, "trigger")

        # top pt reweighting
        if self.is_MC and sample.startswith("TT"):
            def top_pt_weight(pt):
                return op.exp(-2.02274e-01 + 1.09734e-04*pt + -1.30088e-07*pt**2 + (5.83494e+01/(pt+1.96252e+02)))

            def getTopPtWeight(tree):
                lastCopy = op.select(
                    tree.GenPart, lambda p: (op.static_cast("int", p.statusFlags) >> 13) & 1)
                tops = op.select(lastCopy, lambda p: p.pdgId == 6)
                antitops = op.select(lastCopy, lambda p: p.pdgId == -6)
                weight = op.switch(op.AND(op.rng_len(tops) >= 1, op.rng_len(antitops) >= 1),
                                   op.sqrt(top_pt_weight(
                                       tops[0].pt) * top_pt_weight(antitops[0].pt)),
                                   1.)
                return weight

            logger.info(
                "Applying Top Pt reweighting (only for TTbar samples)")

            sel = sel.refine("topPt", weight=op.systematic(
                getTopPtWeight(tree), noTopPt=op.c_float(1.)))
        else:
            sel = sel.refine("topPt", weight=op.c_float(1.))
        self.yields.add(sel, "topPt reweighting")
        return sel

    def btagSF(self, sel):
        # btagging SF
        if self.is_MC:
            from bamboo.scalefactors import get_bTagSF_itFit, makeBtagWeightItFit
            logger.info("Applying btagging SF for "+sel.name)
            def btvSF(flav): return get_bTagSF_itFit(
                BTV_SF_JSONFiles[self.era], "particleNet", "btagPNetB", flav, sel)
            btvWeight = makeBtagWeightItFit(self.ak4Jets, btvSF)
            sel = sel.refine(sel.name+"_btagSF", weight=btvWeight)
        else:
            sel = sel.refine(sel.name+"_btagSF", weight=op.c_float(1.))
        self.yields.add(sel, "btagging SF")
        return sel

    def muonSF(self, sel):
        # Muon SF
        if self.is_MC:
            from bamboo.scalefactors import get_correction
            logger.info("Applying Muon SF for "+sel.name)
            muonIDSF = get_correction(
                MUO_SF_JSONFiles[self.era],
                "NUM_LooseID_DEN_TrackerMuons",
                params={"pt": lambda mu: mu.pt,
                        "abseta": lambda mu: op.abs(mu.eta)},
                systParam="scale_factors",
                systNomName="nominal",
                systName="syst",
                sel=sel
            )
            if sel.name in ['DL_boosted_mumu', 'DL_resolved_mumu']:
                # pt and eta cut here since correction are available only for pt > 15 and |eta| < 2.4
                sel = sel.refine(sel.name+"_muonSF", cut=[
                    op.AND(self.firstOSMuMu[0].pt >= 15,
                           self.firstOSMuMu[1].pt >= 15,
                           op.abs(
                        self.firstOSMuMu[0].eta) < 2.4,
                        op.abs(
                        self.firstOSMuMu[1].eta) < 2.4
                    )],
                    weight=[muonIDSF(self.firstOSMuMu[0]),
                            muonIDSF(self.firstOSMuMu[1])]
                )
            elif sel.name in ['DL_boosted_emu', 'DL_resolved_emu']:
                sel = sel.refine(sel.name+"_muonSF",
                                 cut=[op.AND(self.firstOSElMu[1].pt >= 15,
                                             op.abs(self.firstOSElMu[1].eta) < 2.4)],
                                 weight=muonIDSF(self.firstOSElMu[1])
                                 )
            else:
                sel = sel.refine(sel.name+"_muonSF", weight=op.c_float(1.))
        self.yields.add(sel, "muon SF")
        return sel

    def electronSF(self, sel):
        # Electron SF
        if self.is_MC:
            from bamboo.scalefactors import get_correction
            logger.info("Applying Electron SF for "+sel.name)

            os.mkdir(
                "2022Re-recoBCD") if not os.path.exists("2022Re-recoBCD") else None
            os.mkdir(
                "2022Re-recoE+PromptFG") if not os.path.exists("2022Re-recoE+PromptFG") else None
            os.system("xrdcp root://cms-xrd-global.cern.ch///store/group/phys_egamma/correctionlibJSONs/Run3_2022_recoBCDE_PromptFG/2022Re-recoBCD/electron.json.gz 2022Re-recoBCD/") if not os.path.exists(
                "2022Re-recoBCD/electron.json.gz") else None
            os.system("xrdcp root://cms-xrd-global.cern.ch///store/group/phys_egamma/correctionlibJSONs/Run3_2022_recoBCDE_PromptFG/2022Re-recoE+PromptFG/electron.json.gz 2022Re-recoE+PromptFG/") if not os.path.exists(
                "2022Re-recoE+PromptFG/electron.json.gz") else None

            electronIDSF = get_correction(
                EL_SF_JSONFileDirs[self.era]+"/electron.json.gz",
                "Electron-ID-SF",
                params={"pt": lambda ele: ele.pt, "eta": lambda ele: ele.eta,
                        "year": EL_SF_JSONFileDirs[self.era], "WorkingPoint": "Loose"},
                systParam="ValType",
                systNomName="sf",
                sel=sel
            )
            if sel.name in ['DL_boosted_ee', 'DL_resolved_ee']:
                # pt cut here since correction are available only for pt > 10
                sel = sel.refine(sel.name+"_electronSF", cut=[
                    op.AND(self.firstOSElEl[0].pt >= 10,
                           self.firstOSElEl[1].pt >= 10
                           )],
                    weight=[electronIDSF(self.firstOSElEl[0]),
                            electronIDSF(self.firstOSElEl[1])]
                )
            elif sel.name in ['DL_boosted_emu', 'DL_resolved_emu']:
                sel = sel.refine(sel.name+"_electronSF",
                                 cut=[self.firstOSElMu[0].pt >= 10],
                                 weight=electronIDSF(self.firstOSElMu[0])
                                 )
            else:
                sel = sel.refine(sel.name+"_electronSF", weight=op.c_float(1.))
        self.yields.add(sel, "electron SF")
        return sel
