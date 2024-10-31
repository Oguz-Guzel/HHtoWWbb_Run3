import re
import logging

from bamboo import treefunctions as op
from bamboo.scalefactors import get_correction

logger = logging.getLogger(__name__)

jsonPathBase = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/"

BTV_SF_JSONFiles = {
    "2022": jsonPathBase + "BTV/2022_Summer22/btagging.json.gz",
    "2022EE": jsonPathBase + "BTV/2022_Summer22EE/btagging.json.gz",
    "2023": jsonPathBase + "BTV/2023_Summer23/btagging.json.gz",
    "2023BPix": jsonPathBase + "BTV/2023_Summer23BPix/btagging.json.gz",
}

MUO_SF_JSONFiles = {
    "2022": jsonPathBase + "MUO/2022_Summer22/muon_Z.json.gz",
    "2022EE": jsonPathBase + "MUO/2022_Summer22EE/muon_Z.json.gz",
    "2023": jsonPathBase + "MUO/2023_Summer23/muon_Z.json.gz",
    "2023BPix": jsonPathBase + "MUO/2023_Summer23BPix/muon_Z.json.gz",
}

EL_SF_JSONFiles = {
    "2022": (jsonPathBase + "EGM/2022_Summer22/electron.json.gz", "2022Re-recoBCD"),
    "2022EE": (jsonPathBase + "EGM/2022_Summer22EE/electron.json.gz", "2022Re-recoE+PromptFG"),
    "2023": (jsonPathBase + "EGM/2023_Summer23/electron.json.gz", "2023PromptC"),
    "2023BPix": (jsonPathBase + "EGM/2023_Summer23BPix/electron.json.gz", "2023PromptD"),
}

jetVeto_JSONFiles = {
    "2022": (jsonPathBase + "JME/2022_Summer22/jetvetomaps.json.gz", "Summer22_23Sep2023_RunCD_V1"),
    "2022EE": (jsonPathBase + "JME/2022_Summer22EE/jetvetomaps.json.gz", "Summer22EE_23Sep2023_RunEFG_V1",),
    "2023": (jsonPathBase + "JME/2023_Summer23/jetvetomaps.json.gz", "Summer23Prompt23_RunC_V1"),
    "2023BPix": (jsonPathBase + "JME/2023_Summer23BPix/jetvetomaps.json.gz", "Summer23BPixPrompt23_RunD_V1"),
}


def getRunEra(sample):
    """Return run era (A/B/...) for data sample"""
    result = re.search(r'Run20..([A-Z]?)', sample)
    if result is None:
        return "MC"
    else:
        return result.group(1)


class ScaleFactors():
    """Class to define scale factors"""

    def commonSF(self, tree, sel, sample):
        isMC = self.is_MC
        # top pt reweighting
        if isMC and sample.startswith("TT"):
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
                getTopPtWeight(tree)))
        else:
            sel = sel.refine("topPt", weight=op.c_float(1.))
        self.yields.add(sel, "topPt reweighting")

        # btagging SF
        if isMC:
            from bamboo.scalefactors import get_bTagSF_itFit, makeBtagWeightItFit
            logger.info("Applying btagging SF")

            def btvSF(flav): return get_bTagSF_itFit(
                BTV_SF_JSONFiles[self.era], "particleNet", "btagPNetB", flav, sel, decorr_eras=True, era=self.era)
            btvWeight = makeBtagWeightItFit(self.ak4Jets, btvSF)
        else:
            btvWeight = op.c_float(1.)

        sel = sel.refine("btagSF", weight=btvWeight)

        self.yields.add(sel, "btagging SF")

        # btag reweighting

        if isMC:
            logger.info("Applying btag reweighting")
            btag_corr = get_correction(
                f"{self.git_project_dir}/data/{self.era[:4]}_btagSF_reweight.json.gz",
                "Ratio_btagSF_shape",
                params={
                    "year": self.era,
                    # 1. to make it a float
                    "jet_multiplicity": 1.*op.rng_len(self.ak4Jets)
                },
                sel=sel
            )
        # None since the object is already in the btag_corr i.e. self.ak4Jets
        btag_rescale = btag_corr(None) if isMC else op.c_float(1.)

        sel = sel.refine("btagReweight", weight=btag_rescale)

        self.yields.add(sel, "btag reweight")

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
                        "eta": lambda mu: op.abs(mu.eta)},
                systParam="scale_factors",
                systNomName="nominal",
                systName="syst",
                sel=sel
            )
            if sel.name in ['DL_boosted_mumu', 'DL_resolved_mumu']:
                # pt and eta cut here since correction are available only when pt >= 15 and |eta| < 2.4
                sel = sel.refine(sel.name+"_muonSF", cut=[
                    op.AND(self.firstMuTightPair[0].pt >= 15,
                           self.firstMuTightPair[1].pt >= 15,
                           op.abs(
                        self.firstMuTightPair[0].eta) < 2.4,
                        op.abs(
                        self.firstMuTightPair[1].eta) < 2.4
                    )],
                    weight=[muonIDSF(self.firstMuTightPair[0]),
                            muonIDSF(self.firstMuTightPair[1])]
                )
            elif sel.name in ['DL_boosted_emu', 'DL_resolved_emu']:
                sel = sel.refine(sel.name+"_muonSF",
                                 cut=[op.AND(self.firstEmuTightPair[1].pt >= 15,
                                             op.abs(self.firstEmuTightPair[1].eta) < 2.4)],
                                 weight=muonIDSF(self.firstEmuTightPair[1])
                                 )
            else:
                sel = sel.refine(sel.name+"_muonSF", weight=op.c_float(1.))
        else:
            sel = sel.refine(sel.name+"_muonSF", weight=op.c_float(1.))
        # self.yields.add(sel, "muon SF")
        return sel

    def electronSF(self, sel):
        # Electron SF
        if self.is_MC:
            from bamboo.scalefactors import get_correction
            logger.info("Applying Electron SF for "+sel.name)

            params = {"pt": lambda e: e.pt,
                      "eta": lambda e: e.eta,
                      "year": EL_SF_JSONFiles[self.era][1],
                      "WorkingPoint": "Loose"}
            if self.era in ['2023', '2023BPix']:
                params["phi"] = lambda e: e.phi

            electronIDSF = get_correction(
                EL_SF_JSONFiles[self.era][0],
                "Electron-ID-SF",
                params=params,
                systParam="ValType",
                systNomName="sf",
                sel=sel
            )
            if sel.name in ['DL_boosted_ee', 'DL_resolved_ee']:
                # pt cut here since correction's available only when pt >= 10
                sel = sel.refine(sel.name+"_electronSF", cut=[
                    op.AND(self.firstElTightPair[0].pt >= 10,
                           self.firstElTightPair[1].pt >= 10
                           )],
                    weight=[electronIDSF(self.firstElTightPair[0]),
                            electronIDSF(self.firstElTightPair[1])]
                )
            elif sel.name in ['DL_boosted_emu', 'DL_resolved_emu']:
                sel = sel.refine(sel.name+"_electronSF",
                                 cut=[self.firstEmuTightPair[0].pt >= 10],
                                 weight=electronIDSF(self.firstEmuTightPair[0])
                                 )
            else:
                sel = sel.refine(sel.name+"_electronSF", weight=op.c_float(1.))
        else:
            sel = sel.refine(sel.name+"_electronSF", weight=op.c_float(1.))
        # self.yields.add(sel, "electron SF")
        return sel
