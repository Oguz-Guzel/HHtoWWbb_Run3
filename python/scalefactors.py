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

MUON_SF_JSONFiles = {
    "2022": jsonPathBase + "MUO/2022_Summer22/muon_Z.json.gz",
    "2022EE": jsonPathBase + "MUO/2022_Summer22EE/muon_Z.json.gz",
    "2023": jsonPathBase + "MUO/2023_Summer23/muon_Z.json.gz",
    "2023BPix": jsonPathBase + "MUO/2023_Summer23BPix/muon_Z.json.gz",
}

EGamma_SF_JSONFiles = {
    "2022": (jsonPathBase + "EGM/2022_Summer22/electron.json.gz", "2022Re-recoBCD"),
    "2022EE": (jsonPathBase + "EGM/2022_Summer22EE/electron.json.gz", "2022Re-recoE+PromptFG"),
    "2023": (jsonPathBase + "EGM/2023_Summer23/electron.json.gz", "2023PromptC"),
    "2023BPix": (jsonPathBase + "EGM/2023_Summer23BPix/electron.json.gz", "2023PromptD"),
}


class ScaleFactors():
    """Class to define scale factors"""

    def top_pT_reweight(self, GenPartBranch, sel, sample):
        """ Apply top p_T reweighting."""
        if sample.startswith("TT"):
            def top_pt_weight(pt):
                return op.exp(-2.02274e-01 + 1.09734e-04*pt + -1.30088e-07*pt**2 + (5.83494e+01/(pt+1.96252e+02)))

            def getTopPtWeight(GenPart):
                lastCopy = op.select(
                    GenPart, lambda p: (op.static_cast("int", p.statusFlags) >> 13) & 1)
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
                getTopPtWeight(GenPartBranch)))
        else:
            sel = sel.refine("topPt", weight=op.c_float(1.))
        self.yields.add(sel, "topPt reweighting")

        return sel

    def btagSF(self, sel, jets, json_tagger="particleNet", jet_tagger="btagPNetB", apply_reweighting=True):
        """Apply btagging SF"""
        if self.is_MC:
            from bamboo.scalefactors import get_bTagSF_itFit, makeBtagWeightItFit
            logger.info("Applying btagging SF for "+sel.name)
            def btvSF(flav): return get_bTagSF_itFit(
                BTV_SF_JSONFiles[self.era], json_tagger, jet_tagger, flav, sel=sel, decorr_eras=True, era=self.era)
            btvWeight = makeBtagWeightItFit(jets, btvSF)
            if apply_reweighting:
                btag_corr = get_correction(
                    f"{self.git_project_dir}/data/{self.era[:4]}_btagSF_reweight.json.gz", # the file from the btagReweighting.py
                    "Ratio_btagSF_shape",
                    params={
                        "year": self.era,
                        # 1. to make it a float
                        "jet_multiplicity": 1.*op.rng_len(jets)
                    },
                    sel=sel
                )
                # None since the object is already in the btag_corr i.e. self.ak4Jets
                btag_reweight = btag_corr(None)
            else: 
                btag_reweight = op.c_float(1.)
        else:
            btvWeight = op.c_float(1.)
            btag_reweight = op.c_float(1.)
        
        sel = sel.refine(sel.name+"_btagSF", weight=btvWeight*btag_reweight)

        self.yields.add(sel, sel.name)

        return sel

    def muonSF(self, sel):
        """Apply Muon SF"""
        if self.is_MC:
            from bamboo.scalefactors import get_correction
            logger.info("Applying Muon SF for "+sel.name)

            # Muon SF
            muon_ID_sf = get_correction(
                MUON_SF_JSONFiles[self.era],
                "NUM_MediumID_DEN_TrackerMuons",  # NUM_MediumPromptID_DEN_TrackerMuons, too ?
                params={"pt": lambda mu: mu.pt,
                        "eta": lambda mu: op.abs(mu.eta)},
                systParam="scale_factors",
                systNomName="nominal",
                systName="syst",
                sel=sel
            )

            # Muon Isolation SF
            muon_ISO_sf = get_correction(
                MUON_SF_JSONFiles[self.era],
                # since muon iso is miniPFreliso and id is medium
                "NUM_MediumMiniIso_DEN_MediumID",
                params={"pt": lambda mu: mu.pt,
                        "eta": lambda mu: op.abs(mu.eta),
                        },
                systParam="scale_factors",
                systNomName="nominal",
                systName="syst",
                sel=sel
            )

            # Muon Trigger SF
            muon_trigger_sf = get_correction(
                MUON_SF_JSONFiles[self.era],
                "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium",
                params={
                    "pt": lambda mu: op.max(mu.pt, 26.0),
                    "eta": lambda mu: op.abs(mu.eta)
                },
                systParam="scale_factors",
                systNomName="nominal",
                systName="syst",
                sel=sel
            )

            if sel.name == 'muPairMultiplicitySel':
                # pt and eta cut here since correction are available only when pt >= 15 and |eta| < 2.4
                sel = sel.refine(sel.name+"_muonSF", cut=[
                    op.AND(self.firstMuTightPair[0].pt >= 15,
                           self.firstMuTightPair[1].pt >= 15,
                           op.abs(
                        self.firstMuTightPair[0].eta) < 2.4,
                        op.abs(
                        self.firstMuTightPair[1].eta) < 2.4
                    )],
                    weight=[muon_ID_sf(self.firstMuTightPair[0]),
                            muon_ID_sf(self.firstMuTightPair[1]),
                            muon_trigger_sf(self.firstMuTightPair[0]),
                            muon_trigger_sf(self.firstMuTightPair[1]),
                            muon_ISO_sf(self.firstMuTightPair[0]),
                            muon_ISO_sf(self.firstMuTightPair[1])]
                )
            elif sel.name == 'emuPairMultiplicitySel':
                sel = sel.refine(sel.name+"_muonSF",
                                 cut=[op.AND(self.firstEmuTightPair[1].pt >= 15,
                                             op.abs(self.firstEmuTightPair[1].eta) < 2.4)],
                                 weight=[muon_ID_sf(self.firstEmuTightPair[1]),
                                         muon_trigger_sf(
                                             self.firstEmuTightPair[1]),
                                         muon_ISO_sf(self.firstEmuTightPair[1])]
                                 )
            else:
                sel = sel.refine(sel.name+"_muonSF", weight=op.c_float(1.))
        else:
            sel = sel.refine(sel.name+"_muonSF", weight=op.c_float(1.))
        
        return sel

    def electronSF(self, sel):
        """Apply Electron SF"""
        if self.is_MC:
            from bamboo.scalefactors import get_correction
            logger.info("Applying Electron SF for "+sel.name)

            params = {"pt": lambda e: e.pt,
                      "eta": lambda e: e.eta,
                      "year": EGamma_SF_JSONFiles[self.era][1],
                      "WorkingPoint": "wp90iso"}

            # add phi for 2023 and 2023BPix
            if self.era in ['2023', '2023BPix']:
                params["phi"] = lambda e: e.phi

            # Electron ID SF
            electron_ID_sf = get_correction(
                EGamma_SF_JSONFiles[self.era][0],
                "Electron-ID-SF",
                params=params,
                systParam="ValType",
                systNomName="sf",
                sel=sel
            )

            # Electron Trigger SF
            el_trigger_sf = get_correction(
                (EGamma_SF_JSONFiles[self.era][0]).replace(
                    "electron", "electronHlt"),
                "Electron-HLT-SF",
                params={"pt": lambda el: el.pt,
                        "eta": lambda el: el.eta,
                        "Path": "HLT_SF_Ele30_MVAiso90ID",
                        "year": EGamma_SF_JSONFiles[self.era][1]
                        },
                systParam="ValType",
                systNomName="sf",
                sel=sel
            )

            if sel.name == 'elPairMultiplicitySel':
                # pt cut here since correction's available only for pt >= 10
                sel = sel.refine(sel.name+"_electron_ID_SF", cut=[
                    op.AND(self.firstElTightPair[0].pt >= 10,
                           self.firstElTightPair[1].pt >= 10
                           )],
                    weight=[
                        el_trigger_sf(self.firstElTightPair[0]),
                        electron_ID_sf(self.firstElTightPair[0]),
                        electron_ID_sf(self.firstElTightPair[1])]
                )

            elif sel.name == 'emuPairMultiplicitySel':
                sel = sel.refine(sel.name+"_electron_ID_SF",
                                 cut=[self.firstEmuTightPair[0].pt >= 10],
                                 weight=[
                                     el_trigger_sf(self.firstElTightPair[0]),
                                     electron_ID_sf(self.firstEmuTightPair[0])]
                                 )

            else:
                sel = sel.refine(sel.name+"_electron_ID_SF",
                                 weight=op.c_float(1.))
        else:
            sel = sel.refine(sel.name+"_electron_ID_SF", weight=op.c_float(1.))

        return sel

    def NoiseFilters(self, FlagBranch, sel):
        "https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC"
        flags = [FlagBranch.goodVertices,
                 FlagBranch.globalSuperTightHalo2016Filter,
                 FlagBranch.EcalDeadCellTriggerPrimitiveFilter,
                 FlagBranch.BadPFMuonFilter,
                 FlagBranch.BadPFMuonDzFilter,
                 FlagBranch.hfNoisyHitsFilter,
                 FlagBranch.eeBadScFilter,
                 FlagBranch.ecalBadCalibFilter]
        sel = sel.refine('NoiseFilters', cut=flags)
        self.yields.add(sel, 'Noise filters')
        return sel
