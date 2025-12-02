import os
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
    "2022EE": (
        jsonPathBase + "EGM/2022_Summer22EE/electron.json.gz",
        "2022Re-recoE+PromptFG",
    ),
    "2023": (jsonPathBase + "EGM/2023_Summer23/electron.json.gz", "2023PromptC"),
    "2023BPix": (
        jsonPathBase + "EGM/2023_Summer23BPix/electron.json.gz",
        "2023PromptD",
    ),
}

DY_and_Recoil_JSONFiles = {
    "2022": "DY_pTll_recoil_corrections_2022preEE_v2.json.gz",
    "2022EE": "DY_pTll_recoil_corrections_2022postEE_v2.json.gz",
    "2023": "DY_pTll_recoil_corrections_2023preBPix_v2.json.gz",
    "2023BPix": "DY_pTll_recoil_corrections_2023postBPix_v2.json.gz",
}

JetVeto_JSONFiles = {
    "2022": (jsonPathBase + "JME/2022_Summer22/jetvetomaps.json.gz", "Summer22_23Sep2023_RunCD_V1"),
    "2022EE": (jsonPathBase + "JME/2022_Summer22EE/jetvetomaps.json.gz", "Summer22EE_23Sep2023_RunEFG_V1",),
    "2023": (jsonPathBase + "JME/2023_Summer23/jetvetomaps.json.gz", "Summer23Prompt23_RunC_V1"),
    "2023BPix": (jsonPathBase + "JME/2023_Summer23BPix/jetvetomaps.json.gz", "Summer23BPixPrompt23_V3"),
}

sampleNumDict = {
    "WtoLNu-2Jets": 150,

    "WtoLNu-2Jets_PTLNu-40to100_1J": 151,
    "WtoLNu-2Jets_PTLNu-40to100_2J": 152,
    "WtoLNu-2Jets_PTLNu-100to200_1J": 153,
    "WtoLNu-2Jets_PTLNu-100to200_2J": 154,
    "WtoLNu-2Jets_PTLNu-200to400_1J": 155,
    "WtoLNu-2Jets_PTLNu-200to400_2J": 156,
    "WtoLNu-2Jets_PTLNu-400to600_1J": 157,
    "WtoLNu-2Jets_PTLNu-400to600_2J": 158,

    "WtoLNu-2Jets_PTLNu-600_1J": 159,
    "WtoLNu-2Jets_PTLNu-600_2J": 160,

    "WtoLNu-2Jets_0J": 161,
    "WtoLNu-2Jets_1J": 162,
    "WtoLNu-2Jets_2J": 163,

    "DYto2L-2Jets_MLL-50": 170,

    "DYto2L-2Jets_MLL-50_PTLL-40to100_1J": 171,
    "DYto2L-2Jets_MLL-50_PTLL-40to100_2J": 172,
    "DYto2L-2Jets_MLL-50_PTLL-100to200_1J": 173,
    "DYto2L-2Jets_MLL-50_PTLL-100to200_2J": 174,
    "DYto2L-2Jets_MLL-50_PTLL-200to400_1J": 175,
    "DYto2L-2Jets_MLL-50_PTLL-200to400_2J": 176,
    "DYto2L-2Jets_MLL-50_PTLL-400to600_1J": 177,
    "DYto2L-2Jets_MLL-50_PTLL-400to600_2J": 178,
    "DYto2L-2Jets_MLL-50_PTLL-600_1J": 179,
    "DYto2L-2Jets_MLL-50_PTLL-600_2J": 180,

    "DYto2L-2Jets_MLL-50_0J": 181,
    "DYto2L-2Jets_MLL-50_1J": 182,
    "DYto2L-2Jets_MLL-50_2J": 183,

    "Zto2Nu-2Jets_PTNuNu-40to100_1J": 191,
    "Zto2Nu-2Jets_PTNuNu-40to100_2J": 192,
    "Zto2Nu-2Jets_PTNuNu-100to200_1J": 193,
    "Zto2Nu-2Jets_PTNuNu-100to200_2J": 194,
    "Zto2Nu-2Jets_PTNuNu-200to400_1J": 195,
    "Zto2Nu-2Jets_PTNuNu-200to400_2J": 196,
    "Zto2Nu-2Jets_PTNuNu-400to600_1J": 197,
    "Zto2Nu-2Jets_PTNuNu-400to600_2J": 198,
    "Zto2Nu-2Jets_PTNuNu-600_1J": 199,
    "Zto2Nu-2Jets_PTNuNu-600_2J": 200,

}


class ScaleFactors:
    """Class to define scale factors"""

    def __init__(self, parent):
        self.parent = parent
        self.di_lepton_trigger_JSONFiles = {
            "2022": (
                os.path.join(
                    self.parent.git_project_dir,
                    "data",
                    "2022_di_lepton_trigger_scale_factors.json",
                ),
                "trigger_scale_factors_2d",
            ),
            "2023": (
                os.path.join(
                    self.parent.git_project_dir,
                    "data",
                    "2023_di_lepton_trigger_scale_factors.json",
                ),
                "trigger_scale_factors_2d",
            ),
        }

    def jet_veto_map(self, tree, sel):
        # https://cms-jerc.web.cern.ch/Recommendations/#run-3
        logger.info("Applying JetVeto")
        corr = get_correction(
            JetVeto_JSONFiles[self.parent.era][0],
            JetVeto_JSONFiles[self.parent.era][1],
            params={
                "type": "jetvetomap",
                "eta": lambda j: j.eta,
                "phi": lambda j: j.phi,
            },
            sel=sel
        )
    
        jets_to_veto = op.select(
            tree.Jet, lambda j:
            op.AND(j.pt > 15,
                   (j.jetId & 8) != 0,
                   (j.chEmEF + j.neEmEF) < 0.9)
        )
    
        veto_cuts = op.rng_any(
            jets_to_veto, lambda j: corr(j) != 0
        )
    
        sel = sel.refine("JetVetoMaps", cut=veto_cuts)
    
        self.parent.yields.add(sel, "JetVetoMaps")
    
        return sel

    def NoiseFilters(self, FlagBranch, sel):
        "https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC"
        flags = [
            FlagBranch.goodVertices,
            FlagBranch.globalSuperTightHalo2016Filter,
            FlagBranch.EcalDeadCellTriggerPrimitiveFilter,
            FlagBranch.BadPFMuonFilter,
            FlagBranch.BadPFMuonDzFilter,
            FlagBranch.hfNoisyHitsFilter,
            FlagBranch.eeBadScFilter,
            FlagBranch.ecalBadCalibFilter,
        ]
        sel = sel.refine("NoiseFilters", cut=flags)
        self.parent.yields.add(sel, "Noise filters")
        return sel

    def top_pT_reweight(self, GenPartBranch, sel, sample):
        """Apply top p_T reweighting. Check for more
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the"""
        if sample.startswith("TT"):

            def top_pt_weight(pt):
                return op.exp(
                    -2.02274e-01
                    + 1.09734e-04 * pt
                    + -1.30088e-07 * op.pow(pt, 2)
                    + (5.83494e01 / (pt + 1.96252e02))
                )

            def getTopPtWeight(GenPart):
                lastCopy = op.select(
                    GenPart, lambda p: (op.static_cast(
                        "int", p.statusFlags) >> 13) & 1
                )
                tops = op.select(lastCopy, lambda p: p.pdgId == 6)
                antitops = op.select(lastCopy, lambda p: p.pdgId == -6)
                weight = op.switch(
                    op.AND(op.rng_len(tops) >= 1, op.rng_len(antitops) >= 1),
                    op.sqrt(top_pt_weight(tops[0].pt)
                            * top_pt_weight(antitops[0].pt)),
                    1.0,
                )
                return weight

            logger.info("Applying Top Pt reweighting (only to TTbar samples)")

            w = getTopPtWeight(GenPartBranch)

            # nominal, up, down definitions
            # Systematics are symmetric in log-weight space; also makes the full effect the 1\sigma.
            w_nom = op.c_float(1.0)
            w_up = w
            w_down = 1.0 / w

            sel = sel.refine(
                "topPt",
                weight=op.systematic(
                    w_nom,
                    "topPtRW",
                    up=w_up,
                    down=w_down,
                ),
            )

        else:
            sel = sel.refine("topPt", weight=op.c_float(1.0))
        self.parent.yields.add(sel, "topPt reweighting")

        return sel

    def btagSF(
        self,
        sel,
        jets,
        json_tagger="particleNet",
        jet_tagger="btagPNetB",
        btagReweightStudy=False,
    ):
        """Apply btagging SF"""
        if self.parent.is_MC:
            from bamboo.scalefactors import get_bTagSF_itFit, makeBtagWeightItFit

            logger.info("Applying btagging SF for " + sel.name)

            def btvSF(flav):
                return get_bTagSF_itFit(
                    BTV_SF_JSONFiles[self.parent.era],
                    json_tagger,
                    jet_tagger,
                    flav,
                    sel=sel,
                    decorr_eras=True,
                    era=self.parent.era,
                )

            btvWeight = makeBtagWeightItFit(jets, btvSF)
            if not btagReweightStudy:
                btag_corr = get_correction(
                    f"{self.parent.git_project_dir}/data/{self.parent.era[:4]}_btagSF_reweight.json.gz",
                    "Ratio_btagSF_shape",
                    params={
                        "year": self.parent.era,
                        # 1. to make it a float
                        "jet_multiplicity": 1.0 * op.rng_len(jets),
                    },
                    sel=sel,
                )
                # None since the object is already in the btag_corr i.e. self.ak4Jets
                btag_reweight = btag_corr(None)
            else:
                btag_reweight = op.c_float(1.0)
        else:
            btvWeight = op.c_float(1.0)
            btag_reweight = op.c_float(1.0)

        sel = sel.refine(sel.name + "_btagSF", weight=btvWeight)
        self.parent.yields.add(sel, sel.name)

        sel = sel.refine(sel.name + "_btagRW", weight=btag_reweight)
        self.parent.yields.add(sel, sel.name)

        return sel

    def mumuSF(self, sel):
        """Apply lepton scalefactors for muon pair"""
        if self.parent.is_MC:
            logger.info("Applying Muon SF for " + sel.name)
            # Muon ID SF
            systName = "syst"
            self.muon_ID_sf = get_correction(
                MUON_SF_JSONFiles[self.parent.era],
                "NUM_MediumID_DEN_TrackerMuons",  # NUM_MediumPromptID_DEN_TrackerMuons, too ?
                systVariations={
                    "muonIdSFup": f"{systName}up",
                    "muonIdSFdown": f"{systName}down",
                },
                params={"pt": lambda mu: mu.pt,
                        "eta": lambda mu: op.abs(mu.eta)},
                systParam="scale_factors",
                systNomName="nominal",
                systName=systName,
                defineOnFirstUse=False,
                sel=sel,
            )

            # Muon ISO SF
            self.muon_ISO_sf = get_correction(
                MUON_SF_JSONFiles[self.parent.era],
                # since muon iso is miniPFreliso and id is medium
                "NUM_TightPFIso_DEN_MediumID",
                systVariations={
                    "muonIsoSFup": f"{systName}up",
                    "muonIsoSFdown": f"{systName}down",
                },
                params={
                    "pt": lambda mu: mu.pt,
                    "eta": lambda mu: op.abs(mu.eta),
                },
                systParam="scale_factors",
                systNomName="nominal",
                systName=systName,
                defineOnFirstUse=False,
                sel=sel,
            )

            # single muon Trigger SF
            self.muon_single_TRG_SF = get_correction(
                MUON_SF_JSONFiles[self.parent.era],
                "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                systVariations={
                    "muonTrgSFup": f"{systName}up",
                    "muonTrgSFdown": f"{systName}down",
                },
                params={"pt": lambda mu: mu.pt,
                        "eta": lambda mu: op.abs(mu.eta)},
                systParam="scale_factors",
                systNomName="nominal",
                systName=systName,
                defineOnFirstUse=False,
                sel=sel,
            )
            # pt and eta cut here since correction are available only when pt >= 15 and |eta| < 2.4
            sel = sel.refine(
                "mumu_leading_ID_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[0].pt >= 15.0,
                            op.abs(self.parent.tightMuons[0].eta) < 2.4,
                        ),
                        self.muon_ID_sf(self.parent.tightMuons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "MuMu leading ID SF")
            sel = sel.refine(
                "mumu_subleading_ID_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[1].pt >= 15.0,
                            op.abs(self.parent.tightMuons[1].eta) < 2.4,
                        ),
                        self.muon_ID_sf(self.parent.tightMuons[1]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "MuMu sub-leading ID SF")
            sel = sel.refine(
                "mumu_leading_ISO_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[0].pt >= 15.0,
                            op.abs(self.parent.tightMuons[0].eta) < 2.4,
                        ),
                        self.muon_ISO_sf(self.parent.tightMuons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "MuMu leading ISO SF")
            sel = sel.refine(
                "mumu_subleading_ISO_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[1].pt >= 15.0,
                            op.abs(self.parent.tightMuons[1].eta) < 2.4,
                        ),
                        self.muon_ISO_sf(self.parent.tightMuons[1]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "MuMu sub-leading ISO SF")
            sel = sel.refine(
                "mumu_leading_TRG_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[0].pt >= 26.0,
                            op.abs(self.parent.tightMuons[0].eta) < 2.4,
                        ),
                        self.muon_single_TRG_SF(self.parent.tightMuons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "MuMu single TRG SF")
        else:
            # followings are added to avoid cut-flow breaking because
            # the selection yields are not shown when it's not available
            sel = sel.refine("mumu_leading_ID_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "MuMu leading ID SF")
            sel = sel.refine("mumu_subleading_ID_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "MuMu sub-leading ID SF")
            sel = sel.refine("mumu_leading_ISO_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "MuMu leading ISO SF")
            sel = sel.refine("mumu_subleading_ISO_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "MuMu sub-leading ISO SF")
            sel = sel.refine("mumu_leading_TRG_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "MuMu single TRG SF")

        return sel

    def elelSF(self, sel):
        """Apply lepton scalefactors for electron pair"""
        if self.parent.is_MC:
            logger.info("Applying Electron SF for " + sel.name)

            params = {
                "pt": lambda e: e.pt,
                "eta": lambda e: e.eta,
                "year": EGamma_SF_JSONFiles[self.parent.era][1],
                "WorkingPoint": "wp90iso",
            }

            # add phi for 2023 and 2023BPix
            if self.parent.era in ["2023", "2023BPix"]:
                params["phi"] = lambda e: e.phi

            systNomName = "sf"
            # Electron ID SF
            self.el_ID_sf = get_correction(
                EGamma_SF_JSONFiles[self.parent.era][0],
                "Electron-ID-SF",
                systVariations={
                    "elIdSFup": f"{systNomName}up",
                    "elIdSFdown": f"{systNomName}down",
                },
                params=params,
                systParam="ValType",
                systNomName=systNomName,
                defineOnFirstUse=False,
                sel=sel,
            )

            # single electron trigger SF
            self.elel_single_TRG_SF = get_correction(
                (EGamma_SF_JSONFiles[self.parent.era][0]).replace(
                    "electron", "electronHlt"
                ),
                "Electron-HLT-SF",
                systVariations={
                    "elTrgSFup": f"{systNomName}up",
                    "elTrgSFdown": f"{systNomName}down",
                },
                params={
                    "pt": lambda e: e.pt,
                    "eta": lambda e: e.eta,
                    "Path": "HLT_SF_Ele30_MVAiso90ID",
                    "year": EGamma_SF_JSONFiles[self.parent.era][1],
                },
                systParam="ValType",
                systNomName=systNomName,
                defineOnFirstUse=False,
                sel=sel,
            )
            # pt cut here since corrections are available for certain ranges
            sel = sel.refine(
                "elel_leading_ID_SF",
                weight=[
                    op.switch(
                        self.parent.tightElectrons[0].pt >= 10,
                        self.el_ID_sf(self.parent.tightElectrons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElEl leading ID SF")
            sel = sel.refine(
                "elel_subleading_ID_SF",
                weight=[
                    op.switch(
                        self.parent.tightElectrons[1].pt >= 10,
                        self.el_ID_sf(self.parent.tightElectrons[1]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElEl subleading ID SF")
            sel = sel.refine(
                "elel_leading_TRG_SF",
                weight=[
                    op.switch(
                        self.parent.tightElectrons[0].pt >= 25,
                        self.elel_single_TRG_SF(self.parent.tightElectrons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElEl single TRG SF")
        else:
            sel = sel.refine("elel_leading_ID_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElEl leading ID SF")
            sel = sel.refine("elel_subleading_ID_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElEl subleading ID SF")
            sel = sel.refine("elel_leading_TRG_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElEl single TRG SF")

        return sel

    def elmuSF(self, sel):
        """Apply lepton scalefactors for electron-muon pair."""
        if self.parent.is_MC:
            logger.info("Applying Electron SF for " + sel.name)

            sel = sel.refine(
                "elmu_el_ID_SF",
                weight=[
                    op.switch(
                        self.parent.tightElectrons[0].pt >= 10,
                        self.el_ID_sf(self.parent.tightElectrons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElMu el ID SF")
            sel = sel.refine(
                "elmu_el_TRG_SF",
                weight=[
                    op.switch(
                        self.parent.tightElectrons[0].pt >= 25,
                        self.elel_single_TRG_SF(self.parent.tightElectrons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElMu el single TRG SF")
            sel = sel.refine(
                "elmu_mu_ID_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[0].pt >= 15.0,
                            op.abs(self.parent.tightMuons[0].eta) < 2.4,
                        ),
                        self.muon_ID_sf(self.parent.tightMuons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElMu mu ID SF")
            sel = sel.refine(
                "elmu_mu_ISO_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[0].pt >= 15.0,
                            op.abs(self.parent.tightMuons[0].eta) < 2.4,
                        ),
                        self.muon_ISO_sf(self.parent.tightMuons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElMu mu ISO SF")
            sel = sel.refine(
                "elmu_mu_TRG_SF",
                weight=[
                    op.switch(
                        op.AND(
                            self.parent.tightMuons[0].pt >= 26.0,
                            op.abs(self.parent.tightMuons[0].eta) < 2.4,
                        ),
                        self.muon_single_TRG_SF(self.parent.tightMuons[0]),
                        op.c_float(1.0),
                    )
                ],
            )
            self.parent.yields.add(sel, "ElMu mu single TRG SF")

        else:
            sel = sel.refine("elmu_el_ID_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElMu el ID SF")
            sel = sel.refine("elmu_el_TRG_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElMu el TRG SF")
            sel = sel.refine("elmu_mu_ID_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElMu mu ID SF")
            sel = sel.refine("elmu_mu_ISO_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElMu mu ISO SF")
            sel = sel.refine("elmu_mu_TRG_SF", weight=op.c_float(1.0))
            self.parent.yields.add(sel, "ElMu mu TRG SF")

        return sel

    def Z_pT_reweight(self, sel, sample, GenPartBranch, pdgId):
        """Apply DY ptll and recoil corrections  for given lepton pair."""
        if self.parent.is_MC and sample.startswith("DY"):
            from bamboo.scalefactors import get_correction

            logger.info(
                "Applying DY ptll and recoil corrections for " + sel.name)

            DY_and_Recoil_path = (
                self.parent.git_project_dir + "/data/hleprare/DYandRecoilCorrlib/"
            )

            # for N_unc - consult the json file

            N_unc = 1

            systVariations = {f"ZpTup": f"up{i}" for i in range(1, N_unc + 1)}
            systVariations.update(
                {f"ZpTdown": f"down{i}" for i in range(1, N_unc + 1)}
            )

            get_Z_pT_corr = get_correction(
                DY_and_Recoil_path + DY_and_Recoil_JSONFiles[self.parent.era],
                "DY_pTll_reweighting",
                params={
                    "order": "NLO",
                    "ptll": lambda leptons: op.rng_sum(leptons, lambda l: l.pt),
                },
                systNomName="nom",
                systVariations=systVariations,
                systParam="syst",
                sel=sel,
            )

            def get_gen_parts(GenPartBranch, pdgId):
                gen_leptons = op.sort(
                    op.select(
                        GenPartBranch,
                        lambda p: op.AND(
                            p.status == 1,
                            op.abs(p.pdgId) == pdgId,
                            ((op.static_cast("int", p.statusFlags) << 8) & 1),
                        ),
                    ),
                    lambda p: -p.pt,
                )
                return gen_leptons

            sel = sel.refine(
                sel.name + "_ZpT",
                weight=get_Z_pT_corr(get_gen_parts(GenPartBranch, pdgId)),
            )
        else:
            sel = sel.refine(sel.name + "_ZpT", weight=op.c_float(1.0))
        self.parent.yields.add(sel, sel.name)
        return sel

    def dilepton_trg_sf(self, sel):
        if "mumu" in sel.name:
            leading_lepton_pt = self.parent.tightMuons[0].pt
            subleading_lepton_pt = self.parent.tightMuons[1].pt
        elif "ee" in sel.name:
            leading_lepton_pt = self.parent.tightElectrons[0].pt
            subleading_lepton_pt = self.parent.tightElectrons[1].pt
        elif "emu" in sel.name:
            leading_lepton_pt = op.switch(
                self.parent.tightElectrons[0].pt > self.parent.tightMuons[0].pt,
                self.parent.tightElectrons[0].pt,
                self.parent.tightMuons[0].pt,
            )
            subleading_lepton_pt = op.switch(
                self.parent.tightElectrons[0].pt > self.parent.tightMuons[0].pt,
                self.parent.tightMuons[0].pt,
                self.parent.tightElectrons[0].pt,
            )
        else:
            raise RuntimeError(
                "Selection name must include one of these values: ee, mumu, emu."
            )
        di_lepton_TRG_SF = get_correction(
            self.di_lepton_trigger_JSONFiles[self.parent.era[:4]][0],
            self.di_lepton_trigger_JSONFiles[self.parent.era[:4]][1],
            systParam="systematic",
            systVariations={
                "diTRG_up": "up",
                "diTRG_down": "down",
            },
            systNomName="nominal",
            params={
                "channel": sel.name,
                "pt_leading": leading_lepton_pt,
                "pt_subleading": subleading_lepton_pt,
            },
            defineOnFirstUse=False,
            sel=sel,
        )

        sel = sel.refine(
            sel.name+"_di_lepton_TRG_SF", weight=di_lepton_TRG_SF(None)
        )
        self.parent.yields.add(sel, "di-lepton TRG SF")

        return sel

    def V_Jets_Stitching(self, LHEBranch, sel, sample):
        """Apply V+Jets stitching for MC samples DY MLL > 50. pT binned, Jet multiplicity binned and inclusive."""
        if LHEBranch is not None:
            logger.info("Applying V+Jets stitching for " + sel.name)
            stitch_map_json = f"{self.parent.git_project_dir}/data/Run3NLOStitching.json"

            # Get the base sample name and corresponding sample number
            base_sample = sample.rsplit(f"_{self.parent.era}", 1)[0]
            sampleNum = sampleNumDict.get(base_sample, None)

            _weight = op.c_float(1.)
            if sampleNum is None:
                logger.warning(
                    f"Warning! Sample {base_sample} not found in sampleNumDict for stitching")
            else:
                # Use get_correction for the stitching weights
                era_key_map = {
                    "2022": "NLO_stitch_22",
                    "2022EE": "NLO_stitch_22EE",
                    "2023": "NLO_stitch_23",
                    "2023BPix": "NLO_stitch_23BPix"
                }

                era_key = era_key_map.get(self.parent.era)
                if era_key is None:
                    logger.warning(
                        f"Warning! Era {self.parent.era} is not valid for stitching")
                else:
                    LHE_NpNLO = LHEBranch.NpNLO
                    LHE_VpT = LHEBranch.Vpt

                    # Define VpT bin edges and corresponding bin indices using op.multiSwitch
                    def get_vpt_bin(vpt):
                        return op.multiSwitch(
                            (vpt <= 0, 0.0),
                            (vpt <= 40, 1.0),
                            (vpt <= 100, 2.0),
                            (vpt <= 200, 3.0),
                            (vpt <= 400, 4.0),
                            (vpt <= 600, 5.0),
                            6.0
                        )

                    vpt_bin = get_vpt_bin(LHE_VpT)
                    binVal = vpt_bin + op.c_float(LHE_NpNLO) * 7.0

                    # Use get_correction with the stitching JSON
                    stitching_corr = get_correction(
                        stitch_map_json,
                        era_key,
                        params={
                            "axis0": float(sampleNum),
                            "axis1": lambda obj: binVal
                        },
                        sel=sel
                    )

                    _weight = stitching_corr(None)

            sel = sel.refine(sel.name+"_VJetsStitching", weight=_weight)
        else:
            sel = sel.refine(sel.name+"_VJetsStitching", weight=op.c_float(1.))
        return sel
