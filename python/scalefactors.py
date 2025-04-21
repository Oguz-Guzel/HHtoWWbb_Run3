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

DY_and_Recoil_JSONFiles = {
    "2022": "DY_pTll_recoil_corrections_2022preEE_v2.json.gz",
    "2022EE": "DY_pTll_recoil_corrections_2022postEE_v2.json.gz",
    "2023": "DY_pTll_recoil_corrections_2023preBPix_v2.json.gz",
    "2023BPix": "DY_pTll_recoil_corrections_2023postBPix_v2.json.gz",
}


class ScaleFactors():
    """Class to define scale factors"""

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

    def btagSF(self, sel, jets, json_tagger="particleNet", jet_tagger="btagPNetB", btagReweightStudy=False):
        """Apply btagging SF"""
        if self.is_MC:
            from bamboo.scalefactors import get_bTagSF_itFit, makeBtagWeightItFit
            logger.info("Applying btagging SF for "+sel.name)
            def btvSF(flav): return get_bTagSF_itFit(
                BTV_SF_JSONFiles[self.era], json_tagger, jet_tagger, flav, sel=sel, decorr_eras=True, era=self.era)
            btvWeight = makeBtagWeightItFit(jets, btvSF)
            if not btagReweightStudy:
                btag_corr = get_correction(
                    f"{self.git_project_dir}/data/{self.era[:4]}_btagSF_reweight_v1_2_3.json.gz",
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

        sel = sel.refine(sel.name+"_btagSF", weight=btvWeight)
        self.yields.add(sel, sel.name)
        
        sel = sel.refine(sel.name+"_btagRW", weight=btag_reweight)
        self.yields.add(sel, sel.name)

        return sel

    def mumuSF(self, sel):
        """Apply lepton scalefactors for muon pair"""
        if self.is_MC:
            logger.info("Applying Muon SF for "+sel.name)
            # Muon ID SF
            systName = "syst"
            self.muon_ID_sf = get_correction(
                MUON_SF_JSONFiles[self.era],
                "NUM_MediumID_DEN_TrackerMuons",  # NUM_MediumPromptID_DEN_TrackerMuons, too ?
                systVariations={"muonIdSFup": f"{systName}up",
                                "muonIdSFdown": f"{systName}down"},
                params={"pt": lambda mu: mu.pt,
                        "eta": lambda mu: op.abs(mu.eta)},
                systParam="scale_factors",
                systNomName="nominal",
                systName=systName,
                defineOnFirstUse=False,
                sel=sel
            )

            # Muon ISO SF
            self.muon_ISO_sf = get_correction(
                MUON_SF_JSONFiles[self.era],
                # since muon iso is miniPFreliso and id is medium
                "NUM_TightPFIso_DEN_MediumID",
                systVariations={"muonIsoSFup": f"{systName}up",
                                "muonIsoSFdown": f"{systName}down"},
                params={"pt": lambda mu: mu.pt,
                        "eta": lambda mu: op.abs(mu.eta),
                        },
                systParam="scale_factors",
                systNomName="nominal",
                systName=systName,
                defineOnFirstUse=False,
                sel=sel
            )

            # Muon Trigger SF
            self.muon_TRG_sf = get_correction(
                MUON_SF_JSONFiles[self.era],
                "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight",
                systVariations={"muonTrgSFup": f"{systName}up",
                                "muonTrgSFdown": f"{systName}down"},
                params={
                    "pt": lambda mu: mu.pt,
                    "eta": lambda mu: op.abs(mu.eta)
                },
                systParam="scale_factors",
                systNomName="nominal",
                systName=systName,
                defineOnFirstUse=False,
                sel=sel
            )
            # if 'muPairMultiplicitySel' in sel.name:
            # pt and eta cut here since correction are available only when pt >= 15 and |eta| < 2.4
            sel = sel.refine('mumu_leading_ID_SF',
                             weight=[op.switch(
                                 op.AND(self.tightMuons[0].pt >= 15.,
                                        op.abs(self.tightMuons[0].eta) < 2.4,),
                                     self.muon_ID_sf(self.tightMuons[0]),
                                     op.c_float(1.))]
                             )
            sel = sel.refine('mumu_subleading_ID_SF',
                             weight=[op.switch(
                                 op.AND(self.tightMuons[1].pt >= 15.,
                                        op.abs(self.tightMuons[1].eta) < 2.4,),
                                     self.muon_ID_sf(self.tightMuons[1]),
                                     op.c_float(1.))]
                             )

            sel = sel.refine('mumu_leading_ISO_SF',
                             weight=[op.switch(
                                 op.AND(self.tightMuons[0].pt >= 15.,
                                        op.abs(self.tightMuons[0].eta) < 2.4,),
                                     self.muon_ISO_sf(self.tightMuons[0]),
                                     op.c_float(1.))]
                             )
            sel = sel.refine('mumu_subleading_ISO_SF',
                             weight=[op.switch(
                                 op.AND(self.tightMuons[1].pt >= 15.,
                                        op.abs(self.tightMuons[1].eta) < 2.4,),
                                     self.muon_ISO_sf(self.tightMuons[1]),
                                     op.c_float(1.))]
                             )
            sel = sel.refine('mumu_leading_TRG_SF',
                             weight=[op.switch(
                                 op.AND(self.tightMuons[0].pt >= 26.,
                                        op.abs(self.tightMuons[0].eta) < 2.4,),
                                     self.muon_TRG_sf(self.tightMuons[0]),
                                     op.c_float(1.))]
                             )
            sel = sel.refine('mumu_subleading_TRG_SF',
                             weight=[op.switch(
                                 op.AND(self.tightMuons[1].pt >= 26.,
                                        op.abs(self.tightMuons[1].eta) < 2.4,),
                                     self.muon_TRG_sf(self.tightMuons[1]),
                                     op.c_float(1.))]
                             )
        else:
            # followings are added to avoid cut-flow breaking because
            # the selection yields are not shown when it's not available
            # for any sample
            sel = sel.refine("mumu_leading_ID_SF",
                             weight=op.c_float(1.))
            sel = sel.refine("mumu_subleading_ID_SF",
                             weight=op.c_float(1.))
            sel = sel.refine("mumu_leading_ISO_SF",
                             weight=op.c_float(1.))
            sel = sel.refine("mumu_subleading_ISO_SF",
                             weight=op.c_float(1.))
            sel = sel.refine("mumu_leading_TRG_SF",
                             weight=op.c_float(1.))
            sel = sel.refine("mumu_subleading_TRG_SF",
                             weight=op.c_float(1.))

        self.yields.add(sel, sel.name)

        return sel

    def elelSF(self, sel):
        """Apply lepton scalefactors for electron pair"""
        if self.is_MC:
            logger.info("Applying Electron SF for "+sel.name)

            params = {"pt": lambda e: e.pt,
                      "eta": lambda e: e.eta,
                      "year": EGamma_SF_JSONFiles[self.era][1],
                      "WorkingPoint": "wp90iso"}

            # add phi for 2023 and 2023BPix
            if self.era in ['2023', '2023BPix']:
                params["phi"] = lambda e: e.phi

            systName = "sf"
            # Electron ID SF
            self.el_ID_sf = get_correction(
                EGamma_SF_JSONFiles[self.era][0],
                "Electron-ID-SF",
                systVariations={"elIdSFup": f"{systName}up",
                                "elIdSFdown": f"{systName}down"},
                params=params,
                systParam="ValType",
                systNomName=systName,
                defineOnFirstUse=False,
                sel=sel
            )

            # Electron Trigger SF
            self.el_TRG_sf = get_correction(
                (EGamma_SF_JSONFiles[self.era][0]).replace(
                    "electron", "electronHlt"),
                "Electron-HLT-SF",
                systVariations={"elTrgSFup": f"{systName}up",
                                "elTrgSFdown": f"{systName}down"},
                params={"pt": lambda el: el.pt,
                        "eta": lambda el: el.eta,
                        "Path": "HLT_SF_Ele30_MVAiso90ID",
                        "year": EGamma_SF_JSONFiles[self.era][1]
                        },
                systParam="ValType",
                systNomName=systName,
                defineOnFirstUse=False,
                sel=sel
            )
            # pt cut here since ID correction is available only for pt >= 10
            sel = sel.refine('elel_leading_ID_SF',
                             weight=[op.switch(
                                 self.tightElectrons[0].pt >= 10, self.el_ID_sf(
                                     self.tightElectrons[0]),
                                 op.c_float(1.))]
                             )
            sel = sel.refine('elel_subleading_ID_SF',
                             weight=[op.switch(
                                 self.tightElectrons[1].pt >= 10, self.el_ID_sf(
                                     self.tightElectrons[1]),
                                 op.c_float(1.))]
                             )
            sel = sel.refine('elel_leading_TRG_SF',
                             weight=[op.switch(
                                 self.tightElectrons[0].pt >= 25, self.el_TRG_sf(
                                     self.tightElectrons[0]),
                                 op.c_float(1.))]
                             )
            sel = sel.refine('elel_subleading_TRG_SF',
                             weight=[op.switch(
                                 self.tightElectrons[1].pt >= 25, self.el_TRG_sf(
                                     self.tightElectrons[1]),
                                 op.c_float(1.))]
                             )
        else:
            sel = sel.refine('elel_leading_ID_SF', weight=op.c_float(1.))
            sel = sel.refine('elel_subleading_ID_SF', weight=op.c_float(1.))
            sel = sel.refine('elel_leading_TRG_SF', weight=op.c_float(1.))
            sel = sel.refine('elel_subleading_TRG_SF', weight=op.c_float(1.))

        self.yields.add(sel, sel.name)

        return sel

    def elmuSF(self, sel):
        """Apply lepton scalefactors for electron-muon pair."""
        if self.is_MC:
            logger.info("Applying Electron SF for "+sel.name)
            sel = sel.refine('elmu_el_ID_SF',
                            weight=[op.switch(
                                self.tightElectrons[0].pt >= 10, self.el_ID_sf(
                                    self.tightElectrons[0]),
                                op.c_float(1.))]
                            )
            sel = sel.refine('elmu_mu_ID_SF',
                            weight=[op.switch(
                                op.AND(self.tightMuons[0].pt >= 15.,
                                        op.abs(self.tightMuons[0].eta) < 2.4),
                                self.muon_ID_sf(self.tightMuons[0]),
                                op.c_float(1.))]
                            )
            sel = sel.refine('elmu_mu_ISO_SF',
                            weight=[op.switch(
                                op.AND(self.tightMuons[0].pt >= 15.,
                                        op.abs(self.tightMuons[0].eta) < 2.4),
                                self.muon_ISO_sf(self.tightMuons[0]),
                                op.c_float(1.))]
                            )
            sel = sel.refine('elmu_el_TRG_SF',
                            weight=[op.switch(
                                self.tightElectrons[0].pt >= 25, self.el_TRG_sf(
                                    self.tightElectrons[0]),
                                op.c_float(1.))]
                            )
            sel = sel.refine('elmu_mu_TRG_SF',
                            weight=[op.switch(
                                op.AND(self.tightMuons[0].pt >= 26.,
                                        op.abs(self.tightMuons[0].eta) < 2.4),
                                self.muon_TRG_sf(self.tightMuons[0]),
                                op.c_float(1.))]
                            )
        else:
            sel = sel.refine('elmu_el_ID_SF', weight=op.c_float(1.))
            sel = sel.refine('elmu_mu_ID_SF', weight=op.c_float(1.))
            sel = sel.refine('elmu_mu_ISO_SF', weight=op.c_float(1.))
            sel = sel.refine('elmu_el_TRG_SF', weight=op.c_float(1.))
            sel = sel.refine('elmu_mu_TRG_SF', weight=op.c_float(1.))
        return sel

    def Z_pT_reweight(self, sel, sample, GenPartBranch, pdgId):
        """Apply DY Z pT reweighting for given lepton pair."""
        if self.is_MC and sample.startswith("DY"):
            from bamboo.scalefactors import get_correction
            logger.info("Applying DY Z pT reweighting for " + sel.name)

            DY_and_Recoil_path = self.git_project_dir+"/data/hleprare/DYandRecoilCorrlib/"

            # for N_unc - consult the json file

            N_unc = 10

            systVariations = {
                f"ZpT{i}up": f"up{i}" for i in range(1, N_unc+1)}
            systVariations.update(
                {f"ZpT{i}down": f"down{i}" for i in range(1, N_unc+1)})

            get_Z_pT_corr = get_correction(
                DY_and_Recoil_path + DY_and_Recoil_JSONFiles[self.era],
                "DY_pTll_reweighting",
                params={
                    "order": 'NLO',
                    "ptll": lambda leptons: op.rng_sum(leptons, lambda l: l.pt),
                },
                systNomName="nom",
                systVariations=systVariations,
                systParam="syst",
                sel=sel
            )

            def get_gen_parts(GenPartBranch, pdgId):
                gen_leptons = op.sort(op.select(
                    GenPartBranch, lambda p: op.AND(
                        p.status == 1,
                        op.abs(p.pdgId) == pdgId,
                        ((op.static_cast("int", p.statusFlags) << 8) & 1),
                    )), lambda p: -p.pt)
                return gen_leptons

            sel = sel.refine(
                sel.name+"_ZpT", weight=get_Z_pT_corr(get_gen_parts(GenPartBranch, pdgId)))
        else:
            sel = sel.refine(sel.name+"_ZpT", weight=op.c_float(1.))
        self.yields.add(sel, sel.name)
        return sel
