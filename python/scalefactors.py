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
    "2023BPix": (jsonPathBase + "JME/2023_Summer23BPix/jetvetomaps.json.gz", "Summer23BPixPrompt23_RunD_V1"),
}

AK8_BBCC_SF_FILES = {
    "2022": (
        "ak8_sf_msdtest_Pt-combined_2022_preEE.json",
        "HHbbww_2022_preEE_SF_bb",
        "HHbbww_2022_preEE_SF_cc",
    ),
    "2022EE": (
        "ak8_sf_msdtest_Pt-combined_2022_postEE.json",
        "HHbbww_2022_postEE_SF_bb",
        "HHbbww_2022_postEE_SF_cc",
    ),
    "2023": (
        "ak8_sf_msdtest_Pt-combined_2023_preBPix.json",
        "HHbbww_2023_preBPix_SF_bb",
        "HHbbww_2023_preBPix_SF_cc",
    ),
    "2023BPix": (
        "ak8_sf_msdtest_Pt-combined_2023_postBPix.json",
        "HHbbww_2023_postBPix_SF_bb",
        "HHbbww_2023_postBPix_SF_cc",
    ),
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
        self.di_lepton_TRG_JSONFiles = {
            "ee":
            os.path.join(
                self.parent.git_project_dir,
                "data",
                "sf_ee_trg_lepton0_pt-trg_lepton1_pt-trig_idsV4_syst.json",
            ),
            "mixed":
            os.path.join(
                self.parent.git_project_dir,
                "data",
                "sf_mixed_trg_lepton0_pt-trg_lepton1_pt-trig_idsV4_syst.json",
            ),
            "mm":
            os.path.join(
                self.parent.git_project_dir,
                "data",
                "sf_mm_trg_lepton0_pt-trg_lepton1_pt-trig_idsV4_syst.json",
            ),
        }

    def _muRF_syst_name(self, sample, sampleCfg):
        """Return the muRF systematic name for a given sample.

        Uncertainties are treated as uncorrelated between processes.
        """
        if not sampleCfg or not sample:
            return None

        smp_lower = sample.lower()
        group = sampleCfg.get("group")

        def has_any(*tokens):
            return any(tok in smp_lower for tok in tokens)

        if has_any("ttw", "ttz", "ttv"):
            return "muRF_ttV"
        if group == "TT" or has_any("ttto", "tt_", "ttbar"):
            return "muRF_ttbar"
        if group == "DY" or has_any("dyto", "dy_", "dyjets"):
            return "muRF_DY"
        if group == "ST" or has_any("tw", "tbarw", "tchannel", "schannel"):
            return "muRF_ST"
        if group == "VV" or has_any("ww", "wz", "zz"):
            return "muRF_VV"

        # Single-H backgrounds
        if has_any("ggh", "glugluhto", "vbfh", "wh", "zh", "tth", "thq", "thw"):
            return "muRF_singleH"

        # Signal processes
        if has_any("bbww"):
            return "muRF_bbWW"
        if has_any("bbzz"):
            return "muRF_bbZZ"
        if has_any("bbtautau", "bbtt"):
            return "muRF_bbtautau"

        return None

    def _muRF_norm_factors(self, sample, sampleCfg, envelope_idx=(0, 1, 2, 3, 5, 6, 7, 8)):
        """Return normalisation factors to keep muR/muF variations shape-only.

        Priority:
        1) Explicit numbers in sampleCfg["muRF_norm"] (dict or iterable).
        2) Inclusive LHEScale sum of weights (sampleCfg["LHEScaleSumw"], length >= 9).
           The envelope is built from the sums and inverted to re-normalise.
        Falls back to (1.0, 1.0) if nothing is available.
        """
        norm_up = 1.0
        norm_down = 1.0

        if sampleCfg is None or not sample:
            logger.warning(
                "muRF_norm_factors: sampleCfg missing; falling back to no normalisation (may leak yield)."
            )
            return norm_up, norm_down

        # Load per-sample sums for muRF scale variations from YAML when not provided inline.
        if "LHEScaleSumw" not in sampleCfg:
            import yaml

            sumw_path = os.path.join(
                self.parent.git_project_dir, "data", f"LHEScaleSumw_{self.parent.era[:4]}.yaml"
            )
            try:
                with open(sumw_path) as handle:
                    sumw_map = yaml.safe_load(handle).get("LHEScaleSumw", {})
                if sample in sumw_map:
                    sampleCfg["LHEScaleSumw"] = sumw_map[sample]
            except FileNotFoundError:
                logger.warning("muRF_norm_factors: missing %s; skipping LHEScaleSumw injection", sumw_path)
            except Exception as exc:
                logger.warning("muRF_norm_factors: failed to load %s (%s)", sumw_path, exc)

        norm_cfg = sampleCfg.get("muRF_norm")
        if isinstance(norm_cfg, dict):
            return norm_cfg.get("up", 1.0), norm_cfg.get("down", 1.0)
        if isinstance(norm_cfg, (list, tuple)) and len(norm_cfg) >= 2:
            return norm_cfg[0], norm_cfg[1]

        sumw = sampleCfg.get("LHEScaleSumw")
        if isinstance(sumw, (list, tuple)) and len(sumw) >= 9:
            nom = sumw[4]
            if nom != 0:
                ratios = [sumw[i] / nom for i in envelope_idx]
                max_ratio = max(ratios)
                min_ratio = min(ratios)
                norm_up = 1.0 / max_ratio if max_ratio != 0 else 1.0
                norm_down = 1.0 / min_ratio if min_ratio != 0 else 1.0
            else:
                logger.warning("LHEScaleSumw nominal is zero; skipping muRF normalisation")

        # Warn if normalisation is missing; keep running to avoid job failure.
        if norm_up == 1.0 and norm_down == 1.0:
            logger.warning(
                "muRF shape-only: no muRF_norm or LHEScaleSumw found; variations will include normalisation effects."
            )

        return norm_up, norm_down

    def _pdf_shape_syst_name(self, sample, sampleCfg):
        """Return the PDF shape systematic name for a given sample.

        Uncertainties are treated as uncorrelated between processes.
        """
        if not sampleCfg or not sample:
            return None

        smp_lower = sample.lower()
        group = sampleCfg.get("group")

        def has_any(*tokens):
            return any(tok in smp_lower for tok in tokens)

        if has_any("ttw", "ttz", "ttv"):
            return "pdfShape_ttV"
        if group == "TT" or has_any("ttto", "tt_", "ttbar"):
            return "pdfShape_ttbar"
        if group == "DY" or has_any("dyto", "dy_", "dyjets"):
            return "pdfShape_DY"
        if group == "ST" or has_any("tw", "tbarw", "tchannel", "schannel"):
            return "pdfShape_ST"
        if group == "VV" or has_any("ww", "wz", "zz"):
            return "pdfShape_VV"

        # Single-H backgrounds
        if has_any("ggh", "glugluhto", "vbfh", "wh", "zh", "tth", "thq", "thw"):
            return "pdfShape_singleH"

        # Signal processes
        if has_any("bbww"):
            return "pdfShape_bbWW"
        if has_any("bbzz"):
            return "pdfShape_bbZZ"
        if has_any("bbtautau", "bbtt"):
            return "pdfShape_bbtautau"

        return None

    def _ps_fsr_syst_name(self, sample, sampleCfg):
        """Return the PS FSR systematic name for a given sample.

        Uncertainties are treated as uncorrelated between processes.
        """
        if not sampleCfg or not sample:
            return None

        smp_lower = sample.lower()
        group = sampleCfg.get("group")

        def has_any(*tokens):
            return any(tok in smp_lower for tok in tokens)

        if group == "TT":
            return "psFSR_ttbar"
        if group == "DY":
            return "psFSR_DY"
        if group == "ST":
            return "psFSR_ST"
        if group == "VV":
            return "psFSR_VV"
        if has_any("ttll", "ttz", "ttv"):
            return "psFSR_ttV"
        if has_any("glugluhto", "vbfhto", "wplush", "zh", "tth", "thq", "thw"):
            return "psFSR_singleH"

        return None

    def muRF_scale_weights(self, tree, sel, sample, sampleCfg):
        """Apply muR/muF scale uncertainty from LHEScaleWeight with envelope.

        The variation is kept shape-only by normalising to the inclusive
        phase-space yield. Normalisation factors are taken from
        sampleCfg["muRF_norm"] if provided, otherwise derived from
        sampleCfg["LHEScaleSumw"] (required, len >= 9) to rescale the envelope.
        Missing inputs raise to avoid leaking normalisation effects into the
        shape nuisance.
        """
        if not self.parent.is_MC:
            return sel

        if tree is None or not hasattr(tree, "LHEScaleWeight"):
            return sel

        syst_name = self._muRF_syst_name(sample, sampleCfg)
        if syst_name is None:
            return sel

        # Retrieve normalization inputs; warns (not raises) if missing to keep jobs running.
        norm_up, norm_down = self._muRF_norm_factors(sample, sampleCfg)

        weights = tree.LHEScaleWeight

        def _get_weight(idx):
            return op.switch(
                op.rng_len(weights) > idx,
                weights[idx],
                op.c_float(1.0),
            )

        nominal = _get_weight(4)

        def _safe_div(num, den):
            return op.switch(den != 0, num / den, op.c_float(1.0))

        ratios = [
            _safe_div(_get_weight(i), nominal) for i in [0, 1, 2, 3, 5, 6, 7, 8]
        ]

        def _max_list(vals):
            m = vals[0]
            for v in vals[1:]:
                m = op.max(m, v)
            return m

        def _min_list(vals):
            m = vals[0]
            for v in vals[1:]:
                m = op.min(m, v)
            return m

        up_ratio = _max_list(ratios)
        down_ratio = _min_list(ratios)

        sel = sel.refine(
            f"{syst_name}",
            weight=op.systematic(
                op.c_float(1.0),
                syst_name,
                up=up_ratio * op.c_float(norm_up),
                down=down_ratio * op.c_float(norm_down),
            ),
        )

        return sel

    def pdf_shape_weights(self, tree, sel, sample, sampleCfg):
        """Apply PDF shape uncertainty from LHEPdfWeight with envelope.

        The variation is shape-only; optional normalization factors can be
        provided in sampleCfg["pdf_norm"] as {up: <val>, down: <val>} or
        [up, down].
        """
        if not self.parent.is_MC:
            return sel

        if tree is None or not hasattr(tree, "LHEPdfWeight"):
            return sel

        syst_name = self._pdf_shape_syst_name(sample, sampleCfg)
        if syst_name is None:
            return sel

        weights = tree.LHEPdfWeight

        def _get_weight(idx):
            return op.switch(
                op.rng_len(weights) > idx,
                weights[idx],
                op.c_float(1.0),
            )

        # If weights are already normalized to nominal, the nominal is 1.
        # Otherwise, assume the first weight is the nominal reference.
        nominal = _get_weight(0)

        def _safe_div(num, den):
            return op.switch(den != 0, num / den, op.c_float(1.0))

        # Use a fixed maximum number of PDF replicas (default 100) with guards
        n_pdf_replicas = 100
        ratios = [
            _safe_div(_get_weight(i), nominal) for i in range(1, n_pdf_replicas + 1)
        ]

        up_ratio = ratios[0]
        down_ratio = ratios[0]
        for r in ratios[1:]:
            up_ratio = op.max(up_ratio, r)
            down_ratio = op.min(down_ratio, r)

        norm_up = 1.0
        norm_down = 1.0
        if sampleCfg is not None:
            norm_cfg = sampleCfg.get("pdf_norm")
            if isinstance(norm_cfg, dict):
                norm_up = norm_cfg.get("up", 1.0)
                norm_down = norm_cfg.get("down", 1.0)
            elif isinstance(norm_cfg, (list, tuple)) and len(norm_cfg) >= 2:
                norm_up, norm_down = norm_cfg[0], norm_cfg[1]

        sel = sel.refine(
            f"{syst_name}",
            weight=op.systematic(
                op.c_float(1.0),
                syst_name,
                up=up_ratio * op.c_float(norm_up),
                down=down_ratio * op.c_float(norm_down),
            ),
        )

        return sel

    def ps_isr_fsr_weights(self, tree, sel, sample, sampleCfg):
        """Apply PS ISR/FSR scale uncertainties from PSWeight.

        ISR is treated as correlated among processes (single nuisance).
        FSR is treated as uncorrelated between processes.
        Optional normalization factors can be provided in sampleCfg:
          psISR_norm / psFSR_norm as {up: <val>, down: <val>} or [up, down].
        """
        if not self.parent.is_MC:
            return sel

        if tree is None or not hasattr(tree, "PSWeight"):
            return sel

        weights = tree.PSWeight

        def _get_weight(idx):
            return op.switch(
                op.rng_len(weights) > idx,
                weights[idx],
                op.c_float(1.0),
            )

        # NanoAOD default scheme: 0/1 = ISR down/up, 2/3 = FSR down/up
        isr_down = _get_weight(0)
        isr_up = _get_weight(1)
        fsr_down = _get_weight(2)
        fsr_up = _get_weight(3)

        isr_norm_up = 1.0
        isr_norm_down = 1.0
        fsr_norm_up = 1.0
        fsr_norm_down = 1.0
        if sampleCfg is not None:
            isr_norm_cfg = sampleCfg.get("psISR_norm")
            if isinstance(isr_norm_cfg, dict):
                isr_norm_up = isr_norm_cfg.get("up", 1.0)
                isr_norm_down = isr_norm_cfg.get("down", 1.0)
            elif isinstance(isr_norm_cfg, (list, tuple)) and len(isr_norm_cfg) >= 2:
                isr_norm_up, isr_norm_down = isr_norm_cfg[0], isr_norm_cfg[1]

            fsr_norm_cfg = sampleCfg.get("psFSR_norm")
            if isinstance(fsr_norm_cfg, dict):
                fsr_norm_up = fsr_norm_cfg.get("up", 1.0)
                fsr_norm_down = fsr_norm_cfg.get("down", 1.0)
            elif isinstance(fsr_norm_cfg, (list, tuple)) and len(fsr_norm_cfg) >= 2:
                fsr_norm_up, fsr_norm_down = fsr_norm_cfg[0], fsr_norm_cfg[1]

        # ISR: correlated among processes
        sel = sel.refine(
            "psISR",
            weight=op.systematic(
                op.c_float(1.0),
                "psISR",
                up=isr_up * op.c_float(isr_norm_up),
                down=isr_down * op.c_float(isr_norm_down),
            ),
        )

        # FSR: uncorrelated between processes
        fsr_name = self._ps_fsr_syst_name(sample, sampleCfg)
        if fsr_name is not None:
            sel = sel.refine(
                fsr_name,
                weight=op.systematic(
                    op.c_float(1.0),
                    fsr_name,
                    up=fsr_up * op.c_float(fsr_norm_up),
                    down=fsr_down * op.c_float(fsr_norm_down),
                ),
            )

        return sel

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
                   #    (j.jetId & 8) != 0, # killing all events, investigate !
                   (j.chEmEF + j.neEmEF) < 0.9
                   )
        )

        veto_cuts = op.rng_any(
            jets_to_veto, lambda j: corr(j) == 0
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
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the
        https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2024/019"""
        if sample.startswith("TT"):

            def top_pt_weight(pt):
                return op.product((0.103 * op.exp(-0.0118 * pt) - 0.000134 * pt + 0.973),(0.991 + 0.000075*pt) )

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
            # Recommendation: evaluate uncertainty by comparing with/without reweighting.
            # Do NOT apply the reweighting in the opposite direction.
            w_nom = w
            # the variations represent the no-reweighting case
            w_up = op.c_float(1.0)
            w_down = op.c_float(1.0)

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
                    f"{self.parent.git_project_dir}/data/{self.parent.era[:4]}_btagSF_reweight.json",
                    "Ratio_btagSF_shape",
                    params={
                        "year": self.parent.era,
                        # 1. to make it float
                        "jet_multiplicity": 1.0 * op.rng_len(jets),
                    },
                    sel=sel,
                )
                # None since the object is already in the btag_corr
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

    def ak8_bbcc_sf(self, sel, jets):
        """Apply AK8 bb/cc corrections for boosted jets (with systematics)."""
        if self.parent.is_MC:
            sf_file, bb_name, cc_name = AK8_BBCC_SF_FILES[self.parent.era]
            sf_path = os.path.join(self.parent.git_project_dir, "data", sf_file)

            def _corrs(syst):
                bb_corr = get_correction(
                    sf_path,
                    bb_name,
                    params={"pt": lambda j: j.pt, "systematic": syst},
                    sel=sel,
                )
                cc_corr = get_correction(
                    sf_path,
                    cc_name,
                    params={"pt": lambda j: j.pt, "systematic": syst},
                    sel=sel,
                )
                return bb_corr, cc_corr

            def _weight_for(syst):
                bb_corr, cc_corr = _corrs(syst)

                def jet_sf(j):
                    pnet_pass = j.particleNetWithMass_HbbvsQCD >= 0.92
                    is_bb = op.AND(j.hadronFlavour == 5, pnet_pass)
                    is_cc = op.AND(j.hadronFlavour == 4, pnet_pass)
                    return op.switch(
                        is_bb,
                        bb_corr(j),
                        op.switch(is_cc, cc_corr(j), op.c_float(1.0)),
                    )

                return op.rng_product(jets, jet_sf)

            weight_nom = _weight_for("nominal")
            weight_up = _weight_for("up")
            weight_down = _weight_for("down")
            weight_tau21_up = _weight_for("tau21Up")
            weight_tau21_down = _weight_for("tau21Down")
            weight_msd_up = _weight_for("msdUp")
            weight_msd_down = _weight_for("msdDown")

            sel = sel.refine(
                "ak8BBCCSF",
                weight=op.systematic(
                    weight_nom,
                    "ak8BBCC",
                    up=weight_up,
                    down=weight_down,
                ),
            )
            sel = sel.refine(
                "ak8BBCCSF_tau21",
                weight=op.systematic(
                    weight_nom,
                    "ak8BBCC_tau21",
                    up=weight_tau21_up,
                    down=weight_tau21_down,
                ),
            )
            sel = sel.refine(
                "ak8BBCCSF_msd",
                weight=op.systematic(
                    weight_nom,
                    "ak8BBCC_msd",
                    up=weight_msd_up,
                    down=weight_msd_down,
                ),
            )
        else:
            sel = sel.refine("ak8BBCCSF", weight=op.c_float(1.0))
            sel = sel.refine("ak8BBCCSF_tau21", weight=op.c_float(1.0))
            sel = sel.refine("ak8BBCCSF_msd", weight=op.c_float(1.0))

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
                defineOnFirstUse=False,
                sel=sel,
            )

            # Muon ISO SF
            self.muon_ISO_sf = get_correction(
                MUON_SF_JSONFiles[self.parent.era],
                # since muon iso is miniPFreliso and id is medium
                "NUM_TightPFIso_DEN_MediumID",
                systVariations={
                    f"muonIsoSFup": f"{systName}up",
                    f"muonIsoSFdown": f"{systName}down",
                },
                params={
                    "pt": lambda mu: mu.pt,
                    "eta": lambda mu: op.abs(mu.eta),
                },
                systParam="scale_factors",
                systNomName="nominal",
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
                        self.parent.tightElectrons[0].pt >= 30, # 30 since Ele30_WPTight_Gsf HTL path
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
                        self.parent.tightElectrons[0].pt >= 30, # 30 since Ele30_WPTight_Gsf HTL path
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


    def dilepton_trg_sf(self, sel):
        channel = "mumu" if "mumu" in sel.name else "elel" if "elel" in sel.name else "elmu" if "elmu" in sel.name else None
        if channel is None:
            logger.warning(f"Selection name provided: {sel.name}")
            raise RuntimeError(
                "Final state selection name must include one of these values: ee, mumu, emu."
            )
        if self.parent.is_MC:
            if channel == "mumu":
                ch = "mm"
                leading_lepton_pt = self.parent.tightMuons[0].pt
                subleading_lepton_pt = self.parent.tightMuons[1].pt
            elif channel == "elel":
                ch = "ee"
                leading_lepton_pt = self.parent.tightElectrons[0].pt
                subleading_lepton_pt = self.parent.tightElectrons[1].pt
            elif "elmu" in sel.name:
                ch = "mixed"
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
                logger.warning(f"Selection name provided: {sel.name}")
                raise RuntimeError(
                    "Final state selection name must include one of these values: ee, mumu, emu."
                )
            systVariations = {
                "dileptonTRGSFup": "up",
                "dileptonTRGSFdown": "down",
            }
            di_lepton_TRG_SF = get_correction(
                self.di_lepton_TRG_JSONFiles[ch],
                f"sf_{ch}_trg_lepton0_pt-trg_lepton1_pt-trig_ids",
                params={
                    "trg_lepton0_pt": lambda l: leading_lepton_pt,
                    "trg_lepton1_pt": lambda l: subleading_lepton_pt,
                },
                systNomName="nominal",
                systVariations=systVariations,
                systParam="systematic",
                defineOnFirstUse=False,
                sel=sel,
            )
            sel = sel.refine(
                channel+"_di_lepton_TRG_SF", weight=di_lepton_TRG_SF(None)
            )
        else:
            sel = sel.refine(channel+"_di_lepton_TRG_SF",
                             weight=op.c_float(1.))
        self.parent.yields.add(sel, channel + " di-lepton TRG SF")
        return sel

    def Z_pT_reweight(self, sel, sample, GenPartBranch):
        """Apply DY ptll and recoil corrections for given lepton pair."""
        if "mumu" in sel.name:
            pdgId = 13
            ch = "MuMu"
        elif "elel" in sel.name:
            pdgId = 11
            ch = "ElEl"
        else:
            logger.warning(f"Selection name provided: {sel.name}")
            raise RuntimeError(
                "Final state selection name must include one of these values: elel, mumu."
            )
        if self.parent.is_MC and sample.startswith("DY"):

            logger.info(
                "Applying DY ptll and recoil corrections for " + ch + " channel")

            DY_and_Recoil_path = (
                self.parent.git_project_dir + "/data/hleprare/DYandRecoilCorrlib/"
            )

            # for N_unc - consult the json file

            N_unc = 10

            systVariations = {f"ZpT{i}up": f"up{i}" for i in range(1, N_unc + 1)}
            systVariations.update(
                {f"ZpT{i}down": f"down{i}" for i in range(1, N_unc + 1)}
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
                ch + "_ZpT",
                weight=get_Z_pT_corr(get_gen_parts(GenPartBranch, pdgId)),
            )
        else:
            sel = sel.refine(ch + "_ZpT", weight=op.c_float(1.0))
        self.parent.yields.add(sel, ch + " Z pT reweighting")
        return sel

    def V_Jets_Stitching(self, LHEBranch, sel, sample):
        """Apply V+Jets stitching for MC samples DY MLL > 50. pT binned, Jet multiplicity binned and inclusive."""
        if LHEBranch is not None:
            logger.info("Applying V+Jets stitching")
            stitch_map_json = f"{self.parent.git_project_dir}/data/Run3NLOStitching.json"

            # Get the base sample name and corresponding sample number
            base_sample = sample.rsplit(f"_{self.parent.era}", 1)[0]
            sampleNum = sampleNumDict.get(base_sample, None)

            _weight = op.c_float(1.)
            if sampleNum is None:
                if "10to50" in base_sample:
                    pass
                else:
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

                    # Define VpT bin edges and corresponding bin indices
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
