import os
import re
import logging
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.analysismodules import NanoAODHistoModule

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


class NanoBaseHHWWbb(NanoAODHistoModule):
    """Base module for HH->WWbb analysis"""

    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument(
            "-c",
            "--channel",
            dest="channel",
            type=str,
            default="DL",
            help="Channel to be selected between SL and DL",
        )
        parser.add_argument(
            "--mvaModel",
            dest="mvaModel",
            type=str,
            default="data/model.onnx",
            help="Path to XGBoost model.",
        )
        parser.add_argument(
            "--backend",
            type=str,
            default="dataframe",
            help="Backend to use, 'dataframe' (default), 'lazy', or 'compiled'",
        )
        parser.add_argument(
            "--sync", action="store_true", default=False, help="Run synchronisation"
        )

    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):

        # Define the git project's directory
        self.git_project_dir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..")
        )
        self.era = sampleCfg["era"] if sampleCfg else None
        self.is_MC = self.isMC(sample)
        self.sampleCfg = sampleCfg

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

        # Inclusive cross-section uncertainties (rate-only)
        noSel = self._apply_inclusive_xs_uncertainties(noSel, sample, sampleCfg)

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
        addHLTPath("Muon_", "IsoMu24")
        addHLTPath("EGamma_", "Ele30_WPTight_Gsf")
        addHLTPath("EGamma_", "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
        addHLTPath("EGamma_", "DoubleEle33_CaloIdL_MW")
        addHLTPath("EGamma_", "Ele50_CaloIdVT_GsfTrkIdT_PFJet165")
        addHLTPath("MuonEG_", "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
        addHLTPath("MuonEG_", "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
        addHLTPath("MuonEG_", "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")
        addHLTPath("MuonEG_", "IsoMu24")
        addHLTPath("MuonEG_", "Ele30_WPTight_Gsf")
        addHLTPath("MuonEG_", "Ele50_CaloIdVT_GsfTrkIdT_PFJet165")

        

        if self.is_MC:
            noSel = noSel.refine(
                "triggers",
                cut=(op.OR(*chain.from_iterable(self.triggers_per_PD.values()))),
            )
        else:
            noSel = noSel.refine(
                "triggers",
                cut=makeMultiPrimaryDatasetTriggerSelection(
                    sample, self.triggers_per_PD
                ),
            )

        self.yields.add(noSel, "triggers")

        return tree, noSel, be, lumiArgs

    def _apply_inclusive_xs_uncertainties(self, sel, sample, sampleCfg):
        """Apply inclusive cross-section uncertainties as rate-only systematics.

        These are intended for use in the final fit as correlated rate uncertainties.
        TTbar and DY are excluded (free-floating normalizations).
        """
        if not sampleCfg:
            return sel

        if sampleCfg.get("type") == "data":
            return sel

        if sampleCfg.get("group") in {"TT", "DY"}:
            return sel

        # HH signal samples are not covered by the tables; skip by default
        if sampleCfg.get("type") == "signal":
            return sel

        smp_name = sample or ""
        smp_lower = smp_name.lower()

        def has_any(*tokens):
            return any(tok in smp_lower for tok in tokens)

        def add_rate_uncertainty(sel_in, syst_name, up_frac, down_frac):
            return sel_in.refine(
                f"{syst_name}__xsrate",
                weight=op.systematic(
                    op.c_float(1.0),
                    syst_name,
                    up=op.c_float(1.0 + up_frac),
                    down=op.c_float(1.0 - down_frac),
                ),
            )

        # Process matching helpers
        is_tW = has_any("twminus", "tbarwplus")
        is_st_tchannel = has_any("TBbarto", "TbarBto")
        is_st_schannel = has_any("TBbarQ", "TbarBQ")
        is_ttW = has_any("ttll")
        is_ttZ = has_any("ttz")
        is_tttt = has_any("tttt")
        is_Wjets = sampleCfg.get("group") == "VJets"
        is_WW = has_any("ww_")
        is_WZ = has_any("wz_")
        is_ZZ = has_any("zz_")

        # Single-H backgrounds
        is_Hggf = has_any("glugluhto")
        is_Hvbf = has_any("vbfhto")
        is_WH = has_any("wplush_")
        is_ttH = has_any("tth_")
        is_ZH = has_any("zh_")
        is_tHq = has_any("thq_")
        is_tHW = has_any("thw")

        # QCD scale uncertainties (asymmetric)
        if is_st_tchannel:
            sel = add_rate_uncertainty(sel, "QCDscale_ttbar", 0.011, 0.008)
        if is_st_schannel:
            sel = add_rate_uncertainty(sel, "QCDscale_ttbar", 0.005, 0.004)
        if is_tW:
            sel = add_rate_uncertainty(sel, "QCDscale_ttbar", 0.023, 0.022)
        if is_ttW:
            sel = add_rate_uncertainty(sel, "QCDscale_ttbar", 0.261, 0.162)
        if is_ttZ:
            sel = add_rate_uncertainty(sel, "QCDscale_ttbar", 0.086, 0.095)
        if is_tttt:
            sel = add_rate_uncertainty(sel, "QCDscale_ttbar", 0.082, 0.175)

        if is_Wjets:
            sel = add_rate_uncertainty(sel, "QCDscale_V", 0.012, 0.013)

        if is_WW:
            sel = add_rate_uncertainty(sel, "QCDscale_VV", 0.025, 0.022)
        if is_WZ:
            sel = add_rate_uncertainty(sel, "QCDscale_VV", 0.041, 0.032)
        if is_ZZ:
            sel = add_rate_uncertainty(sel, "QCDscale_VV", 0.029, 0.027)

        if is_Hggf:
            sel = add_rate_uncertainty(sel, "QCDscale_ggH", 0.046, 0.067)
        if is_Hvbf:
            sel = add_rate_uncertainty(sel, "QCDscale_qqH", 0.005, 0.003)
        if is_WH:
            sel = add_rate_uncertainty(sel, "QCDscale_VH", 0.004, 0.007)
        if is_ZH:
            sel = add_rate_uncertainty(sel, "QCDscale_VH", 0.037, 0.032)
        if is_ttH:
            sel = add_rate_uncertainty(sel, "QCDscale_ttH", 0.060, 0.093)
        if is_tHq:
            sel = add_rate_uncertainty(sel, "QCDscale_ttH", 0.065, 0.148)
        if is_tHW:
            sel = add_rate_uncertainty(sel, "QCDscale_ttH", 0.050, 0.068)

        # PDF+alphaS uncertainties (symmetric unless noted)
        if is_tW:
            sel = add_rate_uncertainty(sel, "PDFalphaS_gq", 0.027, 0.027)
        if is_tttt:
            sel = add_rate_uncertainty(sel, "PDFalphaS_gg", 0.067, 0.067)
        # NOTE: ttZ appears in multiple initial-state groupings in the table.
        # Here we keep it in the qq-group (adjust if you prefer gg).
        if is_ttZ:
            sel = add_rate_uncertainty(sel, "PDFalphaS_qq", 0.023, 0.023)
        if is_st_tchannel:
            sel = add_rate_uncertainty(sel, "PDFalphaS_qq", 0.015, 0.009)
        if is_Wjets:
            sel = add_rate_uncertainty(sel, "PDFalphaS_qq", 0.007, 0.007)
        if is_ttW:
            sel = add_rate_uncertainty(sel, "PDFalphaS_qq", 0.021, 0.021)

        if is_Hggf:
            sel = add_rate_uncertainty(sel, "PDFalphaS_ggH", 0.032, 0.032)
        if is_Hvbf:
            sel = add_rate_uncertainty(sel, "PDFalphaS_qqH", 0.021, 0.021)
        if is_WH:
            sel = add_rate_uncertainty(sel, "PDFalphaS_VH", 0.018, 0.018)
        if is_ZH:
            sel = add_rate_uncertainty(sel, "PDFalphaS_VH", 0.016, 0.016)
        if is_ttH:
            sel = add_rate_uncertainty(sel, "PDFalphaS_ttH", 0.035, 0.035)
        if is_tHq:
            sel = add_rate_uncertainty(sel, "PDFalphaS_ttH", 0.037, 0.037)
        if is_tHW:
            sel = add_rate_uncertainty(sel, "PDFalphaS_ttH", 0.063, 0.063)

        return sel

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """Postprocess: run plotIt

        The list of plots is created if needed (from a representative file,
        this enables rerunning the postprocessing step on the results files),
        and then plotIt is executed
        """
        import os

        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir, config=config)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport

        plotList_cutflowreport = [
            ap for ap in self.plotList if isinstance(ap, CutFlowReport)
        ]
        plotList_plotIt = [
            ap
            for ap in self.plotList
            if (isinstance(ap, Plot) or isinstance(ap, DerivedPlot))
            and len(ap.binnings) == 1
        ]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            from bamboo.analysisutils import printCutFlowReports

            printCutFlowReports(
                config,
                plotList_cutflowreport,
                workdir=workdir,
                resultsdir=resultsdir,
                readCounters=self.readCounters,
                eras=(eraMode, eras),
                verbose=self.args.verbose,
            )
        if plotList_plotIt and not self.args.sync:
            from bamboo.analysisutils import writePlotIt, runPlotIt

            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(
                config,
                plotList_plotIt,
                cfgName,
                eras=eras,
                workdir=workdir,
                resultsdir=resultsdir,
                readCounters=self.readCounters,
                plotDefaults=self.plotDefaults,
                vetoFileAttributes=self.__class__.CustomSampleAttributes,
            )
            # runPlotIt(
            #     cfgName,
            #     workdir=workdir,
            #     plotIt=self.args.plotIt,
            #     eras=(eraMode, eras),
            #     verbose=self.args.verbose,
            # )
            # hadd signal files and create another plots.yml called plots_full.yml
            plotsDir = "plots_full"
            import os
            import shutil

            outDir = os.path.join(resultsdir, "normalizedSummedSignal")
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            utils.custom_Plotit(
                cfgName,
                workdir,
                resultsdir,
                outDir,
                self.readCounters,
                config,
                plotIt=self.args.plotIt,
                verbose=True,
            )
            # end of merging signal samples and plotting

        # create pdf presentation
        if not self.mvaModel and not self.args.sync:
            try:
                for era in eras:
                    os.system(
                        utils.runPDF(
                            workdir=workdir, channel=self.args.channel, era=era
                        )
                    )
                    logger.info(f"PDF presentation created for era {era}.\n")
                os.system(
                    utils.runPDF(
                        workdir=workdir, channel=self.args.channel, plotsDir=plotsDir
                    )
                )
                logger.info(f"PDF presentation created for all eras combined.")
            except Exception as e:
                logger.info(e)

        from bamboo.plots import Skim

        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]

        from bamboo.analysisutils import loadPlotIt

        _, samples, _, _, _ = loadPlotIt(
            config,
            [],
            eras=self.args.eras[1],
            workdir=workdir,
            resultsdir=resultsdir,
            readCounters=self.readCounters,
            vetoFileAttributes=self.__class__.CustomSampleAttributes,
        )

        # create sync skims if asked
        if self.args.sync:
            if skims:
                from bamboo.root import gbl
                import pandas as pd

                sync_dfs = []
                for skim in skims:
                    frames = []
                    for smp in samples:
                        for cb in smp.files if hasattr(smp, "files") else [smp]:
                            tree = cb.tFile.Get(skim.treeName)
                            if not tree:
                                logger.info(
                                    "WARNING: skim tree %s not found in file %s"
                                    % (skim.treeName, cb.tFile.GetName())
                                )
                                logger.info("         skipping...")
                            else:
                                N = tree.GetEntries()
                                cols = gbl.ROOT.RDataFrame(tree).AsNumpy()
                                if "sync" not in skim.name:
                                    cols["weight"] *= cb.scale
                                    cols["process"] = [smp.name] * len(cols["weight"])
                                frames.append(pd.DataFrame(cols))
                    df = pd.concat(frames)
                    df = df[self.order]
                    sync_dfs.append(df)
                df = pd.concat(sync_dfs)
                df = df.sort_values(by="event_no")
                syncFileName = f"{self.channel}_sync.csv"
                df.to_csv(os.path.join(resultsdir, syncFileName))
                logger.info(f"Saved dataframe for sync to {syncFileName}")
            else:
                logger.warning("No skims are found, hence sync file is not produced.")
