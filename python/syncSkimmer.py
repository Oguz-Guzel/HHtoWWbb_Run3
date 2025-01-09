
from bamboo.plots import Plot, Skim
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection
from scalefactors import ScaleFactors as sf
import definitions as defs


class syncSkimmer(NanoBaseHHWWbb):
    """ Class to create control plots, cutflow reports and skims"""

    def __init__(self, args):
        super(syncSkimmer, self).__init__(args)
        self.channel = self.args.channel
        self.mvaModels = None

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # cutflow report
        plots.append(self.yields)

        # define objects
        defs.defineObjects(self, tree)

        # top pT reweighting
        noSel = sf.top_pT_reweight(self, tree, noSel, sample)

        # btag scale factors
        noSel = sf.btagSF(self, noSel)

        # btag rescaling
        noSel = sf.btagRescale(self, noSel)

        # Noise filters
        noSel = sf.NoiseFilters(self, tree, noSel)

        # get DL selections
        [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
            DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
            DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu], _ = makeDLSelection(
                self, noSel)

        # example code
        DL_1_preElectron = noSel.refine(
            "DL_1_preElectron", cut=op.rng_len(self.preElectrons) == 1)
        self.yields.add(DL_1_preElectron, '1 pre-Electron')

        DL_1_preMuon = noSel.refine(
            "DL_1_preMuon", cut=op.rng_len(self.preMuons) == 1)
        self.yields.add(DL_1_preMuon, '1 pre-Muon')
        # end of example code

        # cutflow report for DL channel
        self.yields.add(DL_boosted_ee, 'DL boosted ee')
        self.yields.add(DL_boosted_mumu, 'DL boosted mumu')
        self.yields.add(DL_boosted_emu, 'DL boosted emu')
        self.yields.add(DL_resolved_1b_ee, 'DL resolved 1b ee')
        self.yields.add(DL_resolved_2b_ee, 'DL resolved 2b ee')
        self.yields.add(DL_resolved_1b_mumu, 'DL resolved 1b mumu')
        self.yields.add(DL_resolved_2b_mumu, 'DL resolved 2b mumu')
        self.yields.add(DL_resolved_1b_emu, 'DL resolved 1b emu')
        self.yields.add(DL_resolved_2b_emu, 'DL resolved 2b emu')

        ### Syncronisaton ###
        syncVars_DL = {
            "event_no": tree.event,
            "lumi_block": tree.luminosityBlock,
            "is_dl_ee": op.switch(op.rng_len(self.tightElectrons) == 2, 1, 0),
            "is_dl_mumu": op.switch(op.rng_len(self.tightMuons) == 2, 1, 0),
            "is_dl_emu": op.switch(op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons)) == 1, 1, 0),
            "nLooseElectron": op.rng_len(self.preElectrons),
            "nFakeElectron": op.rng_len(self.fakeElectrons),
            "nTightElectron": op.rng_len(self.tightElectrons),
            "nLooseMuon": op.rng_len(self.preMuons),
            "nFakeMuon": op.rng_len(self.fakeMuons),
            "nTightMuon": op.rng_len(self.tightMuons),
            "lepton0pt": op.multiSwitch(
                (op.rng_len(self.tightElectrons)
                    == 2, self.tightElectrons[0].pt),
                (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].pt),
                (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                    self.tightElectrons[0].pt >= self.tightMuons[0].pt, self.tightElectrons[0].pt, self.tightMuons[0].pt)),
                op.c_float(-9999.)
            ),
            "lepton1pt": op.multiSwitch(
                (op.rng_len(self.tightElectrons)
                    == 2, self.tightElectrons[1].pt),
                (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].pt),
                (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                    self.tightElectrons[0].pt >= self.tightMuons[0].pt, self.tightMuons[0].pt, self.tightElectrons[0].pt)),
                op.c_float(-9999.)
            ),
            "nAK4": op.rng_len(self.ak4Jets),
            "nAK4bJet": op.rng_len(self.ak4BJets),
            "nAK8bJet": op.rng_len(self.ak8BJets),
            "ak4jet0pt": op.switch(op.rng_len(self.ak4Jets) > 0, self.ak4Jets[0].pt, op.c_float(-9999.)),
            "ak4jet0eta": op.switch(op.rng_len(self.ak4Jets) > 0, self.ak4Jets[0].eta, op.c_float(-9999.)),
            "ak4jet1pt": op.switch(op.rng_len(self.ak4Jets) > 1, self.ak4Jets[1].pt, op.c_float(-9999.)),
            "ak4jet1eta": op.switch(op.rng_len(self.ak4Jets) > 1, self.ak4Jets[1].eta, op.c_float(-9999.)),
            "ak8jetpt": op.switch(op.rng_len(self.ak8Jets) > 0, self.ak8Jets[0].pt, op.c_float(-9999.)),
            "ak8jeteta": op.switch(op.rng_len(self.ak8Jets) > 0, self.ak8Jets[0].eta, op.c_float(-9999.)),
        }
        # to order the columns
        self.order = [key for key in syncVars_DL.keys()]
        plots.extend([
            Skim("DL_resolved_1b_ee_sync", syncVars_DL, DL_resolved_1b_ee),
            Skim("DL_resolved_2b_ee_sync", syncVars_DL, DL_resolved_2b_ee),
            Skim("DL_resolved_1b_mumu_sync",
                 syncVars_DL, DL_resolved_1b_mumu),
            Skim("DL_resolved_2b_mumu_sync",
                 syncVars_DL, DL_resolved_2b_mumu),
            Skim("DL_resolved_1b_emu_sync",
                 syncVars_DL, DL_resolved_1b_emu),
            Skim("DL_resolved_2b_emu_sync",
                 syncVars_DL, DL_resolved_2b_emu),
            Skim("DL_boosted_ee_sync", syncVars_DL, DL_boosted_ee),
            Skim("DL_boosted_mumu_sync", syncVars_DL, DL_boosted_mumu),
            Skim("DL_boosted_emu_sync", syncVars_DL, DL_boosted_emu),
            Plot.make1D("DL_boosted_nfatJet_ee", op.rng_len(self.ak8Jets), DL_boosted_ee, EqBin(
                10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet")
        ])

        return plots
