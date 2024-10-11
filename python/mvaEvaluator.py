import os

from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection
from scalefactors import ScaleFactors as sf
import definitions as defs
from utils import labeler


class mvaEvaluator(NanoBaseHHWWbb):
    """ Class to create control plots, cutflow reports and skims"""

    def __init__(self, args):
        super(mvaEvaluator, self).__init__(args)
        self.channel = self.args.channel
        self.mvaModels = self.args.mvaModels

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # cutflow report
        plots.append(self.yields)

        # define objects
        defs.defineObjects(self, tree)

        # common scale factors
        noSel = sf.commonSF(self, tree, noSel, sample)

        if self.channel == 'DL':
            # get DL selections
            DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu, \
                DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu, \
                DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu = makeDLSelection(
                    self, noSel)

            # muonSF
            DL_boosted_mumu = sf.muonSF(self, DL_boosted_mumu)
            DL_boosted_emu = sf.muonSF(self, DL_boosted_emu)
            DL_resolved_1b_mumu = sf.muonSF(self, DL_resolved_1b_mumu)
            DL_resolved_1b_emu = sf.muonSF(self, DL_resolved_1b_emu)
            DL_resolved_2b_emu = sf.muonSF(self, DL_resolved_2b_emu)
            DL_resolved_2b_mumu = sf.muonSF(self, DL_resolved_2b_mumu)

            # electronSF
            DL_boosted_ee = sf.electronSF(self, DL_boosted_ee)
            DL_boosted_emu = sf.electronSF(self, DL_boosted_emu)
            DL_resolved_1b_ee = sf.electronSF(self, DL_resolved_1b_ee)
            DL_resolved_1b_emu = sf.electronSF(self, DL_resolved_1b_emu)
            DL_resolved_2b_ee = sf.electronSF(self, DL_resolved_2b_ee)
            DL_resolved_2b_emu = sf.electronSF(self, DL_resolved_2b_emu)

        # mva variables

        l1_Px = op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2,
             self.tightElectrons[0].p4.Px()),  # if nElectrons = 2
            (op.rng_len(self.tightMuons) == 2,
             self.tightMuons[0].p4.Px()),  # elif nMuons = 2
            (op.switch(  # else meaning nElectrons = nMuons = 1 since no other case in the DL channel
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Px(), self.tightMuons[0].p4.Px()))
        )
        l2_Px = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.Px()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Px()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Px(), self.tightMuons[0].p4.Px()))
        )
        l1_Py = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].p4.Py()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.Py()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Py(), self.tightMuons[0].p4.Py()))
        )
        l2_Py = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.Py()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Py()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Py(), self.tightMuons[0].p4.Py()))
        )
        l1_Pz = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].p4.Pz()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.Pz()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Pz(), self.tightMuons[0].p4.Pz()))
        )
        l2_Pz = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.Pz()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Pz()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Pz(), self.tightMuons[0].p4.Pz()))
        )
        l1_E = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].p4.E()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.E()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.E(), self.tightMuons[0].p4.E()))
        )
        l2_E = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.E()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.E()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.E(), self.tightMuons[0].p4.E()))
        )
        l1_pdgId = op.multiSwitch(  # static_cast is used to convert the pdgId to float
            (op.rng_len(self.tightElectrons) == 2, op.static_cast(
                'float', self.tightElectrons[0].pdgId)),
            (op.rng_len(self.tightMuons) == 2, op.static_cast(
                'float', self.tightMuons[0].pdgId)),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, op.static_cast('float', self.tightElectrons[0].pdgId), op.static_cast('float', self.tightMuons[0].pdgId)))
        )
        l2_pdgId = op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.static_cast(
                'float', self.tightElectrons[1].pdgId)),
            (op.rng_len(self.tightMuons) == 2, op.static_cast(
                'float', self.tightMuons[1].pdgId)),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, op.static_cast('float', self.tightElectrons[0].pdgId), op.static_cast('float', self.tightMuons[0].pdgId)))
        )
        l1_charge = op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.static_cast(
                'float', self.tightElectrons[0].charge)),
            (op.rng_len(self.tightMuons) == 2, op.static_cast(
                'float', self.tightMuons[0].charge)),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, op.static_cast('float', self.tightElectrons[0].charge), op.static_cast('float', self.tightMuons[0].charge)))
        )
        l2_charge = op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.static_cast(
                'float', self.tightElectrons[1].charge)),
            (op.rng_len(self.tightMuons) == 2, op.static_cast(
                'float', self.tightMuons[1].charge)),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, op.static_cast('float', self.tightElectrons[0].charge), op.static_cast('float', self.tightMuons[0].charge)))
        )

        j1_Px = self.ak4Jets[0].p4.Px()
        j1_Py = self.ak4Jets[0].p4.Py()
        j1_Pz = self.ak4Jets[0].p4.Pz()
        j1_E = self.ak4Jets[0].p4.E()
        j1_btag = self.ak4Jets[0].btagPNetB
        j2_Px = self.ak4Jets[1].p4.Px()
        j2_Py = self.ak4Jets[1].p4.Py()
        j2_Pz = self.ak4Jets[1].p4.Pz()
        j2_E = self.ak4Jets[1].p4.E()
        j2_btag = self.ak4Jets[1].btagPNetB

        met_Px = op.product(tree.MET.pt, op.cos(tree.MET.phi))
        met_Py = op.product(tree.MET.pt, op.sin(tree.MET.phi))
        met_E = tree.MET.pt

        #############################################################################
        #                            MVA evaluation                                 #
        #############################################################################
        if self.mvaModels and self.channel == 'DL':
            DL_DNN = {**labeler('DL DNN score - blinded'),
                      'blinded-range': [0.25, 0.999]}
            DL_DNN_EE = {**labeler('DL DNN score EE - blinded'),
                         'blinded-range': [0.25, 0.999]}
            DL_DNN_MuMu = {
                **labeler('DL DNN score MuMu - blinded'), 'blinded-range': [0.25, 0.999]}
            DL_DNN_EMu = {
                **labeler('DL DNN score EMu - blinded'), 'blinded-range': [0.25, 0.999]}

            # prepare the input for the model
            l1 = op.array('float', *[l1_Px, l1_Py, l1_Pz, l1_E, l1_pdgId, l1_charge,])
            l2 = op.array('float', *[l2_Px, l2_Py, l2_Pz, l2_E, l2_pdgId, l2_charge,])
            j1 = op.array('float', *[j1_Px, j1_Py, j1_Pz, j1_E, j1_btag])
            j2 = op.array('float', *[j2_Px, j2_Py, j2_Pz, j2_E, j2_btag])
            met = op.array('float', *[met_Px, met_Py, met_E])
            
            # load the model
            split_var = 'even' if tree.event % 2 == 1 else 'odd'
            model = os.path.join(self.mvaModels, f"{split_var}_model/model.onnx")
            # evaluate the model
            dnn = op.mvaEvaluator(model, otherArgs='output')
            DNN_output = dnn(l1, l2, j1, j2, met)

            dnn_score_1b_ee = Plot.make1D("dnn_score_1b_ee", DNN_output[0], DL_resolved_1b_ee, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EE
            )
            dnn_score_2b_ee = Plot.make1D("dnn_score_2b_ee", DNN_output[0], DL_resolved_2b_ee, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EE
            )
            dnn_score_1b_emu = Plot.make1D("dnn_score_1b_emu", DNN_output[0], DL_resolved_1b_emu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EMu
            )
            dnn_score_2b_emu = Plot.make1D("dnn_score_2b_emu", DNN_output[0], DL_resolved_2b_emu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EMu
            )
            dnn_score_1b_mumu = Plot.make1D("dnn_score_1b_mumu", DNN_output[0], DL_resolved_1b_mumu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_MuMu
            )
            dnn_score_2b_mumu = Plot.make1D("dnn_score_2b_mumu", DNN_output[0], DL_resolved_2b_mumu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_MuMu
            )
            dnn_score_boosted_ee = Plot.make1D("dnn_score_boosted_ee", DNN_output[0], DL_boosted_ee, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN
            )
            dnn_score_boosted_emu = Plot.make1D("dnn_score_boosted_emu", DNN_output[0], DL_boosted_emu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN
            )
            dnn_score_boosted_mumu = Plot.make1D("dnn_score_boosted_mumu", DNN_output[0], DL_boosted_mumu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN
            )
            plots.extend([
                dnn_score_1b_ee,
                dnn_score_2b_ee,
                dnn_score_1b_emu,
                dnn_score_2b_emu,
                dnn_score_1b_mumu,
                dnn_score_2b_mumu,
                dnn_score_boosted_ee,
                dnn_score_boosted_emu,
                dnn_score_boosted_mumu,
                SummedPlot("DL_dnn_score", [
                           dnn_score_1b_ee, dnn_score_2b_ee, dnn_score_1b_emu, dnn_score_2b_emu, dnn_score_1b_mumu, dnn_score_2b_mumu, dnn_score_boosted_ee, dnn_score_boosted_emu, dnn_score_boosted_mumu], title="DL DNN score", plotopts=DL_DNN
                           )
            ])

        return plots
