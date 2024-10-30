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

            l1, l2, j1, j2, met = defs.ml_input_features(self, tree)
            # prepare the input for the model
            l1 = op.array('float', *l1.values())
            l2 = op.array('float', *l2.values())
            j1 = op.array('float', *j1.values())
            j2 = op.array('float', *j2.values())
            met = op.array('float', *met.values())

            # load the model
            split_var = 'even' if tree.event % 2 == 1 else 'odd'
            model = os.path.join(
                self.mvaModels, f"{split_var}_model/model_simplified.onnx")
            # evaluate the model
            dnn = op.mvaEvaluator(model, otherArgs='output')
            DNN_output = dnn(l1, l2, j1, j2, met)

            signal_node = DNN_output[0]

            dnn_score_1b_ee = Plot.make1D("dnn_score_1b_ee", signal_node, DL_resolved_1b_ee, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EE
            )
            dnn_score_2b_ee = Plot.make1D("dnn_score_2b_ee", signal_node, DL_resolved_2b_ee, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EE
            )
            dnn_score_1b_emu = Plot.make1D("dnn_score_1b_emu", signal_node, DL_resolved_1b_emu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EMu
            )
            dnn_score_2b_emu = Plot.make1D("dnn_score_2b_emu", signal_node, DL_resolved_2b_emu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EMu
            )
            dnn_score_1b_mumu = Plot.make1D("dnn_score_1b_mumu", signal_node, DL_resolved_1b_mumu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_MuMu
            )
            dnn_score_2b_mumu = Plot.make1D("dnn_score_2b_mumu", signal_node, DL_resolved_2b_mumu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_MuMu
            )
            dnn_score_boosted_ee = Plot.make1D("dnn_score_boosted_ee", signal_node, DL_boosted_ee, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN
            )
            dnn_score_boosted_emu = Plot.make1D("dnn_score_boosted_emu", signal_node, DL_boosted_emu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN
            )
            dnn_score_boosted_mumu = Plot.make1D("dnn_score_boosted_mumu", signal_node, DL_boosted_mumu, EqBin(
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
