import os
import logging

from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection
from scalefactors import ScaleFactors as sf
import definitions as defs
from utils import labeler

logger = logging.getLogger(__name__)


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

        #############################################################################
        #                            MVA evaluation                                 #
        #############################################################################
        DL_label = {**labeler('DL TNN score - blinded'),
                    'blinded-range': [0.25, 0.999]}
        DL_ee_label = {**labeler('DL TNN score EE - blinded'),
                       'blinded-range': [0.25, 0.999]}
        DL_mumu_label = {
            **labeler('DL TNN score MuMu - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_emu_label = {
            **labeler('DL TNN score EMu - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_boosted_label = {
            **labeler('DL Boosted TNN score - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_resolved1b_label = {
            **labeler('DL Resolved 1b TNN score - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_resolved2b_label = {
            **labeler('DL Resolved 2b TNN score - blinded'), 'blinded-range': [0.25, 0.999]}
        l1, l2, j1, j2, met = defs.ml_input_features(self, tree)
        # prepare the input for the model
        l1 = op.array('float', *l1.values())
        l2 = op.array('float', *l2.values())
        j1 = op.array('float', *j1.values())
        j2 = op.array('float', *j2.values())
        met = op.array('float', *met.values())

        # load the model file
        split_var = 'even' if tree.event % 2 == 1 else 'odd'
        modelFile = os.path.join(
            self.mvaModels, f"{split_var}_model/model_simplified.onnx")
        logger.info(f"Using model file: {modelFile}")

        # evaluate the model
        tnn_model = op.mvaEvaluator(modelFile, otherArgs='output')
        tnn_output = tnn_model(l1, l2, j1, j2, met)

        signal_node = tnn_output[0]

        tnn_score_1b_ee = Plot.make1D("tnn_score_1b_ee", signal_node, DL_resolved_1b_ee, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score 1b ee", plotopts=DL_ee_label
        )
        tnn_score_2b_ee = Plot.make1D("tnn_score_2b_ee", signal_node, DL_resolved_2b_ee, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score 2b ee", plotopts=DL_ee_label
        )
        tnn_score_1b_emu = Plot.make1D("tnn_score_1b_emu", signal_node, DL_resolved_1b_emu, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score 1b emu", plotopts=DL_emu_label
        )
        tnn_score_2b_emu = Plot.make1D("tnn_score_2b_emu", signal_node, DL_resolved_2b_emu, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score 1b emu", plotopts=DL_emu_label
        )
        tnn_score_1b_mumu = Plot.make1D("tnn_score_1b_mumu", signal_node, DL_resolved_1b_mumu, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score 1b mumu", plotopts=DL_mumu_label
        )
        tnn_score_2b_mumu = Plot.make1D("tnn_score_2b_mumu", signal_node, DL_resolved_2b_mumu, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score 2b mumu", plotopts=DL_mumu_label
        )
        tnn_score_boosted_ee = Plot.make1D("tnn_score_boosted_ee", signal_node, DL_boosted_ee, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score boosted ee", plotopts=DL_label
        )
        tnn_score_boosted_emu = Plot.make1D("tnn_score_boosted_emu", signal_node, DL_boosted_emu, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score bossted emu", plotopts=DL_label
        )
        tnn_score_boosted_mumu = Plot.make1D("tnn_score_boosted_mumu", signal_node, DL_boosted_mumu, EqBin(
            100, 0, 1.), title='TNN', xTitle="TNN Score bossted mumu", plotopts=DL_label
        )
        boosted_tnn_score_plots = [
            tnn_score_boosted_ee, tnn_score_boosted_emu, tnn_score_boosted_mumu]

        resolved1b_tnn_score_plots = [
            tnn_score_1b_ee, tnn_score_1b_emu, tnn_score_1b_mumu]

        resolved2b_tnn_score_plots = [
            tnn_score_2b_ee, tnn_score_2b_emu, tnn_score_2b_mumu]

        tnn_score_booosted = SummedPlot("DL_boosted_tnn_score",
                                        boosted_tnn_score_plots,
                                        title="TNN",
                                        xTitle="DL Boosted TNN score",
                                        plotopts=DL_boosted_label,
                                        )
        tnn_score_resolved1b = SummedPlot("DL_resolved1b_tnn_score",
                                          resolved1b_tnn_score_plots,
                                          title="TNN",
                                          xTitle="DL Resolved 1b TNN score",
                                          plotopts=DL_resolved1b_label,
                                          )
        tnn_score_resolved2b = SummedPlot("DL_resolved2b_tnn_score",
                                          resolved2b_tnn_score_plots,
                                          title="TNN",
                                          xTitle="DL Resolved 2b TNN score",
                                          plotopts=DL_resolved2b_label,
                                          )
        tnn_score_DL = SummedPlot("DL_tnn_score",
                                  [*boosted_tnn_score_plots,
                                   *resolved1b_tnn_score_plots,
                                   *resolved2b_tnn_score_plots],
                                  title="DL TNN score",
                                  plotopts=DL_label,
                                  )
        plots.extend([
            *resolved1b_tnn_score_plots,
            *resolved2b_tnn_score_plots,
            *boosted_tnn_score_plots,
            tnn_score_booosted,
            tnn_score_resolved1b,
            tnn_score_resolved2b,
            tnn_score_DL
        ])

        return plots
