import os
import logging

from bamboo.plots import Plot, SummedPlot, Skim
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

import definitions as defs
from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection
from utils import labeler, ml_input_var_binning

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

        # get DL selections
        [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
            DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
            DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu], _ = makeDLSelection(
                self, noSel, tree, sample)

        # fetch and prepare the input for the model evaluation
        l1, l2, j1, j2, met = defs.ml_input_features(self)

        ml_vars = {
            "event_no": tree.event,
            "weight": noSel.weight,
        }
        ml_vars = ml_vars | l1 | l2 | j1 | j2 | met
        # line above is equivalent to the following (concetanation of dictionaries in dim=0)
        # ml_vars =  {**ml_vars, **l1, **l2, **j1, **j2, **met}

        l1 = op.array('float', *l1.values())
        l2 = op.array('float', *l2.values())
        j1 = op.array('float', *j1.values())
        j2 = op.array('float', *j2.values())
        met = op.array('float', *met.values())

        # load the model file
        split_var = 'even' if tree.event % 2 == 1 else 'odd'
        models_dir = os.path.join(self.git_project_dir, self.mvaModels)
        modelFile = os.path.join(
            models_dir, f"{split_var}_model/model_simplified.onnx")
        logger.info(
            f"Using the following directory for ML models: {models_dir}"
        )

        # evaluate the model
        ml_model = op.mvaEvaluator(modelFile, otherArgs='output')
        ml_output = ml_model(l1, l2, j1, j2, met)

        signal_node = ml_output[0]

        # prepare the labels
        DL_label = {**labeler('DL ML score - blinded'),
                    'blinded-range': [0.25, 0.999]}
        DL_ee_label = {**labeler('DL ML score EE - blinded'),
                       'blinded-range': [0.25, 0.999]}
        DL_mumu_label = {
            **labeler('DL ML score MuMu - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_emu_label = {
            **labeler('DL ML score EMu - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_boosted_label = {
            **labeler('DL Boosted ML score - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_resolved1b_label = {
            **labeler('DL Resolved 1b ML score - blinded'), 'blinded-range': [0.25, 0.999]}
        DL_resolved2b_label = {
            **labeler('DL Resolved 2b ML score - blinded'), 'blinded-range': [0.25, 0.999]}

        # create the plots
        ml_score_1b_ee = Plot.make1D("ml_score_1b_ee", signal_node, DL_resolved_1b_ee, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score 1b ee", plotopts=DL_ee_label
        )
        ml_score_2b_ee = Plot.make1D("ml_score_2b_ee", signal_node, DL_resolved_2b_ee, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score 2b ee", plotopts=DL_ee_label
        )
        ml_score_1b_emu = Plot.make1D("ml_score_1b_emu", signal_node, DL_resolved_1b_emu, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score 1b emu", plotopts=DL_emu_label
        )
        ml_score_2b_emu = Plot.make1D("ml_score_2b_emu", signal_node, DL_resolved_2b_emu, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score 1b emu", plotopts=DL_emu_label
        )
        ml_score_1b_mumu = Plot.make1D("ml_score_1b_mumu", signal_node, DL_resolved_1b_mumu, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score 1b mumu", plotopts=DL_mumu_label
        )
        ml_score_2b_mumu = Plot.make1D("ml_score_2b_mumu", signal_node, DL_resolved_2b_mumu, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score 2b mumu", plotopts=DL_mumu_label
        )
        ml_score_boosted_ee = Plot.make1D("ml_score_boosted_ee", signal_node, DL_boosted_ee, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score boosted ee", plotopts=DL_label
        )
        ml_score_boosted_emu = Plot.make1D("ml_score_boosted_emu", signal_node, DL_boosted_emu, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score bossted emu", plotopts=DL_label
        )
        ml_score_boosted_mumu = Plot.make1D("ml_score_boosted_mumu", signal_node, DL_boosted_mumu, EqBin(
            100, 0, 1.), title='ML', xTitle="ML Score bossted mumu", plotopts=DL_label
        )
        boosted_ml_score_plots = [
            ml_score_boosted_ee, ml_score_boosted_emu, ml_score_boosted_mumu]

        resolved1b_ml_score_plots = [
            ml_score_1b_ee, ml_score_1b_emu, ml_score_1b_mumu]

        resolved2b_ml_score_plots = [
            ml_score_2b_ee, ml_score_2b_emu, ml_score_2b_mumu]

        ml_score_booosted = SummedPlot("DL_boosted_ml_score",
                                       boosted_ml_score_plots,
                                       title="ML",
                                       xTitle="DL Boosted ML score",
                                       plotopts=DL_boosted_label,
                                       )
        ml_score_resolved1b = SummedPlot("DL_resolved1b_ml_score",
                                         resolved1b_ml_score_plots,
                                         title="ML",
                                         xTitle="DL Resolved 1b ML score",
                                         plotopts=DL_resolved1b_label,
                                         )
        ml_score_resolved2b = SummedPlot("DL_resolved2b_ml_score",
                                         resolved2b_ml_score_plots,
                                         title="ML",
                                         xTitle="DL Resolved 2b ML score",
                                         plotopts=DL_resolved2b_label,
                                         )
        ml_score_DL = SummedPlot("DL_ml_score",
                                 [*boosted_ml_score_plots,
                                  *resolved1b_ml_score_plots,
                                  *resolved2b_ml_score_plots],
                                 title="DL ML score",
                                 plotopts=DL_label,
                                 )
        plots.extend([
            *resolved1b_ml_score_plots,
            *resolved2b_ml_score_plots,
            *boosted_ml_score_plots,
            ml_score_booosted,
            ml_score_resolved1b,
            ml_score_resolved2b,
            ml_score_DL
        ])

        event_selections = [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
                            DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
                            DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu]

        # add skims that hold variables for the ML
        for sel in event_selections:
            plots.append(
                Skim(sel.name+"_ml_vars", ml_vars, sel)
            )

        # We're not interested in the following two variables' match between data and MC.
        # Hence they're not included in the input feature plots.
        ml_vars.pop('event_no')
        ml_vars.pop('weight')

        for name, var in ml_vars.items():
            ml_plots = []
            for selection in event_selections:
                ml_plots.append(
                    Plot.make1D(name+"_"+selection.name, var, selection,
                                ml_input_var_binning(name), title=name, xTitle=name)
                )
            plots.extend(ml_plots)
            plots.append(
                SummedPlot(
                    name+"_summed", [plt for plt in ml_plots if plt.name.startswith(name)], title=name+"_summed")
            )

        # labels on plots
        DLboostedEE_label = labeler('DL boosted EE')
        DLboostedMuMu_label = labeler('DL boosted MuMu')
        DLboostedEMU_label = labeler('DL boosted EMu')

        DLresolved_1b_EE_label = labeler('DL resolved 1b EE')
        DLresolved_1b_MuMu_label = labeler('DL resolved 1b MuMu')
        DLresolved_1b_EMu_label = labeler('DL resolved 1b EMu')

        DLresolved_2b_EE_label = labeler('DL resolved 2b EE')
        DLresolved_2b_MuMu_label = labeler('DL resolved 2b MuMu')
        DLresolved_2b_EMu_label = labeler('DL resolved 2b EMu')

        plots.extend([
            # Boosted - fatjet eta
            Plot.make1D("DL_boosted_fatJet_eta_ee", self.ak8BJets[0].eta, DL_boosted_ee, EqBin(
                30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEE_label),
            Plot.make1D("DL_boosted_fatJet_eta_mumu", self.ak8BJets[0].eta, DL_boosted_mumu, EqBin(
                30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedMuMu_label),
            Plot.make1D("DL_boosted_fatJet_eta_emu", self.ak8BJets[0].eta, DL_boosted_emu, EqBin(
                30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEMU_label),

            # Resolved 1b - leading b-jet eta
            Plot.make1D("DL_resolved_1b_leadingJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_1b_ee, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_1b_EE_label),
            Plot.make1D("DL_resolved_1b_leadingJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_1b_mumu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_1b_MuMu_label),
            Plot.make1D("DL_resolved_1b_leadingJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_1b_emu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_1b_EMu_label),

            # Resolved 2b - leading b-jet eta
            Plot.make1D("DL_resolved_2b_leadingJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_2b_ee, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_2b_EE_label),
            Plot.make1D("DL_resolved_2b_leadingJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_2b_mumu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_2b_MuMu_label),
            Plot.make1D("DL_resolved_2b_leadingJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_2b_emu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_2b_EMu_label),
        ])

        return plots
