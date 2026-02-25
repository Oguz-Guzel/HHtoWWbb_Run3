import logging

from bamboo.plots import Plot, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from definitions import ml_input_features
from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection
from utils import labeler

logger = logging.getLogger(__name__)


class mvaEvaluator(NanoBaseHHWWbb):
    """Class to create MVA distribution plots"""

    def __init__(self, args):
        super().__init__(args)
        self.mvaModel = self.args.mvaModel

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # cutflow report
        plots.append(self.yields)

        # get DL selections
        [
            DL_boosted_ee,
            DL_boosted_mumu,
            DL_boosted_emu,
            DL_resolved_1b_ee,
            DL_resolved_1b_mumu,
            DL_resolved_1b_emu,
            DL_resolved_2b_ee,
            DL_resolved_2b_mumu,
            DL_resolved_2b_emu,
            DL_VBF_resolved_ee,
            DL_VBF_resolved_mumu,
            DL_VBF_resolved_emu,
            DL_VBF_boosted_ee,
            DL_VBF_boosted_mumu,
            DL_VBF_boosted_emu,
        ] = makeDLSelection(self, noSel, tree, sample)

        ml_vars = ml_input_features(self)
        # ml_vars["event_no"] = tree.event
        # ml_vars["weight"] = noSel.weight

        # prepare the input variables for the ML model
        inputs = ml_vars.values()

        cast_inputs = []
        for var in inputs:
            cast_inputs.append(op.static_cast("float", var))

        inputs = op.array("float", *cast_inputs)

        # load the model file
        model_file = self.mvaModel
        logger.info(f"Using the ML model: {model_file}")

        ml_evaluator = op.mvaEvaluator(model_file, otherArgs="probabilities")

        ml_output = ml_evaluator(inputs)

        # get the ML scores
        signal_node = ml_output[1]

        # prepare the labels
        DL_label = {"blinded-range": [0.25, 0.999]}
        DL_ee_label = {
            **labeler("DL ML score EE"),
            # "blinded-range": [0.25, 0.999],
        }
        DL_mumu_label = {
            **labeler("DL ML score MuMu"),
            # "blinded-range": [0.25, 0.999],
        }
        DL_emu_label = {
            **labeler("DL ML score EMu"),
            # "blinded-range": [0.25, 0.999],
        }
        DL_boosted_label = {
            **labeler("DL Boosted ML score"),
            # "blinded-range": [0.25, 0.999],
        }
        DL_resolved1b_label = {
            **labeler("DL Resolved 1b ML score"),
            # "blinded-range": [0.25, 0.999],
        }
        DL_resolved2b_label = {
            **labeler("DL Resolved 2b ML score"),
            # "blinded-range": [0.25, 0.999],
        }

        # create the plots
        ml_score_1b_ee = Plot.make1D(
            "ml_score_1b_ee",
            signal_node,
            DL_resolved_1b_ee,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score 1b ee",
            plotopts=DL_ee_label,
        )
        ml_score_2b_ee = Plot.make1D(
            "ml_score_2b_ee",
            signal_node,
            DL_resolved_2b_ee,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score 2b ee",
            plotopts=DL_ee_label,
        )
        ml_score_1b_emu = Plot.make1D(
            "ml_score_1b_emu",
            signal_node,
            DL_resolved_1b_emu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score 1b emu",
            plotopts=DL_emu_label,
        )
        ml_score_2b_emu = Plot.make1D(
            "ml_score_2b_emu",
            signal_node,
            DL_resolved_2b_emu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score 1b emu",
            plotopts=DL_emu_label,
        )
        ml_score_1b_mumu = Plot.make1D(
            "ml_score_1b_mumu",
            signal_node,
            DL_resolved_1b_mumu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score 1b mumu",
            plotopts=DL_mumu_label,
        )
        ml_score_2b_mumu = Plot.make1D(
            "ml_score_2b_mumu",
            signal_node,
            DL_resolved_2b_mumu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score 2b mumu",
            plotopts=DL_mumu_label,
        )
        ml_score_boosted_ee = Plot.make1D(
            "ml_score_boosted_ee",
            signal_node,
            DL_boosted_ee,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score boosted ee",
            # plotopts=DL_label,
        )
        ml_score_boosted_emu = Plot.make1D(
            "ml_score_boosted_emu",
            signal_node,
            DL_boosted_emu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score bossted emu",
            # plotopts=DL_label,
        )
        ml_score_boosted_mumu = Plot.make1D(
            "ml_score_boosted_mumu",
            signal_node,
            DL_boosted_mumu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score bossted mumu",
            # plotopts=DL_label,
        )
        ml_score_VBF_resolved_ee = Plot.make1D(
            "ml_score_VBF_resolved_ee",
            signal_node,
            DL_VBF_resolved_ee,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score VBF resolved ee",
            plotopts=DL_ee_label,
        )
        ml_score_VBF_resolved_emu = Plot.make1D(
            "ml_score_VBF_resolved_emu",
            signal_node,
            DL_VBF_resolved_emu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score VBF resolved emu",
            plotopts=DL_emu_label,
        )
        ml_score_VBF_resolved_mumu = Plot.make1D(
            "ml_score_VBF_resolved_mumu",
            signal_node,
            DL_VBF_resolved_mumu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score VBF resolved mumu",
            plotopts=DL_mumu_label,
        )
        ml_score_VBF_boosted_ee = Plot.make1D(
            "ml_score_VBF_boosted_ee",
            signal_node,
            DL_VBF_boosted_ee,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score VBF boosted ee",
            plotopts=DL_ee_label,
        )
        ml_score_VBF_boosted_emu = Plot.make1D(
            "ml_score_VBF_boosted_emu",
            signal_node,
            DL_VBF_boosted_emu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score VBF boosted emu",
            plotopts=DL_emu_label,
        )
        ml_score_VBF_boosted_mumu = Plot.make1D(
            "ml_score_VBF_boosted_mumu",
            signal_node,
            DL_VBF_boosted_mumu,
            EqBin(50, 0, 1.0),
            title="ML",
            xTitle="ML Score VBF boosted mumu",
            plotopts=DL_mumu_label,
        )
        boosted_ml_score_plots = [
            ml_score_boosted_ee,
            ml_score_boosted_emu,
            ml_score_boosted_mumu,
        ]

        resolved1b_ml_score_plots = [ml_score_1b_ee, ml_score_1b_emu, ml_score_1b_mumu]

        resolved2b_ml_score_plots = [ml_score_2b_ee, ml_score_2b_emu, ml_score_2b_mumu]

        VBF_resolved_ml_score_plots = [
            ml_score_VBF_resolved_ee,
            ml_score_VBF_resolved_emu,
            ml_score_VBF_resolved_mumu,
        ]

        VBF_boosted_ml_score_plots = [
            ml_score_VBF_boosted_ee,
            ml_score_VBF_boosted_emu,
            ml_score_VBF_boosted_mumu,
        ]

        ml_score_booosted = SummedPlot(
            "DL_boosted_ml_score",
            boosted_ml_score_plots,
            title="ML",
            xTitle="DL Boosted ML score",
            plotopts=DL_boosted_label,
        )
        ml_score_resolved1b = SummedPlot(
            "DL_resolved1b_ml_score",
            resolved1b_ml_score_plots,
            title="ML",
            xTitle="DL Resolved 1b ML score",
            plotopts=DL_resolved1b_label,
        )
        ml_score_resolved2b = SummedPlot(
            "DL_resolved2b_ml_score",
            resolved2b_ml_score_plots,
            title="ML",
            xTitle="DL Resolved 2b ML score",
            plotopts=DL_resolved2b_label,
        )
        ml_score_VBF_resolved = SummedPlot(
            "DL_VBF_resolved_ml_score",
            [*VBF_resolved_ml_score_plots],
            title="DL ML score",
            # plotopts=DL_label,
        )
        ml_score_VBF_boosted = SummedPlot(
            "DL_VBF_boosted_ml_score",
            [*VBF_boosted_ml_score_plots],
            title="DL ML score",
            # plotopts=DL_label,
        )
        ml_score_DL = SummedPlot(
            "DL_ml_score",
            [
                *boosted_ml_score_plots,
                *resolved1b_ml_score_plots,
                *resolved2b_ml_score_plots,
                *VBF_resolved_ml_score_plots,
                *VBF_boosted_ml_score_plots,
            ],
            title="DL ML score",
            xTitle="ggF/HH Binary Score Distribution",
            # plotopts=DL_label,
        )
        plots.extend(
            [
                *resolved1b_ml_score_plots,
                *resolved2b_ml_score_plots,
                *boosted_ml_score_plots,
                *VBF_resolved_ml_score_plots,
                *VBF_boosted_ml_score_plots,
                ml_score_booosted,
                ml_score_resolved1b,
                ml_score_resolved2b,
                ml_score_VBF_resolved,
                ml_score_VBF_boosted,
                ml_score_DL,
            ]
        )

        return plots
