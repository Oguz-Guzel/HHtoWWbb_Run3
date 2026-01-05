
from bamboo.plots import Plot, SummedPlot, CategorizedSelection

from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection
from definitions import ml_input_features
from utils import ml_input_var_binning

class mvaInputPlotter(NanoBaseHHWWbb):
    """ Class to create MVA input variable distribution plots """

    def __init__(self, args):
        super().__init__(args)
        self.mvaModel = "/tmp/" # 

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # cutflow report
        plots.append(self.yields)

        # get DL selections
        [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
         DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
         DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu,
         DL_VBF_resolved_ee, DL_VBF_resolved_mumu, DL_VBF_resolved_emu,
         DL_VBF_boosted_ee, DL_VBF_boosted_mumu, DL_VBF_boosted_emu] = makeDLSelection(
            self, noSel, tree, sample)
        
        event_selections = [
            DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
            DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
            DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu,
            DL_VBF_resolved_ee, DL_VBF_resolved_mumu, DL_VBF_resolved_emu,
            DL_VBF_boosted_ee, DL_VBF_boosted_mumu, DL_VBF_boosted_emu
        ]

        # Following is the code to plot the input features for the ML model.
        # It takes some time to run, so it's commented out.
        ml_vars = ml_input_features(self)
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

        # # The following code is an example of how to print stuff to the terminal
        # # It is not necessary for the functionality of the code, but it can be useful for
        # # debugging or logging purposes.
        # # It uses the bamboo ROOT library to declare a function that prints the entry number,
        # # event number, and the ML output to the terminal.

        # from bamboo.root import gbl
        # from bamboo.analysisutils import addPrintout

        # gbl.gInterpreter.Declare("""
        #     bool bamboo_printEntry(long entry, long event, const std::vector<float>& ml_output) {
        #         std::cout << "Processing entry #" << entry << ": event " << event << " ml score [";
        #         for (size_t i = 0; i < ml_output.size(); ++i) {
        #             std::cout << ml_output[i];
        #             if (i < ml_output.size() - 1) std::cout << ", ";
        #         }
        #         std::cout << "]" << std::endl;
        #         return true;
        #     }
        # """)
        # addPrintout(DL_resolved_1b_ee, "bamboo_printEntry", op.extVar("ULong_t", "rdfentry_"), tree.event, ml_output)

        return plots
