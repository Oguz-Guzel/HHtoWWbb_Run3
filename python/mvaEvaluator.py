
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
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].p4.Px()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.Px()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Px(), self.tightMuons[0].p4.Px())),
            op.c_float(-9999.)
        )
        l2_Px = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.Px()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Px()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightMuons[0].p4.Px(), self.tightElectrons[0].p4.Px())),
            op.c_float(-9999.)
        )
        l1_Py = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].p4.Py()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.Py()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Py(), self.tightMuons[0].p4.Py())),
            op.c_float(-9999.)
        )
        l2_Py = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.Py()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Py()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightMuons[0].p4.Py(), self.tightElectrons[0].p4.Py())),
            op.c_float(-9999.)
        )
        l1_Pz = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].p4.Pz()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.Pz()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.Pz(), self.tightMuons[0].p4.Pz())),
            op.c_float(-9999.)
        )
        l2_Pz = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.Pz()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Pz()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightMuons[0].p4.Pz(), self.tightElectrons[0].p4.Pz())),
            op.c_float(-9999.)
        )
        l1_E = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].p4.E()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.E()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.E(), self.tightMuons[0].p4.E())),
            op.c_float(-9999.)
        )
        l2_E = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].p4.E()),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.E()),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].p4.E(), self.tightMuons[0].p4.E())),
            op.c_float(-9999.)
        )
        l1_pdgId = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].pdgId),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].pdgId),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].pdgId, self.tightMuons[0].pdgId)),
            op.c_int(-9999.)
        )
        l2_pdgId = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].pdgId),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].pdgId),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].pdgId, self.tightMuons[0].pdgId)),
            op.c_int(-9999.)
        )
        l1_charge = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[0].charge),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].charge),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].charge, self.tightMuons[0].charge)),
            op.c_int(-9999.)
        )
        l2_charge = op.multiSwitch(
            (op.rng_len(self.tightElectrons) ==
             2, self.tightElectrons[1].charge),
            (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].charge),
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt, self.tightElectrons[0].charge, self.tightMuons[0].charge)),
            op.c_int(-9999.)
        )

        j1_Px = op.switch(op.rng_len(self.ak4Jets) > 0,
                          self.ak4Jets[0].p4.Px(), op.c_float(-9999.))
        j1_Py = op.switch(op.rng_len(self.ak4Jets) > 0,
                          self.ak4Jets[0].p4.Py(), op.c_float(-9999.))
        j1_Pz = op.switch(op.rng_len(self.ak4Jets) > 0,
                          self.ak4Jets[0].p4.Pz(), op.c_float(-9999.))
        j1_E = op.switch(op.rng_len(self.ak4Jets) > 0,
                         self.ak4Jets[0].p4.E(), op.c_float(-9999.))
        j1_btag = op.switch(op.rng_len(self.ak4Jets) > 0,
                            self.ak4Jets[0].btagPNetB, op.c_float(-9999.))
        j2_Px = op.switch(op.rng_len(self.ak4Jets) > 1,
                          self.ak4Jets[1].p4.Px(), op.c_float(-9999.))
        j2_Py = op.switch(op.rng_len(self.ak4Jets) > 1,
                          self.ak4Jets[1].p4.Py(), op.c_float(-9999.))
        j2_Pz = op.switch(op.rng_len(self.ak4Jets) > 1,
                          self.ak4Jets[1].p4.Pz(), op.c_float(-9999.))
        j2_E = op.switch(op.rng_len(self.ak4Jets) > 1,
                         self.ak4Jets[1].p4.E(), op.c_float(-9999.))
        j2_btag = op.switch(op.rng_len(self.ak4Jets) > 1,
                            self.ak4Jets[1].btagPNetB, op.c_float(-9999.))
        j3_Px = op.switch(op.rng_len(self.ak4Jets) > 2,
                          self.ak4Jets[2].p4.Px(), op.c_float(-9999.))
        j3_Py = op.switch(op.rng_len(self.ak4Jets) > 2,
                          self.ak4Jets[2].p4.Py(), op.c_float(-9999.))
        j3_Pz = op.switch(op.rng_len(self.ak4Jets) > 2,
                          self.ak4Jets[2].p4.Pz(), op.c_float(-9999.))
        j3_E = op.switch(op.rng_len(self.ak4Jets) > 2,
                         self.ak4Jets[2].p4.E(), op.c_float(-9999.))
        j3_btag = op.switch(op.rng_len(self.ak4Jets) > 2,
                            self.ak4Jets[2].btagPNetB, op.c_float(-9999.))
        j4_Px = op.switch(op.rng_len(self.ak4Jets) > 3,
                          self.ak4Jets[3].p4.Px(), op.c_float(-9999.))
        j4_Py = op.switch(op.rng_len(self.ak4Jets) > 3,
                          self.ak4Jets[3].p4.Py(), op.c_float(-9999.))
        j4_Pz = op.switch(op.rng_len(self.ak4Jets) > 3,
                          self.ak4Jets[3].p4.Pz(), op.c_float(-9999.))
        j4_E = op.switch(op.rng_len(self.ak4Jets) > 3,
                         self.ak4Jets[3].p4.E(), op.c_float(-9999.))
        j4_btag = op.switch(op.rng_len(self.ak4Jets) > 3,
                            self.ak4Jets[3].btagPNetB, op.c_float(-9999.))

        met_Px = op.product(tree.MET.pt, op.cos(tree.MET.phi))
        met_Py = op.product(tree.MET.pt, op.sin(tree.MET.phi))
        met_E = tree.MET.pt

        mvaVars_DL = {
            "event_no": tree.event,
            "weight": noSel.weight,
            "l1_Px": l1_Px, "l1_Py": l1_Py, "l1_Pz": l1_Pz, "l1_E": l1_E, "l1_pdgId": l1_pdgId, "l1_charge": l1_charge,
            "l2_Px": l2_Px, "l2_Py": l2_Py, "l2_Pz": l2_Pz, "l2_E": l2_E, "l2_pdgId": l2_pdgId, "l2_charge": l2_charge,
            "j1_Px": j1_Px, "j1_Py": j1_Py, "j1_Pz": j1_Pz, "j1_E": j1_E, "j1_btag": j1_btag,
            "j2_Px": j2_Px, "j2_Py": j2_Py, "j2_Pz": j2_Pz, "j2_E": j2_E, "j2_btag": j2_btag,
            "j3_Px": j3_Px, "j3_Py": j3_Py, "j3_Pz": j3_Pz, "j3_E": j3_E, "j3_btag": j3_btag,
            "j4_Px": j4_Px, "j4_Py": j4_Py, "j4_Pz": j4_Pz, "j4_E": j4_E, "j4_btag": j4_btag,
            "met_Px": met_Px, "met_Py": met_Py, "met_E": met_E,
        }

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
            DL_DNN_InvM_cat1_label = {
                **labeler('DL DNN cat. 1 - blinded'), 'blinded-range': [0., 299.99]}
            DL_DNN_InvM_cat2_label = {
                **labeler('DL DNN cat. 2 - blinded'), 'blinded-range': [0., 299.99]}
            DL_DNN_InvM_cat3_label = {
                **labeler('DL DNN cat. 3 - blinded'), 'blinded-range': [0., 299.99]}
            DL_DNN_InvM_cat4_label = {
                **labeler('DL DNN cat. 4 - blinded'), 'blinded-range': [0., 299.99]}
            mvaVars_DL.pop("weight", None)

            split_var = 'even' if tree.event % 2 == 0 else 'odd'
            if split_var == 'odd':
                model = self.mvaModels + "/even_model.pth"
            elif split_var == 'even':
                model = self.mvaModels + "/odd_model.pth"
            else:
                print("Please provide a valid split variable !")

            dnn = op.mvaEvaluator(model, mvaType="Torch", nameHint="DL_DNN_nameHint")
            input_vars = [op.static_cast('float', v)
                          for v in mvaVars_DL.values()]
            DNN_inputs = op.array('float', *input_vars)
            DNN_output = dnn(DNN_inputs, defineOnFirstUse=False)

            # DNN cuts
            DL_resolved_1b_ee_DNNcat1 = DL_resolved_1b_ee.refine(
                "DL_resolved_1b_eeDNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_1b_ee_DNNcat2 = DL_resolved_1b_ee.refine(
                "DL_resolved_1b_eeDNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_1b_ee_DNNcat3 = DL_resolved_1b_ee.refine(
                "DL_resolved_1b_eeDNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_1b_ee_DNNcat4 = DL_resolved_1b_ee.refine(
                "DL_resolved_1b_eeDNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_resolved_2b_ee_DNNcat1 = DL_resolved_2b_ee.refine(
                "DL_resolved_2b_eeDNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_2b_ee_DNNcat2 = DL_resolved_2b_ee.refine(
                "DL_resolved_2b_eeDNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_2b_ee_DNNcat3 = DL_resolved_2b_ee.refine(
                "DL_resolved_2b_eeDNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_2b_ee_DNNcat4 = DL_resolved_2b_ee.refine(
                "DL_resolved_2b_eeDNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_boosted_ee_DNNcat1 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_boosted_ee_DNNcat2 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_boosted_ee_DNNcat3 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_boosted_ee_DNNcat4 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_resolved_1b_mumu_DNNcat1 = DL_resolved_1b_mumu.refine(
                "DL_resolved_1b_mumu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_1b_mumu_DNNcat2 = DL_resolved_1b_mumu.refine(
                "DL_resolved_1b_mumu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_1b_mumu_DNNcat3 = DL_resolved_1b_mumu.refine(
                "DL_resolved_1b_mumu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_1b_mumu_DNNcat4 = DL_resolved_1b_mumu.refine(
                "DL_resolved_1b_mumu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_resolved_2b_mumu_DNNcat1 = DL_resolved_2b_mumu.refine(
                "DL_resolved_2b_mumu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_2b_mumu_DNNcat2 = DL_resolved_2b_mumu.refine(
                "DL_resolved_2b_mumu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_2b_mumu_DNNcat3 = DL_resolved_2b_mumu.refine(
                "DL_resolved_2b_mumu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_2b_mumu_DNNcat4 = DL_resolved_2b_mumu.refine(
                "DL_resolved_2b_mumu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_boosted_mumu_DNNcat1 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_boosted_mumu_DNNcat2 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_boosted_mumu_DNNcat3 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_boosted_mumu_DNNcat4 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_resolved_1b_emu_DNNcat1 = DL_resolved_1b_emu.refine(
                "DL_resolved_1b_emu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_1b_emu_DNNcat2 = DL_resolved_1b_emu.refine(
                "DL_resolved_1b_emu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_1b_emu_DNNcat3 = DL_resolved_1b_emu.refine(
                "DL_resolved_1b_emu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_1b_emu_DNNcat4 = DL_resolved_1b_emu.refine(
                "DL_resolved_1b_emu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_resolved_2b_emu_DNNcat1 = DL_resolved_2b_emu.refine(
                "DL_resolved_2b_emu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_2b_emu_DNNcat2 = DL_resolved_2b_emu.refine(
                "DL_resolved_2b_emu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_2b_emu_DNNcat3 = DL_resolved_2b_emu.refine(
                "DL_resolved_2b_emu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_2b_emu_DNNcat4 = DL_resolved_2b_emu.refine(
                "DL_resolved_2b_emu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_boosted_emu_DNNcat1 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_boosted_emu_DNNcat2 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_boosted_emu_DNNcat3 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_boosted_emu_DNNcat4 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            # yields
            self.yields.add(DL_resolved_1b_ee_DNNcat4,
                            'DL_resolved_1b_ee_DNNcat4')
            self.yields.add(DL_resolved_2b_ee_DNNcat4,
                            'DL_resolved_2b_ee_DNNcat4')
            self.yields.add(DL_resolved_1b_mumu_DNNcat4,
                            'DL_resolved_1b_mumu_DNNcat4')
            self.yields.add(DL_resolved_2b_mumu_DNNcat4,
                            'DL_resolved_2b_mumu_DNNcat4')
            self.yields.add(DL_resolved_1b_emu_DNNcat4,
                            'DL_resolved_1b_emu_DNNcat4')
            self.yields.add(DL_resolved_2b_emu_DNNcat4,
                            'DL_resolved_2b_emu_DNNcat4')

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

            mElEl = op.invariant_mass(
                self.firstElTightPair[0].p4, self.firstElTightPair[1].p4)

            # DL ee DNN cat 1
            DL_resolved_1b_InvM_ee_DNNcat1 = Plot.make1D("DL_resolved_1b_InvM_ee_DNNcat1", mElEl, DL_resolved_1b_ee_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            DL_resolved_2b_InvM_ee_DNNcat1 = Plot.make1D("DL_resolved_2b_InvM_ee_DNNcat1", mElEl, DL_resolved_2b_ee_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            DL_boosted_InvM_ee_DNNcat1 = Plot.make1D("DL_boosted_InvM_ee_DNNcat1", mElEl, DL_boosted_ee_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            plots.extend(
                [DL_resolved_1b_InvM_ee_DNNcat1,
                 DL_resolved_2b_InvM_ee_DNNcat1,
                 DL_boosted_InvM_ee_DNNcat1,
                 SummedPlot("DL_InvM_ee_DNNcat1", [
                            DL_resolved_1b_InvM_ee_DNNcat1, DL_resolved_2b_InvM_ee_DNNcat1, DL_boosted_InvM_ee_DNNcat1], title="DL m(ee) DNN cat1")
                 ])

            # DL ee DNN cat 2
            DL_resolved_1b_InvM_ee_DNNcat2 = Plot.make1D("DL_resolved_1b_InvM_ee_DNNcat2", mElEl, DL_resolved_1b_ee_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            DL_resolved_2b_InvM_ee_DNNcat2 = Plot.make1D("DL_resolved_2b_InvM_ee_DNNcat2", mElEl, DL_resolved_2b_ee_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            DL_boosted_InvM_ee_DNNcat2 = Plot.make1D("DL_boosted_InvM_ee_DNNcat2", mElEl, DL_boosted_ee_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            plots.extend(
                [DL_resolved_1b_InvM_ee_DNNcat2,
                 DL_resolved_2b_InvM_ee_DNNcat2,
                 DL_boosted_InvM_ee_DNNcat2,
                 SummedPlot("DL_InvM_ee_DNNcat2", [
                            DL_resolved_1b_InvM_ee_DNNcat2, DL_resolved_2b_InvM_ee_DNNcat2, DL_boosted_InvM_ee_DNNcat2], title="DL m(ee) DNN cat2")
                 ])

            # DL ee DNN cat 3
            DL_resolved_1b_InvM_ee_DNNcat3 = Plot.make1D("DL_resolved_1b_InvM_ee_DNNcat3", mElEl, DL_resolved_1b_ee_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            DL_resolved_2b_InvM_ee_DNNcat3 = Plot.make1D("DL_resolved_2b_InvM_ee_DNNcat3", mElEl, DL_resolved_2b_ee_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            DL_boosted_InvM_ee_DNNcat3 = Plot.make1D("DL_boosted_InvM_ee_DNNcat3", mElEl, DL_boosted_ee_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            plots.extend(
                [DL_resolved_1b_InvM_ee_DNNcat3,
                 DL_resolved_2b_InvM_ee_DNNcat3,
                 DL_boosted_InvM_ee_DNNcat3,
                 SummedPlot("DL_InvM_ee_DNNcat3", [
                            DL_resolved_1b_InvM_ee_DNNcat3, DL_resolved_2b_InvM_ee_DNNcat3, DL_boosted_InvM_ee_DNNcat3], title="DL m(ee) DNN cat3")
                 ])

            # DL ee DNN cat 4
            DL_resolved_1b_InvM_ee_DNNcat4 = Plot.make1D("DL_resolved_1b_InvM_ee_DNNcat4", mElEl, DL_resolved_1b_ee_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            DL_resolved_2b_InvM_ee_DNNcat4 = Plot.make1D("DL_resolved_2b_InvM_ee_DNNcat4", mElEl, DL_resolved_2b_ee_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            DL_boosted_InvM_ee_DNNcat4 = Plot.make1D("DL_boosted_InvM_ee_DNNcat4", mElEl, DL_boosted_ee_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            plots.extend(
                [DL_resolved_1b_InvM_ee_DNNcat4,
                 DL_resolved_2b_InvM_ee_DNNcat4,
                 DL_boosted_InvM_ee_DNNcat4,
                 SummedPlot("DL_InvM_ee_DNNcat4", [
                            DL_resolved_1b_InvM_ee_DNNcat4, DL_resolved_2b_InvM_ee_DNNcat4, DL_boosted_InvM_ee_DNNcat4], title="DL m(ee) DNN cat4")
                 ])

            mMuMu = op.invariant_mass(
                self.firstMuTightPair[0].p4, self.firstMuTightPair[1].p4)

            # DL mumu DNN cat 1
            DL_resolved_1b_InvM_mumu_DNNcat1 = Plot.make1D("DL_resolved_1b_InvM_mumu_DNNcat1", mMuMu, DL_resolved_1b_mumu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            DL_resolved_2b_InvM_mumu_DNNcat1 = Plot.make1D("DL_resolved_2b_InvM_mumu_DNNcat1", mMuMu, DL_resolved_2b_mumu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            DL_boosted_InvM_mumu_DNNcat1 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat1", mMuMu, DL_boosted_mumu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            plots.extend(
                [DL_resolved_1b_InvM_mumu_DNNcat1,
                 DL_resolved_2b_InvM_mumu_DNNcat1,
                 DL_boosted_InvM_mumu_DNNcat1,
                 SummedPlot("DL_InvM_mumu_DNNcat1", [
                            DL_resolved_1b_InvM_mumu_DNNcat1, DL_resolved_2b_InvM_mumu_DNNcat1, DL_boosted_InvM_mumu_DNNcat1], title="DL m(mumu) DNN cat1")
                 ])

            # DL mumu DNN cat 2
            DL_resolved_1b_InvM_mumu_DNNcat2 = Plot.make1D("DL_resolved_1b_InvM_mumu_DNNcat2", mMuMu, DL_resolved_1b_mumu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            DL_resolved_2b_InvM_mumu_DNNcat2 = Plot.make1D("DL_resolved_2b_InvM_mumu_DNNcat2", mMuMu, DL_resolved_2b_mumu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            DL_boosted_InvM_mumu_DNNcat2 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat2", mMuMu, DL_boosted_mumu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            plots.extend(
                [DL_resolved_1b_InvM_mumu_DNNcat2,
                 DL_resolved_2b_InvM_mumu_DNNcat2,
                 DL_boosted_InvM_mumu_DNNcat2,
                 SummedPlot("DL_InvM_mumu_DNNcat2", [
                            DL_resolved_1b_InvM_mumu_DNNcat2, DL_resolved_2b_InvM_mumu_DNNcat2, DL_boosted_InvM_mumu_DNNcat2], title="DL m(mumu) DNN cat2")
                 ])

            # DL mumu DNN cat 3
            DL_resolved_1b_InvM_mumu_DNNcat3 = Plot.make1D("DL_resolved_1b_InvM_mumu_DNNcat3", mMuMu, DL_resolved_1b_mumu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            DL_resolved_2b_InvM_mumu_DNNcat3 = Plot.make1D("DL_resolved_2b_InvM_mumu_DNNcat3", mMuMu, DL_resolved_2b_mumu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            DL_boosted_InvM_mumu_DNNcat3 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat3", mMuMu, DL_boosted_mumu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            plots.extend(
                [DL_resolved_1b_InvM_mumu_DNNcat3,
                 DL_resolved_2b_InvM_mumu_DNNcat3,
                 DL_boosted_InvM_mumu_DNNcat3,
                 SummedPlot("DL_InvM_mumu_DNNcat3", [
                            DL_resolved_1b_InvM_mumu_DNNcat3, DL_resolved_2b_InvM_mumu_DNNcat3, DL_boosted_InvM_mumu_DNNcat3], title="DL m(mumu) DNN cat3")
                 ])

            # DL mumu DNN cat 4
            DL_resolved_1b_InvM_mumu_DNNcat4 = Plot.make1D("DL_resolved_1b_InvM_mumu_DNNcat4", mMuMu, DL_resolved_1b_mumu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            DL_resolved_2b_InvM_mumu_DNNcat4 = Plot.make1D("DL_resolved_2b_InvM_mumu_DNNcat4", mMuMu, DL_resolved_2b_mumu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            DL_boosted_InvM_mumu_DNNcat4 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat4", mMuMu, DL_boosted_mumu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            plots.extend(
                [DL_resolved_1b_InvM_mumu_DNNcat4,
                 DL_resolved_2b_InvM_mumu_DNNcat4,
                 DL_boosted_InvM_mumu_DNNcat4,
                 SummedPlot("DL_InvM_mumu_DNNcat4", [
                            DL_resolved_1b_InvM_mumu_DNNcat4, DL_resolved_2b_InvM_mumu_DNNcat4, DL_boosted_InvM_mumu_DNNcat4], title="DL m(mumu) DNN cat4")
                 ])

            mElMu = op.invariant_mass(
                self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4
            )

            # DL emu DNN cat 1
            DL_resolved_1b_InvM_emu_DNNcat1 = Plot.make1D("DL_resolved_1b_InvM_emu_DNNcat1", mElMu, DL_resolved_1b_emu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            DL_resolved_2b_InvM_emu_DNNcat1 = Plot.make1D("DL_resolved_2b_InvM_emu_DNNcat1", mElMu, DL_resolved_2b_emu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            DL_boosted_InvM_emu_DNNcat1 = Plot.make1D("DL_boosted_InvM_emu_DNNcat1", mElMu, DL_boosted_emu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat1_label)
            plots.extend(
                [DL_resolved_1b_InvM_emu_DNNcat1,
                 DL_resolved_1b_InvM_emu_DNNcat1,
                 DL_boosted_InvM_emu_DNNcat1,
                 SummedPlot("DL_InvM_emu_DNNcat1", [
                            DL_resolved_1b_InvM_emu_DNNcat1, DL_resolved_2b_InvM_emu_DNNcat1, DL_boosted_InvM_emu_DNNcat1], title="DL m(elmu) DNN cat1")
                 ])

            # DL emu DNN cat 2
            DL_resolved_1b_InvM_emu_DNNcat2 = Plot.make1D("DL_resolved_1b_InvM_emu_DNNcat2", mElMu, DL_resolved_1b_emu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            DL_resolved_2b_InvM_emu_DNNcat2 = Plot.make1D("DL_resolved_2b_InvM_emu_DNNcat2", mElMu, DL_resolved_2b_emu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            DL_boosted_InvM_emu_DNNcat2 = Plot.make1D("DL_boosted_InvM_emu_DNNcat2", mElMu, DL_boosted_emu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat2_label)
            plots.extend(
                [DL_resolved_1b_InvM_emu_DNNcat2,
                 DL_resolved_1b_InvM_emu_DNNcat2,
                 DL_boosted_InvM_emu_DNNcat2,
                 SummedPlot("DL_InvM_emu_DNNcat2", [
                            DL_resolved_1b_InvM_emu_DNNcat2, DL_resolved_2b_InvM_emu_DNNcat2, DL_boosted_InvM_emu_DNNcat2], title="DL m(elmu) DNN cat2")
                 ])

            # DL emu DNN cat 3
            DL_resolved_1b_InvM_emu_DNNcat3 = Plot.make1D("DL_resolved_1b_InvM_emu_DNNcat3", mElMu, DL_resolved_1b_emu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            DL_resolved_2b_InvM_emu_DNNcat3 = Plot.make1D("DL_resolved_2b_InvM_emu_DNNcat3", mElMu, DL_resolved_2b_emu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            DL_boosted_InvM_emu_DNNcat3 = Plot.make1D("DL_boosted_InvM_emu_DNNcat3", mElMu, DL_boosted_emu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat3_label)
            plots.extend(
                [DL_resolved_1b_InvM_emu_DNNcat3,
                 DL_resolved_1b_InvM_emu_DNNcat3,
                 DL_boosted_InvM_emu_DNNcat3,
                 SummedPlot("DL_InvM_emu_DNNcat3", [
                            DL_resolved_1b_InvM_emu_DNNcat3, DL_resolved_2b_InvM_emu_DNNcat3, DL_boosted_InvM_emu_DNNcat3], title="DL m(elmu) DNN cat3")
                 ])

            # DL emu DNN cat 4
            DL_resolved_1b_InvM_emu_DNNcat4 = Plot.make1D("DL_resolved_1b_InvM_emu_DNNcat4", mElMu, DL_resolved_1b_emu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            DL_resolved_2b_InvM_emu_DNNcat4 = Plot.make1D("DL_resolved_2b_InvM_emu_DNNcat4", mElMu, DL_resolved_2b_emu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            DL_boosted_InvM_emu_DNNcat4 = Plot.make1D("DL_boosted_InvM_emu_DNNcat4", mElMu, DL_boosted_emu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_cat4_label)
            plots.extend(
                [DL_resolved_1b_InvM_emu_DNNcat4,
                 DL_resolved_1b_InvM_emu_DNNcat4,
                 DL_boosted_InvM_emu_DNNcat4,
                 SummedPlot("DL_InvM_emu_DNNcat4", [
                            DL_resolved_1b_InvM_emu_DNNcat4, DL_resolved_2b_InvM_emu_DNNcat4, DL_boosted_InvM_emu_DNNcat4], title="DL m(emu) DNN cat4")
                 ])

        return plots
