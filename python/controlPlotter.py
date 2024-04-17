
from bamboo.plots import Plot, Skim, SummedPlot
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection, makeSLSelection
from scalefactors import ScaleFactors as sf
import definitions as defs
from utils import labeler


class controlPlotter(NanoBaseHHWWbb):
    """ Class to create control plots, cutflow reports and skims"""

    def __init__(self, args):
        super(controlPlotter, self).__init__(args)
        self.channel = self.args.channel
        self.mvaModels = self.args.mvaModels
        self.sync = self.args.sync

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
            DL_boosted_ee, DL_boosted_mumu, \
                DL_boosted_emu, DL_resolved_ee, \
                DL_resolved_mumu, DL_resolved_emu = makeDLSelection(
                    self, noSel)
            
            if self.sync:
                DL_1_preElectron = noSel.refine("DL_1_preElectron", cut=op.rng_len(self.preElectrons) == 1)
                self.yields.add(DL_1_preElectron, '1 pre-Electron')
                
                DL_1_preMuon = noSel.refine("DL_1_preMuon", cut=op.rng_len(self.preMuons) == 1)
                self.yields.add(DL_1_preMuon, '1 pre-Muon')

            # muonSF
            DL_boosted_mumu = sf.muonSF(self, DL_boosted_mumu)
            DL_resolved_mumu = sf.muonSF(self, DL_resolved_mumu)
            DL_boosted_emu = sf.muonSF(self, DL_boosted_emu)
            DL_resolved_emu = sf.muonSF(self, DL_resolved_emu)

            # electronSF
            DL_boosted_ee = sf.electronSF(self, DL_boosted_ee)
            DL_resolved_ee = sf.electronSF(self, DL_resolved_ee)
            DL_boosted_emu = sf.electronSF(self, DL_boosted_emu)
            DL_resolved_emu = sf.electronSF(self, DL_resolved_emu)

            # btagging SF, only for resolved since the btagSF function takes only ak4 jets for now, will update once ak8 corrections are available
            DL_resolved_ee = sf.btagSF(self, DL_resolved_ee)
            DL_resolved_mumu = sf.btagSF(self, DL_resolved_mumu)
            DL_resolved_emu = sf.btagSF(self, DL_resolved_emu)

            # cutflow report for DL channel
            self.yields.add(DL_boosted_ee, 'DL boosted ee')
            self.yields.add(DL_boosted_mumu, 'DL boosted mumu')
            self.yields.add(DL_boosted_emu, 'DL boosted emu')
            self.yields.add(DL_resolved_ee, 'DL resolved ee')
            self.yields.add(DL_resolved_mumu, 'DL resolved mumu')
            self.yields.add(DL_resolved_emu, 'DL resolved emu')

            # labels on plots
            DLboostedEE_label = labeler('DL boosted EE')
            DLboostedMuMu_label = labeler('DL boosted MuMu')
            DLboostedEMU_label = labeler('DL boosted EMu')
            DLresolvedEE_label = labeler('DL resolved EE')
            DLresolvedMuMu_label = labeler('DL resolved MuMu')
            DLresolvedEMu_label = labeler('DL resolved EMu')

        if self.channel == 'SL':
            # get SL selections
            SL_resolved_pre, SL_resolved_1b_e, SL_resolved_2b_e, \
                SL_resolved_1b_mu, SL_resolved_2b_mu, SL_boosted, \
                SL_boosted_e, SL_boosted_mu, SL_e = makeSLSelection(
                    self, noSel)

            # cutflow report for SL channel
            self.yields.add(SL_boosted_e, 'SL boosted e')
            self.yields.add(SL_boosted_mu, 'SL boosted mu')
            self.yields.add(SL_resolved_1b_e, 'SL resolved e')
            self.yields.add(SL_resolved_1b_mu, 'SL resolved mu')

            # labels on plots
            SLboostedE_label = labeler('SL boosted E')
            SLboostedMu_label = labeler('SL boosted Mu')
            SLresolvedE_label = labeler('SL resolved E')
            SLresolvedMu_label = labeler('SL resolved Mu')

        # mva variables
        lepton0conept = op.multiSwitch(
            # if there are 2 electrons
            (op.rng_len(self.tightElectrons) == 2,
             self.electron_conept[self.tightElectrons[0].idx]),
            # elif there are 2 muons
            (op.rng_len(self.tightMuons) == 2,
             self.muon_conept[self.tightMuons[0].idx]),
            # elif there is 1 electron and 1 muon
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1),
             op.switch(
                # if electron pt is greater than muon pt
                self.tightElectrons[0].pt >= self.tightMuons[0].pt,
                self.electron_conept[
                    self.tightElectrons[0].idx],
                # else
                self.muon_conept[self.tightMuons[0].idx])),
            # else
            op.c_float(0.)
        )

        lepton1conept = op.multiSwitch(
            # if there are 2 electrons
            (op.rng_len(self.tightElectrons) == 2,
             self.electron_conept[self.tightElectrons[1].idx]),
            # elif there are 2 muons
            (op.rng_len(self.tightMuons) == 2,
             self.muon_conept[self.tightMuons[0].idx]),
            # elif there is 1 electron and 1 muon
            (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                # if electron pt is greater than muon pt
                self.tightElectrons[0].pt >= self.tightMuons[0].pt, self.muon_conept[self.tightMuons[0].idx],
                # else
                self.electron_conept[self.tightElectrons[1].idx])),
            op.c_float(0.)
        )

        invmLL = op.multiSwitch(
            # if there are 2 electrons
            (op.rng_len(self.tightElectrons) == 2,
             op.invariant_mass(self.tightElectrons[0].p4, self.tightElectrons[1].p4)),
            # elif there are 2 muons
            (op.rng_len(self.tightMuons) == 2,
             op.invariant_mass(self.tightMuons[0].p4, self.tightMuons[1].p4)),
            # elif there is 1 electron and 1 muon
            (op.AND(op.rng_len(self.tightMuons) == 1, op.rng_len(self.tightElectrons) == 1),
             op.invariant_mass(
                self.tightElectrons[0].p4, self.tightMuons[0].p4)),
            op.c_float(0.)
        )

        drLL = op.multiSwitch(
            # if there are 2 electrons
            (op.rng_len(self.tightElectrons) == 2,
             op.deltaR(self.tightElectrons[0].p4, self.tightElectrons[1].p4)),
            # elif there are 2 muons
            (op.rng_len(self.tightMuons) == 2,
             op.deltaR(self.tightMuons[0].p4, self.tightMuons[1].p4)),
            # elif there is 1 electron and 1 muon
            (op.AND(op.rng_len(self.tightMuons) == 1, op.rng_len(self.tightElectrons) == 1),
             op.deltaR(
                self.tightElectrons[0].p4, self.tightMuons[0].p4)),
            op.c_float(0.)
        )
        DRTwoak4bjets = op.multiSwitch(
            # if there are 2 b-tagged ak4 jets
            (op.rng_len(self.ak4BJets) == 2,
             op.deltaR(self.ak4BJets[0].p4, self.ak4BJets[1].p4)),
            op.c_float(0.)
        )
        inmvAK4bJJ = op.multiSwitch(
            # if there are 2 b-tagged ak4 jets
            (op.rng_len(self.ak4BJets) == 2,
             op.invariant_mass(op.sum(self.ak4BJets[0].p4, self.ak4BJets[1].p4))),  # double check
            op.c_float(0.)
        )

        mvaVars_DL = {
            "weight": noSel.weight,
            "DRll": drLL,
            "DRtwoAK4bjets": DRTwoak4bjets,
            "InvMak4bJJ": inmvAK4bJJ,
            "InvMll": invmLL,
            "MET": tree.MET.pt,
            "ak4jet0eta": op.switch(op.rng_len(self.ak4Jets) > 0, self.ak4Jets[0].eta, op.c_float(0.)),
            "ak4jet0pt": op.switch(op.rng_len(self.ak4Jets) > 0, self.ak4Jets[0].pt, op.c_float(0.)),
            "ak4jet1pt": op.switch(op.rng_len(self.ak4Jets) > 1, self.ak4Jets[1].pt, op.c_float(0.)),
            "ak4jet1eta": op.switch(op.rng_len(self.ak4Jets) > 1, self.ak4Jets[1].eta, op.c_float(0.)),
            "ak8bjetpt": op.switch(op.rng_len(self.ak8BJets) > 0, self.ak8BJets[0].pt, op.c_float(0.)),
            "ak8bjeteta": op.switch(op.rng_len(self.ak8BJets) > 0, self.ak8BJets[0].eta, op.c_float(0.)),
            "ak8jetpt": op.switch(op.rng_len(self.ak8Jets) > 0, self.ak8Jets[0].pt, op.c_float(0.)),
            "ak8jeteta": op.switch(op.rng_len(self.ak8Jets) > 0, self.ak8Jets[0].eta, op.c_float(0.)),
            "lepton0Cone_pt": lepton0conept,
            "lepton1Conept": lepton1conept,
            "nAK4": op.rng_len(self.ak4Jets),
            "nAK4bJet": op.rng_len(self.ak4BJets),
            "nak8": op.rng_len(self.ak8Jets),
            "nAK8bJet": op.rng_len(self.ak8BJets)
        }

        ### Syncronisaton ###
        if self.args.sync and self.channel == 'DL' and not self.mvaModels:
            syncVars_DL = {
                "event_no": tree.event,
                "lumi_block": tree.luminosityBlock,
                "is_dl_ee": op.switch(op.rng_len(self.tightElectrons) == 2, 1, op.c_float(0.)),
                "is_dl_mumu": op.switch(op.rng_len(self.tightMuons) == 2, 1, op.c_float(0.)),
                "is_dl_emu": op.switch(op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons)) == 1, 1, op.c_float(0.)),
                "nLooseElectron": op.rng_len(self.preElectrons),
                "nFakeElectron": op.rng_len(self.fakeElectrons),
                "nTightElectron": op.rng_len(self.tightElectrons),
                "nLooseMuon": op.rng_len(self.preMuons),
                "nFakeMuon": op.rng_len(self.fakeMuons),
                "nTightMuon": op.rng_len(self.tightMuons),
                "lepton0pt": op.multiSwitch(
                    (op.rng_len(self.tightElectrons) == 2, self.tightElectrons[0].pt),
                    (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].pt),
                    (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                        self.tightElectrons[0].pt >= self.tightMuons[0].pt, self.tightElectrons[0].pt, self.tightMuons[0].pt)),
                 op.c_float(0.)
                ),
                "lepton1pt": op.multiSwitch(
                    (op.rng_len(self.tightElectrons) == 2, self.tightElectrons[1].pt),
                    (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].pt),
                    (op.AND(op.rng_len(self.tightElectrons) == 1, op.rng_len(self.tightMuons) == 1), op.switch(
                        self.tightElectrons[0].pt >= self.tightMuons[0].pt, self.tightMuons[0].pt, self.tightElectrons[0].pt)),
                op.c_float(0.)
                ),
                "nAK4": op.rng_len(self.ak4Jets),
                "nAK4bJet": op.rng_len(self.ak4BJets),
                "nAK8bJet": op.rng_len(self.ak8BJets),
                "ak4jet0pt": op.switch(op.rng_len(self.ak4Jets) > 0, self.ak4Jets[0].pt, op.c_float(0.)),
                "ak4jet0eta": op.switch(op.rng_len(self.ak4Jets) > 0, self.ak4Jets[0].eta, op.c_float(0.)),
                "ak4jet1pt": op.switch(op.rng_len(self.ak4Jets) > 1, self.ak4Jets[1].pt, op.c_float(0.)),
                "ak4jet1eta": op.switch(op.rng_len(self.ak4Jets) > 1, self.ak4Jets[1].eta, op.c_float(0.)),
                "ak8jetpt": op.switch(op.rng_len(self.ak8Jets) > 0, self.ak8Jets[0].pt, op.c_float(0.)),
                "ak8jeteta": op.switch(op.rng_len(self.ak8Jets) > 0, self.ak8Jets[0].eta, op.c_float(0.)),
            }
            # to order the columns
            self.order = [key for key in syncVars_DL.keys()]
            plots.extend([
                Skim("DL_resolved_ee_sync", syncVars_DL, DL_resolved_ee),
                Skim("DL_resolved_mumu_sync", syncVars_DL, DL_resolved_mumu),
                Skim("DL_resolved_emu_sync", syncVars_DL, DL_resolved_emu),
                Skim("DL_boosted_ee_sync", syncVars_DL, DL_boosted_ee),
                Skim("DL_boosted_mumu_sync", syncVars_DL, DL_boosted_mumu),
                Skim("DL_boosted_emu_sync", syncVars_DL, DL_boosted_emu),
                Plot.make1D("DL_boosted_nfatJet_ee", op.rng_len(self.ak8Jets), DL_boosted_ee, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedEE_label)
            ])
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
            DL_DNN_InvM_Cat1_label = {
                **labeler('DL DNN cat. 1 - blinded'), 'blinded-range': [0., 299.99]}
            DL_DNN_InvM_Cat2_label = {
                **labeler('DL DNN cat. 2 - blinded'), 'blinded-range': [0., 299.99]}
            DL_DNN_InvM_Cat3_label = {
                **labeler('DL DNN cat. 3 - blinded'), 'blinded-range': [0., 299.99]}
            DL_DNN_InvM_Cat4_label = {
                **labeler('DL DNN cat. 4 - blinded'), 'blinded-range': [0., 299.99]}
            mvaVars_DL.pop("weight", None)

            import random
            split_var = random.randint(0, 1)
            if split_var == 0:
                model = self.mvaModels + "/model_even/model.onnx"
            elif split_var == 1:
                model = self.mvaModels + "/model_odd/model.onnx"
            else:
                print("ERROR: split_var is not 0 or 1")
                exit(1)

            dnn = op.mvaEvaluator(model, otherArgs=("predictions"))
            input_vars = [op.static_cast('float', v)
                          for v in mvaVars_DL.values()]
            DNN_inputs = op.array('float', *input_vars)
            DNN_output = dnn(DNN_inputs)

            # DNN cuts
            DL_resolved_ee_DNNcat1 = DL_resolved_ee.refine(
                "DL_resolved_eeDNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_ee_DNNcat2 = DL_resolved_ee.refine(
                "DL_resolved_eeDNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_ee_DNNcat3 = DL_resolved_ee.refine(
                "DL_resolved_eeDNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_ee_DNNcat4 = DL_resolved_ee.refine(
                "DL_resolved_eeDNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_boosted_ee_DNNcat1 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_boosted_ee_DNNcat2 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_boosted_ee_DNNcat3 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_boosted_ee_DNNcat4 = DL_boosted_ee.refine(
                "DL_boosted_ee_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_resolved_mumu_DNNcat1 = DL_resolved_mumu.refine(
                "DL_resolved_mumu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_mumu_DNNcat2 = DL_resolved_mumu.refine(
                "DL_resolved_mumu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_mumu_DNNcat3 = DL_resolved_mumu.refine(
                "DL_resolved_mumu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_mumu_DNNcat4 = DL_resolved_mumu.refine(
                "DL_resolved_mumu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_boosted_mumu_DNNcat1 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_boosted_mumu_DNNcat2 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_boosted_mumu_DNNcat3 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_boosted_mumu_DNNcat4 = DL_boosted_mumu.refine(
                "DL_boosted_mumu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_resolved_elmu_DNNcat1 = DL_resolved_emu.refine(
                "DL_resolved_emu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_resolved_elmu_DNNcat2 = DL_resolved_emu.refine(
                "DL_resolved_emu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_resolved_elmu_DNNcat3 = DL_resolved_emu.refine(
                "DL_resolved_emu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_resolved_elmu_DNNcat4 = DL_resolved_emu.refine(
                "DL_resolved_emu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            DL_boosted_elmu_DNNcat1 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat1", cut=op.in_range(0.1, DNN_output[0], 0.6))
            DL_boosted_elmu_DNNcat2 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat2", cut=op.in_range(0.6, DNN_output[0], 0.8))
            DL_boosted_elmu_DNNcat3 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat3", cut=op.in_range(0.8, DNN_output[0], 0.92))
            DL_boosted_elmu_DNNcat4 = DL_boosted_emu.refine(
                "DL_boosted_emu_DNNcat4", cut=op.in_range(0.92, DNN_output[0], 1.0))

            # to make sure yields are working fine
            self.yields.add(DL_resolved_ee_DNNcat4, 'DL_resolved_ee_DNNcat4')

            dnn_score_ee = Plot.make1D("dnn_score_ee", DNN_output[0], DL_resolved_ee, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EE
            )
            dnn_score_emu = Plot.make1D("dnn_score_emu", DNN_output[0], DL_resolved_emu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_EMu
            )
            dnn_score_mumu = Plot.make1D("dnn_score_mumu", DNN_output[0], DL_resolved_mumu, EqBin(
                100, 0, 1.), title='DNN', xTitle="DNN Score", plotopts=DL_DNN_MuMu
            )
            plots.extend([
                dnn_score_ee,
                dnn_score_emu,
                dnn_score_mumu,
                SummedPlot("DL_dnn_score", [
                           dnn_score_ee, dnn_score_emu, dnn_score_mumu], title="DL DNN score", plotopts=DL_DNN
                           )
            ])

            mElEl = op.invariant_mass(
                self.firstOSElEl[0].p4, self.firstOSElEl[1].p4)

            # DL ee DNN cat 1
            DL_resolved_InvM_ee_DNNcat1 = Plot.make1D("DL_resolved_InvM_ee_DNNcat1", mElEl, DL_resolved_ee_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat1_label)
            DL_boosted_InvM_ee_DNNcat1 = Plot.make1D("DL_boosted_InvM_ee_DNNcat1", mElEl, DL_boosted_ee_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat1_label)
            plots.extend(
                [DL_resolved_InvM_ee_DNNcat1,
                 DL_boosted_InvM_ee_DNNcat1,
                 SummedPlot("DL_InvM_ee_DNNcat1", [
                            DL_resolved_InvM_ee_DNNcat1, DL_boosted_InvM_ee_DNNcat1], title="DL m(ee) DNN cat1")
                 ])

            # DL ee DNN cat 2
            DL_resolved_InvM_ee_DNNcat2 = Plot.make1D("DL_resolved_InvM_ee_DNNcat2", mElEl, DL_resolved_ee_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat2_label)
            DL_boosted_InvM_ee_DNNcat2 = Plot.make1D("DL_boosted_InvM_ee_DNNcat2", mElEl, DL_boosted_ee_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat2_label)
            plots.extend(
                [DL_resolved_InvM_ee_DNNcat2,
                 DL_boosted_InvM_ee_DNNcat2,
                 SummedPlot("DL_InvM_ee_DNNcat2", [
                            DL_resolved_InvM_ee_DNNcat2, DL_boosted_InvM_ee_DNNcat2], title="DL m(ee) DNN cat2")
                 ])

            # DL ee DNN cat 3
            DL_resolved_InvM_ee_DNNcat3 = Plot.make1D("DL_resolved_InvM_ee_DNNcat3", mElEl, DL_resolved_ee_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat3_label)
            DL_boosted_InvM_ee_DNNcat3 = Plot.make1D("DL_boosted_InvM_ee_DNNcat3", mElEl, DL_boosted_ee_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat3_label)
            plots.extend(
                [DL_resolved_InvM_ee_DNNcat3,
                 DL_boosted_InvM_ee_DNNcat3,
                 SummedPlot("DL_InvM_ee_DNNcat3", [
                            DL_resolved_InvM_ee_DNNcat3, DL_boosted_InvM_ee_DNNcat3], title="DL m(ee) DNN cat3")
                 ])

            # DL ee DNN cat 4
            DL_resolved_InvM_ee_DNNcat4 = Plot.make1D("DL_resolved_InvM_ee_DNNcat4", mElEl, DL_resolved_ee_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat4_label)
            DL_boosted_InvM_ee_DNNcat4 = Plot.make1D("DL_boosted_InvM_ee_DNNcat4", mElEl, DL_boosted_ee_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat4_label)
            plots.extend(
                [DL_resolved_InvM_ee_DNNcat4,
                 DL_boosted_InvM_ee_DNNcat4,
                 SummedPlot("DL_InvM_ee_DNNcat4", [
                            DL_resolved_InvM_ee_DNNcat4, DL_boosted_InvM_ee_DNNcat4], title="DL m(ee) DNN cat4")
                 ])

            mMuMu = op.invariant_mass(
                self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4)

            # DL mumu DNN cat 1
            DL_resolved_InvM_mumu_DNNcat1 = Plot.make1D("DL_resolved_InvM_mumu_DNNcat1", mMuMu, DL_resolved_mumu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat1_label)
            DL_boosted_InvM_mumu_DNNcat1 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat1", mMuMu, DL_boosted_mumu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat1_label)
            plots.extend(
                [DL_resolved_InvM_mumu_DNNcat1,
                 DL_boosted_InvM_mumu_DNNcat1,
                 SummedPlot("DL_InvM_mumu_DNNcat1", [
                            DL_resolved_InvM_mumu_DNNcat1, DL_boosted_InvM_mumu_DNNcat1], title="DL m(mumu) DNN cat1")
                 ])

            # DL mumu DNN cat 2
            DL_resolved_InvM_mumu_DNNcat2 = Plot.make1D("DL_resolved_InvM_mumu_DNNcat2", mMuMu, DL_resolved_mumu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat2_label)
            DL_boosted_InvM_mumu_DNNcat2 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat2", mMuMu, DL_boosted_mumu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat2_label)
            plots.extend(
                [DL_resolved_InvM_mumu_DNNcat2,
                 DL_boosted_InvM_mumu_DNNcat2,
                 SummedPlot("DL_InvM_mumu_DNNcat2", [
                            DL_resolved_InvM_mumu_DNNcat2, DL_boosted_InvM_mumu_DNNcat2], title="DL m(mumu) DNN cat2")
                 ])

            # DL mumu DNN cat 3
            DL_resolved_InvM_mumu_DNNcat3 = Plot.make1D("DL_resolved_InvM_mumu_DNNcat3", mMuMu, DL_resolved_mumu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat3_label)
            DL_boosted_InvM_mumu_DNNcat3 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat3", mMuMu, DL_boosted_mumu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat3_label)
            plots.extend(
                [DL_resolved_InvM_mumu_DNNcat3,
                 DL_boosted_InvM_mumu_DNNcat3,
                 SummedPlot("DL_InvM_mumu_DNNcat3", [
                            DL_resolved_InvM_mumu_DNNcat3, DL_boosted_InvM_mumu_DNNcat3], title="DL m(mumu) DNN cat3")
                 ])

            # DL mumu DNN cat 4
            DL_resolved_InvM_mumu_DNNcat4 = Plot.make1D("DL_resolved_InvM_mumu_DNNcat4", mMuMu, DL_resolved_mumu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat4_label)
            DL_boosted_InvM_mumu_DNNcat4 = Plot.make1D("DL_boosted_InvM_mumu_DNNcat4", mMuMu, DL_boosted_mumu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat4_label)
            plots.extend(
                [DL_resolved_InvM_mumu_DNNcat4,
                 DL_boosted_InvM_mumu_DNNcat4,
                 SummedPlot("DL_InvM_mumu_DNNcat4", [
                            DL_resolved_InvM_mumu_DNNcat4, DL_boosted_InvM_mumu_DNNcat4], title="DL m(mumu) DNN cat4")
                 ])

            mElMu = op.invariant_mass(
                self.firstOSElMu[0].p4, self.firstOSElMu[1].p4
            )

            # DL elmu DNN cat 1
            DL_resolved_InvM_elmu_DNNcat1 = Plot.make1D("DL_resolved_InvM_elmu_DNNcat1", mElMu, DL_resolved_elmu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat1_label)
            DL_boosted_InvM_elmu_DNNcat1 = Plot.make1D("DL_boosted_InvM_elmu_DNNcat1", mElMu, DL_boosted_elmu_DNNcat1, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat1_label)
            plots.extend(
                [DL_resolved_InvM_elmu_DNNcat1,
                 DL_boosted_InvM_elmu_DNNcat1,
                 SummedPlot("DL_InvM_elmu_DNNcat1", [
                            DL_resolved_InvM_elmu_DNNcat1, DL_boosted_InvM_elmu_DNNcat1], title="DL m(elmu) DNN cat1")
                 ])

            # DL elmu DNN cat 2
            DL_resolved_InvM_elmu_DNNcat2 = Plot.make1D("DL_resolved_InvM_elmu_DNNcat2", mElMu, DL_resolved_elmu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat2_label)
            DL_boosted_InvM_elmu_DNNcat2 = Plot.make1D("DL_boosted_InvM_elmu_DNNcat2", mElMu, DL_boosted_elmu_DNNcat2, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat2_label)
            plots.extend(
                [DL_resolved_InvM_elmu_DNNcat2,
                 DL_boosted_InvM_elmu_DNNcat2,
                 SummedPlot("DL_InvM_elmu_DNNcat2", [
                            DL_resolved_InvM_elmu_DNNcat2, DL_boosted_InvM_elmu_DNNcat2], title="DL m(elmu) DNN cat2")
                 ])

            # DL elmu DNN cat 3
            DL_resolved_InvM_elmu_DNNcat3 = Plot.make1D("DL_resolved_InvM_elmu_DNNcat3", mElMu, DL_resolved_elmu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat3_label)
            DL_boosted_InvM_elmu_DNNcat3 = Plot.make1D("DL_boosted_InvM_elmu_DNNcat3", mElMu, DL_boosted_elmu_DNNcat3, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat3_label)
            plots.extend(
                [DL_resolved_InvM_elmu_DNNcat3,
                 DL_boosted_InvM_elmu_DNNcat3,
                 SummedPlot("DL_InvM_elmu_DNNcat3", [
                            DL_resolved_InvM_elmu_DNNcat3, DL_boosted_InvM_elmu_DNNcat3], title="DL m(elmu) DNN cat3")
                 ])

            # DL elmu DNN cat 4
            DL_resolved_InvM_elmu_DNNcat4 = Plot.make1D("DL_resolved_InvM_elmu_DNNcat4", mElMu, DL_resolved_elmu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat4_label)
            DL_boosted_InvM_elmu_DNNcat4 = Plot.make1D("DL_boosted_InvM_elmu_DNNcat4", mElMu, DL_boosted_elmu_DNNcat4, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of el-mu pair (GeV/c^{2})", plotopts=DL_DNN_InvM_Cat4_label)
            plots.extend(
                [DL_resolved_InvM_elmu_DNNcat4,
                 DL_boosted_InvM_elmu_DNNcat4,
                 SummedPlot("DL_InvM_elmu_DNNcat4", [
                            DL_resolved_InvM_elmu_DNNcat4, DL_boosted_InvM_elmu_DNNcat4], title="DL m(elmu) DNN cat4")
                 ])

        #############################################################################
        #                                 Plots                                     #
        #############################################################################

        if self.channel == 'DL' and not self.mvaModels and not self.sync:
            plots.extend([
                #########################################
                #                 Skims                 #
                #########################################

                Skim("DL_resolved_ee_mva", mvaVars_DL, DL_resolved_ee),
                Skim("DL_resolved_mumu_mva", mvaVars_DL, DL_resolved_mumu),
                Skim("DL_resolved_emu_mva", mvaVars_DL, DL_resolved_emu),
                Skim("DL_boosted_ee_mva", mvaVars_DL, DL_boosted_ee),
                Skim("DL_boosted_mumu_mva", mvaVars_DL, DL_boosted_mumu),
                Skim("DL_boosted_emu_mva", mvaVars_DL, DL_boosted_emu),

                #########################################
                ######                             ######
                ######       DL boosted plots      ######
                ######                             ######
                #########################################

                # number of ak8 b-jets
                Plot.make1D("DL_boosted_nfatJet_ee", op.rng_len(self.ak8Jets), DL_boosted_ee, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nfatJet_mumu", op.rng_len(self.ak8Jets), DL_boosted_mumu, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nfatJet_emu", op.rng_len(self.ak8Jets), DL_boosted_emu, EqBin(
                    10, 0, 10), title="N(ak8jet)", xTitle="Number of fatjet", plotopts=DLboostedEMU_label),

                # fatjet pt
                Plot.make1D("DL_boosted_fatJet_pt_ee", self.ak8Jets[0].pt, DL_boosted_ee, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_fatJet_pt_mumu", self.ak8Jets[0].pt, DL_boosted_mumu, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_fatJet_pt_emu", self.ak8Jets[0].pt, DL_boosted_emu, EqBin(
                    100, 200, 800), title="pT(ak8jet)", xTitle="Fatjet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # subjet1 pt
                Plot.make1D("DL_boosted_subjet1_pt_ee", self.ak8Jets[0].subJet1.pt, DL_boosted_ee, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subjet1_pt_mumu", self.ak8Jets[0].subJet1.pt, DL_boosted_mumu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subjet1_pt_emu", self.ak8Jets[0].subJet1.pt, DL_boosted_emu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="First sub-jet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # subjet2 pt
                Plot.make1D("DL_boosted_subjet2_pt_ee", self.ak8Jets[0].subJet2.pt, DL_boosted_ee, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subjet2_pt_mumu", self.ak8Jets[0].subJet2.pt, DL_boosted_mumu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subjet2_pt_emu", self.ak8Jets[0].subJet2.pt, DL_boosted_emu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Second sub-jet p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # fatjet eta
                Plot.make1D("DL_boosted_fatJet_eta_ee", self.ak8Jets[0].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_fatJet_eta_mumu", self.ak8Jets[0].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_fatJet_eta_emu", self.ak8Jets[0].eta, DL_boosted_emu, EqBin(
                    30, -3, 3), title="eta(ak8jet)", xTitle="Fatjet \eta", plotopts=DLboostedEMU_label),

                # Invariant mass of leptons
                Plot.make1D("DL_boosted_InvM_ee", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_boosted_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_mumu", op.invariant_mass(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_boosted_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons  (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_emu", op.invariant_mass(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_boosted_emu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair  (GeV/c^{2})", plotopts=DLboostedEMU_label),

                # pt of the di-lepton
                Plot.make1D("DL_boosted_dileptonPt_ee", op.sum(self.firstOSElEl[0].pt, self.firstOSElEl[1].pt), DL_boosted_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electrons (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_dileptonPt_mumu", op.sum(self.firstOSMuMu[0].pt, self.firstOSMuMu[1].pt), DL_boosted_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of muons  (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_dileptonPt_emu", op.sum(self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_boosted_emu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electron-muon pair  (GeV/c^{2})", plotopts=DLboostedEMU_label),

                # MET pt
                Plot.make1D("DL_boosted_MET_pt_ee", tree.MET.pt, DL_boosted_ee, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_MET_pt_mumu", tree.MET.pt, DL_boosted_mumu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_MET_pt_emu", tree.MET.pt, DL_boosted_emu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLboostedEMU_label),

                # MET phi
                Plot.make1D("DL_boosted_MET_phi_ee", tree.MET.phi, DL_boosted_ee, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_MET_phi_mumu", tree.MET.phi, DL_boosted_mumu, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_MET_phi_emu", tree.MET.phi, DL_boosted_emu, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLboostedEMU_label),

                # total charge of leptons
                Plot.make1D("DL_boosted_totalCharge_ee", op.sum(self.firstOSElEl[0].charge, self.firstOSElEl[1].charge), DL_boosted_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_totalCharge_mumu", op.sum(self.firstOSMuMu[0].charge, self.firstOSMuMu[1].charge), DL_boosted_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons ", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_totalCharge_emu", op.sum(self.firstOSElMu[0].charge, self.firstOSElMu[1].charge), DL_boosted_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair ", plotopts=DLboostedEMU_label),

                # invariant mass of subjets
                Plot.make1D("DL_boosted_InvM_jj_ee", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_ee, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_jj_mumu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_mumu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_jj_emu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_emu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEMU_label),

                # leading lepton pt
                Plot.make1D("DL_boosted_leadingLepton_pt_ee", self.firstOSElEl[0].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_mumu", self.firstOSMuMu[0].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEMU_label),
                Plot.make1D("DL_boosted_electron_pt_emu", self.firstOSElMu[0].pt, DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading electron (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_muon_pt_emu", self.firstOSElMu[1].pt, DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading muon (GeV/c)", plotopts=DLboostedMuMu_label),
                # sub-leading lepton pt
                Plot.make1D("DL_boosted_subleadingLepton_pt_ee", self.firstOSElEl[1].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_mumu", self.firstOSMuMu[1].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].pt, self.firstOSElMu[0].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEMU_label),

                # leading lepton eta
                Plot.make1D("DL_boosted_leadingLepton_eta_ee", self.firstOSElEl[0].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_mumu", self.firstOSMuMu[0].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].eta, self.firstOSElMu[1].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEMU_label),
                Plot.make1D("DL_boosted_electron_eta_emu", self.firstOSElMu[0].eta, DL_boosted_emu, EqBin(
                    30, -3, 3), title="leadingLeptonEta", xTitle="eta of the electron (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_muon_eta_emu", self.firstOSElMu[1].eta, DL_boosted_emu, EqBin(
                    30, -3, 3), title="leadingLeptonEta", xTitle="eta of the muon (GeV/c)", plotopts=DLboostedMuMu_label),

                # sub-leading lepton eta
                Plot.make1D("DL_boosted_subleadingLepton_eta_ee", self.firstOSElEl[1].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_mumu", self.firstOSMuMu[1].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].eta, self.firstOSElMu[0].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEMU_label),

                # DR between leading and sub-leading lepton
                Plot.make1D("DL_boosted_DR_leptons_ee", op.deltaR(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leptons_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leptons_emu", op.deltaR(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEMU_label),

                # DR between leading lepton and ak8 jet
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_ee", op.deltaR(self.firstOSElEl[0].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEMU_label),

                # DR between subleading lepton and ak8 jet
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_ee", op.deltaR(self.firstOSElEl[1].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_mumu", op.deltaR(self.firstOSMuMu[1].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].p4, self.firstOSElMu[0].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedEMU_label),

                # number of electrons
                Plot.make1D("DL_boosted_nElectrons_ee", op.rng_len(self.tightElectrons), DL_boosted_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_boosted_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nElectrons_emu", op.rng_len(self.tightElectrons), DL_boosted_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEMU_label),

                # number of muons
                Plot.make1D("DL_boosted_nMuons_ee", op.rng_len(self.tightMuons), DL_boosted_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_nMuons_mumu", op.rng_len(self.tightMuons), DL_boosted_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_nMuons_emu", op.rng_len(self.tightMuons), DL_boosted_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLboostedEMU_label),

                #########################################
                ######                             ######
                ######      DL resolved plots      ######
                ######                             ######
                #########################################

                # number of ak4 bjets
                Plot.make1D("DL_resolved_nAK4bJets_ee", op.rng_len(self.ak4BJets), DL_resolved_ee, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nAK4bJets_mumu", op.rng_len(self.ak4BJets), DL_resolved_mumu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nAK4bJets_emu", op.rng_len(self.ak4BJets), DL_resolved_emu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolvedEMu_label),

                # ak4 bjet pt
                Plot.make1D("DL_resolved_ak4BJet_pt_ee", self.ak4BJets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_ak4BJet_pt_mumu", self.ak4BJets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_ak4BJet_pt_emu", self.ak4BJets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # ak4 bjet eta
                Plot.make1D("DL_resolved_ak4BJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_ak4BJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_ak4BJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolvedEMu_label),

                # number of ak4 jets
                Plot.make1D("DL_resolved_nak4Jets_ee", op.rng_len(self.ak4Jets), DL_resolved_ee, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nak4Jets_mumu", op.rng_len(self.ak4Jets), DL_resolved_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nak4Jets_emu", op.rng_len(self.ak4Jets), DL_resolved_emu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolvedEMu_label),

                # leading jet pt
                Plot.make1D("DL_resolved_leadingJet_pt_ee", self.ak4Jets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingJet_pt_mumu", self.ak4Jets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingJet_pt_emu", self.ak4Jets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # leading jet eta
                Plot.make1D("DL_resolved_leadingJet_eta_ee", self.ak4Jets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingJet_eta_mumu", self.ak4Jets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingJet_eta_emu", self.ak4Jets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolvedEMu_label),

                # btagging score of the jet
                Plot.make1D("DL_resolved_jet_btagScore_ee", self.ak4BJets[0].btagPNetB, DL_resolved_ee, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_jet_btagScore_mumu", self.ak4BJets[0].btagPNetB, DL_resolved_mumu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_jet_btagScore_emu", self.ak4BJets[0].btagPNetB, DL_resolved_emu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet", plotopts=DLresolvedEMu_label),

                # sub-leading jet pt
                Plot.make1D("DL_resolved_subleadingJet_pt_ee", self.ak4Jets[1].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingJet_pt_mumu", self.ak4Jets[1].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingJet_pt_emu", self.ak4Jets[1].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # sub-leading jet eta
                Plot.make1D("DL_resolved_subleadingJet_eta_ee", self.ak4Jets[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingJet_eta_mumu", self.ak4Jets[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingJet_eta_emu", self.ak4Jets[1].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolvedEMu_label),

                # btagging score of the jet
                Plot.make1D("DL_resolved_subleadingJet_btagScore_ee", self.ak4BJets[1].btagPNetB, DL_resolved_ee, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingJet_btagScore_mumu", self.ak4BJets[1].btagPNetB, DL_resolved_mumu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingJet_btagScore_emu", self.ak4BJets[1].btagPNetB, DL_resolved_emu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet", plotopts=DLresolvedEMu_label),

                # DR between leading and sub-leading jet
                Plot.make1D("DL_resolved_DR_jets_ee", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_jets_mumu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_jets_emu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolvedEMu_label),

                # Invariant mass of leptons
                Plot.make1D("DL_resolved_InvM_ee", op.invariant_mass(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_InvM_mumu", op.invariant_mass(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_InvM_emu", op.invariant_mass(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (GeV/c^{2})", plotopts=DLresolvedEMu_label),

                # pt of the di-lepton
                Plot.make1D("DL_resolved_dileptonPt_ee", op.sum(self.firstOSElEl[0].pt, self.firstOSElEl[1].pt), DL_resolved_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electrons (GeV/c^{2})", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_dileptonPt_mumu", op.sum(self.firstOSMuMu[0].pt, self.firstOSMuMu[1].pt), DL_resolved_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of muons  (GeV/c^{2})", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_dileptonPt_emu", op.sum(self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_resolved_emu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electron-muon pair  (GeV/c^{2})", plotopts=DLresolvedEMu_label),

                # MET pt
                Plot.make1D("DL_resolved_MET_pt_ee", tree.MET.pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_MET_pt_mumu", tree.MET.pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_MET_pt_emu", tree.MET.pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLresolvedEMu_label),

                # MET phi
                Plot.make1D("DL_resolved_MET_phi_ee", tree.MET.phi, DL_resolved_ee, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_MET_phi_mumu", tree.MET.phi, DL_resolved_mumu, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_MET_phi_emu", tree.MET.phi, DL_resolved_emu, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLresolvedEMu_label),

                # total charge of leptons
                Plot.make1D("DL_resolved_totalCharge_ee", op.sum(self.firstOSElEl[0].charge, self.firstOSElEl[1].charge), DL_resolved_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_totalCharge_mumu", op.sum(self.firstOSMuMu[0].charge, self.firstOSMuMu[1].charge), DL_resolved_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_totalCharge_emu", op.sum(self.firstOSElMu[0].charge, self.firstOSElMu[1].charge), DL_resolved_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair", plotopts=DLresolvedEMu_label),

                # leading lepton pt
                Plot.make1D("DL_resolved_leadingLepton_pt_ee", self.firstOSElEl[0].pt, DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingLepton_pt_mumu", self.firstOSMuMu[0].pt, DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].pt, self.firstOSElMu[1].pt), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolvedEMu_label),

                # sub-leading lepton pt
                Plot.make1D("DL_resolved_subleadingLepton_pt_ee", self.firstOSElEl[1].pt, DL_resolved_ee, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingLepton_pt_mumu", self.firstOSMuMu[1].pt, DL_resolved_mumu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingLepton_pt_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].pt, self.firstOSElMu[0].pt), DL_resolved_emu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolvedEMu_label),

                # leading lepton eta
                Plot.make1D("DL_resolved_leadingLepton_eta_ee", self.firstOSElEl[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_leadingLepton_eta_mumu", self.firstOSMuMu[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_leadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].eta, self.firstOSElMu[1].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolvedEMu_label),

                # sub-leading lepton eta
                Plot.make1D("DL_resolved_subleadingLepton_eta_ee", self.firstOSElEl[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_subleadingLepton_eta_mumu", self.firstOSMuMu[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_subleadingLepton_eta_emu", op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].eta, self.firstOSElMu[0].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolvedEMu_label),

                # DR between leading and sub-leading lepton
                Plot.make1D("DL_resolved_DR_leptons_ee", op.deltaR(self.firstOSElEl[0].p4, self.firstOSElEl[1].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_leptons_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.firstOSMuMu[1].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_leptons_emu", op.deltaR(self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolvedEMu_label),

                # DR between leading lepton and ak4 b jet
                Plot.make1D("DL_resolved_DR_leadingleptonANDak4bjet_ee", op.deltaR(self.firstOSElEl[0].p4, self.ak4BJets[0].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_leadingleptonANDak4bjet_mumu", op.deltaR(self.firstOSMuMu[0].p4, self.ak4BJets[0].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_leadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[0].p4, self.firstOSElMu[1].p4), self.ak4BJets[0].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolvedEMu_label),

                # DR between sub-leading lepton and ak4 b jet
                Plot.make1D("DL_resolved_DR_subleadingleptonANDak4bjet_ee", op.deltaR(self.firstOSElEl[1].p4, self.ak4BJets[0].p4), DL_resolved_ee, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_DR_subleadingleptonANDak4bjet_mumu", op.deltaR(self.firstOSMuMu[1].p4, self.ak4BJets[0].p4), DL_resolved_mumu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_DR_subleadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstOSElMu[0].pt >= self.firstOSElMu[1].pt), self.firstOSElMu[1].p4, self.firstOSElMu[0].p4), self.ak4BJets[0].p4), DL_resolved_emu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolvedEMu_label),

                # number of electrons
                Plot.make1D("DL_resolved_nElectrons_ee", op.rng_len(self.tightElectrons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nElectrons_emu", op.rng_len(self.tightElectrons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEMu_label),

                # number of muons
                Plot.make1D("DL_resolved_nMuons_ee", op.rng_len(self.tightMuons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEE_label),
                Plot.make1D("DL_resolved_nMuons_mumu", op.rng_len(self.tightMuons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedMuMu_label),
                Plot.make1D("DL_resolved_nMuons_emu", op.rng_len(self.tightMuons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolvedEMu_label),
            ])
        if self.channel == "SL":
            # if self.args.sync:
            #     plots.extend([Skim("SL_resolved_1b_e_sync",
            #                        syncVars_SL, SL_resolved_1b_e),  # update sync_vars
            #                   Skim("SL_resolved_2b_e_sync",
            #                        syncVars_SL, SL_resolved_2b_e),
            #                   Skim("SL_resolved_1b_mu_sync",
            #                        syncVars_SL, SL_resolved_1b_mu),
            #                   Skim("SL_resolved_2b_mu_sync",
            #                        syncVars_SL, SL_resolved_2b_mu)])
            plots.extend([

                #########################################
                #                 Skims                 #
                #########################################

                Skim("SL_resolved_1b_e_mva",
                     mvaVars_SL_resolved, SL_resolved_1b_e),
                Skim("SL_resolved_2b_e_mva",
                     mvaVars_SL_resolved, SL_resolved_2b_e),
                Skim("SL_resolved_1b_mu_mva",
                     mvaVars_SL_resolved, SL_resolved_1b_mu),
                Skim("SL_resolved_2b_mu_mva",
                     mvaVars_SL_resolved, SL_resolved_2b_mu),

                #########################################
                ######                             ######
                ######       SL boosted plots      ######
                ######                             ######
                #########################################

                # number of fat b-jets
                Plot.make1D("SL_boosted_nfatJet_e", op.rng_len(self.ak8BJets), SL_boosted_e, EqBin(
                    10, 0, 10), title="N(ak8bjet)", xTitle="Number of fat b-jet", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_nfatJet_mu", op.rng_len(self.ak8BJets), SL_boosted_mu, EqBin(
                    10, 0, 10), title="N(ak8bjet)", xTitle="Number of fat b-jet", plotopts=SLboostedMu_label),

                # fatjet pt
                Plot.make1D("SL_boosted_fatJet_pt_e", self.ak8BJets[0].pt, SL_boosted_e, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="Fat b-jet p_{T} (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_fatJet_pt_mu", self.ak8BJets[0].pt, SL_boosted_mu, EqBin(
                    400, 200, 1000), title="pT(j)", xTitle="Fat b-jet p_{T} (GeV/c)", plotopts=SLboostedMu_label),

                # subjet1 pt
                Plot.make1D("SL_boosted_subjet1_pt_e", self.ak8BJets[0].subJet1.pt, SL_boosted_e, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="Subjet 1 p_{T} (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet1_pt_mu", self.ak8BJets[0].subJet1.pt, SL_boosted_mu, EqBin(
                    50, 0, 500), title=" pT(subjet1)", xTitle="Subjet 1 p_{T} (GeV/c)", plotopts=SLboostedMu_label),

                # subjet2 pt
                Plot.make1D("SL_boosted_subjet2_pt_e", self.ak8BJets[0].subJet2.pt, SL_boosted_e, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Subjet 2 p_{T} (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet2_pt_mu", self.ak8BJets[0].subJet2.pt, SL_boosted_mu, EqBin(
                    50, 0, 500), title=" pT(subjet2)", xTitle="Subjet 2 p_{T} (GeV/c)", plotopts=SLboostedMu_label),

                # ak8jet eta
                Plot.make1D("SL_boosted_fatJet_eta_e", self.ak8BJets[0].eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(j)", xTitle="eta(j)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_fatJet_eta_mu", self.ak8BJets[0].eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(j)", xTitle="eta(j)", plotopts=SLboostedMu_label),

                # subjet1 eta
                Plot.make1D("SL_boosted_subjet1_eta_e", self.ak8BJets[0].subJet1.eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(subjet1)", xTitle="Subjet 1 \eta", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet1_eta_mu", self.ak8BJets[0].subJet1.eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(subjet1)", xTitle="Subjet 1 \eta", plotopts=SLboostedMu_label),

                # subjet2 eta
                Plot.make1D("SL_boosted_subjet2_eta_e", self.ak8BJets[0].subJet2.eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="eta(subjet2)", xTitle="Subjet 2 \eta", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_subjet2_eta_mu", self.ak8BJets[0].subJet2.eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="eta(subjet2)", xTitle="Subjet 2 \eta", plotopts=SLboostedMu_label),

                # Invariant mass of subjets
                Plot.make1D("SL_boosted_InvM_jj_e", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_e, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_InvM_jj_mu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), SL_boosted_mu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=SLboostedMu_label),

                # lepton pt
                Plot.make1D("SL_boosted_lepton_pt_e", self.tightElectrons[0].pt, SL_boosted_e, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_lepton_pt_mu", self.tightMuons[0].pt, SL_boosted_mu, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLboostedMu_label),

                # lepton eta
                Plot.make1D("SL_boosted_lepton_eta_e", self.tightElectrons[0].eta, SL_boosted_e, EqBin(
                    30, -3, 3), title="lepton eta", xTitle="lepton \eta", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_lepton_eta_mu", self.tightMuons[0].eta, SL_boosted_mu, EqBin(
                    30, -3, 3), title="lepton eta", xTitle="lepton \eta", plotopts=SLboostedMu_label),

                # DR between lepton and ak8 jet
                Plot.make1D("SL_boosted_DR_leptonANDak8bjet_e", op.deltaR(self.tightElectrons[0].p4, self.ak8Jets[0].p4), SL_boosted_e, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_DR_leptonANDak8bjet_mu", op.deltaR(self.tightMuons[0].p4, self.ak8Jets[0].p4), SL_boosted_mu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=SLboostedMu_label),

                # MET pt
                Plot.make1D("SL_boosted_MET_pt_e", tree.MET.pt, SL_boosted_e, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=SLboostedE_label),
                Plot.make1D("SL_boosted_MET_pt_mu", tree.MET.pt, SL_boosted_mu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=SLboostedMu_label),


                #########################################
                ######                             ######
                ######      SL resolved plots      ######
                ######                             ######
                #########################################

                # number of ak4 jets
                Plot.make1D("SL_resolved_nJets_e", op.rng_len(self.ak4BJets), SL_resolved_1b_e, EqBin(
                    15, 0., 15.), xTitle="Number of b-jets", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_nJets_mu", op.rng_len(self.ak4BJets), SL_resolved_1b_mu, EqBin(
                    15, 0., 15.), xTitle="Number of b-jets", plotopts=SLresolvedMu_label),

                # leading jet pt
                Plot.make1D("SL_resolved_leadingJet_pt_e", self.ak4BJets[0].pt, SL_resolved_1b_e, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_leadingJet_pt_mu", self.ak4BJets[0].pt, SL_resolved_1b_mu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="pT(j1) (GeV/c)", plotopts=SLresolvedMu_label),

                # leading jet eta
                Plot.make1D("SL_resolved_leadingJet_eta_e", self.ak4BJets[0].eta, SL_resolved_1b_e, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="B-jet \eta", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_leadingJet_eta_mu", self.ak4BJets[0].eta, SL_resolved_1b_mu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="eta(j1)", plotopts=SLresolvedMu_label),

                # sub-leading jet pt
                Plot.make1D("SL_resolved_subleadingJet_pt_e", self.ak4Jets[1].pt, SL_resolved_1b_e, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_subleadingJet_pt_mu", self.ak4Jets[1].pt, SL_resolved_1b_mu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="pT(j2) (GeV/c)", plotopts=SLresolvedMu_label),

                # sub-leading jet eta
                Plot.make1D("SL_resolved_subleadingJet_eta_e", self.ak4Jets[1].eta, SL_resolved_1b_e, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="eta(j2)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_subleadingJet_eta_mu", self.ak4Jets[1].eta, SL_resolved_1b_mu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="eta(j2)", plotopts=SLresolvedMu_label),

                # DR between leading and sub-leading jet
                Plot.make1D("SL_resolved_DR_jets_e", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), SL_resolved_1b_e, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_DR_jets_mu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), SL_resolved_1b_mu, EqBin(
                    100, 0, 10), title="DR(j1,j2)", xTitle="DR(j1,j2)", plotopts=SLresolvedMu_label),

                # DR between  lepton and ak4 b jet
                Plot.make1D("SL_resolved_DR_leptonANDak4bjet_e", op.deltaR(self.tightElectrons[0].p4, self.ak4BJets[0].p4), SL_resolved_1b_e, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(lepton, ak8bjet)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_DR_leptonANDak4bjet_mu", op.deltaR(self.tightMuons[0].p4, self.ak4BJets[0].p4), SL_resolved_1b_mu, EqBin(
                    35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(lepton, ak8bjet)", plotopts=SLresolvedMu_label),

                # lepton pt
                Plot.make1D("SL_resolved_lepton_pt_e", self.tightElectrons[0].pt, SL_resolved_1b_e, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_lepton_pt_mu", self.tightMuons[0].pt, SL_resolved_1b_mu, EqBin(
                    100, 0., 300.), title="lepton pT", xTitle="p_{T} of the lepton (GeV/c)", plotopts=SLresolvedMu_label),

                # lepton eta
                Plot.make1D("SL_resolved_lepton_eta_e", self.tightElectrons[0].eta, SL_resolved_1b_e, EqBin(
                    30, -3, 3), title="lepton pT", xTitle="\eta of the lepton (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_lepton_eta_mu", self.tightMuons[0].eta, SL_resolved_1b_mu, EqBin(
                    30, -3, 3), title="lepton pT", xTitle="\eta of the lepton (GeV/c)", plotopts=SLresolvedMu_label),

                # MET pt
                Plot.make1D("SL_resolved_MET_pt_e", tree.MET.pt, SL_resolved_1b_e, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=SLresolvedE_label),
                Plot.make1D("SL_resolved_MET_pt_mu", tree.MET.pt, SL_resolved_1b_mu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=SLresolvedMu_label),
            ])

        return plots
