
from bamboo.plots import Plot, Skim
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection, makeSLSelection
from scalefactors import ScaleFactors as sf
import definitions as defs
from utils import labeler


class cutflowAnalysis(NanoBaseHHWWbb):
    """ Class to create control plots, cutflow reports and skims"""

    def __init__(self, args):
        super(cutflowAnalysis, self).__init__(args)
        self.channel = self.args.channel

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # add cutflow report
        plots.append(self.yields)

        # define objects
        defs.defineObjects(self, tree)

        # common scale factors
        noSel = sf.commonSF(self, tree, noSel, sample)

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

        event_selections = [
            DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
            DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
            DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu]

        # cutflow report for DL channel
        for sel in event_selections:
            self.yields.add(sel, sel.name)

        # initiate TNN input features' dictionary
        tnn_input_vars = {
            "event_no": tree.event,
            "weight": noSel.weight,
        }

        # get the input features
        l1, l2, j1, j2, met = defs.DNN_variables(self, tree)

        # concatenate with the initial dictionary
        tnn_input_vars = tnn_input_vars | l1 | l2 | j1 | j2 | met
        # line above is equivalent to the following (concetanation of dictionaries)
        # tnn_input_vars =  {**tnn_input_vars, **l1, **l2, **j1, **j2, **met}

        # create skims which will hold DNN input features
        for selection in event_selections:
            plots.append(
                Skim(selection.name+'_tnn_input_vars',
                     tnn_input_vars, selection)
            )

        def tnn_input_var_binning(var_name):
            "Function to return binning, min and max values for the TNN input feature plots."
            if "_Px" in var_name or "_Py" in var_name:
                N, mn, mx = 200, -1000, 1000
            elif "_Pz" in var_name:
                N, mn, mx = 200, -4000, 4000
            elif "_E" in var_name:
                N, mn, mx = 250, 0, 2500
            elif "_charge" in var_name:
                N, mn, mx = 5, -2.5, 2.5
            elif "_btag" in var_name:
                N, mn, mx = 100, 0, 1
            elif "_pdgId" in var_name:
                N, mn, mx = 30, -15, 15

            return EqBin(N, mn, mx)

        # We're not interested in the following two variables' match between data and MC.
        # Hence they're not included in the input feature plots.
        tnn_input_vars.pop('event_no')
        tnn_input_vars.pop('weight')

        for selection in event_selections:
            for name, var in tnn_input_vars.items():
                plots.append(
                    Plot.make1D(name+"_"+selection.name, var, selection,
                                tnn_input_var_binning(name), title=name, xTitle=name)
                )

        boosted_categories = [
            sel for sel in event_selections if 'boosted' in sel.name]
        resolved_categories = [
            sel for sel in event_selections if 'resolved' in sel.name]

        vars_to_plot_for_boosted = {
            'n_ak8': [op.rng_len(self.ak8Jets), EqBin(10, 0, 10), 'Number of AK8 Jets'],
            'ak8_pT': [self.ak8Jets[0].pt, EqBin(100, 200, 800), 'AK8 Jet p_T'],
            'ak8_eta': [self.ak8Jets[0].eta, EqBin(30, -3, 3), 'AK8 Jet \eta'],
            'ak8_subjet1_pT': [self.ak8Jets[0].subJet1.pt, EqBin(50, 0, 500), 'AK8 leading sub-jet p_T'],
            'ak8_subjet2_pT': [self.ak8Jets[0].subJet2.pt, EqBin(50, 0, 500), 'AK8 sub-leading sub-jet p_T'],
            'met_pT': [tree.MET.pt, EqBin(100, 0, 500), 'MET p_{T} (GeV/c)'],
            'met_phi': [tree.MET.phi, EqBin(100, 0, 500), 'MET \phi'],
            'invM_ak8_subjets': [op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4),
                                 EqBin(100, 0., 200.), 'Invariant Mass of sub-jets (GeV/c^{2})'],
            'n_electrons': [op.rng_len(self.tightElectrons), EqBin(3, 0, 3), 'Number of electrons'],
            'n_muons': [op.rng_len(self.tightMuons), EqBin(3, 0, 3), 'Number of muons'],
        }

        for sel in boosted_categories:
            for var_name, var in vars_to_plot_for_boosted.items():
                plots.append(
                    Plot.make1D(
                        sel.name+"_"+var_name, var[0], sel, var[1], title=var_name, xTitle=var[2], plotopts=labeler(sel.name))
                )

        plots.extend([
            # Invariant mass of leptons
            Plot.make1D("DL_boosted_InvM_ee", op.invariant_mass(self.firstElTightPair[0].p4, self.firstElTightPair[1].p4), DL_boosted_ee, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})"),
            Plot.make1D("DL_boosted_InvM_mumu", op.invariant_mass(self.firstMuTightPair[0].p4, self.firstMuTightPair[1].p4), DL_boosted_mumu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons  (GeV/c^{2})"),
            Plot.make1D("DL_boosted_InvM_emu", op.invariant_mass(self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), DL_boosted_emu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair  (GeV/c^{2})"),

            # pt of the di-lepton
            Plot.make1D("DL_boosted_dileptonPt_ee", op.sum(self.firstElTightPair[0].pt, self.firstElTightPair[1].pt), DL_boosted_ee, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electrons (GeV/c^{2})"),
            Plot.make1D("DL_boosted_dileptonPt_mumu", op.sum(self.firstMuTightPair[0].pt, self.firstMuTightPair[1].pt), DL_boosted_mumu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of muons  (GeV/c^{2})"),
            Plot.make1D("DL_boosted_dileptonPt_emu", op.sum(self.firstEmuTightPair[0].pt, self.firstEmuTightPair[1].pt), DL_boosted_emu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electron-muon pair  (GeV/c^{2})"),

            # total charge of leptons
            Plot.make1D("DL_boosted_totalCharge_ee", op.sum(self.firstElTightPair[0].charge, self.firstElTightPair[1].charge), DL_boosted_ee, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons"),
            Plot.make1D("DL_boosted_totalCharge_mumu", op.sum(self.firstMuTightPair[0].charge, self.firstMuTightPair[1].charge), DL_boosted_mumu, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons "),
            Plot.make1D("DL_boosted_totalCharge_emu", op.sum(self.firstEmuTightPair[0].charge, self.firstEmuTightPair[1].charge), DL_boosted_emu, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair "),

            # leading lepton pt
            Plot.make1D("DL_boosted_leadingLepton_pt_ee", self.firstElTightPair[0].pt, DL_boosted_ee, EqBin(
                100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_boosted_leadingLepton_pt_mumu", self.firstMuTightPair[0].pt, DL_boosted_mumu, EqBin(
                100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_boosted_leadingLepton_pt_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].pt, self.firstEmuTightPair[1].pt), DL_boosted_emu, EqBin(
                100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_boosted_electron_pt_emu", self.firstEmuTightPair[0].pt, DL_boosted_emu, EqBin(
                100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading electron (GeV/c)"),
            Plot.make1D("DL_boosted_muon_pt_emu", self.firstEmuTightPair[1].pt, DL_boosted_emu, EqBin(
                100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading muon (GeV/c)"),

            # sub-leading lepton pt
            Plot.make1D("DL_boosted_subleadingLepton_pt_ee", self.firstElTightPair[1].pt, DL_boosted_ee, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),
            Plot.make1D("DL_boosted_subleadingLepton_pt_mumu", self.firstMuTightPair[1].pt, DL_boosted_mumu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),
            Plot.make1D("DL_boosted_subleadingLepton_pt_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].pt, self.firstEmuTightPair[0].pt), DL_boosted_emu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),

            # leading lepton eta
            Plot.make1D("DL_boosted_leadingLepton_eta_ee", self.firstElTightPair[0].eta, DL_boosted_ee, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_boosted_leadingLepton_eta_mumu", self.firstMuTightPair[0].eta, DL_boosted_mumu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_boosted_leadingLepton_eta_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].eta, self.firstEmuTightPair[1].eta), DL_boosted_emu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_boosted_electron_eta_emu", self.firstEmuTightPair[0].eta, DL_boosted_emu, EqBin(
                30, -3, 3), title="leadingLeptonEta", xTitle="eta of the electron (GeV/c)"),
            Plot.make1D("DL_boosted_muon_eta_emu", self.firstEmuTightPair[1].eta, DL_boosted_emu, EqBin(
                30, -3, 3), title="leadingLeptonEta", xTitle="eta of the muon (GeV/c)"),

            # sub-leading lepton eta
            Plot.make1D("DL_boosted_subleadingLepton_eta_ee", self.firstElTightPair[1].eta, DL_boosted_ee, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),
            Plot.make1D("DL_boosted_subleadingLepton_eta_mumu", self.firstMuTightPair[1].eta, DL_boosted_mumu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),
            Plot.make1D("DL_boosted_subleadingLepton_eta_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].eta, self.firstEmuTightPair[0].eta), DL_boosted_emu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),

            # DR between leading and sub-leading lepton
            Plot.make1D("DL_boosted_DR_leptons_ee", op.deltaR(self.firstElTightPair[0].p4, self.firstElTightPair[1].p4), DL_boosted_ee, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),
            Plot.make1D("DL_boosted_DR_leptons_mumu", op.deltaR(self.firstMuTightPair[0].p4, self.firstMuTightPair[1].p4), DL_boosted_mumu, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),
            Plot.make1D("DL_boosted_DR_leptons_emu", op.deltaR(self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), DL_boosted_emu, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),

            # DR between leading lepton and ak8 jet
            Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_ee", op.deltaR(self.firstElTightPair[0].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),
            Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_mumu", op.deltaR(self.firstMuTightPair[0].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),
            Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),

            # DR between subleading lepton and ak8 jet
            Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_ee", op.deltaR(self.firstElTightPair[1].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),
            Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_mumu", op.deltaR(self.firstMuTightPair[1].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),
            Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].p4, self.firstEmuTightPair[0].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),


            #########################################
            ######                             ######
            ######      DL resolved plots      ######
            ######                             ######
            #########################################

            # number of ak4 bjets
            Plot.make1D("DL_resolved_1b_nAK4bJets_ee", op.rng_len(self.ak4BJets), DL_resolved_1b_ee, EqBin(
                10, 0., 10), xTitle="Number of AK4 B-jets"),
            Plot.make1D("DL_resolved_1b_nAK4bJets_mumu", op.rng_len(self.ak4BJets), DL_resolved_1b_mumu, EqBin(
                10, 0., 10), xTitle="Number of AK4 B-jets"),
            Plot.make1D("DL_resolved_1b_nAK4bJets_emu", op.rng_len(self.ak4BJets), DL_resolved_1b_emu, EqBin(
                10, 0., 10), xTitle="Number of AK4 B-jets"),
            Plot.make1D("DL_resolved_2b_nAK4bJets_ee", op.rng_len(self.ak4BJets), DL_resolved_2b_ee, EqBin(
                10, 0., 10), xTitle="Number of AK4 B-jets"),
            Plot.make1D("DL_resolved_2b_nAK4bJets_mumu", op.rng_len(self.ak4BJets), DL_resolved_2b_mumu, EqBin(
                10, 0., 10), xTitle="Number of AK4 B-jets"),
            Plot.make1D("DL_resolved_2b_nAK4bJets_emu", op.rng_len(self.ak4BJets), DL_resolved_2b_emu, EqBin(
                10, 0., 10), xTitle="Number of AK4 B-jets"),

            # ak4 bjet pt
            Plot.make1D("DL_resolved_1b_ak4BJet_pt_ee", self.ak4BJets[0].pt, DL_resolved_1b_ee, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_ak4BJet_pt_mumu", self.ak4BJets[0].pt, DL_resolved_1b_mumu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_ak4BJet_pt_emu", self.ak4BJets[0].pt, DL_resolved_1b_emu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_ak4BJet_pt_ee", self.ak4BJets[0].pt, DL_resolved_2b_ee, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_ak4BJet_pt_mumu", self.ak4BJets[0].pt, DL_resolved_2b_mumu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_ak4BJet_pt_emu", self.ak4BJets[0].pt, DL_resolved_2b_emu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)"),

            # ak4 bjet eta
            Plot.make1D("DL_resolved_1b_ak4BJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_1b_ee, EqBin(
                30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta"),
            Plot.make1D("DL_resolved_1b_ak4BJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_1b_mumu, EqBin(
                30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta"),
            Plot.make1D("DL_resolved_1b_ak4BJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_1b_emu, EqBin(
                30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta"),
            Plot.make1D("DL_resolved_2b_ak4BJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_2b_ee, EqBin(
                30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta"),
            Plot.make1D("DL_resolved_2b_ak4BJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_2b_mumu, EqBin(
                30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta"),
            Plot.make1D("DL_resolved_2b_ak4BJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_2b_emu, EqBin(
                30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta"),

            # number of ak4 jets
            Plot.make1D("DL_resolved_1b_nak4Jets_ee", op.rng_len(self.ak4Jets), DL_resolved_1b_ee, EqBin(
                15, 0., 15.), xTitle="Number of AK4 jets"),
            Plot.make1D("DL_resolved_1b_nak4Jets_mumu", op.rng_len(self.ak4Jets), DL_resolved_1b_mumu, EqBin(
                15, 0., 15.), xTitle="Number of AK4 jets"),
            Plot.make1D("DL_resolved_1b_nak4Jets_emu", op.rng_len(self.ak4Jets), DL_resolved_1b_emu, EqBin(
                15, 0., 15.), xTitle="Number of AK4 jets"),
            Plot.make1D("DL_resolved_2b_nak4Jets_ee", op.rng_len(self.ak4Jets), DL_resolved_2b_ee, EqBin(
                15, 0., 15.), xTitle="Number of AK4 jets"),
            Plot.make1D("DL_resolved_2b_nak4Jets_mumu", op.rng_len(self.ak4Jets), DL_resolved_2b_mumu, EqBin(
                15, 0., 15.), xTitle="Number of AK4 jets"),
            Plot.make1D("DL_resolved_2b_nak4Jets_emu", op.rng_len(self.ak4Jets), DL_resolved_2b_emu, EqBin(
                15, 0., 15.), xTitle="Number of AK4 jets"),

            # leading jet pt
            Plot.make1D("DL_resolved_1b_leadingJet_pt_ee", self.ak4Jets[0].pt, DL_resolved_1b_ee, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_leadingJet_pt_mumu", self.ak4Jets[0].pt, DL_resolved_1b_mumu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_leadingJet_pt_emu", self.ak4Jets[0].pt, DL_resolved_1b_emu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_leadingJet_pt_ee", self.ak4Jets[0].pt, DL_resolved_2b_ee, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_leadingJet_pt_mumu", self.ak4Jets[0].pt, DL_resolved_2b_mumu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_leadingJet_pt_emu", self.ak4Jets[0].pt, DL_resolved_2b_emu, EqBin(
                100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)"),

            # leading jet eta
            Plot.make1D("DL_resolved_1b_leadingJet_eta_ee", self.ak4Jets[0].eta, DL_resolved_1b_ee, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta"),
            Plot.make1D("DL_resolved_1b_leadingJet_eta_mumu", self.ak4Jets[0].eta, DL_resolved_1b_mumu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta"),
            Plot.make1D("DL_resolved_1b_leadingJet_eta_emu", self.ak4Jets[0].eta, DL_resolved_1b_emu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta"),
            Plot.make1D("DL_resolved_2b_leadingJet_eta_ee", self.ak4Jets[0].eta, DL_resolved_2b_ee, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta"),
            Plot.make1D("DL_resolved_2b_leadingJet_eta_mumu", self.ak4Jets[0].eta, DL_resolved_2b_mumu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta"),
            Plot.make1D("DL_resolved_2b_leadingJet_eta_emu", self.ak4Jets[0].eta, DL_resolved_2b_emu, EqBin(
                30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta"),

            # btagging score of the jet
            Plot.make1D("DL_resolved_1b_jet_btagScore_ee", self.ak4BJets[0].btagPNetB, DL_resolved_1b_ee, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet"),
            Plot.make1D("DL_resolved_1b_jet_btagScore_mumu", self.ak4BJets[0].btagPNetB, DL_resolved_1b_mumu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet"),
            Plot.make1D("DL_resolved_1b_jet_btagScore_emu", self.ak4BJets[0].btagPNetB, DL_resolved_1b_emu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet"),
            Plot.make1D("DL_resolved_2b_jet_btagScore_ee", self.ak4BJets[0].btagPNetB, DL_resolved_2b_ee, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet"),
            Plot.make1D("DL_resolved_2b_jet_btagScore_mumu", self.ak4BJets[0].btagPNetB, DL_resolved_2b_mumu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet"),
            Plot.make1D("DL_resolved_2b_jet_btagScore_emu", self.ak4BJets[0].btagPNetB, DL_resolved_2b_emu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet"),

            # sub-leading jet pt
            Plot.make1D("DL_resolved_1b_subleadingJet_pt_ee", self.ak4Jets[1].pt, DL_resolved_1b_ee, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_subleadingJet_pt_mumu", self.ak4Jets[1].pt, DL_resolved_1b_mumu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_subleadingJet_pt_emu", self.ak4Jets[1].pt, DL_resolved_1b_emu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_subleadingJet_pt_ee", self.ak4Jets[1].pt, DL_resolved_2b_ee, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_subleadingJet_pt_mumu", self.ak4Jets[1].pt, DL_resolved_2b_mumu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_subleadingJet_pt_emu", self.ak4Jets[1].pt, DL_resolved_2b_emu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)"),

            # sub-leading jet eta
            Plot.make1D("DL_resolved_1b_subleadingJet_eta_ee", self.ak4Jets[1].eta, DL_resolved_1b_ee, EqBin(
                30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta"),
            Plot.make1D("DL_resolved_1b_subleadingJet_eta_mumu", self.ak4Jets[1].eta, DL_resolved_1b_mumu, EqBin(
                30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta"),
            Plot.make1D("DL_resolved_1b_subleadingJet_eta_emu", self.ak4Jets[1].eta, DL_resolved_1b_emu, EqBin(
                30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta"),
            Plot.make1D("DL_resolved_2b_subleadingJet_eta_ee", self.ak4Jets[1].eta, DL_resolved_2b_ee, EqBin(
                30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta"),
            Plot.make1D("DL_resolved_2b_subleadingJet_eta_mumu", self.ak4Jets[1].eta, DL_resolved_2b_mumu, EqBin(
                30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta"),
            Plot.make1D("DL_resolved_2b_subleadingJet_eta_emu", self.ak4Jets[1].eta, DL_resolved_2b_emu, EqBin(
                30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta"),

            # btagging score of the jet
            Plot.make1D("DL_resolved_1b_subleadingJet_btagScore_ee", self.ak4BJets[1].btagPNetB, DL_resolved_1b_ee, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet"),
            Plot.make1D("DL_resolved_1b_subleadingJet_btagScore_mumu", self.ak4BJets[1].btagPNetB, DL_resolved_1b_mumu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet"),
            Plot.make1D("DL_resolved_1b_subleadingJet_btagScore_emu", self.ak4BJets[1].btagPNetB, DL_resolved_1b_emu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet"),
            Plot.make1D("DL_resolved_2b_subleadingJet_btagScore_ee", self.ak4BJets[1].btagPNetB, DL_resolved_2b_ee, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet"),
            Plot.make1D("DL_resolved_2b_subleadingJet_btagScore_mumu", self.ak4BJets[1].btagPNetB, DL_resolved_2b_mumu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet"),
            Plot.make1D("DL_resolved_2b_subleadingJet_btagScore_emu", self.ak4BJets[1].btagPNetB, DL_resolved_2b_emu, EqBin(
                100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet"),

            # DR between leading and sub-leading jet
            Plot.make1D("DL_resolved_1b_DR_jets_ee", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_1b_ee, EqBin(
                35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets"),
            Plot.make1D("DL_resolved_1b_DR_jets_mumu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_1b_mumu, EqBin(
                35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets"),
            Plot.make1D("DL_resolved_1b_DR_jets_emu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_1b_emu, EqBin(
                35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets"),
            Plot.make1D("DL_resolved_2b_DR_jets_ee", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_2b_ee, EqBin(
                35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets"),
            Plot.make1D("DL_resolved_2b_DR_jets_mumu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_2b_mumu, EqBin(
                35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets"),
            Plot.make1D("DL_resolved_2b_DR_jets_emu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_2b_emu, EqBin(
                35, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets"),

            # Invariant mass of leptons
            Plot.make1D("DL_resolved_1b_InvM_ee", op.invariant_mass(self.firstElTightPair[0].p4, self.firstElTightPair[1].p4), DL_resolved_1b_ee, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})"),
            Plot.make1D("DL_resolved_1b_InvM_mumu", op.invariant_mass(self.firstMuTightPair[0].p4, self.firstMuTightPair[1].p4), DL_resolved_1b_mumu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})"),
            Plot.make1D("DL_resolved_1b_InvM_emu", op.invariant_mass(self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), DL_resolved_1b_emu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (GeV/c^{2})"),
            Plot.make1D("DL_resolved_2b_InvM_ee", op.invariant_mass(self.firstElTightPair[0].p4, self.firstElTightPair[1].p4), DL_resolved_2b_ee, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})"),
            Plot.make1D("DL_resolved_2b_InvM_mumu", op.invariant_mass(self.firstMuTightPair[0].p4, self.firstMuTightPair[1].p4), DL_resolved_2b_mumu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})"),
            Plot.make1D("DL_resolved_2b_InvM_emu", op.invariant_mass(self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), DL_resolved_2b_emu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (GeV/c^{2})"),

            # pt of the di-lepton
            Plot.make1D("DL_resolved_1b_dileptonPt_ee", op.sum(self.firstElTightPair[0].pt, self.firstElTightPair[1].pt), DL_resolved_1b_ee, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electrons (GeV/c^{2})"),
            Plot.make1D("DL_resolved_1b_dileptonPt_mumu", op.sum(self.firstMuTightPair[0].pt, self.firstMuTightPair[1].pt), DL_resolved_1b_mumu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of muons  (GeV/c^{2})"),
            Plot.make1D("DL_resolved_1b_dileptonPt_emu", op.sum(self.firstEmuTightPair[0].pt, self.firstEmuTightPair[1].pt), DL_resolved_1b_emu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electron-muon pair  (GeV/c^{2})"),
            Plot.make1D("DL_resolved_2b_dileptonPt_ee", op.sum(self.firstElTightPair[0].pt, self.firstElTightPair[1].pt), DL_resolved_2b_ee, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electrons (GeV/c^{2})"),
            Plot.make1D("DL_resolved_2b_dileptonPt_mumu", op.sum(self.firstMuTightPair[0].pt, self.firstMuTightPair[1].pt), DL_resolved_2b_mumu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of muons  (GeV/c^{2})"),
            Plot.make1D("DL_resolved_2b_dileptonPt_emu", op.sum(self.firstEmuTightPair[0].pt, self.firstEmuTightPair[1].pt), DL_resolved_2b_emu, EqBin(
                60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electron-muon pair  (GeV/c^{2})"),

            # MET pt
            Plot.make1D("DL_resolved_1b_MET_pt_ee", tree.MET.pt, DL_resolved_1b_ee, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_MET_pt_mumu", tree.MET.pt, DL_resolved_1b_mumu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_1b_MET_pt_emu", tree.MET.pt, DL_resolved_1b_emu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_MET_pt_ee", tree.MET.pt, DL_resolved_2b_ee, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_MET_pt_mumu", tree.MET.pt, DL_resolved_2b_mumu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)"),
            Plot.make1D("DL_resolved_2b_MET_pt_emu", tree.MET.pt, DL_resolved_2b_emu, EqBin(
                100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)"),

            # MET phi
            Plot.make1D("DL_resolved_1b_MET_phi_ee", tree.MET.phi, DL_resolved_1b_ee, EqBin(
                7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)"),
            Plot.make1D("DL_resolved_1b_MET_phi_mumu", tree.MET.phi, DL_resolved_1b_mumu, EqBin(
                7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)"),
            Plot.make1D("DL_resolved_1b_MET_phi_emu", tree.MET.phi, DL_resolved_1b_emu, EqBin(
                7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)"),
            Plot.make1D("DL_resolved_2b_MET_phi_ee", tree.MET.phi, DL_resolved_2b_ee, EqBin(
                7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)"),
            Plot.make1D("DL_resolved_2b_MET_phi_mumu", tree.MET.phi, DL_resolved_2b_mumu, EqBin(
                7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)"),
            Plot.make1D("DL_resolved_2b_MET_phi_emu", tree.MET.phi, DL_resolved_2b_emu, EqBin(
                7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)"),

            # total charge of leptons
            Plot.make1D("DL_resolved_1b_totalCharge_ee", op.sum(self.firstElTightPair[0].charge, self.firstElTightPair[1].charge), DL_resolved_1b_ee, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons"),
            Plot.make1D("DL_resolved_1b_totalCharge_mumu", op.sum(self.firstMuTightPair[0].charge, self.firstMuTightPair[1].charge), DL_resolved_1b_mumu, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons"),
            Plot.make1D("DL_resolved_1b_totalCharge_emu", op.sum(self.firstEmuTightPair[0].charge, self.firstEmuTightPair[1].charge), DL_resolved_1b_emu, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair"),
            Plot.make1D("DL_resolved_2b_totalCharge_ee", op.sum(self.firstElTightPair[0].charge, self.firstElTightPair[1].charge), DL_resolved_2b_ee, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons"),
            Plot.make1D("DL_resolved_2b_totalCharge_mumu", op.sum(self.firstMuTightPair[0].charge, self.firstMuTightPair[1].charge), DL_resolved_2b_mumu, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons"),
            Plot.make1D("DL_resolved_2b_totalCharge_emu", op.sum(self.firstEmuTightPair[0].charge, self.firstEmuTightPair[1].charge), DL_resolved_2b_emu, EqBin(
                5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair"),

            # leading lepton pt
            Plot.make1D("DL_resolved_1b_leadingLepton_pt_ee", self.firstElTightPair[0].pt, DL_resolved_1b_ee, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_1b_leadingLepton_pt_mumu", self.firstMuTightPair[0].pt, DL_resolved_1b_mumu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_1b_leadingLepton_pt_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].pt, self.firstEmuTightPair[1].pt), DL_resolved_1b_emu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_2b_leadingLepton_pt_ee", self.firstElTightPair[0].pt, DL_resolved_2b_ee, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_2b_leadingLepton_pt_mumu", self.firstMuTightPair[0].pt, DL_resolved_2b_mumu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_2b_leadingLepton_pt_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].pt, self.firstEmuTightPair[1].pt), DL_resolved_2b_emu, EqBin(
                100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)"),

            # sub-leading lepton pt
            Plot.make1D("DL_resolved_1b_subleadingLepton_pt_ee", self.firstElTightPair[1].pt, DL_resolved_1b_ee, EqBin(
                50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_1b_subleadingLepton_pt_mumu", self.firstMuTightPair[1].pt, DL_resolved_1b_mumu, EqBin(
                50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_1b_subleadingLepton_pt_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].pt, self.firstEmuTightPair[0].pt), DL_resolved_1b_emu, EqBin(
                50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_2b_subleadingLepton_pt_ee", self.firstElTightPair[1].pt, DL_resolved_2b_ee, EqBin(
                50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_2b_subleadingLepton_pt_mumu", self.firstMuTightPair[1].pt, DL_resolved_2b_mumu, EqBin(
                50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),
            Plot.make1D("DL_resolved_2b_subleadingLepton_pt_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].pt, self.firstEmuTightPair[0].pt), DL_resolved_2b_emu, EqBin(
                50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)"),

            # leading lepton eta
            Plot.make1D("DL_resolved_1b_leadingLepton_eta_ee", self.firstElTightPair[0].eta, DL_resolved_1b_ee, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_resolved_1b_leadingLepton_eta_mumu", self.firstMuTightPair[0].eta, DL_resolved_1b_mumu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_resolved_1b_leadingLepton_eta_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].eta, self.firstEmuTightPair[1].eta), DL_resolved_1b_emu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_resolved_2b_leadingLepton_eta_ee", self.firstElTightPair[0].eta, DL_resolved_2b_ee, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_resolved_2b_leadingLepton_eta_mumu", self.firstMuTightPair[0].eta, DL_resolved_2b_mumu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),
            Plot.make1D("DL_resolved_2b_leadingLepton_eta_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].eta, self.firstEmuTightPair[1].eta), DL_resolved_2b_emu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton"),

            # sub-leading lepton eta
            Plot.make1D("DL_resolved_1b_subleadingLepton_eta_ee", self.firstElTightPair[1].eta, DL_resolved_1b_ee, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),
            Plot.make1D("DL_resolved_1b_subleadingLepton_eta_mumu", self.firstMuTightPair[1].eta, DL_resolved_1b_mumu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),
            Plot.make1D("DL_resolved_1b_subleadingLepton_eta_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].eta, self.firstEmuTightPair[0].eta), DL_resolved_1b_emu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),
            Plot.make1D("DL_resolved_2b_subleadingLepton_eta_ee", self.firstElTightPair[1].eta, DL_resolved_2b_ee, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),
            Plot.make1D("DL_resolved_2b_subleadingLepton_eta_mumu", self.firstMuTightPair[1].eta, DL_resolved_2b_mumu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),
            Plot.make1D("DL_resolved_2b_subleadingLepton_eta_emu", op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].eta, self.firstEmuTightPair[0].eta), DL_resolved_2b_emu, EqBin(
                30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton"),

            # DR between leading and sub-leading lepton
            Plot.make1D("DL_resolved_1b_DR_leptons_ee", op.deltaR(self.firstElTightPair[0].p4, self.firstElTightPair[1].p4), DL_resolved_1b_ee, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),
            Plot.make1D("DL_resolved_1b_DR_leptons_mumu", op.deltaR(self.firstMuTightPair[0].p4, self.firstMuTightPair[1].p4), DL_resolved_1b_mumu, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),
            Plot.make1D("DL_resolved_1b_DR_leptons_emu", op.deltaR(self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), DL_resolved_1b_emu, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),
            Plot.make1D("DL_resolved_2b_DR_leptons_ee", op.deltaR(self.firstElTightPair[0].p4, self.firstElTightPair[1].p4), DL_resolved_2b_ee, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),
            Plot.make1D("DL_resolved_2b_DR_leptons_mumu", op.deltaR(self.firstMuTightPair[0].p4, self.firstMuTightPair[1].p4), DL_resolved_2b_mumu, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),
            Plot.make1D("DL_resolved_2b_DR_leptons_emu", op.deltaR(self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), DL_resolved_2b_emu, EqBin(
                35, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons"),

            # DR between leading lepton and ak4 b jet
            Plot.make1D("DL_resolved_1b_DR_leadingleptonANDak4bjet_ee", op.deltaR(self.firstElTightPair[0].p4, self.ak4BJets[0].p4), DL_resolved_1b_ee, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_1b_DR_leadingleptonANDak4bjet_mumu", op.deltaR(self.firstMuTightPair[0].p4, self.ak4BJets[0].p4), DL_resolved_1b_mumu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_1b_DR_leadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), self.ak4BJets[0].p4), DL_resolved_1b_emu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_2b_DR_leadingleptonANDak4bjet_ee", op.deltaR(self.firstElTightPair[0].p4, self.ak4BJets[0].p4), DL_resolved_2b_ee, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_2b_DR_leadingleptonANDak4bjet_mumu", op.deltaR(self.firstMuTightPair[0].p4, self.ak4BJets[0].p4), DL_resolved_2b_mumu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_2b_DR_leadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[0].p4, self.firstEmuTightPair[1].p4), self.ak4BJets[0].p4), DL_resolved_2b_emu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)"),

            # DR between sub-leading lepton and ak4 b jet
            Plot.make1D("DL_resolved_1b_DR_subleadingleptonANDak4bjet_ee", op.deltaR(self.firstElTightPair[1].p4, self.ak4BJets[0].p4), DL_resolved_1b_ee, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_1b_DR_subleadingleptonANDak4bjet_mumu", op.deltaR(self.firstMuTightPair[1].p4, self.ak4BJets[0].p4), DL_resolved_1b_mumu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_1b_DR_subleadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].p4, self.firstEmuTightPair[0].p4), self.ak4BJets[0].p4), DL_resolved_1b_emu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_2b_DR_subleadingleptonANDak4bjet_ee", op.deltaR(self.firstElTightPair[1].p4, self.ak4BJets[0].p4), DL_resolved_2b_ee, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_2b_DR_subleadingleptonANDak4bjet_mumu", op.deltaR(self.firstMuTightPair[1].p4, self.ak4BJets[0].p4), DL_resolved_2b_mumu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),
            Plot.make1D("DL_resolved_2b_DR_subleadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.firstEmuTightPair[0].pt >= self.firstEmuTightPair[1].pt), self.firstEmuTightPair[1].p4, self.firstEmuTightPair[0].p4), self.ak4BJets[0].p4), DL_resolved_2b_emu, EqBin(
                35, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)"),

            # number of electrons
            Plot.make1D("DL_resolved_1b_nElectrons_ee", op.rng_len(self.tightElectrons), DL_resolved_1b_ee, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_1b_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_resolved_1b_mumu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_1b_nElectrons_emu", op.rng_len(self.tightElectrons), DL_resolved_1b_emu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_2b_nElectrons_ee", op.rng_len(self.tightElectrons), DL_resolved_2b_ee, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_2b_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_resolved_2b_mumu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_2b_nElectrons_emu", op.rng_len(self.tightElectrons), DL_resolved_2b_emu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),

            # number of muons
            Plot.make1D("DL_resolved_1b_nMuons_ee", op.rng_len(self.tightMuons), DL_resolved_1b_ee, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_1b_nMuons_mumu", op.rng_len(self.tightMuons), DL_resolved_1b_mumu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_1b_nMuons_emu", op.rng_len(self.tightMuons), DL_resolved_1b_emu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_2b_nMuons_ee", op.rng_len(self.tightMuons), DL_resolved_2b_ee, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_2b_nMuons_mumu", op.rng_len(self.tightMuons), DL_resolved_2b_mumu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
            Plot.make1D("DL_resolved_2b_nMuons_emu", op.rng_len(self.tightMuons), DL_resolved_2b_emu, EqBin(
                3, 0, 3), title="N(el)", xTitle="Number of electrons"),
        ])

        return plots
