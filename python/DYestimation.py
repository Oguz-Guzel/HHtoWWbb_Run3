
from bamboo.plots import Plot, Skim
from bamboo.plots import EquidistantBinning as EqBin
from bamboo import treefunctions as op

from baseAnalysis import NanoBaseHHWWbb
from selections import makeDLSelection

from definitions import ml_input_features
from utils import labeler


class DYestimation(NanoBaseHHWWbb):
    """ Class to create control plots, cutflow reports and skims"""

    def __init__(self, args):
        super().__init__(args)
        self.channel = self.args.channel
        self.mvaModels = None

    def definePlots(self, tree, noSel, sample=None, sampleCfg=None):
        plots = []

        # cutflow report
        plots.append(self.yields)

        # get DL selections
        [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu, DL_resolved_ee, DL_resolved_mumu, DL_resolved_emu] = makeDLSelection(
            self, noSel, tree, sample, DYControlRegion=True)

        Z_peak_selections = [DL_boosted_ee, DL_boosted_mumu,
                             DL_boosted_emu, DL_resolved_ee, DL_resolved_mumu, DL_resolved_emu]

        # labels on plots
        DLboostedEE_label = labeler('DL boosted EE')
        DLboostedMuMu_label = labeler('DL boosted MuMu')
        DLboostedEMU_label = labeler('DL boosted EMu')

        DLresolved_1b_EE_label = labeler('DL resolved 1b EE')
        DLresolved_1b_MuMu_label = labeler('DL resolved 1b MuMu')
        DLresolved_1b_EMu_label = labeler('DL resolved 1b EMu')

        ml_vars = {
            "event_no": tree.event,
            "weight": noSel.weight,
        }
        l1, l2, j1, j2, met, _ = ml_input_features(self)

        ml_vars = ml_vars | l1 | l2 | j1 | j2 | met

        # add skims that hold variables for the ML model
        for sel in Z_peak_selections:
            plots.append(
                Skim(sel.name+"_ml_vars", ml_vars, sel)
            )

        if self.channel == 'DL':
            plots.extend([
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
                Plot.make1D("DL_boosted_InvM_ee", op.invariant_mass(self.tightElectrons[0].p4, self.tightElectrons[1].p4), DL_boosted_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_mumu", op.invariant_mass(self.tightMuons[0].p4, self.tightMuons[1].p4), DL_boosted_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons  (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_emu", op.invariant_mass(self.tightElectrons[0].p4, self.tightMuons[0].p4), DL_boosted_emu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair  (GeV/c^{2})", plotopts=DLboostedEMU_label),

                # pt of the di-lepton
                Plot.make1D("DL_boosted_dileptonPt_ee", op.sum(self.tightElectrons[0].pt, self.tightElectrons[1].pt), DL_boosted_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electrons (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_dileptonPt_mumu", op.sum(self.tightMuons[0].pt, self.tightMuons[1].pt), DL_boosted_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of muons  (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_dileptonPt_emu", op.sum(self.tightElectrons[0].pt, self.tightMuons[0].pt), DL_boosted_emu, EqBin(
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
                Plot.make1D("DL_boosted_totalCharge_ee", op.sum(self.tightElectrons[0].charge, self.tightElectrons[1].charge), DL_boosted_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_totalCharge_mumu", op.sum(self.tightMuons[0].charge, self.tightMuons[1].charge), DL_boosted_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons ", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_totalCharge_emu", op.sum(self.tightElectrons[0].charge, self.tightMuons[0].charge), DL_boosted_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair ", plotopts=DLboostedEMU_label),

                # invariant mass of subjets
                Plot.make1D("DL_boosted_InvM_jj_ee", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_ee, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_InvM_jj_mumu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_mumu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_InvM_jj_emu", op.invariant_mass(self.ak8BJets[0].subJet1.p4, self.ak8BJets[0].subJet2.p4), DL_boosted_emu, EqBin(
                    100, 0., 200.), title="InvM(jj)", xTitle="Invariant Mass of sub-jets (GeV/c^{2})", plotopts=DLboostedEMU_label),

                # leading lepton pt
                Plot.make1D("DL_boosted_leadingLepton_pt_ee", self.tightElectrons[0].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_mumu", self.tightMuons[0].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_pt_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightElectrons[0].pt, self.tightMuons[0].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLboostedEMU_label),
                Plot.make1D("DL_boosted_electron_pt_emu", self.tightElectrons[0].pt, DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading electron (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_muon_pt_emu", self.tightMuons[0].pt, DL_boosted_emu, EqBin(
                    100, 0., 300.), title="leadingLeptonPt", xTitle="p_{T} of the leading muon (GeV/c)", plotopts=DLboostedMuMu_label),
                # sub-leading lepton pt
                Plot.make1D("DL_boosted_subleadingLepton_pt_ee", self.tightElectrons[1].pt, DL_boosted_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_mumu", self.tightMuons[1].pt, DL_boosted_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_pt_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightMuons[0].pt, self.tightElectrons[0].pt), DL_boosted_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLboostedEMU_label),

                # leading lepton eta
                Plot.make1D("DL_boosted_leadingLepton_eta_ee", self.tightElectrons[0].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_mumu", self.tightMuons[0].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_leadingLepton_eta_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightElectrons[0].eta, self.tightMuons[0].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLboostedEMU_label),
                Plot.make1D("DL_boosted_electron_eta_emu", self.tightElectrons[0].eta, DL_boosted_emu, EqBin(
                    30, -3, 3), title="leadingLeptonEta", xTitle="eta of the electron (GeV/c)", plotopts=DLboostedEMU_label),
                Plot.make1D("DL_boosted_muon_eta_emu", self.tightMuons[0].eta, DL_boosted_emu, EqBin(
                    30, -3, 3), title="leadingLeptonEta", xTitle="eta of the muon (GeV/c)", plotopts=DLboostedEMU_label),

                # sub-leading lepton eta
                Plot.make1D("DL_boosted_subleadingLepton_eta_ee", self.tightElectrons[1].eta, DL_boosted_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_mumu", self.tightMuons[1].eta, DL_boosted_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_subleadingLepton_eta_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightMuons[0].eta, self.tightElectrons[0].eta), DL_boosted_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLboostedEMU_label),

                # DR between leading and sub-leading lepton
                Plot.make1D("DL_boosted_DR_leptons_ee", op.deltaR(self.tightElectrons[0].p4, self.tightElectrons[1].p4), DL_boosted_ee, EqBin(
                    70, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leptons_mumu", op.deltaR(self.tightMuons[0].p4, self.tightMuons[1].p4), DL_boosted_mumu, EqBin(
                    70, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leptons_emu", op.deltaR(self.tightElectrons[0].p4, self.tightMuons[0].p4), DL_boosted_emu, EqBin(
                    70, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLboostedEMU_label),

                # DR between leading lepton and ak8 jet
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_ee", op.deltaR(self.tightElectrons[0].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_mumu", op.deltaR(self.tightMuons[0].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_leadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightElectrons[0].p4, self.tightMuons[0].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLboostedEMU_label),

                # DR between subleading lepton and ak8 jet
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_ee", op.deltaR(self.tightElectrons[1].p4, self.ak8Jets[0].p4), DL_boosted_ee, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedEE_label),
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_mumu", op.deltaR(self.tightMuons[1].p4, self.ak8Jets[0].p4), DL_boosted_mumu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedMuMu_label),
                Plot.make1D("DL_boosted_DR_subleadingleptonANDak8bjet_emu", op.deltaR(op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightMuons[0].p4, self.tightElectrons[0].p4), self.ak8Jets[0].p4), DL_boosted_emu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLboostedEMU_label),

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
                Plot.make1D("DL_resolved_1b_nAK4bJets_ee", op.rng_len(self.ak4BJets), DL_resolved_ee, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_nAK4bJets_mumu", op.rng_len(self.ak4BJets), DL_resolved_mumu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_nAK4bJets_emu", op.rng_len(self.ak4BJets), DL_resolved_emu, EqBin(
                    10, 0., 10), xTitle="Number of AK4 B-jets", plotopts=DLresolved_1b_EMu_label),
                # ak4 bjet pt
                Plot.make1D("DL_resolved_1b_ak4BJet_pt_ee", self.ak4BJets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_ak4BJet_pt_mumu", self.ak4BJets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_ak4BJet_pt_emu", self.ak4BJets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="AK4 B-jet jet p_{T} (GeV/c)", plotopts=DLresolved_1b_EMu_label),
                # ak4 bjet eta
                Plot.make1D("DL_resolved_1b_ak4BJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_ak4BJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_ak4BJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="pT(j1)", xTitle="AK4 B-jet \eta", plotopts=DLresolved_1b_EMu_label),
                # number of ak4 jets
                Plot.make1D("DL_resolved_1b_nak4Jets_ee", op.rng_len(self.ak4Jets), DL_resolved_ee, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_nak4Jets_mumu", op.rng_len(self.ak4Jets), DL_resolved_mumu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_nak4Jets_emu", op.rng_len(self.ak4Jets), DL_resolved_emu, EqBin(
                    15, 0., 15.), xTitle="Number of AK4 jets", plotopts=DLresolved_1b_EMu_label),
                # leading jet pt
                Plot.make1D("DL_resolved_1b_leadingJet_pt_ee", self.ak4Jets[0].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_leadingJet_pt_mumu", self.ak4Jets[0].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_leadingJet_pt_emu", self.ak4Jets[0].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j1)", xTitle="Leading jet p_{T} (GeV/c)", plotopts=DLresolved_1b_EMu_label),
                # leading jet eta
                Plot.make1D("DL_resolved_1b_leadingJet_eta_ee", self.ak4BJets[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_leadingJet_eta_mumu", self.ak4BJets[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_leadingJet_eta_emu", self.ak4BJets[0].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j1)", xTitle="Leading jet \eta", plotopts=DLresolved_1b_EMu_label),
                # btagging score of the jet
                Plot.make1D("DL_resolved_1b_jet_btagScore_ee", self.ak4BJets[0].btagPNetB, DL_resolved_ee, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_jet_btagScore_mumu", self.ak4BJets[0].btagPNetB, DL_resolved_mumu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_jet_btagScore_emu", self.ak4BJets[0].btagPNetB, DL_resolved_emu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the leading jet", plotopts=DLresolved_1b_EMu_label),
                # sub-leading jet pt
                Plot.make1D("DL_resolved_1b_subleadingJet_pt_ee", self.ak4Jets[1].pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_subleadingJet_pt_mumu", self.ak4Jets[1].pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_subleadingJet_pt_emu", self.ak4Jets[1].pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="Sub-leading jet p_{T} (GeV/c)", plotopts=DLresolved_1b_EMu_label),
                # sub-leading jet eta
                Plot.make1D("DL_resolved_1b_subleadingJet_eta_ee", self.ak4Jets[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_subleadingJet_eta_mumu", self.ak4Jets[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_subleadingJet_eta_emu", self.ak4Jets[1].eta, DL_resolved_emu, EqBin(
                    30, -3, 3), title="eta(j2)", xTitle="Sub-leading jet \eta", plotopts=DLresolved_1b_EMu_label),
                # btagging score of the jet
                Plot.make1D("DL_resolved_1b_subleadingJet_btagScore_ee", self.ak4Jets[1].btagPNetB, DL_resolved_ee, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_subleadingJet_btagScore_mumu", self.ak4Jets[1].btagPNetB, DL_resolved_mumu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_subleadingJet_btagScore_emu", self.ak4Jets[1].btagPNetB, DL_resolved_emu, EqBin(
                    100, 0, 1), title="btagScore(ak4jet)", xTitle="B-tagging score of the subleadingJet jet", plotopts=DLresolved_1b_EMu_label),
                # DR between leading and sub-leading jet
                Plot.make1D("DL_resolved_1b_DR_jets_ee", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_ee, EqBin(
                    70, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_DR_jets_mumu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_mumu, EqBin(
                    70, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_DR_jets_emu", op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4), DL_resolved_emu, EqBin(
                    70, 0, 7), title="DR(j1,j2)", xTitle="Angular distance between jets", plotopts=DLresolved_1b_EMu_label),
                # Invariant mass of leptons
                Plot.make1D("DL_resolved_1b_InvM_ee", op.invariant_mass(self.tightElectrons[0].p4, self.tightElectrons[1].p4), DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electrons (GeV/c^{2})", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_InvM_mumu", op.invariant_mass(self.tightMuons[0].p4, self.tightMuons[1].p4), DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of muons (GeV/c^{2})", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_InvM_emu", op.invariant_mass(self.tightElectrons[0].p4, self.tightMuons[0].p4), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="Invariant Mass of electron-muon pair (GeV/c^{2})", plotopts=DLresolved_1b_EMu_label),
                # pt of the di-lepton
                Plot.make1D("DL_resolved_1b_dileptonPt_ee", op.sum(self.tightElectrons[0].pt, self.tightElectrons[1].pt), DL_resolved_ee, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electrons (GeV/c^{2})", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_dileptonPt_mumu", op.sum(self.tightMuons[0].pt, self.tightMuons[1].pt), DL_resolved_mumu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of muons  (GeV/c^{2})", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_dileptonPt_emu", op.sum(self.tightElectrons[0].pt, self.tightMuons[0].pt), DL_resolved_emu, EqBin(
                    60, 0., 300.), title="InvM(ll)", xTitle="P_{T} of electron-muon pair  (GeV/c^{2})", plotopts=DLresolved_1b_EMu_label),
                # MET pt
                Plot.make1D("DL_resolved_1b_MET_pt_ee", tree.MET.pt, DL_resolved_ee, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_MET_pt_mumu", tree.MET.pt, DL_resolved_mumu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_MET_pt_emu", tree.MET.pt, DL_resolved_emu, EqBin(
                    100, 0, 500), title="pT(j2)", xTitle="MET p_{T} (GeV/c)", plotopts=DLresolved_1b_EMu_label),
                # MET phi
                Plot.make1D("DL_resolved_1b_MET_phi_ee", tree.MET.phi, DL_resolved_ee, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_MET_phi_mumu", tree.MET.phi, DL_resolved_mumu, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_MET_phi_emu", tree.MET.phi, DL_resolved_emu, EqBin(
                    7, -3.5, 3.5), title="pT(j2)", xTitle="MET phi (GeV/c)", plotopts=DLresolved_1b_EMu_label),
                # total charge of leptons
                Plot.make1D("DL_resolved_1b_totalCharge_ee", op.sum(self.tightElectrons[0].charge, self.tightElectrons[1].charge), DL_resolved_ee, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electrons", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_totalCharge_mumu", op.sum(self.tightMuons[0].charge, self.tightMuons[1].charge), DL_resolved_mumu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of muons", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_totalCharge_emu", op.sum(self.tightElectrons[0].charge, self.tightMuons[0].charge), DL_resolved_emu, EqBin(
                    5, -2.5, 2.5), title="total charge", xTitle="Total charge of electron-muon pair", plotopts=DLresolved_1b_EMu_label),
                # leading lepton pt
                Plot.make1D("DL_resolved_1b_leadingLepton_pt_ee", self.tightElectrons[0].pt, DL_resolved_ee, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_leadingLepton_pt_mumu", self.tightMuons[0].pt, DL_resolved_mumu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_leadingLepton_pt_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightElectrons[0].pt, self.tightMuons[0].pt), DL_resolved_emu, EqBin(
                    100, 0., 300.), title="InvM(ll)", xTitle="p_{T} of the leading lepton (GeV/c)", plotopts=DLresolved_1b_EMu_label),
                # sub-leading lepton pt
                Plot.make1D("DL_resolved_1b_subleadingLepton_pt_ee", self.tightElectrons[1].pt, DL_resolved_ee, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_subleadingLepton_pt_mumu", self.tightMuons[1].pt, DL_resolved_mumu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_subleadingLepton_pt_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightMuons[0].pt, self.tightElectrons[0].pt), DL_resolved_emu, EqBin(
                    50, 0., 200.), title="InvM(ll)", xTitle="p_{T} of the sub-leading lepton (GeV/c)", plotopts=DLresolved_1b_EMu_label),
                # leading lepton eta
                Plot.make1D("DL_resolved_1b_leadingLepton_eta_ee", self.tightElectrons[0].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_leadingLepton_eta_mumu", self.tightMuons[0].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_leadingLepton_eta_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightElectrons[0].eta, self.tightMuons[0].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the leading lepton", plotopts=DLresolved_1b_EMu_label),
                # sub-leading lepton eta
                Plot.make1D("DL_resolved_1b_subleadingLepton_eta_ee", self.tightElectrons[1].eta, DL_resolved_ee, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_subleadingLepton_eta_mumu", self.tightMuons[1].eta, DL_resolved_mumu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_subleadingLepton_eta_emu", op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightMuons[0].eta, self.tightElectrons[0].eta), DL_resolved_emu, EqBin(
                    30, -3, 3), title="InvM(ll)", xTitle="\eta of the sub-leading lepton", plotopts=DLresolved_1b_EMu_label),
                # DR between leading and sub-leading lepton
                Plot.make1D("DL_resolved_1b_DR_leptons_ee", op.deltaR(self.tightElectrons[0].p4, self.tightElectrons[1].p4), DL_resolved_ee, EqBin(
                    70, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_DR_leptons_mumu", op.deltaR(self.tightMuons[0].p4, self.tightMuons[1].p4), DL_resolved_mumu, EqBin(
                    70, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_DR_leptons_emu", op.deltaR(self.tightElectrons[0].p4, self.tightMuons[0].p4), DL_resolved_emu, EqBin(
                    70, 0, 7), title="DR(l1,l2)", xTitle="Angular distance between leptons", plotopts=DLresolved_1b_EMu_label),
                # DR between leading lepton and ak4 b jet
                Plot.make1D("DL_resolved_1b_DR_leadingleptonANDak4bjet_ee", op.deltaR(self.tightElectrons[0].p4, self.ak4BJets[0].p4), DL_resolved_ee, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_DR_leadingleptonANDak4bjet_mumu", op.deltaR(self.tightMuons[0].p4, self.ak4BJets[0].p4), DL_resolved_mumu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_DR_leadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightElectrons[0].p4, self.tightMuons[0].p4), self.ak4BJets[0].p4), DL_resolved_emu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(leading-lepton, ak8bjet)", plotopts=DLresolved_1b_EMu_label),
                # DR between sub-leading lepton and ak4 b jet
                Plot.make1D("DL_resolved_1b_DR_subleadingleptonANDak4bjet_ee", op.deltaR(self.tightElectrons[1].p4, self.ak4BJets[0].p4), DL_resolved_ee, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_DR_subleadingleptonANDak4bjet_mumu", op.deltaR(self.tightMuons[1].p4, self.ak4BJets[0].p4), DL_resolved_mumu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_DR_subleadingleptonANDak4bjet_emu", op.deltaR(op.switch((self.tightElectrons[0].pt >= self.tightMuons[0].pt), self.tightMuons[0].p4, self.tightElectrons[0].p4), self.ak4BJets[0].p4), DL_resolved_emu, EqBin(
                    70, 0, 7), title="DR(l1,ak8)", xTitle="\Delta R(subleading-lepton, ak8bjet)", plotopts=DLresolved_1b_EMu_label),
                # number of electrons
                Plot.make1D("DL_resolved_1b_nElectrons_ee", op.rng_len(self.tightElectrons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_nElectrons_mumu", op.rng_len(self.tightElectrons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_nElectrons_emu", op.rng_len(self.tightElectrons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolved_1b_EMu_label),
                # number of muons
                Plot.make1D("DL_resolved_1b_nMuons_ee", op.rng_len(self.tightMuons), DL_resolved_ee, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolved_1b_EE_label),
                Plot.make1D("DL_resolved_1b_nMuons_mumu", op.rng_len(self.tightMuons), DL_resolved_mumu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolved_1b_MuMu_label),
                Plot.make1D("DL_resolved_1b_nMuons_emu", op.rng_len(self.tightMuons), DL_resolved_emu, EqBin(
                    3, 0, 3), title="N(el)", xTitle="Number of electrons", plotopts=DLresolved_1b_EMu_label),])

        return plots
