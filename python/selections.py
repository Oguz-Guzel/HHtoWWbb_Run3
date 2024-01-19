from bamboo import treefunctions as op

# common variables for DL and SL channels


def OSDilepton(dilep): dilep[0].charge != dilep[1].charge


Zmass = 91.1876


def lowMllCut(dileptons): return op.NOT(op.rng_any(
    dileptons, lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))


def outZ(dileptons): return op.NOT(op.rng_any(
    dileptons, lambda dilep: op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - Zmass) <= 10.))

# end common variables


def makeDLSelection(self, noSel):
    """Selections for the DL channel
    return the following selections:
    - DL_boosted_ee: boosted selection for ee channel
    - DL_boosted_mumu: boosted selection for mumu channel
    - DL_boosted_emu: boosted selection for emu channel
    - DL_resolved_1b_ee: resolved selection for ee channel with at least one b-tagged ak4 jet
    - DL_resolved_1b_mumu: resolved selection for mumu channel with at least one b-tagged ak4 jet
    - DL_resolved_1b_emu: resolved selection for emu channel with at least one b-tagged ak4 jet
    """

    # Pt cuts : subleading above 15 GeV and leading above 25 GeV
    def ptCutElEl(dilep): return op.AND(
        self.electron_conept[dilep[0].idx] > 15,
        self.electron_conept[dilep[1].idx] > 15,
        op.OR(self.electron_conept[dilep[0].idx] >
              25, self.electron_conept[dilep[1].idx] > 25)
    )

    def ptCutMuMu(dilep): return op.AND(
        self.muon_conept[dilep[0].idx] > 15,
        self.muon_conept[dilep[1].idx] > 15,
        op.OR(self.muon_conept[dilep[0].idx] > 25,
              self.muon_conept[dilep[1].idx] > 25)
    )

    def ptCutElMu(dilep): return op.AND(
        self.electron_conept[dilep[0].idx] > 15,
        self.muon_conept[dilep[1].idx] > 15,
        op.OR(self.electron_conept[dilep[0].idx] >
              25, self.muon_conept[dilep[1].idx] > 25)
    )

    ElElSel = noSel.refine("ElElSel", cut=[
        op.rng_len(self.ElElFakePair) >= 1,
        op.OR(op.rng_len(self.fakeMuons) == 0,
              op.AND(op.rng_len(self.fakeMuons) == 1,
                     self.electron_conept[self.ElElFakePair[0]
                                          [0].idx] > self.muon_conept[self.fakeMuons[0].idx],
                     self.electron_conept[self.ElElFakePair[0][1].idx] > self.muon_conept[self.fakeMuons[0].idx]),
              op.AND(op.rng_len(self.fakeMuons) >= 2,
                     self.electron_conept[self.ElElFakePair[0]
                                          [0].idx] > self.muon_conept[self.fakeMuons[0].idx],
                     self.electron_conept[self.ElElFakePair[0]
                                          [1].idx] > self.muon_conept[self.fakeMuons[0].idx],
                     self.electron_conept[self.ElElFakePair[0]
                                          [0].idx] > self.muon_conept[self.fakeMuons[1].idx],
                     self.electron_conept[self.ElElFakePair[0][1].idx] > self.muon_conept[self.fakeMuons[1].idx]))])

    MuMuSel = noSel.refine("MuMuSel", cut=[
        op.rng_len(self.MuMuFakePair) >= 1,
        op.OR(op.rng_len(self.fakeElectrons) == 0,
              op.AND(op.rng_len(self.fakeElectrons) == 1,
                     self.muon_conept[self.MuMuFakePair[0]
                                      [0].idx] > self.electron_conept[self.fakeElectrons[0].idx],
                     self.muon_conept[self.MuMuFakePair[0][1].idx] > self.electron_conept[self.fakeElectrons[0].idx]),
              op.AND(op.rng_len(self.fakeElectrons) >= 2,
                     self.muon_conept[self.MuMuFakePair[0]
                                      [0].idx] > self.electron_conept[self.fakeElectrons[0].idx],
                     self.muon_conept[self.MuMuFakePair[0]
                                      [1].idx] > self.electron_conept[self.fakeElectrons[0].idx],
                     self.muon_conept[self.MuMuFakePair[0]
                                      [0].idx] > self.electron_conept[self.fakeElectrons[1].idx],
                     self.muon_conept[self.MuMuFakePair[0][1].idx] > self.electron_conept[self.fakeElectrons[1].idx]))])

    ElMuSel = noSel.refine("ElMuSel", cut=[
        op.rng_len(self.ElMuFakePair) >= 1,
        op.OR(op.AND(op.rng_len(self.fakeElectrons) == 1,
                     op.rng_len(self.fakeMuons) == 1),
              op.AND(op.rng_len(self.fakeElectrons) >= 2,
                     op.rng_len(self.fakeMuons) == 1,
                     self.muon_conept[self.ElMuFakePair[0][1].idx] > self.electron_conept[self.fakeElectrons[1].idx]),
              op.AND(op.rng_len(self.fakeMuons) >= 2,
                     op.rng_len(self.fakeElectrons) == 1,
                     self.electron_conept[self.ElMuFakePair[0][0].idx] > self.muon_conept[self.fakeMuons[1].idx]),
              op.AND(op.rng_len(self.fakeElectrons) >= 2,
                     op.rng_len(self.fakeMuons) >= 2,
                     self.muon_conept[self.ElMuFakePair[0]
                                      [1].idx] > self.electron_conept[self.fakeElectrons[1].idx],
                     self.electron_conept[self.ElMuFakePair[0][0].idx] > self.muon_conept[self.fakeMuons[1].idx]))])

    # OS dilepton selections #
    ElElSel = ElElSel.refine(
        'ElElSelOS', cut=[OSDilepton(self.ElElFakePair[0])])
    MuMuSel = MuMuSel.refine(
        'MuMuSelOS', cut=[OSDilepton(self.MuMuFakePair[0])])
    ElMuSel = ElMuSel.refine(
        'ElMuSelOS', cut=[OSDilepton(self.ElMuFakePair[0])])

    # pt cuts #
    ElElSel = ElElSel.refine(
        'ElElSelOSPtCuts', cut=[ptCutElEl(self.ElElFakePair[0])])
    MuMuSel = MuMuSel.refine(
        'MuMuSelOSPtCuts', cut=[ptCutMuMu(self.MuMuFakePair[0])])
    ElMuSel = ElMuSel.refine(
        'ElMuSelOSPtCuts', cut=[ptCutElMu(self.ElMuFakePair[0])])

    # mll cut #
    mllCut = [lowMllCut(self.ElElDileptonPreSel), lowMllCut(
        self.MuMuDileptonPreSel), lowMllCut(self.ElMuDileptonPreSel)]
    ElElSel = ElElSel.refine(
        "ElElSelOSPtCutsPreMllCut", cut=mllCut)
    MuMuSel = MuMuSel.refine(
        "MuMuSelOSPtCutsPreMllCut", cut=mllCut)
    ElMuSel = ElMuSel.refine(
        "ElMuSelOSPtCutsPreMllCut", cut=mllCut)

    # Zveto #
    outZCut = [outZ(self.OSElElDileptonPreSel),
               outZ(self.OSMuMuDileptonPreSel)]
    ElElSel = ElElSel.refine(
        "ElElSelOSPtCutsPreMllCutOutZ", cut=outZCut)
    MuMuSel = MuMuSel.refine(
        "MuMuSelOSPtCutsPreMllCutOutZ", cut=outZCut)
    ElMUSel = ElMuSel.refine(
        "ElMuSelOSPtCutsPreMllCutOutZ", cut=outZCut)

    # tight selection #
    ElElSel = ElElSel.refine("ElElSelOSPtCutsPreMllCutOutZTightSelected", cut=[
        self.tightpair_ElEl(self.ElElFakePair[0]),
        op.rng_len(self.tightElectrons) == 2,
        op.rng_len(self.tightMuons) == 0,
        self.ElElTightPair[0][0].idx == self.ElElFakePair[0][0].idx,
        self.ElElTightPair[0][1].idx == self.ElElFakePair[0][1].idx])

    # MuMuSel = MuMuSel.refine("MuMuSelOSPtCutsPreMllCutOutZTightSelected", cut=[
    #     self.tightpair_MuMu(self.MuMuFakePair[0]),
    #     op.rng_len(self.tightElectrons) == 0,
    #     op.rng_len(self.tightMuons) == 2,
    #     self.MuMuTightPair[0][0].idx == self.MuMuFakePair[0][0].idx,
    #     self.MuMuTightPair[0][1].idx == self.MuMuFakePair[0][1].idx])

    # ElMUSel = ElMuSel.refine("ElMuSelOSPtCutsPreMllCutOutZTightSelected", cut=[
    #     self.tightpair_ElMu(self.ElMuFakePair[0]),
    #     op.rng_len(self.tightElectrons) == 1,
    #     op.rng_len(self.tightMuons) == 1,
    #     self.ElMuTightPair[0][0].idx == self.ElMuFakePair[0][0].idx,
    #     self.ElMuTightPair[0][1].idx == self.ElMuFakePair[0][1].idx])

    self.firstOSElEl = self.ElElTightPair[0]
    self.firstOSMuMu = self.MuMuTightPair[0]
    self.firstOSElMu = self.ElMuTightPair[0]

    # boosted -> at least one b-tagged ak8 jet
    DL_boosted_ee = ElElSel.refine(
        'DL_boosted_ee', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_mumu = MuMuSel.refine(
        'DL_boosted_mumu', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_emu = ElMUSel.refine(
        'DL_boosted_emu', cut=(op.rng_len(self.ak8BJets) >= 1))

    # resolved -> and at least two ak4 jets with at least one b-tagged and no ak8 jets
    DL_resolved_ee = ElElSel.refine('DL_resolved_ee',
                                    cut=(op.AND(op.rng_len(self.ak4Jets) >= 2, op.rng_len(self.ak4BJets) >= 1, op.rng_len(self.ak8Jets) == 0)))

    DL_resolved_mumu = MuMuSel.refine('DL_resolved_mumu',
                                      cut=(op.AND(op.rng_len(self.ak4Jets) >= 2, op.rng_len(self.ak4BJets) >= 1, op.rng_len(self.ak8Jets) == 0)))

    DL_resolved_emu = ElMUSel.refine('DL_resolved_emu',
                                     cut=(op.AND(op.rng_len(self.ak4Jets) >= 2, op.rng_len(self.ak4BJets) >= 1, op.rng_len(self.ak8Jets) == 0)))

    DL_selections = [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
                     DL_resolved_ee, DL_resolved_mumu, DL_resolved_emu]

    return DL_selections


def makeSLSelection(self, noSel):
    """ Selections for the SL channel
    return the following selections:
    - SL_resolved: resolved selection
    - SL_resolved_e: resolved selection for e channel
    - SL_resolved_mu: resolved selection for mu channel
    - SL_boosted: boosted selection
    - SL_boosted_e: boosted selection for e channel
    - SL_boosted_mu: boosted selection for mu channel"""

    def elPtCut(lep): return self.electron_conept[lep[0].idx] > 32.0

    def muPtCut(lep): return self.muon_conept[lep[0].idx] > 25.0

    def tau_h_veto(taus): return op.rng_len(taus) == 0

    # OS loose lepton pairs of same type to be vetoed around Z peak
    ElElLoosePairs = op.combine(
        self.clElectrons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    MuMuLoosePairs = op.combine(
        self.muons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    ElMuLoosePairs = op.combine(
        (self.clElectrons, self.muons), N=2, pred=lambda el, mu: el.charge != mu.charge)

    # Z-veto : reject events with dileptons of same type with mass around Z peak
    outZCut = op.AND(outZ(ElElLoosePairs), outZ(MuMuLoosePairs))

    # low Mll cut : reject events with dilepton mass below 12 GeV
    mllCut = op.AND(lowMllCut(ElElLoosePairs), lowMllCut(
        MuMuLoosePairs), lowMllCut(ElMuLoosePairs))

    OSoutZelelSel = noSel.refine('OSoutZsel', cut=op.AND(
        mllCut, outZCut, tau_h_veto(self.cleanedTaus)))

    SL_resolved = OSoutZelelSel.refine('SL_resolved', cut=[
        op.rng_len(self.ak4Jets) >= 3,
        op.rng_len(self.ak4BJets) >= 1,
        op.rng_len(self.ak8BJets) == 0])

    SL_resolved_e = SL_resolved.refine('SL_resolved_e', cut=[
        elPtCut(self.tightElectrons),
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 0])

    SL_resolved_mu = SL_resolved.refine('SL_resolved_mu', cut=[
        muPtCut(self.tightMuons),
        op.rng_len(self.tightElectrons) == 0,
        op.rng_len(self.tightMuons) == 1])

    SL_boosted = OSoutZelelSel.refine('SL_boosted', cut=[
        op.rng_len(self.ak8BJets) >= 1,
        op.rng_len(self.ak4JetsCleanedFromAk8b) >= 1])

    SL_boosted_e = SL_boosted.refine('SL_boosted_e', cut=[
        elPtCut(self.tightElectrons),
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 0])

    SL_boosted_mu = SL_boosted.refine('SL_boosted_mu', cut=[
        muPtCut(self.tightMuons),
        op.rng_len(self.tightMuons) == 1,
        op.rng_len(self.tightElectrons) == 0])

    SL_selections = [SL_resolved, SL_resolved_e,
                     SL_resolved_mu, SL_boosted, SL_boosted_e, SL_boosted_mu]

    return SL_selections
