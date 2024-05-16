from bamboo import treefunctions as op
from bamboo.plots import CategorizedSelection

# common variables for DL and SL channels

Zmass = 91.1876

def tau_h_veto(taus) -> bool:
    "Veto events with hadronic taus."
    return op.rng_len(taus) == 0

def lowMllCut(dileptons) -> bool:
    " Minimum dilepton invariant mass cut, that's 12 GeV."
    return op.NOT(op.rng_any(
        dileptons, lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))


def outZ(dileptons) -> bool:
    "Reject events with same-flavoured dilepton mass around Z peak."
    return op.NOT(op.rng_any(
        dileptons, lambda dilep: op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - Zmass) <= 10.))

# end common variables


def makeDLSelection(self, sel):
    """Creates a list of selection objects for the Dilepton final state.

    Args:
        sel: the selection to be used as a base for the Dilepton selections.

    Returns:
        A list of selection objects for the Dilepton analysis.
    """

    # Pt cuts : leading above 25 GeV and sub-leading above 15 GeV
    def ptCutSameFlavourPair(dilep) -> bool:
        """Minimum pT cut for the same flavour leptons.
        Leading lepton pT > 25 GeV and sub-leading lepton pT > 15 GeV.
        There is no need to check for the opposite order since we're taking
        the objects from the same collection that is already sorted by pt.
        """
        return op.AND(
            dilep[0].pt > 25,
            dilep[1].pt > 15,
        )

    def ptCutDifferentFlavourPair(dilep) -> bool:
        """Minimum pT cut for the different flavour leptons.
        Leading lepton pT > 25 GeV and sub-leading lepton pT > 15 GeV.
        Here is necessary to check for the opposite order since we're taking
        the objects from different collections hence we don't know which one has a higher pT."""
        return op.AND(
            dilep[0].pt > 15,
            dilep[1].pt > 15,
            op.OR(dilep[0].pt > 25, dilep[1].pt > 25)
        )

    # OS loose lepton pairs of same type to be vetoed around Z peak
    elLoosePair = op.combine(
        self.clElectrons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    muLoosePair = op.combine(
        self.preMuons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    emuLoosePair = op.combine(
        (self.clElectrons, self.preMuons), N=2, pred=lambda el, mu: el.charge != mu.charge)

    # OS tight dilepton collections
    elTightPair = op.combine(self.tightElectrons, N=2,
                             pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    muTightPair = op.combine(self.tightMuons, N=2,
                             pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    emuTightPair = op.combine((self.tightElectrons, self.tightMuons),
                              N=2, pred=lambda el, mu: el.charge != mu.charge)

    # the actual leptons that will be used for the analysis
    self.firstElTightPair = elTightPair[0]
    self.firstMuTightPair = muTightPair[0]
    self.firstEmuTightPair = emuTightPair[0]
    
    # tau veto
    tauVetoSel = sel.refine('tauVetoSel', cut=tau_h_veto(self.cleanedTaus))

    # minimum pT cut : at least one lepton pair with leading lepton above 25 GeV
    elPairMinPtSel = tauVetoSel.refine(
        'elPairMinPtSel', cut=[ptCutSameFlavourPair(self.firstElTightPair)])
    muPairMinPtSel = tauVetoSel.refine(
        'muPairMinPtSel', cut=[ptCutSameFlavourPair(self.firstMuTightPair)])
    emuPairMinPtSel = tauVetoSel.refine(
        'emuPairMinPtSel', cut=[ptCutDifferentFlavourPair(self.firstEmuTightPair)])

    # low Mll cut : reject events with dilepton mass below 12 GeV
    mllCut = op.AND(lowMllCut(elLoosePair), lowMllCut(
        muLoosePair), lowMllCut(emuLoosePair))

    # Z-veto : reject events with same flavour-lepton pair with mass around Z peak.
    outZCut = op.AND(outZ(elLoosePair), outZ(muLoosePair))

    outZelPairSel = elPairMinPtSel.refine(
        'outZelPairSel', cut=op.AND(mllCut, outZCut))
    outZmuPairSel = muPairMinPtSel.refine(
        'outZmuPairSel', cut=op.AND(mllCut, outZCut))
    outZemuPairSel = emuPairMinPtSel.refine(
        'outZemuPairSel', cut=op.AND(mllCut, outZCut))

    # di-lepton multiplicity cut
    elPairMultiplicitySel = outZelPairSel.refine('elPairMultiplicitySel', cut=[op.AND(
        op.rng_len(elTightPair) == 1,
        op.rng_len(muTightPair) == 0,
        op.rng_len(emuTightPair) == 0
    )])
    mumuPairMultiplicitySel = outZmuPairSel.refine('mumuPairMultiplicitySel', cut=[op.AND(
        op.rng_len(elTightPair) == 0,
        op.rng_len(muTightPair) == 1,
        op.rng_len(emuTightPair) == 0
    )])
    emuPairMultiplicitySel = outZemuPairSel.refine('emuPairMultiplicitySel', cut=[op.AND(
        op.rng_len(elTightPair) == 0,
        op.rng_len(muTightPair) == 0,
        op.rng_len(emuTightPair) == 1,
    )])

    # boosted -> at least one b-tagged ak8 jet
    DL_boosted_ee = elPairMultiplicitySel.refine(
        'DL_boosted_ee', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_mumu = mumuPairMultiplicitySel.refine(
        'DL_boosted_mumu', cut=(op.rng_len(self.ak8BJets) >= 1))
    DL_boosted_emu = emuPairMultiplicitySel.refine(
        'DL_boosted_emu', cut=(op.rng_len(self.ak8BJets) >= 1))

    # resolved -> and at least two ak4 jets with 1 or at least 2 b-tagged jets and no ak8 jets
    DL_resolved_1b_ee = elPairMultiplicitySel.refine(
        'DL_resolved_1b_ee',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) == 1,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_resolved_2b_ee = elPairMultiplicitySel.refine(
        'DL_resolved_2b_ee',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) >= 2,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_resolved_1b_mumu = mumuPairMultiplicitySel.refine(
        'DL_resolved_1b_mumu',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) == 1,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_resolved_2b_mumu = mumuPairMultiplicitySel.refine(
        'DL_resolved_2b_mumu',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) >= 2,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_resolved_1b_emu = emuPairMultiplicitySel.refine(
        'DL_resolved_1b_emu',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) == 1,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_resolved_2b_emu = emuPairMultiplicitySel.refine(
        'DL_resolved_2b_emu',
        cut=(op.AND(
            op.rng_len(self.ak4Jets) >= 2,
            op.rng_len(self.ak4BJets) >= 2,
            op.rng_len(self.ak8Jets) == 0))
    )

    DL_selections = [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
                     DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
                     DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu]

    return DL_selections


def makeSLSelection(self, sel):
    """Creates a list of selection objects for the Single Lepton final state.

    Args:
        sel: the selection to be used as a base for the Single Lepton selections.

    Returns:
        A list of selection objects for the Single Lepton analysis.
    """

    def elPtCut(lep): return lep[0].pt > 15.0

    def muPtCut(lep): return lep[0].pt > 15.0

    # OS loose lepton pairs of same type to be vetoed around Z peak
    ElElLooseSel = op.combine(
        self.clElectrons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    MuMuLooseSel = op.combine(
        self.preMuons, N=2, pred=lambda lep1, lep2: lep1.charge != lep2.charge)
    ElMuLooseSel = op.combine(
        (self.clElectrons, self.preMuons), N=2, pred=lambda el, mu: el.charge != mu.charge)

    # Z-veto : reject events with dileptons of same type with mass around Z peak
    outZCut = op.AND(outZ(ElElLooseSel), outZ(MuMuLooseSel))

    # low Mll cut : reject events with dilepton mass below 12 GeV
    mllCut = op.AND(lowMllCut(ElElLooseSel), lowMllCut(
        MuMuLooseSel), lowMllCut(ElMuLooseSel))

    OSoutZelelSel = noSel.refine('OSoutZsel', cut=op.AND(
        mllCut, outZCut, tau_h_veto(self.cleanedTaus)))

    SL_resolved_pre = OSoutZelelSel.refine('SL_resolved_pre', cut=[
        op.rng_len(self.ak4Jets) >= 3,
        op.rng_len(self.ak4BJets) >= 1,
        op.rng_len(self.ak8BJets) == 0])

    SL_resolved_1b_pre = SL_resolved_pre.refine('SL_resolved_1b_pre', cut=[
        op.rng_len(self.ak4BJets) == 1])

    SL_resolved_2b_pre = SL_resolved_pre.refine('SL_resolved_2b_pre', cut=[
        op.rng_len(self.ak4BJets) > 1])

    SL_resolved_1b_e = SL_resolved_1b_pre.refine('SL_resolved_1b_e', cut=[
        elPtCut(self.tightElectrons),
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 0])

    SL_resolved_2b_e = SL_resolved_2b_pre.refine('SL_resolved_2b_e', cut=[
        elPtCut(self.tightElectrons),
        op.rng_len(self.tightElectrons) == 1,
        op.rng_len(self.tightMuons) == 0])

    SL_resolved_1b_mu = SL_resolved_1b_pre.refine('SL_resolved_1b_mu', cut=[
        muPtCut(self.tightMuons),
        op.rng_len(self.tightElectrons) == 0,
        op.rng_len(self.tightMuons) == 1])

    SL_resolved_2b_mu = SL_resolved_2b_pre.refine('SL_resolved_2b_mu', cut=[
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
    
    SL_e = CategorizedSelection("SL_e", categories={
        "resolved_1b": (SL_resolved_1b_e, lambda event: True),
        "resolved_2b": (SL_resolved_2b_e, lambda event: True),
        "boosted": (SL_boosted_e, lambda event: True)
    })

    SL_selections = [SL_resolved_pre, SL_resolved_1b_e, SL_resolved_2b_e,
                     SL_resolved_1b_mu, SL_resolved_2b_mu, SL_boosted, SL_boosted_e, SL_boosted_mu, SL_e]

    return SL_selections
