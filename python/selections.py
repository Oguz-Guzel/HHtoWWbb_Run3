from bamboo import treefunctions as op

import definitions as defs
from scalefactors import ScaleFactors as sf

Zmass = 91.1876  # GeV


def lowMllCut(dileptons) -> bool:
    " Minimum dilepton invariant mass cut of 12 GeV."
    return op.NOT(op.rng_any(
        dileptons, lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))


def outZ(dileptons) -> bool:
    "Reject events with same-flavoured dilepton mass around Z peak."
    return op.NOT(op.rng_any(
        dileptons, lambda dilep: op.abs(op.invariant_mass(dilep[0].p4, dilep[1].p4) - Zmass) < 10.))


# lepton Pt cuts : leading above 25 GeV and sub-leading above 15 GeV


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


def makeDLSelection(self, sel, tree, sample, apply_btagReweight=True):
    """Creates a list of selection objects for the Dilepton final state.

    Args:
        sel: the selection to be used as a base for the Dilepton selections.

    Returns:
        A list of selection objects for the Dilepton analysis.
    """

    # call defined objects
    defs.defineObjects(self, tree)

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

    # the actual lepton pairs that will be used for the analysis
    self.firstElTightPair = elTightPair[0]
    self.firstMuTightPair = muTightPair[0]
    self.firstEmuTightPair = emuTightPair[0]

    # top pT reweighting
    # defining the genPart branch here since the it's only available for MC samples
    # and we don't want to pass the tree object to the top pt reweight method
    genPartBranch = tree.GenPart if self.is_MC else None
    sel = sf.top_pT_reweight(self, genPartBranch, sel, sample)

    # Noise filters
    sel = sf.NoiseFilters(self, tree.Flag, sel)

    # minimum pT cut : at least one lepton pair with leading lepton above 25 GeV
    elPairMinPtSel = sel.refine(
        'elPairMinPtSel', cut=[ptCutSameFlavourPair(self.firstElTightPair)])
    muPairMinPtSel = sel.refine(
        'muPairMinPtSel', cut=[ptCutSameFlavourPair(self.firstMuTightPair)])
    emuPairMinPtSel = sel.refine(
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
    muPairMultiplicitySel = outZmuPairSel.refine('muPairMultiplicitySel', cut=[op.AND(
        op.rng_len(elTightPair) == 0,
        op.rng_len(muTightPair) == 1,
        op.rng_len(emuTightPair) == 0
    )])
    emuPairMultiplicitySel = outZemuPairSel.refine('emuPairMultiplicitySel', cut=[op.AND(
        op.rng_len(elTightPair) == 0,
        op.rng_len(muTightPair) == 0,
        op.rng_len(emuTightPair) == 1,
    )])

    # lepton scale factors
    muPairMultiplicitySel = sf.muonSF(self, muPairMultiplicitySel)
    elPairMultiplicitySel = sf.electronSF(self, elPairMultiplicitySel)
    emuPairMultiplicitySel = sf.muonSF(self, emuPairMultiplicitySel)
    emuPairMultiplicitySel = sf.electronSF(self, emuPairMultiplicitySel)

    # boosted pre-final state selections for btag reweighting
    DL_boosted_pre_ee = elPairMultiplicitySel.refine(
        'DL_boosted_pre_ee', cut=op.rng_len(self.ak8Jets) >= 1)
    DL_boosted_pre_mumu = muPairMultiplicitySel.refine(
        'DL_boosted_pre_mumu', cut=op.rng_len(self.ak8Jets) >= 1)
    DL_boosted_pre_emu = emuPairMultiplicitySel.refine(
        'DL_boosted_pre_emu', cut=op.rng_len(self.ak8Jets) >= 1)

    self.yields.add(DL_boosted_pre_ee, 'DL boosted pre ee')
    self.yields.add(DL_boosted_pre_mumu, 'DL boosted pre mumu')
    self.yields.add(DL_boosted_pre_emu, 'DL boosted pre emu')

    # btagging sf and reweighting for boosted to be done before
    # any b-tagged jet selection
    DL_boosted_pre_ee_btagSF = sf.btagSF(
        self, DL_boosted_pre_ee, self.ak8Jets, jet_tagger="particleNet_XbbVsQCD", apply_reweighting=apply_btagReweight)
    DL_boosted_pre_mumu_btagSF = sf.btagSF(
        self, DL_boosted_pre_mumu, self.ak8Jets, jet_tagger="particleNet_XbbVsQCD", apply_reweighting=apply_btagReweight)
    DL_boosted_pre_emu_btagSF = sf.btagSF(
        self, DL_boosted_pre_emu, self.ak8Jets, jet_tagger="particleNet_XbbVsQCD", apply_reweighting=apply_btagReweight)

    # boosted -> at least one b-tagged ak8 jet
    DL_boosted_ee = DL_boosted_pre_ee_btagSF.refine(
        'DL_boosted_ee', cut=op.rng_len(self.ak8BJets) >= 1)
    DL_boosted_mumu = DL_boosted_pre_mumu_btagSF.refine(
        'DL_boosted_mumu', cut=op.rng_len(self.ak8BJets) >= 1)
    DL_boosted_emu = DL_boosted_pre_emu_btagSF.refine(
        'DL_boosted_emu', cut=op.rng_len(self.ak8BJets) >= 1)

    # resolved pre-final state selections for btag reweighting
    DL_resolved_pre_ee = elPairMultiplicitySel.refine(
        'DL_resolved_pre_ee', cut=[op.rng_len(self.ak4Jets) >= 2])
    DL_resolved_pre_mumu = muPairMultiplicitySel.refine(
        'DL_resolved_pre_mumu', cut=[op.rng_len(self.ak4Jets) >= 2])
    DL_resolved_pre_emu = emuPairMultiplicitySel.refine(
        'DL_resolved_pre_emu', cut=[op.rng_len(self.ak4Jets) >= 2])

    self.yields.add(DL_resolved_pre_ee, 'DL resolved pre ee')
    self.yields.add(DL_resolved_pre_mumu, 'DL resolved pre mumu')
    self.yields.add(DL_resolved_pre_emu, 'DL resolved pre emu')

    # btagging sf for resolved.
    # default jet_tagger is btagPNetB, hence we're using that here
    DL_resolved_pre_ee_btagSF = sf.btagSF(
        self, DL_resolved_pre_ee, self.ak4Jets, apply_reweighting=apply_btagReweight)
    DL_resolved_pre_mumu_btagSF = sf.btagSF(
        self, DL_resolved_pre_mumu, self.ak4Jets, apply_reweighting=apply_btagReweight)
    DL_resolved_pre_emu_btagSF = sf.btagSF(
        self, DL_resolved_pre_emu, self.ak4Jets, apply_reweighting=apply_btagReweight)

    # resolved -> and at least two ak4 jets with 1 or at least 2 b-tagged jets and no ak8 jets
    DL_resolved_1b_ee = DL_resolved_pre_ee_btagSF.refine(
        'DL_resolved_1b_ee',
        cut=(op.AND(
            op.rng_len(self.ak4BJets) == 1,
            op.rng_len(self.ak8BJets) == 0))
    )

    DL_resolved_2b_ee = DL_resolved_pre_ee_btagSF.refine(
        'DL_resolved_2b_ee',
        cut=(op.AND(
            op.rng_len(self.ak4BJets) >= 2,
            op.rng_len(self.ak8BJets) == 0))
    )

    DL_resolved_1b_mumu = DL_resolved_pre_mumu_btagSF.refine(
        'DL_resolved_1b_mumu',
        cut=(op.AND(
            op.rng_len(self.ak4BJets) == 1,
            op.rng_len(self.ak8BJets) == 0))
    )

    DL_resolved_2b_mumu = DL_resolved_pre_mumu_btagSF.refine(
        'DL_resolved_2b_mumu',
        cut=(op.AND(
            op.rng_len(self.ak4BJets) >= 2,
            op.rng_len(self.ak8BJets) == 0))
    )

    DL_resolved_1b_emu = DL_resolved_pre_emu_btagSF.refine(
        'DL_resolved_1b_emu',
        cut=(op.AND(
            op.rng_len(self.ak4BJets) == 1,
            op.rng_len(self.ak8BJets) == 0))
    )

    DL_resolved_2b_emu = DL_resolved_pre_emu_btagSF.refine(
        'DL_resolved_2b_emu',
        cut=(op.AND(
            op.rng_len(self.ak4BJets) >= 2,
            op.rng_len(self.ak8BJets) == 0))
    )

    pre_final_state_sels = {
        'DL_boosted_pre_ee': (DL_boosted_pre_ee, DL_boosted_pre_ee_btagSF),
        'DL_boosted_pre_mumu': (DL_boosted_pre_mumu, DL_boosted_pre_mumu_btagSF),
        'DL_boosted_pre_emu': (DL_boosted_pre_emu, DL_boosted_pre_emu_btagSF),
        'DL_resolved_pre_ee': (DL_resolved_pre_ee, DL_resolved_pre_ee_btagSF),
        'DL_resolved_pre_mumu': (DL_resolved_pre_mumu, DL_resolved_pre_mumu_btagSF),
        'DL_resolved_pre_emu': (DL_resolved_pre_emu, DL_resolved_pre_emu_btagSF),
    }

    DL_selections = [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu,
                     DL_resolved_1b_ee, DL_resolved_1b_mumu, DL_resolved_1b_emu,
                     DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu]

    # cutflow reports for the final states
    self.yields.add(DL_boosted_ee, 'DL boosted ee')
    self.yields.add(DL_boosted_mumu, 'DL boosted mumu')
    self.yields.add(DL_boosted_emu, 'DL boosted emu')
    self.yields.add(DL_resolved_1b_ee, 'DL resolved 1b ee')
    self.yields.add(DL_resolved_2b_ee, 'DL resolved 2b ee')
    self.yields.add(DL_resolved_1b_mumu, 'DL resolved 1b mumu')
    self.yields.add(DL_resolved_2b_mumu, 'DL resolved 2b mumu')
    self.yields.add(DL_resolved_1b_emu, 'DL resolved 1b emu')
    self.yields.add(DL_resolved_2b_emu, 'DL resolved 2b emu')

    return DL_selections, pre_final_state_sels
