from bamboo import treefunctions as op

import definitions as defs
from scalefactors import ScaleFactors as sf

Zmass = 91.1876  # GeV


def lowMllCut(lepton_collection) -> bool:
    " Minimum dilepton invariant mass cut of 12 GeV."
    return op.NOT(op.invariant_mass(lepton_collection[0].p4, lepton_collection[1].p4) < 12.)


def outZ(lepton_collection) -> bool:
    "Reject events with same-flavoured dilepton mass around Z peak."
    return op.NOT(op.abs(op.invariant_mass(lepton_collection[0].p4, lepton_collection[1].p4) - Zmass) < 10.)


# lepton Pt cuts : leading above 25 GeV and sub-leading above 15 GeV


def ptCutSameFlavourPair(lepton_collection) -> bool:
    """Minimum pT cut for the same flavour leptons.
    Leading lepton pT > 25 GeV and sub-leading lepton pT > 15 GeV.
    There is no need to check for the opposite order since we're taking
    the objects from the same collection that is already sorted by pt.
    """
    return op.AND(
        lepton_collection[0].pt > 25,
        lepton_collection[1].pt > 15,
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

    # top pT reweighting
    genPartBranch = tree.GenPart if self.is_MC else None
    sel = sf.top_pT_reweight(self, genPartBranch, sel, sample)

    # Noise filters
    sel = sf.NoiseFilters(self, tree.Flag, sel)

    # final states
    elel_sel = sel.refine(
        'eePairSel', cut=[
            ptCutSameFlavourPair(self.tightElectrons),
            op.rng_len(self.tightElectrons) == 2,
            op.rng_len(self.tightMuons) == 0,
            self.tightElectrons[0].charge != self.tightElectrons[1].charge,
            op.NOT(op.invariant_mass(self.tightElectrons[0].p4, self.tightElectrons[1].p4) < 12.),
            op.NOT(op.abs(op.invariant_mass(self.tightElectrons[0].p4, self.tightElectrons[1].p4) - Zmass) < 10.)
            ])
    mumu_sel = sel.refine(
        'mumuPairSel', cut=[
            ptCutSameFlavourPair(self.tightMuons),
            op.rng_len(self.tightMuons) == 2,
            op.rng_len(self.tightElectrons) == 0,
            self.tightMuons[0].charge != self.tightMuons[1].charge,
            op.NOT(op.invariant_mass(self.tightMuons[0].p4, self.tightMuons[1].p4) < 12.),
            op.NOT(op.abs(op.invariant_mass(self.tightMuons[0].p4, self.tightMuons[1].p4) - Zmass) < 10.)
            ])
    elmu_sel = sel.refine(
        'emuPairSel', cut=[
            op.rng_len(self.tightElectrons) == 1,
            op.rng_len(self.tightMuons) == 1,
            self.tightElectrons[0].charge != self.tightMuons[0].charge,
            op.NOT(op.invariant_mass(self.tightElectrons[0].p4, self.tightMuons[0].p4) < 12.),
            ])

    # lepton scale factors
    elel_SF_sel = sf.elelSF(self, elel_sel)
    mumu_SF_sel = sf.mumuSF(self, mumu_sel)
    elmu_SF_sel = sf.elmuSF(self, elmu_sel)

    # DY Z pT reweighting
    elel_SF_sel = sf.Z_pT_reweight_elel(
        self, elel_SF_sel, sample, genPartBranch)
    mumu_SF_sel = sf.Z_pT_reweight_mumu(
        self, mumu_SF_sel, sample, genPartBranch)

    # boosted pre-final state selections for btag reweighting
    DL_boosted_pre_ee = elel_SF_sel.refine(
        'DL_boosted_pre_ee', cut=op.rng_len(self.ak8Jets) >= 1)
    DL_boosted_pre_mumu = mumu_SF_sel.refine(
        'DL_boosted_pre_mumu', cut=op.rng_len(self.ak8Jets) >= 1)
    DL_boosted_pre_emu = elmu_SF_sel.refine(
        'DL_boosted_pre_emu', cut=op.rng_len(self.ak8Jets) >= 1)

    self.yields.add(DL_boosted_pre_ee, 'DL boosted pre ee')
    self.yields.add(DL_boosted_pre_mumu, 'DL boosted pre mumu')
    self.yields.add(DL_boosted_pre_emu, 'DL boosted pre emu')

    # btagging sf and reweighting - to be done before
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
    DL_resolved_pre_ee = elel_SF_sel.refine(
        'DL_resolved_pre_ee', cut=[op.rng_len(self.ak4Jets) >= 2])
    DL_resolved_pre_mumu = mumu_SF_sel.refine(
        'DL_resolved_pre_mumu', cut=[op.rng_len(self.ak4Jets) >= 2])
    DL_resolved_pre_emu = elmu_SF_sel.refine(
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

    DL_VBF_resolved_ee = DL_resolved_pre_ee.refine(
        'DL_VBF_resolved_ee',
        cut=op.rng_len(self.VBFjetPairsResolved) >= 1
    )

    DL_VBF_boosted_ee = DL_boosted_pre_ee.refine(
        'DL_VBF_boosted_ee',
        cut=op.rng_len(self.VBFjetPairsBoosted) >= 1
    )

    DL_VBF_resolved_mumu = DL_resolved_pre_mumu.refine(
        'DL_VBF_resolved_mumu',
        cut=op.rng_len(self.VBFjetPairsResolved) >= 1
    )

    DL_VBF_boosted_mumu = DL_boosted_pre_mumu.refine(
        'DL_VBF_boosted_mumu',
        cut=op.rng_len(self.VBFjetPairsBoosted) >= 1
    )

    DL_VBF_resolved_emu = DL_resolved_pre_emu.refine(
        'DL_VBF_resolved_emu',
        cut=op.rng_len(self.VBFjetPairsResolved) >= 1
    )

    DL_VBF_boosted_emu = DL_boosted_pre_emu.refine(
        'DL_VBF_boosted_emu',
        cut=op.rng_len(self.VBFjetPairsBoosted) >= 1
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
                     DL_resolved_2b_ee, DL_resolved_2b_mumu, DL_resolved_2b_emu,
                     DL_VBF_resolved_ee, DL_VBF_resolved_mumu, DL_VBF_resolved_emu,
                     DL_VBF_boosted_ee, DL_VBF_boosted_mumu, DL_VBF_boosted_emu]

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
    self.yields.add(DL_VBF_resolved_ee, 'DL VBF resolved ee')
    self.yields.add(DL_VBF_boosted_ee, 'DL VBF boosted ee')
    self.yields.add(DL_VBF_resolved_mumu, 'DL VBF resolved mumu')
    self.yields.add(DL_VBF_boosted_mumu, 'DL VBF boosted mumu')
    self.yields.add(DL_VBF_resolved_emu, 'DL VBF resolved emu')
    self.yields.add(DL_VBF_boosted_emu, 'DL VBF boosted emu')

    return pre_final_state_sels if not apply_btagReweight else DL_selections
