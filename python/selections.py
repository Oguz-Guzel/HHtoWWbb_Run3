from bamboo import treefunctions as op

import definitions as defs
from scalefactors import ScaleFactors

Zmass = 90 # it's actually 91.1876 GeV but 90 GeV is used 

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
    the objects from different collections hence we don't know which one has a higher pT.
    """
    return op.AND(
        dilep[0].pt > 15, dilep[1].pt > 15, op.OR(
            dilep[0].pt > 25, dilep[1].pt > 25)
    )


def makeDLSelection(
    analysis,
    sel,
    tree,
    sample,
    btagReweightStudy=False,
    trigger_study=False,
    DYControlRegion=False,
    TTbarControlRegion=False,
):
    """Creates a list of selection objects for the Dilepton final state.

    Args:
        sel: the selection to be used as a base for the Dilepton selections.

    Returns:
        A list of selection objects for the Dilepton analysis.
    """

    # call defined objects
    defs.defineObjects(analysis, tree)

    # call scale factors
    scale_factors = ScaleFactors(analysis)

    # jet veto maps correction
    sel = scale_factors.jet_veto_map(tree, sel)

    # top pT reweighting
    genPartBranch = tree.GenPart if analysis.is_MC else None
    sel = scale_factors.top_pT_reweight(genPartBranch, sel, sample)

    # V+Jets sample stitching, uncomment the following two lines to enable it
    LHEBranch = tree.LHE if (analysis.is_MC and sample.startswith(("DYto2L", "WtoLNu", "Zto2Nu"))) else None
    sel = scale_factors.V_Jets_Stitching(LHEBranch, sel, sample)

    # muR/muF scale uncertainty (envelope from generator weights)
    sel = scale_factors.muRF_scale_weights(tree, sel, sample, analysis.sampleCfg)

    # Noise filters
    sel = scale_factors.NoiseFilters(tree.Flag, sel)

    # lepton selections
    elel_sel = sel.refine(
        "eePairSel",
        cut=[
            ptCutSameFlavourPair(analysis.tightElectrons),
            op.rng_len(analysis.tightElectrons) == 2,
            op.rng_len(analysis.tightMuons) == 0,
            analysis.tightElectrons[0].charge != analysis.tightElectrons[1].charge,
            op.NOT(
                op.invariant_mass(
                    analysis.tightElectrons[0].p4, analysis.tightElectrons[1].p4
                )
                < 12.0
            ),
        ],
    )
    mumu_sel = sel.refine(
        "mumuPairSel",
        cut=[
            ptCutSameFlavourPair(analysis.tightMuons),
            op.rng_len(analysis.tightMuons) == 2,
            op.rng_len(analysis.tightElectrons) == 0,
            analysis.tightMuons[0].charge != analysis.tightMuons[1].charge,
            op.NOT(
                op.invariant_mass(
                    analysis.tightMuons[0].p4, analysis.tightMuons[1].p4)
                < 12.0
            ),
        ],
    )
    elmu_sel = sel.refine(
        "emuPairSel",
        cut=[
            op.rng_len(analysis.tightElectrons) == 1,
            op.rng_len(analysis.tightMuons) == 1,
            analysis.tightElectrons[0].charge != analysis.tightMuons[0].charge,
            op.NOT(
                op.invariant_mass(
                    analysis.tightElectrons[0].p4, analysis.tightMuons[0].p4
                )
                < 12.0
            ),
        ],
    )

    analysis.yields.add(elel_sel, "EE lepton sel")
    analysis.yields.add(mumu_sel, "MuMu lepton sel")
    analysis.yields.add(elmu_sel, "EMu lepton sel")

    if DYControlRegion:
        elel_sel = elel_sel.refine(
            "eePairZpeakSel",
            cut=[
                op.abs(
                    op.invariant_mass(
                        analysis.tightElectrons[0].p4, analysis.tightElectrons[1].p4
                    )
                    - Zmass
                )
                <= 10.0
            ],
        )
        mumu_sel = mumu_sel.refine(
            "mumuPairZpeakSel",
            cut=[
                op.abs(
                    op.invariant_mass(
                        analysis.tightMuons[0].p4, analysis.tightMuons[1].p4
                    )
                    - Zmass
                )
                <= 10.0
            ],
        )
        analysis.yields.add(elel_sel, "EE DY peak sel")
        analysis.yields.add(mumu_sel, "MuMu DY peak sel")
    elif TTbarControlRegion:
        elel_sel = elel_sel.refine(
            "eePairZpeakSel",
            cut=[
                op.invariant_mass(
                    analysis.tightElectrons[0].p4, analysis.tightElectrons[1].p4
                )
                > (Zmass + 10.0)
            ],
        )
        mumu_sel = mumu_sel.refine(
            "mumuPairZpeakSel",
            cut=[
                op.invariant_mass(
                    analysis.tightMuons[0].p4, analysis.tightMuons[1].p4)
                > (Zmass + 10.0)
            ],
        )
        analysis.yields.add(elel_sel, "EE above DY peak sel")
        analysis.yields.add(mumu_sel, "MuMu above DY peak sel")
    else:
        elel_sel = elel_sel.refine(
            "eePairZpeakSel",
            cut=[
                op.invariant_mass(
                    analysis.tightElectrons[0].p4, analysis.tightElectrons[1].p4
                )
                < (Zmass - 10.0)
            ],
        )
        mumu_sel = mumu_sel.refine(
            "mumuPairZpeakSel",
            cut=[
                op.invariant_mass(
                    analysis.tightMuons[0].p4, analysis.tightMuons[1].p4)
                < (Zmass - 10.0)
            ],
        )
        elmu_sel = elmu_sel.refine(
            "emuHigMllSel",
            cut=[
                op.invariant_mass(
                    analysis.tightElectrons[0].p4, analysis.tightMuons[0].p4
                )
                < 100.0 # educated cut
            ],
        )
        analysis.yields.add(elel_sel, "EE below DY peak sel")
        analysis.yields.add(mumu_sel, "MuMu below DY peak sel")
        analysis.yields.add(elmu_sel, "EMu mll below 100 sel")

    # lepton scale factors
    elel_SF_sel = scale_factors.elelSF(elel_sel)
    mumu_SF_sel = scale_factors.mumuSF(mumu_sel)
    elmu_SF_sel = scale_factors.elmuSF(elmu_sel)

    # dilepton trigger scale factors
    elel_SF_sel = scale_factors.dilepton_trg_sf(elel_SF_sel)
    mumu_SF_sel = scale_factors.dilepton_trg_sf(mumu_SF_sel)
    elmu_SF_sel = scale_factors.dilepton_trg_sf(elmu_SF_sel)

    # DY Z pT and recoil correction
    elel_SF_sel = scale_factors.Z_pT_reweight(
        elel_SF_sel, sample, genPartBranch
    )
    mumu_SF_sel = scale_factors.Z_pT_reweight(
        mumu_SF_sel, sample, genPartBranch
    )

    # boosted pre-final state selections for btag reweighting
    DL_boosted_pre_ee = elel_SF_sel.refine(
        "DL_boosted_pre_ee", cut=op.rng_len(analysis.ak8Jets) >= 1
    )
    DL_boosted_pre_mumu = mumu_SF_sel.refine(
        "DL_boosted_pre_mumu", cut=op.rng_len(analysis.ak8Jets) >= 1
    )
    DL_boosted_pre_emu = elmu_SF_sel.refine(
        "DL_boosted_pre_emu", cut=op.rng_len(analysis.ak8Jets) >= 1
    )

    analysis.yields.add(DL_boosted_pre_ee, "DL boosted pre ee")
    analysis.yields.add(DL_boosted_pre_mumu, "DL boosted pre mumu")
    analysis.yields.add(DL_boosted_pre_emu, "DL boosted pre emu")

    # btagging SF and reweighting - to be done before
    # any b-tagged jet selection
    DL_boosted_pre_ee_btagSF = scale_factors.btagSF(
        DL_boosted_pre_ee,
        analysis.ak8Jets,
        jet_tagger="particleNet_XbbVsQCD",
        btagReweightStudy=btagReweightStudy,
    )
    DL_boosted_pre_mumu_btagSF = scale_factors.btagSF(
        DL_boosted_pre_mumu,
        analysis.ak8Jets,
        jet_tagger="particleNet_XbbVsQCD",
        btagReweightStudy=btagReweightStudy,
    )
    DL_boosted_pre_emu_btagSF = scale_factors.btagSF(
        DL_boosted_pre_emu,
        analysis.ak8Jets,
        jet_tagger="particleNet_XbbVsQCD",
        btagReweightStudy=btagReweightStudy,
    )

    if DYControlRegion:
        Z_peak_selections = []
        # no b-jets
        DL_boosted_ee = DL_boosted_pre_ee_btagSF.refine(
            "DL_boosted_ee", cut=op.rng_len(analysis.ak8BJets) == 0
        )
        DL_boosted_mumu = DL_boosted_pre_mumu_btagSF.refine(
            "DL_boosted_mumu", cut=op.rng_len(analysis.ak8BJets) == 0
        )
        DL_boosted_emu = DL_boosted_pre_emu_btagSF.refine(
            "DL_boosted_emu", cut=op.rng_len(analysis.ak8BJets) == 0
        )
        Z_peak_selections.extend(
            [DL_boosted_ee, DL_boosted_mumu, DL_boosted_emu])
    else:
        # boosted -> at least one b-tagged ak8 jet
        DL_boosted_ee = DL_boosted_pre_ee_btagSF.refine(
            "DL_boosted_ee", cut=op.rng_len(analysis.ak8BJets) >= 1
        )
        DL_boosted_mumu = DL_boosted_pre_mumu_btagSF.refine(
            "DL_boosted_mumu", cut=op.rng_len(analysis.ak8BJets) >= 1
        )
        DL_boosted_emu = DL_boosted_pre_emu_btagSF.refine(
            "DL_boosted_emu", cut=op.rng_len(analysis.ak8BJets) >= 1
        )

    # resolved pre-final state selections for btag reweighting
    DL_resolved_pre_ee = elel_SF_sel.refine(
        "DL_resolved_pre_ee", cut=[op.rng_len(analysis.ak4Jets) >= 2]
    )
    DL_resolved_pre_mumu = mumu_SF_sel.refine(
        "DL_resolved_pre_mumu", cut=[op.rng_len(analysis.ak4Jets) >= 2]
    )
    DL_resolved_pre_emu = elmu_SF_sel.refine(
        "DL_resolved_pre_emu", cut=[op.rng_len(analysis.ak4Jets) >= 2]
    )

    analysis.yields.add(DL_resolved_pre_ee, "DL resolved pre ee")
    analysis.yields.add(DL_resolved_pre_mumu, "DL resolved pre mumu")
    analysis.yields.add(DL_resolved_pre_emu, "DL resolved pre emu")

    # btagging scale_factors for resolved.
    # default jet_tagger is btagPNetB, hence we're using that here
    DL_resolved_pre_ee_btagSF = scale_factors.btagSF(
        DL_resolved_pre_ee, analysis.ak4Jets, btagReweightStudy=btagReweightStudy
    )
    DL_resolved_pre_mumu_btagSF = scale_factors.btagSF(
        DL_resolved_pre_mumu, analysis.ak4Jets, btagReweightStudy=btagReweightStudy
    )
    DL_resolved_pre_emu_btagSF = scale_factors.btagSF(
        DL_resolved_pre_emu, analysis.ak4Jets, btagReweightStudy=btagReweightStudy
    )

    if DYControlRegion:
        # no b-jets
        DL_resolved_ee = DL_resolved_pre_ee_btagSF.refine(
            "DL_resolved_ee",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) == 0,
                    op.rng_len(analysis.ak8Jets) == 0,
                )
            ),  # hence no ak8 b-jet
        )
        DL_resolved_mumu = DL_resolved_pre_mumu_btagSF.refine(
            "DL_resolved_mumu",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) == 0,
                    op.rng_len(analysis.ak8Jets) == 0,
                )
            ),
        )
        DL_resolved_emu = DL_resolved_pre_emu_btagSF.refine(
            "DL_resolved_emu",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) == 0,
                    op.rng_len(analysis.ak8Jets) == 0,
                )
            ),
        )
        Z_peak_selections.extend(
            [DL_resolved_ee, DL_resolved_mumu, DL_resolved_emu])
    else:
        # resolved -> and at least two ak4 jets with 1 or at least 2 b-tagged jets and no ak8 jets
        DL_resolved_1b_ee = DL_resolved_pre_ee_btagSF.refine(
            "DL_resolved_1b_ee",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) == 1,
                    op.rng_len(analysis.ak8BJets) == 0,
                )
            ),
        )

        DL_resolved_2b_ee = DL_resolved_pre_ee_btagSF.refine(
            "DL_resolved_2b_ee",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) >= 2,
                    op.rng_len(analysis.ak8BJets) == 0,
                )
            ),
        )

        DL_resolved_1b_mumu = DL_resolved_pre_mumu_btagSF.refine(
            "DL_resolved_1b_mumu",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) == 1,
                    op.rng_len(analysis.ak8BJets) == 0,
                )
            ),
        )

        DL_resolved_2b_mumu = DL_resolved_pre_mumu_btagSF.refine(
            "DL_resolved_2b_mumu",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) >= 2,
                    op.rng_len(analysis.ak8BJets) == 0,
                )
            ),
        )

        DL_resolved_1b_emu = DL_resolved_pre_emu_btagSF.refine(
            "DL_resolved_1b_emu",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) == 1,
                    op.rng_len(analysis.ak8BJets) == 0,
                )
            ),
        )

        DL_resolved_2b_emu = DL_resolved_pre_emu_btagSF.refine(
            "DL_resolved_2b_emu",
            cut=(
                op.AND(
                    op.rng_len(analysis.ak4BJets) >= 2,
                    op.rng_len(analysis.ak8BJets) == 0,
                )
            ),
        )

        DL_VBF_resolved_ee = DL_resolved_pre_ee.refine(
            "DL_VBF_resolved_ee", cut=op.rng_len(analysis.VBFjetPairsResolved) >= 1
        )

        DL_VBF_boosted_ee = DL_boosted_pre_ee.refine(
            "DL_VBF_boosted_ee", cut=op.rng_len(analysis.VBFjetPairsBoosted) >= 1
        )

        DL_VBF_resolved_mumu = DL_resolved_pre_mumu.refine(
            "DL_VBF_resolved_mumu", cut=op.rng_len(analysis.VBFjetPairsResolved) >= 1
        )

        DL_VBF_boosted_mumu = DL_boosted_pre_mumu.refine(
            "DL_VBF_boosted_mumu", cut=op.rng_len(analysis.VBFjetPairsBoosted) >= 1
        )

        DL_VBF_resolved_emu = DL_resolved_pre_emu.refine(
            "DL_VBF_resolved_emu", cut=op.rng_len(analysis.VBFjetPairsResolved) >= 1
        )

        DL_VBF_boosted_emu = DL_boosted_pre_emu.refine(
            "DL_VBF_boosted_emu", cut=op.rng_len(analysis.VBFjetPairsBoosted) >= 1
        )

        pre_final_state_sels = {
            "DL_boosted_pre_ee": (DL_boosted_pre_ee, DL_boosted_pre_ee_btagSF),
            "DL_boosted_pre_mumu": (DL_boosted_pre_mumu, DL_boosted_pre_mumu_btagSF),
            "DL_boosted_pre_emu": (DL_boosted_pre_emu, DL_boosted_pre_emu_btagSF),
            "DL_resolved_pre_ee": (DL_resolved_pre_ee, DL_resolved_pre_ee_btagSF),
            "DL_resolved_pre_mumu": (DL_resolved_pre_mumu, DL_resolved_pre_mumu_btagSF),
            "DL_resolved_pre_emu": (DL_resolved_pre_emu, DL_resolved_pre_emu_btagSF),
        }

        DL_signal_region_selections = [
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
        ]

        # cutflow reports for the final states
        analysis.yields.add(DL_boosted_ee, "DL boosted ee")
        analysis.yields.add(DL_boosted_mumu, "DL boosted mumu")
        analysis.yields.add(DL_boosted_emu, "DL boosted emu")
        analysis.yields.add(DL_resolved_1b_ee, "DL resolved 1b ee")
        analysis.yields.add(DL_resolved_2b_ee, "DL resolved 2b ee")
        analysis.yields.add(DL_resolved_1b_mumu, "DL resolved 1b mumu")
        analysis.yields.add(DL_resolved_2b_mumu, "DL resolved 2b mumu")
        analysis.yields.add(DL_resolved_1b_emu, "DL resolved 1b emu")
        analysis.yields.add(DL_resolved_2b_emu, "DL resolved 2b emu")
        analysis.yields.add(DL_VBF_resolved_ee, "DL VBF resolved ee")
        analysis.yields.add(DL_VBF_boosted_ee, "DL VBF boosted ee")
        analysis.yields.add(DL_VBF_resolved_mumu, "DL VBF resolved mumu")
        analysis.yields.add(DL_VBF_boosted_mumu, "DL VBF boosted mumu")
        analysis.yields.add(DL_VBF_resolved_emu, "DL VBF resolved emu")
        analysis.yields.add(DL_VBF_boosted_emu, "DL VBF boosted emu")

    if btagReweightStudy:
        return pre_final_state_sels
    elif DYControlRegion:
        return Z_peak_selections
    else:
        return DL_signal_region_selections
