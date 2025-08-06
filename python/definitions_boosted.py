from bamboo import treefunctions as op
from bamboo.treeproxies import BoolProxy


# Lepton functions

def hasAssociatedJet(lep):
    """Check if the lepton has an associated jet"""
    return lep.jet.idx != -1


def muon_x(mu):
    """Muon x variable for btag interpolation"""
    return op.min(op.max(0., (0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)


def muon_btagInterpolation(mu):
    """Muon btag interpolation"""
    return muon_x(mu)*0.047 + (1-muon_x(mu))*0.245


def muon_pNetInterpIfMvaFailed(mu):
    """Muon pNet interpolation if MVA failed"""
    return op.OR(op.NOT(
        hasAssociatedJet(mu)), mu.jet.btagPNetB < muon_btagInterpolation(mu))


def lepton_associatedJetLessThanMediumBtag(lep, era):
    """Check if the lepton's associated jet has btag score less than medium WP"""
    return op.OR(
        op.NOT(hasAssociatedJet(lep)), lep.jet.btagPNetB <= ak4MediumBtagWP(era))


def lepton_associatedJetLessThanTightBtag(lep, era):
    """Check if the lepton's associated jet has btag score less than tight WP"""
    return op.OR(
        op.NOT(hasAssociatedJet(lep)), lep.jet.btagPNetB <= ak4TightBtagWP(era))


def jetIdCorrection(j):
    """Correction for jet tight Id in nanoAOD v12.
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13p6TeV#nanoAOD_Flags"""
    passJetIdTight: BoolProxy = (j.jetId & op.c_int(1 << 1)) > 0
    return op.multiSwitch(
        (op.abs(j.eta) <= 2.7,
            passJetIdTight),
        (op.AND(op.abs(j.eta) > 2.7, op.abs(j.eta) <= 3.0),
            op.AND(passJetIdTight, j.neHEF < 0.99)),
        op.AND(passJetIdTight, j.neEmEF < 0.4)
    )

# Object definitions


def muonPreSel(muons):
    """Muon preselection"""
    return op.select(muons, lambda mu: op.AND(
        mu.pt > 5.,
        op.abs(mu.eta) < 2.4,
        op.abs(mu.dxy) < 0.05,
        op.abs(mu.dz) < 0.1,
        mu.miniPFRelIso_all <= 0.4,
        mu.sip3d < 8,
        mu.looseId
    ))


def muonFakeSel(muons, era):
    """Muon fakeable selection"""
    return op.select(muons, lambda mu: op.AND(
        mu.pt > 10.,
        op.OR(lepton_associatedJetLessThanMediumBtag(mu, era), op.AND(mu.jetRelIso < 0.8, muon_pNetInterpIfMvaFailed(mu))))
    )


def muonTightSel(muons):
    """Muon tight selection"""
    return op.select(muons, lambda mu: op.AND(
        mu.mediumPromptId,  # this run3 replacement along with mediumId, for mu.mvaTTH >= 0.50
        mu.mediumId
    ))


def elePreSel(electrons):
    """Electron preselection"""
    return op.select(electrons, lambda el: op.AND(
        el.pt > 5.,
        op.abs(el.eta) < 2.5,
        op.abs(el.dxy) < 0.05,
        op.abs(el.dz) < 0.1,
        el.sip3d < 8,
        el.miniPFRelIso_all <= 0.4,
        op.OR(el.mvaNoIso_WP80, el.mvaIso_WP80),
        el.lostHits <= 1
    ))


def cleanElectrons(electrons, muons):
    """Remove electrons within a cone of DR<0.3 of muons"""
    cleanedElectrons = op.select(electrons, lambda el: op.NOT(
        op.rng_any(
            muons, lambda mu: op.deltaR(el.p4, mu.p4) <= 0.3))
    )
    return cleanedElectrons


def elFakeSel(electrons, era):
    """Electron fakeable selection"""
    return op.select(electrons, lambda el: op.AND(
        el.pt >= 10,
        op.OR(
            op.AND(op.abs(el.eta+el.deltaEtaSC) <= 1.479, el.sieie <= 0.011),
            op.AND(op.abs(el.eta+el.deltaEtaSC) > 1.479, el.sieie <= 0.030)
        ),
        el.hoe <= 0.10,
        el.eInvMinusPInv >= -0.04,
        op.OR(el.mvaTTH >= 0.30, op.AND(el.jetRelIso < 0.7, el.mvaNoIso_WP90)),
        op.switch(
            el.mvaTTH < 0.30,
            lepton_associatedJetLessThanTightBtag(el, era),
            lepton_associatedJetLessThanMediumBtag(el, era)),
        el.lostHits == 0,
        el.convVeto
    ))


def elTightSel(electrons):
    """Electron tight selection"""
    return op.select(electrons, lambda el: el.mvaIso_WP90)


def ak4jetDef(jets):
    """AK4 jet selection"""
    return op.select(jets, lambda jet: op.AND(
        jetIdCorrection(jet),
        jet.pt >= 25.,
        op.abs(jet.eta) <= 2.4,
        jet.btagPNetB >= 0,  # due to some events having negative value for this
    ))


def ak8jetDef(jets):
    """AK8 jet selection"""
    return op.select(jets, lambda jet: op.AND(
        jet.pt >= 200.,
        op.abs(jet.eta) <= 2.4,
        (jet.jetId >> 1 & 0x1) == 1,
        (jet.jetId >> 2 & 0x1) == 1,
        jet.subJet1.isValid,
        jet.subJet2.isValid,
        jet.subJet1.pt >= 20.,
        jet.subJet2.pt >= 20.,
        op.abs(jet.subJet1.eta) <= 2.4,
        op.abs(jet.subJet2.eta) <= 2.4,
        op.AND(jet.msoftdrop >= 30., jet.msoftdrop <= 210.),
        jet.tau2 / jet.tau1 <= 0.75,
        jet.particleNet_XbbVsQCD >= 0,
    ))

# bTagging for jets


def ak4MediumBtagWP(era):
    """Medium btag WP for AK4 jets with particleNet tagger from
    https://btv-wiki.docs.cern.ch/ScaleFactors/"""
    WPs = {
        '2022': 0.245,
        '2022EE': 0.2605,
        '2023': 0.1917,
        '2023BPix': 0.1919
    }
    return WPs[era]


def ak4TightBtagWP(era):
    """Tight btag WP for AK4 jets with particleNet tagger from
    https://btv-wiki.docs.cern.ch/ScaleFactors/"""
    WPs = {
        '2022': 0.6734,
        '2022EE': 0.6915,
        '2023': 0.6172,
        '2023BPix': 0.6133
    }
    return WPs[era]


def ak8Btag(fatjet):
    """Btagging for AK8 jets"""
    return op.AND(
        fatjet.particleNet_XbbVsQCD > 0.4,
        op.OR(fatjet.subJet1.pt >= 30,
              fatjet.subJet2.pt >= 30)
    )


def tauDef(taus):
    """Tau selection"""
    return op.select(taus, lambda tau: op.AND(
        tau.pt > 20.,
        op.abs(tau.eta) < 2.3,
        op.abs(tau.dxy) <= 1000.0,
        op.abs(tau.dz) <= 0.2,
        tau.idDecayModeOldDMs,
        op.OR(tau.decayMode == 0,
              tau.decayMode == 1,
              tau.decayMode == 2,
              tau.decayMode == 10,
              tau.decayMode == 11),
        (tau.idDeepTau2017v2p1VSjet >> 4 & 0x1) == 1,
        (tau.idDeepTau2017v2p1VSe >> 0 & 0x1) == 1,
        (tau.idDeepTau2017v2p1VSmu >> 0 & 0x1) == 1
    ))


def cleanTaus(taus, electrons, muons):
    """Remove taus within a cone of DR<0.3 of electrons and muons"""
    return op.select(taus, lambda tau: op.AND(
        op.NOT(op.rng_any(
            electrons, lambda el: op.deltaR(tau.p4, el.p4) <= 0.3)),
        op.NOT(op.rng_any(
            muons, lambda mu: op.deltaR(tau.p4, mu.p4) <= 0.3))
    ))

# remove jets within cone of DR<0.4 of leading leptons at each channel


def cleaningWithRespectToLeadingLeptons(electrons, muons, DR):
    """Remove jets within a cone of DR<0.4 of leading leptons at each channel"""
    return lambda j: op.multiSwitch(
        # Only electrons
        (op.AND(op.rng_len(electrons) >= 2, op.rng_len(muons) == 0),
            op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, electrons[1].p4) >= DR)),
        # Only muons
        (op.AND(op.rng_len(electrons) == 0, op.rng_len(muons) >= 2),
            op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, muons[1].p4) >= DR)),
        # One electron + one muon
        (op.AND(op.rng_len(electrons) == 1, op.rng_len(muons) == 1),
            op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, muons[0].p4) >= DR)),
        # At least one electron + at least one muon
        (op.AND(op.rng_len(electrons) >= 1, op.rng_len(muons) >= 1),
            op.switch(
            # Electron is the leading lepton
            electrons[0].pt > muons[0].pt,
            op.switch(op.rng_len(electrons) == 1,
                      op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(
                          j.p4, muons[0].p4) >= DR),
                      op.switch(electrons[1].pt > muons[0].pt,
                                op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(
                                    j.p4, electrons[1].p4) >= DR),
                                op.AND(op.deltaR(j.p4, electrons[0].p4) >= DR, op.deltaR(j.p4, muons[0].p4) >= DR))),
            # Muon is the leading lepton
            op.switch(op.rng_len(muons) == 1,
                      op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(
                          j.p4, electrons[0].p4) >= DR),
                      op.switch(muons[1].pt > electrons[0].pt,
                                op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(
                                    j.p4, muons[1].p4) >= DR),
                                op.AND(op.deltaR(j.p4, muons[0].p4) >= DR, op.deltaR(j.p4, electrons[0].p4) >= DR))))),
        op.c_bool(True)
    )


def VBFjetDef(jets):
    """VBF jet selection"""
    return op.select(jets, lambda jet: op.AND(
        jetIdCorrection(jet),
        jet.pt >= 30.,
        op.abs(jet.eta) <= 4.7,
        op.OR(
            jet.pt >= 60.,
            op.abs(jet.eta) < 2.7,
            op.abs(jet.eta) > 3.0
        ),
        jet.btagPNetB >= 0,
    ))


def cleanVBFAk4_lambda(ak4jetsbybtag):
    """Remove jets within a cone of DR<0.8 of the two leading btagged jets"""
    return lambda j: op.multiSwitch(
        (op.rng_len(ak4jetsbybtag) > 1, op.AND(op.deltaR(
            j.p4, ak4jetsbybtag[0].p4) > 0.8, op.deltaR(j.p4, ak4jetsbybtag[1].p4) > 0.8)),
        (op.rng_len(ak4jetsbybtag) == 1, op.deltaR(
            j.p4, ak4jetsbybtag[0].p4) > 0.8),
        op.c_bool(True)
    )


def cleanVBFAk8_lambda(ak8jets, ak8bjets):
    """Remove fat jets within a cone of DR<1.2 of the leading AK8 jet or the leading btagged AK8 jet"""
    return lambda j: op.multiSwitch(
        (op.rng_len(ak8bjets) > 0, op.deltaR(j.p4, ak8bjets[0].p4) > 1.2),
        (op.rng_len(ak8jets) == 1, op.deltaR(j.p4, ak8jets[0].p4) > 1.2),
        op.c_bool(True)
    )


def VBFpair_lambda(j1, j2):
    """VBF pair selection"""
    return op.AND(
        op.invariant_mass(j1.p4, j2.p4) > 500.,
        op.abs(j1.eta - j2.eta) > 3.0,
    )


def defineObjects(self, tree):
    """Define objects for the analysis"""
    # lepton definitions sorted by their pt
    self.preMuons = op.sort(muonPreSel(tree.Muon), lambda mu: -mu.pt)

    self.preElectrons = op.sort(elePreSel(tree.Electron), lambda el: -el.pt)

    # cleaning electrons wrt muons
    self.clElectrons = cleanElectrons(self.preElectrons, self.preMuons)

    # Fakeable leptons
    self.fakeMuons = muonFakeSel(self.preMuons, era=self.era)
    self.fakeElectrons = elFakeSel(self.clElectrons, era=self.era)

    # tight leptons
    self.tightMuons = muonTightSel(self.fakeMuons)
    self.tightElectrons = elTightSel(self.fakeElectrons)

    # Taus
    taus = tauDef(tree.Tau)
    self.cleanedTaus = cleanTaus(taus, self.fakeElectrons, self.fakeMuons)

    # clean jets wrt leptons
    cleanAk4Jets_lambda = cleaningWithRespectToLeadingLeptons(
        self.fakeElectrons, self.fakeMuons, 0.4)
    cleanAk8Jets_lambda = cleaningWithRespectToLeadingLeptons(
        self.fakeElectrons, self.fakeMuons, 0.8)

    # AK4 jets
    self.ak4JetsPreSel = op.sort(ak4jetDef(tree.Jet), lambda jet: -jet.pt)

    self.ak4Jets = op.select(self.ak4JetsPreSel, cleanAk4Jets_lambda)

    ak4jetsbybtag = op.sort(self.ak4Jets, lambda j: -j.btagPNetB)

    self.ak4BJets = op.select(
        self.ak4Jets, lambda j: j.btagPNetB >= ak4MediumBtagWP(self.era))

    # AK8 Jets
    self.ak8JetsDef = ak8jetDef(tree.FatJet)

    ak8JetsPreSel = op.sort(self.ak8JetsDef, lambda j: -j.pt)

    self.ak8Jets = op.select(ak8JetsPreSel, cleanAk8Jets_lambda)

    self.ak8BJets = op.select(self.ak8Jets, ak8Btag)

    # MET

    self.met = tree.MET

    # VBF jets

    VBFjetsPreSel = op.sort(VBFjetDef(tree.Jet), lambda jet: -jet.pt)

    VBFjets = op.select(VBFjetsPreSel, cleanAk4Jets_lambda)

    self.VBFjetsResolved = op.select(
        VBFjets, cleanVBFAk4_lambda(ak4jetsbybtag))

    self.VBFjetsBoosted = op.select(
        VBFjets, cleanVBFAk8_lambda(self.ak8Jets, self.ak8BJets))

    # VBFjetPairs = op.sort(op.combine(
    #     VBFjets, N=2, pred=VBFpair_lambda), lambda pair: -op.invariant_mass(pair[0].p4, pair[1].p4))

    self.VBFjetPairsResolved = op.sort(op.combine(
        self.VBFjetsResolved, N=2, pred=VBFpair_lambda), lambda pair: -op.invariant_mass(pair[0].p4, pair[1].p4))

    self.VBFjetPairsBoosted = op.sort(op.combine(
        self.VBFjetsBoosted, N=2, pred=VBFpair_lambda), lambda pair: -op.invariant_mass(pair[0].p4, pair[1].p4))


def ml_input_features(self):
    """Define variables to be used to create the skims containing them,
        also to use them later on in the DNN evaluation."""
    def get_lepton_callable(var, idx):
        """ Get lepton callable based on the topology."""
        return op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2,
             getattr(self.tightElectrons[idx].p4, var)()),
            (op.rng_len(self.tightMuons) == 2,
             getattr(self.tightMuons[idx].p4, var)()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                getattr(self.tightElectrons[idx].p4, var)(),
                getattr(self.tightMuons[idx].p4, var)()))
        )

    def get_lepton_tree_var(var, idx):
        """Get lepton tree variable based on the topology.
        Unlike get_lepton_callable, this function implements
        conversion to float since charge and pdgId are integers."""
        return op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2,
             op.static_cast('float', getattr(self.tightElectrons[idx], var))),
            (op.rng_len(self.tightMuons) == 2,
             op.static_cast('float', getattr(self.tightMuons[idx], var))),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.static_cast('float', getattr(
                    self.tightElectrons[idx], var)),
                op.static_cast('float', getattr(self.tightMuons[idx], var)))
             ))

    def get_jet_callable(var, idx, jet):
        """Get jet callable based on the topology.
            Priority is given to AK8 jets."""
        return op.switch(
            op.rng_len(jet) > idx,
            getattr(jet[idx].p4, var)(),
            op.c_float(0)
        )

    def get_jet_tree_var(var, idx, jet):
        "Similar to get_jet_callable but returns the tree variable"
        return op.switch(
            op.rng_len(jet) > idx,
            getattr(jet[idx], var),
            op.c_float(0)
        )

    ml_vars = {
        # leptons
        "l1_Px": get_lepton_callable('Px', 0),
        "l1_Py": get_lepton_callable('Py', 0),
        "l1_Pz": get_lepton_callable('Pz', 0),
        "l1_E": get_lepton_callable('E', 0),
        "l1_pdgId": get_lepton_tree_var('pdgId', 0),
        "l1_charge": get_lepton_tree_var('charge', 0),
        "l2_Px": get_lepton_callable('Px', 1),
        "l2_Py": get_lepton_callable('Py', 1),
        "l2_Pz": get_lepton_callable('Pz', 1),
        "l2_E": get_lepton_callable('E', 1),
        "l2_pdgId": get_lepton_tree_var('pdgId', 1),
        "l2_charge": get_lepton_tree_var('charge', 1),
        # jets
        "j1_Px": get_jet_callable('Px', 0, self.ak4Jets),
        "j1_Py": get_jet_callable('Py', 0, self.ak4Jets),
        "j1_Pz": get_jet_callable('Pz', 0, self.ak4Jets),
        "j1_E": get_jet_callable('E', 0, self.ak4Jets),
        "j1_btag": get_jet_tree_var('btagPNetB', 0, self.ak4Jets),
        "j2_Px": get_jet_callable('Px', 1, self.ak4Jets),
        "j2_Py": get_jet_callable('Py', 1, self.ak4Jets),
        "j2_Pz": get_jet_callable('Pz', 1, self.ak4Jets),
        "j2_E": get_jet_callable('E', 1, self.ak4Jets),
        "j2_btag": get_jet_tree_var('btagPNetB', 1, self.ak4Jets),
        "j3_Px": get_jet_callable('Px', 2, self.ak4Jets),
        "j3_Py": get_jet_callable('Py', 2, self.ak4Jets),
        "j3_Pz": get_jet_callable('Pz', 2, self.ak4Jets),
        "j3_E": get_jet_callable('E', 2, self.ak4Jets),
        "j3_btag": get_jet_tree_var('btagPNetB', 2, self.ak4Jets),
        "j4_Px": get_jet_callable('Px', 3, self.ak4Jets),
        "j4_Py": get_jet_callable('Py', 3, self.ak4Jets),
        "j4_Pz": get_jet_callable('Pz', 3, self.ak4Jets),
        "j4_E": get_jet_callable('E', 3, self.ak4Jets),
        "j4_btag": get_jet_tree_var('btagPNetB', 3, self.ak4Jets),
        "j8_Px": get_jet_callable('Px', 0, self.ak8Jets),
        "j8_Py": get_jet_callable('Py', 0, self.ak8Jets),
        "j8_Pz": get_jet_callable('Pz', 0, self.ak8Jets),
        "j8_E": get_jet_callable('E', 0, self.ak8Jets),
        "j8_btag": get_jet_tree_var('particleNet_XbbVsQCD', 0, self.ak8Jets),
        # met
        "met_Px": op.product(self.met.pt, op.cos(self.met.phi)),
        "met_Py": op.product(self.met.pt, op.sin(self.met.phi)),
        "met_E": self.met.pt,

        "dR_l1_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.deltaR(
                self.tightElectrons[0].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.deltaR(
                self.tightMuons[0].p4, self.tightMuons[1].p4)),
            op.deltaR(self.tightElectrons[0].p4, self.tightMuons[0].p4)
        ),

        "dR_j1_j2": op.deltaR(self.ak4Jets[0].p4, self.ak4Jets[1].p4),

        "dR_j1_l1": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.deltaR(
                self.ak4Jets[0].p4, self.tightElectrons[0].p4)),
            (op.rng_len(self.tightMuons) == 2, op.deltaR(
                self.ak4Jets[0].p4, self.tightMuons[0].p4)),
            op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.deltaR(self.ak4Jets[0].p4, self.tightElectrons[0].p4),
                op.deltaR(self.ak4Jets[0].p4, self.tightMuons[0].p4))
        ),

        "dR_j1_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.deltaR(
                self.ak4Jets[0].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.deltaR(
                self.ak4Jets[0].p4, self.tightMuons[1].p4)),
            op.switch(
                self.tightElectrons[1].pt > self.tightMuons[1].pt,
                op.deltaR(self.ak4Jets[0].p4, self.tightMuons[1].p4),
                op.deltaR(self.ak4Jets[0].p4, self.tightElectrons[1].p4))
        ),

        "dR_j2_l1": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.deltaR(
                self.ak4Jets[1].p4, self.tightElectrons[0].p4)),
            (op.rng_len(self.tightMuons) == 2, op.deltaR(
                self.ak4Jets[1].p4, self.tightMuons[0].p4)),
            op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.deltaR(self.ak4Jets[1].p4, self.tightElectrons[0].p4),
                op.deltaR(self.ak4Jets[1].p4, self.tightMuons[0].p4))
        ),

        "dR_j2_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.deltaR(
                self.ak4Jets[1].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.deltaR(
                self.ak4Jets[1].p4, self.tightMuons[1].p4)),
            op.switch(
                self.tightElectrons[1].pt > self.tightMuons[1].pt,
                op.deltaR(self.ak4Jets[1].p4, self.tightMuons[1].p4),
                op.deltaR(self.ak4Jets[1].p4, self.tightElectrons[1].p4))
        ),

        "dR_j8_l1": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.deltaR(
                self.ak8Jets[0].p4, self.tightElectrons[0].p4)),
            (op.rng_len(self.tightMuons) == 2, op.deltaR(
                self.ak8Jets[0].p4, self.tightMuons[0].p4)),
            op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.deltaR(self.ak8Jets[0].p4, self.tightElectrons[0].p4),
                op.deltaR(self.ak8Jets[0].p4, self.tightMuons[0].p4))
        ),

        "dR_j8_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.deltaR(
                self.ak8Jets[0].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.deltaR(
                self.ak8Jets[0].p4, self.tightMuons[1].p4)),
            op.switch(
                self.tightElectrons[1].pt > self.tightMuons[1].pt,
                op.deltaR(self.ak8Jets[0].p4, self.tightMuons[1].p4),
                op.deltaR(self.ak8Jets[0].p4, self.tightElectrons[1].p4))
        ),

        "InvM_j1_j2": op.invariant_mass(self.ak4Jets[0].p4, self.ak4Jets[1].p4),

        "InvM_l1_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.invariant_mass(
                self.tightElectrons[0].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.invariant_mass(
                self.tightMuons[0].p4, self.tightMuons[1].p4)),
            op.invariant_mass(self.tightElectrons[0].p4, self.tightMuons[0].p4)
        ),

        "InvM_j1_l1": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.invariant_mass(
                self.ak4Jets[0].p4, self.tightElectrons[0].p4)),
            (op.rng_len(self.tightMuons) == 2, op.invariant_mass(
                self.ak4Jets[0].p4, self.tightMuons[0].p4)),
            op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.invariant_mass(self.ak4Jets[0].p4, self.tightElectrons[0].p4),
                op.invariant_mass(self.ak4Jets[0].p4, self.tightMuons[0].p4))
        ),

        "InvM_j1_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.invariant_mass(
                self.ak4Jets[0].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.invariant_mass(
                self.ak4Jets[0].p4, self.tightMuons[1].p4)),
            op.switch(
                self.tightElectrons[1].pt > self.tightMuons[1].pt,
                op.invariant_mass(self.ak4Jets[0].p4, self.tightMuons[1].p4),
                op.invariant_mass(self.ak4Jets[0].p4, self.tightElectrons[1].p4))
        ),

        "InvM_j2_l1": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.invariant_mass(
                self.ak4Jets[1].p4, self.tightElectrons[0].p4)),
            (op.rng_len(self.tightMuons) == 2, op.invariant_mass(
                self.ak4Jets[1].p4, self.tightMuons[0].p4)),
            op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.invariant_mass(self.ak4Jets[1].p4, self.tightElectrons[0].p4),
                op.invariant_mass(self.ak4Jets[1].p4, self.tightMuons[0].p4))
        ),

        "InvM_j2_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.invariant_mass(
                self.ak4Jets[1].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.invariant_mass(
                self.ak4Jets[1].p4, self.tightMuons[1].p4)),
            op.switch(
                self.tightElectrons[1].pt > self.tightMuons[1].pt,
                op.invariant_mass(self.ak4Jets[1].p4, self.tightMuons[1].p4),
                op.invariant_mass(self.ak4Jets[1].p4, self.tightElectrons[1].p4))
        ),

        "InvM_j8_l1": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.invariant_mass(
                self.ak8Jets[0].p4, self.tightElectrons[0].p4)),
            (op.rng_len(self.tightMuons) == 2, op.invariant_mass(
                self.ak8Jets[0].p4, self.tightMuons[0].p4)),
            op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.invariant_mass(self.ak8Jets[0].p4, self.tightElectrons[0].p4),
                op.invariant_mass(self.ak8Jets[0].p4, self.tightMuons[0].p4))
        ),

        "InvM_j8_l2": op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.invariant_mass(
                self.ak8Jets[0].p4, self.tightElectrons[1].p4)),
            (op.rng_len(self.tightMuons) == 2, op.invariant_mass(
                self.ak8Jets[0].p4, self.tightMuons[1].p4)),
            op.switch(
                self.tightElectrons[1].pt > self.tightMuons[1].pt,
                op.invariant_mass(self.ak8Jets[0].p4, self.tightMuons[1].p4),
                op.invariant_mass(self.ak8Jets[0].p4, self.tightElectrons[1].p4))
        )

    }

    return ml_vars
