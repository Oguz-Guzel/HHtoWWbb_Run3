from bamboo import treefunctions as op


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

# Object definitions


def muonPreSel(muons):
    """Muon preselection"""
    return op.select(muons, lambda mu: op.AND(
        mu.pt >= 7.,
        op.abs(mu.eta) <= 2.4,
        op.abs(mu.dxy) <= 0.05,
        op.abs(mu.dz) <= 0.1,
        mu.miniPFRelIso_all <= 0.4,
        mu.sip3d <= 8,
        mu.looseId
    ))


def muonFakeSel(muons, era):
    """Muon fakeable selection"""
    return op.select(muons, lambda mu: op.AND(
        mu.pt >= 10.,
        op.OR(lepton_associatedJetLessThanMediumBtag(mu), op.AND(mu.jetRelIso < 0.8, muon_pNetInterpIfMvaFailed(mu))))
    )


def muonTightSel(muons): return op.select(muons, lambda mu: op.AND(
    mu.mediumPromptId,  # this run3 replacement along with mediumId, for mu.mvaTTH >= 0.50
    mu.mediumId
))


def elePreSel(electrons):
    """Electron preselection"""
    return op.select(electrons, lambda el: op.AND(
        el.pt >= 7.,
        op.abs(el.eta) <= 2.5,
        op.abs(el.dxy) <= 0.05,
        op.abs(el.dz) <= 0.1,
        el.sip3d <= 8,
        el.miniPFRelIso_all <= 0.4,
        el.mvaIso_WP90,  # no mvaNoIso_WPL for run3 signal, using this instead
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
            lepton_associatedJetLessThanTightBtag(el),
            lepton_associatedJetLessThanMediumBtag(el)),
        el.lostHits == 0,
        el.convVeto
    ))


def elTightSel(electrons):
    """Electron tight selection"""
    return op.select(electrons, lambda el: el.mvaIso_WP90)


def ak4jetDef(jets):
    """AK4 jet selection"""
    return op.select(jets, lambda jet: op.AND(
        jet.jetId & 2,  # tight
        jet.pt >= 25.,
        op.abs(jet.eta) <= 2.4,
        jet.btagPNetB >= 0,  # due to some events having negative value for this
        # op.OR(((jet.puId >> 2) & 1), jet.pt > 50.) # Jet PU ID bit1 is loose # no puId in Run3 so far
    ))


def ak8jetDef(jets):
    """AK8 jet selection"""
    return op.select(jets, lambda jet: op.AND(
        jet.pt >= 200.,
        op.abs(jet.eta) <= 2.4,
        jet.jetId & 2,  # tight, change to the following, for ak4 too
        #       (jet.jetId>>1 & 0x1) == 1, #pass tigh
        #       (jet.jetId>>2 & 0x1) == 1, #pass tightleptveto
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


def ak8Btag(fatjet): return op.AND(
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


def cleaningWithRespectToLeadingLepton(electrons, muons, DR):
    return lambda jet: op.multiSwitch(
        (op.AND(op.rng_len(electrons) >= 1, op.rng_len(
            muons) == 0), op.deltaR(jet.p4, electrons[0].p4) >= DR),
        (op.AND(op.rng_len(electrons) == 0, op.rng_len(
            muons) >= 1), op.deltaR(jet.p4, muons[0].p4) >= DR),
        (op.AND(op.rng_len(muons) >= 1, op.rng_len(electrons) >= 1), op.switch(
            electrons[0].pt >= muons[0].pt,
            op.deltaR(jet.p4, electrons[0].p4) >= DR,
            op.deltaR(jet.p4, muons[0].p4) >= DR)),
        op.c_bool(True)
    )


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


def defineObjects(self, tree):
    """Define objects for the analysis"""
    # lepton definitions sorted by their pt
    self.preMuons = op.sort(muonPreSel(tree.Muon), lambda mu: -mu.pt)

    self.preElectrons = op.sort(elePreSel(tree.Electron), lambda el: -el.pt)

    # cleaning electrons wrt muons
    self.clElectrons = cleanElectrons(self.preElectrons, self.preMuons)

    # Fakeable leptons
    self.fakeMuons = muonFakeSel(self.preMuons)
    self.fakeElectrons = elFakeSel(self.clElectrons)

    # tight leptons
    self.tightMuons = muonTightSel(self.fakeMuons)
    self.tightElectrons = elTightSel(self.fakeElectrons)

    # Taus
    taus = tauDef(tree.Tau)
    self.cleanedTaus = cleanTaus(taus, self.fakeElectrons, self.fakeMuons)

    # clean jets wrt leptons
    if self.channel == 'DL':
        cleanAk4Jets_lambda = cleaningWithRespectToLeadingLeptons(
            self.fakeElectrons, self.fakeMuons, 0.4)
        cleanAk8Jets_lambda = cleaningWithRespectToLeadingLeptons(
            self.fakeElectrons, self.fakeMuons, 0.8)

    if self.channel == 'SL':
        cleanAk4Jets_lambda = cleaningWithRespectToLeadingLepton(
            self.fakeElectrons, self.fakeMuons, 0.4)
        cleanAk8Jets_lambda = cleaningWithRespectToLeadingLepton(
            self.fakeElectrons, self.fakeMuons, 0.8)

    # AK4 jets
    self.ak4JetsPreSel = op.sort(ak4jetDef(tree.Jet), lambda jet: -jet.pt)

    self.ak4Jets = op.select(self.ak4JetsPreSel, cleanAk4Jets_lambda)

    self.ak4BJets = op.select(
        self.ak4Jets, lambda j: j.btagPNetB >= ak4MediumBtagWP(self.era))

    # AK8 Jets
    self.ak8JetsDef = ak8jetDef(tree.FatJet)

    ak8JetsPreSel = op.sort(self.ak8JetsDef, lambda j: -j.pt)

    self.ak8Jets = op.select(ak8JetsPreSel, cleanAk8Jets_lambda)

    self.ak8BJets = op.select(self.ak8Jets, ak8Btag)

    # MET

    self.met = tree.MET


def ml_input_features(self):
    """Define variables to be used to create the skims containing them, also to use them later on in the DNN evaluation."""
    def get_lepton_callable(var, idx):
        return op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, getattr(
                self.tightElectrons[idx].p4, var)()),
            (op.rng_len(self.tightMuons) == 2, getattr(
                self.tightMuons[idx].p4, var)()),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                getattr(self.tightElectrons[idx].p4, var)(
                ) if idx == 0 else getattr(self.tightMuons[0].p4, var)(),
                getattr(self.tightMuons[idx].p4, var)() if idx == 0 else getattr(self.tightElectrons[0].p4, var)()))
        )

    def get_lepton_tree_var(var, idx):
        return op.multiSwitch(
            (op.rng_len(self.tightElectrons) == 2, op.static_cast('float', getattr(
                self.tightElectrons[idx], var))),
            (op.rng_len(self.tightMuons) == 2, op.static_cast('float', getattr(
                self.tightMuons[idx], var))),
            (op.switch(
                self.tightElectrons[0].pt > self.tightMuons[0].pt,
                op.static_cast('float', getattr(self.tightElectrons[idx], var)) if idx == 0 else op.static_cast('float', getattr(
                    self.tightMuons[0], var)),
                op.static_cast('float', getattr(self.tightMuons[idx], var)) if idx == 0 else op.static_cast('float', getattr(self.tightElectrons[0], var)))
             ))

    # if you add new items in the follwoing dictionaries, make sure to
    # add correspoding binnings in `ml_input_var_binning` function in the utils.py file
    # and arrange inputs to the ML model if you're running with a previously trained model
    lepton1_vars = {"l1_Px": get_lepton_callable('Px', 0),
                    "l1_Py": get_lepton_callable('Py', 0),
                    "l1_Pz": get_lepton_callable('Pz', 0),
                    "l1_E": get_lepton_callable('E', 0),
                    "l1_pdgId": get_lepton_tree_var('pdgId', 0),
                    "l1_charge": get_lepton_tree_var('charge', 0),
                    "l1_pt": get_lepton_callable('pt', 0),
                    "l1_eta": get_lepton_callable('eta', 0)}

    lepton2_vars = {"l2_Px": get_lepton_callable('Px', 1),
                    "l2_Py": get_lepton_callable('Py', 1),
                    "l2_Pz": get_lepton_callable('Pz', 1),
                    "l2_E": get_lepton_callable('E', 1),
                    "l2_pdgId": get_lepton_tree_var('pdgId', 1),
                    "l2_charge": get_lepton_tree_var('charge', 1),
                    "l2_pt": get_lepton_callable('pt', 1),
                    "l2_eta": get_lepton_callable('eta', 1)}

    jet1_vars = {
        "j1_Px": self.ak4Jets[0].p4.Px(),
        "j1_Py": self.ak4Jets[0].p4.Py(),
        "j1_Pz": self.ak4Jets[0].p4.Pz(),
        "j1_E": self.ak4Jets[0].p4.E(),
        "j1_btag": self.ak4Jets[0].btagPNetB,
        "j1_pt": self.ak4Jets[0].pt,
        "j1_eta": self.ak4Jets[0].eta,
        "j1_N": op.static_cast('float', op.rng_len(self.ak4Jets))
    }

    jet2_vars = {
        "j2_Px": self.ak4Jets[1].p4.Px(),
        "j2_Py": self.ak4Jets[1].p4.Py(),
        "j2_Pz": self.ak4Jets[1].p4.Pz(),
        "j2_E": self.ak4Jets[1].p4.E(),
        "j2_btag": self.ak4Jets[1].btagPNetB,
        "j2_pt": self.ak4Jets[1].pt,
        "j2_eta": self.ak4Jets[1].eta,
        "j2_N": op.static_cast('float', op.rng_len(self.ak4Jets))
    }

    met_vars = {
        "met_Px": op.product(self.met.pt, op.cos(self.met.phi)),
        "met_Py": op.product(self.met.pt, op.sin(self.met.phi)),
        "met_E": self.met.pt
    }

    return lepton1_vars, lepton2_vars, jet1_vars, jet2_vars, met_vars
