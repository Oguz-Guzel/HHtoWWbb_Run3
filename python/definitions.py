from bamboo import treefunctions as op


# Lepton functions

def hasAssociatedJet(lep): return lep.jet.idx != -1


def muon_x(mu): return op.min(
    op.max(0., (0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)


def muon_btagInterpolation(mu): return muon_x(
    mu)*0.047 + (1-muon_x(mu))*0.245


def muon_pNetInterpIfMvaFailed(mu): return op.OR(op.NOT(
    hasAssociatedJet(mu)), mu.jet.btagPNetB < muon_btagInterpolation(mu))


def lepton_associatedJetLessThanMediumBtag(lep): return op.OR(
    op.NOT(hasAssociatedJet(lep)), lep.jet.btagPNetB <= 0.245)


def lepton_associatedJetLessThanTightBtag(lep): return op.OR(
    op.NOT(hasAssociatedJet(lep)), lep.jet.btagPNetB <= 0.6734)

# Object definitions


def muonPreSel(muons):
    return op.select(muons, lambda mu: op.AND(
        mu.pt >= 5.,
        op.abs(mu.eta) <= 2.4,
        op.abs(mu.dxy) <= 0.05,
        op.abs(mu.dz) <= 0.1,
        mu.miniPFRelIso_all <= 0.4,
        mu.sip3d <= 8,
        mu.looseId
    ))


def muonFakeSel(muons):
    return op.select(muons, lambda mu: op.AND(
        mu.pt >= 10.,
        op.OR(lepton_associatedJetLessThanMediumBtag(mu), op.AND(mu.jetRelIso < 0.8, muon_pNetInterpIfMvaFailed(mu))))
    )


def muonTightSel(muons): return op.select(muons, lambda mu: op.AND(
    mu.mediumPromptId,  # this run3 replacement along with mediumId, for mu.mvaTTH >= 0.50
    mu.mediumId
))


def elePreSel(electrons):
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
    cleanedElectrons = op.select(electrons, lambda el: op.NOT(
        op.rng_any(
            muons, lambda mu: op.deltaR(el.p4, mu.p4) <= 0.3))
    )
    return cleanedElectrons


def elFakeSel(electrons):
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


def elTightSel(electrons): return op.select(electrons, lambda el:
                                            el.mvaIso_WP90
                                            )


def ak4jetDef(jets):
    return op.select(jets, lambda jet: op.AND(
        jet.jetId & 2,  # tight
        jet.pt >= 25.,
        op.abs(jet.eta) <= 2.4,
        # op.OR(((jet.puId >> 2) & 1), jet.pt > 50.) # Jet PU ID bit1 is loose # no puId in Run3 so far
    ))


def ak8jetDef(jets):
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
        jet.tau2 / jet.tau1 <= 0.75
    ))

# bTagging for ak4 jets


def ak4BtagSel(jet): return jet.btagPNetB > 0.245


def ak8Btag(fatjet): return op.AND(
    fatjet.particleNet_XbbVsQCD > 0.4,
    op.OR(fatjet.subJet1.pt >= 30,
          fatjet.subJet2.pt >= 30)
)


def tauDef(taus):
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

    # AK4 Jets sorted by their pt
    self.ak4JetsPreSel = op.sort(ak4jetDef(tree.Jet), lambda jet: -jet.pt)

    # clean jets wrt leptons
    if self.channel == 'DL':
        self.cleanAk4Jets = cleaningWithRespectToLeadingLeptons(
            self.fakeElectrons, self.fakeMuons, 0.4)
        self.cleanAk8Jets = cleaningWithRespectToLeadingLeptons(
            self.fakeElectrons, self.fakeMuons, 0.8)

    if self.channel == 'SL':
        self.cleanAk4Jets = cleaningWithRespectToLeadingLepton(
            self.fakeElectrons, self.fakeMuons, 0.4)
        self.cleanAk8Jets = cleaningWithRespectToLeadingLepton(
            self.fakeElectrons, self.fakeMuons, 0.8)

    self.ak4Jets = op.select(self.ak4JetsPreSel, self.cleanAk4Jets)
    self.ak4JetsByBtagScore = op.sort(self.ak4Jets, lambda j: -j.btagPNetB)

    self.ak4BJets = op.select(self.ak4Jets, ak4BtagSel)

    # AK8 Jets
    self.ak8JetsDef = ak8jetDef(tree.FatJet)

    if self.channel == 'SL':  # sorted by btag score
        ak8JetsPreSel = op.sort(
            self.ak8JetsDef, lambda j: -j.particleNet_XbbVsQCD)
    if self.channel == 'DL':  # sorted by pt
        ak8JetsPreSel = op.sort(self.ak8JetsDef, lambda j: -j.pt)

    self.ak8Jets = op.select(ak8JetsPreSel, self.cleanAk8Jets)

    self.ak8BJets = op.select(self.ak8Jets, ak8Btag)

    # Ak4 Jet Collection cleaned from Ak8b #
    def cleanAk4FromAk8b(ak4j): return op.AND(op.rng_len(
        self.ak8BJets) > 0, op.deltaR(ak4j.p4, self.ak8BJets[0].p4) > 1.2)

    self.ak4JetsCleanedFromAk8b = op.select(self.ak4Jets, cleanAk4FromAk8b)


def ml_input_features(self, tree):
    """Define variables to be used to create the skims containing them, also to use them later on in the DNN evaluation."""
    l1_Px = op.multiSwitch(
        (op.rng_len(self.tightElectrons) == 2,
            self.tightElectrons[0].p4.Px()),  # if nElectrons = 2
        (op.rng_len(self.tightMuons) == 2,
            self.tightMuons[0].p4.Px()),  # elif nMuons = 2
        (op.switch(  # else meaning nElectrons = nMuons = 1 since no other case in the DL channel
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.Px(), self.tightMuons[0].p4.Px()))
    )
    l2_Px = op.multiSwitch(
        (op.rng_len(self.tightElectrons) ==
            2, self.tightElectrons[1].p4.Px()),
        (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Px()),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.Px(), self.tightMuons[0].p4.Px()))
    )
    l1_Py = op.multiSwitch(
        (op.rng_len(self.tightElectrons) ==
            2, self.tightElectrons[0].p4.Py()),
        (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.Py()),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.Py(), self.tightMuons[0].p4.Py()))
    )
    l2_Py = op.multiSwitch(
        (op.rng_len(self.tightElectrons) ==
            2, self.tightElectrons[1].p4.Py()),
        (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Py()),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.Py(), self.tightMuons[0].p4.Py()))
    )
    l1_Pz = op.multiSwitch(
        (op.rng_len(self.tightElectrons) ==
            2, self.tightElectrons[0].p4.Pz()),
        (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.Pz()),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.Pz(), self.tightMuons[0].p4.Pz()))
    )
    l2_Pz = op.multiSwitch(
        (op.rng_len(self.tightElectrons) ==
            2, self.tightElectrons[1].p4.Pz()),
        (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.Pz()),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.Pz(), self.tightMuons[0].p4.Pz()))
    )
    l1_E = op.multiSwitch(
        (op.rng_len(self.tightElectrons) ==
            2, self.tightElectrons[0].p4.E()),
        (op.rng_len(self.tightMuons) == 2, self.tightMuons[0].p4.E()),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.E(), self.tightMuons[0].p4.E()))
    )
    l2_E = op.multiSwitch(
        (op.rng_len(self.tightElectrons) ==
            2, self.tightElectrons[1].p4.E()),
        (op.rng_len(self.tightMuons) == 2, self.tightMuons[1].p4.E()),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            self.tightElectrons[0].p4.E(), self.tightMuons[0].p4.E()))
    )
    l1_pdgId = op.multiSwitch(  # static_cast is used to convert the pdgId to float
        (op.rng_len(self.tightElectrons) == 2, op.static_cast(
            'float', self.tightElectrons[0].pdgId)),
        (op.rng_len(self.tightMuons) == 2, op.static_cast(
            'float', self.tightMuons[0].pdgId)),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            op.static_cast('float', self.tightElectrons[0].pdgId), op.static_cast('float', self.tightMuons[0].pdgId)))
    )
    l2_pdgId = op.multiSwitch(
        (op.rng_len(self.tightElectrons) == 2, op.static_cast(
            'float', self.tightElectrons[1].pdgId)),
        (op.rng_len(self.tightMuons) == 2, op.static_cast(
            'float', self.tightMuons[1].pdgId)),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            op.static_cast('float', self.tightElectrons[0].pdgId), op.static_cast('float', self.tightMuons[0].pdgId)))
    )
    l1_charge = op.multiSwitch(
        (op.rng_len(self.tightElectrons) == 2, op.static_cast(
            'float', self.tightElectrons[0].charge)),
        (op.rng_len(self.tightMuons) == 2, op.static_cast(
            'float', self.tightMuons[0].charge)),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            op.static_cast('float', self.tightElectrons[0].charge), op.static_cast('float', self.tightMuons[0].charge)))
    )
    l2_charge = op.multiSwitch(
        (op.rng_len(self.tightElectrons) == 2, op.static_cast(
            'float', self.tightElectrons[1].charge)),
        (op.rng_len(self.tightMuons) == 2, op.static_cast(
            'float', self.tightMuons[1].charge)),
        (op.switch(
            self.tightElectrons[0].pt > self.tightMuons[0].pt,
            op.static_cast('float', self.tightElectrons[0].charge), op.static_cast('float', self.tightMuons[0].charge)))
    )

    j1_Px = self.ak4Jets[0].p4.Px()
    j1_Py = self.ak4Jets[0].p4.Py()
    j1_Pz = self.ak4Jets[0].p4.Pz()
    j1_E = self.ak4Jets[0].p4.E()
    j1_btag = self.ak4Jets[0].btagPNetB
    j2_Px = self.ak4Jets[1].p4.Px()
    j2_Py = self.ak4Jets[1].p4.Py()
    j2_Pz = self.ak4Jets[1].p4.Pz()
    j2_E = self.ak4Jets[1].p4.E()
    j2_btag = self.ak4Jets[1].btagPNetB

    met_Px = op.product(tree.MET.pt, op.cos(tree.MET.phi))
    met_Py = op.product(tree.MET.pt, op.sin(tree.MET.phi))
    met_E = tree.MET.pt

    lepton1_vars = {"l1_Px": l1_Px, "l1_Py": l1_Py, "l1_Pz": l1_Pz,
                    "l1_E": l1_E, "l1_pdgId": l1_pdgId, "l1_charge": l1_charge}
    lepton2_vars = {"l2_Px": l2_Px, "l2_Py": l2_Py, "l2_Pz": l2_Pz,
                    "l2_E": l2_E, "l2_pdgId": l2_pdgId, "l2_charge": l2_charge}

    jet1_vars = {"j1_Px": j1_Px, "j1_Py": j1_Py, "j1_Pz": j1_Pz,
                 "j1_E": j1_E, "j1_btag": j1_btag}
    jet2_vars = {"j2_Px": j2_Px, "j2_Py": j2_Py, "j2_Pz": j2_Pz,
                 "j2_E": j2_E, "j2_btag": j2_btag}

    met_vars = {"met_Px": met_Px, "met_Py": met_Py, 'met_E': met_E}

    return lepton1_vars, lepton2_vars, jet1_vars, jet2_vars, met_vars
