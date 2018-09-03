import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepTauNTuplizer',
                                jets       = cms.InputTag("slimmedJetsPuppi"),
                                taus       = cms.InputTag("slimmedTaus"),
                                genparticles = cms.InputTag("prunedGenParticles"),
                                vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                
                                minGenPt  = cms.double(5),
                                maxGenEta = cms.double(3),
                                
                                pupInfo = cms.InputTag("slimmedAddPileupInfo"),
                                rhoInfo = cms.InputTag("fixedGridRhoFastjetAll"),	
                                
                                onlyRecTaus = cms.bool(True),
                                overSample = cms.int32(1),
                                onlyTaus = cms.bool(False),
                                
                                
                                )
