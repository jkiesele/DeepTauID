import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                                jets       = cms.InputTag("slimmedJets"),
                                taus       = cms.InputTag("slimmedTaus"),
                                genparticles = cms.InputTag("prunedGenParticles"),
                                
                                minGenPt  = cms.double(5),
                                maxGenEta = cms.double(3),
                                
                                pupInfo = cms.InputTag("slimmedAddPileupInfo"),
                                rhoInfo = cms.InputTag("fixedGridRhoFastjetAll"),	
                                
                                )
