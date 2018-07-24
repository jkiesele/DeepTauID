import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/70000/B275DA67-7045-E811-BEA8-3417EBE644DA.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/70000/AE730D32-7945-E811-AC36-7845C4FC3C23.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/70000/DEA93B38-8445-E811-97B9-008CFAFBF0BA.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/70000/8A239656-9F45-E811-9D12-7CD30AD08F02.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/70000/748A7E0F-D545-E811-B606-7CD30AB049DA.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/90000/EA4AA6B5-5345-E811-810B-008CFAFBDBA8.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/90000/366330AF-7045-E811-908E-3417EBE644DA.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/90000/4AB0A087-CB45-E811-B813-509A4C74816D.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/20000/02AB313D-AA45-E811-A1B6-7CD30AB15C58.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/20000/FC645839-DD45-E811-99D7-008CFAFFEA48.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/20000/C8A3BC3A-E045-E811-801C-509A4C7480B2.root",
"/store/mc/PhaseIISpr18AODMiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v5-v1/20000/6EC3A6A6-1446-E811-B1F3-3417EBE34CAB.root"

])

secFiles.extend( [
               ] )
