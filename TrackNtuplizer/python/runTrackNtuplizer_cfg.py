import FWCore.ParameterSet.Config as cms
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
        'file:../../../reco_tenmu_fixformarcofix.root'
        ] );

# secFiles.extend( [
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/DC902058-0FD7-E111-A3B9-001731EF61B4.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/C85EBE55-0FD7-E111-9340-001731EF61B4.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/B0B7BB0A-0DD7-E111-A82B-003048679000.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/A26EA106-0FD7-E111-96E7-001731EF61B4.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/86E50D58-0FD7-E111-8C8B-001731EF61B4.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/7E0E18C7-0BD7-E111-8D65-002618943809.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/4ECA1F58-0FD7-E111-A71F-001731EF61B4.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/34CF9F06-0FD7-E111-AD25-001731EF61B4.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/22BCD400-0FD7-E111-848B-001731EF61B4.root',
#     '/store/relval/CMSSW_6_0_0_pre10-START60_V4/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/1E7ECF4F-20D7-E111-B9BF-00261894392B.root'
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/C8A762A4-5BAE-E111-A5F8-003048678F84.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/B8A7A3AA-5BAE-E111-82DC-003048678E6E.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/8013E551-62AE-E111-9D9F-001A92810AE4.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/6A2C81AA-5BAE-E111-AFE2-003048678E6E.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/4C36419F-5BAE-E111-B482-00304867918E.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/4A6564A4-5BAE-E111-8B1E-003048678F84.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/329D89AA-5BAE-E111-A17C-003048678E6E.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/1EAC9A9E-5BAE-E111-877E-002618943829.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/0A8EE99F-5BAE-E111-9218-002618943933.root',
# #    '/store/relval/CMSSW_6_0_0_pre6-START53_V6/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/v1/0000/0A3B1DA1-5BAE-E111-AFC6-00304867904E.root'
#        ] );


#import FWCore.ParameterSet.Config as cms

process = cms.Process("TRKNTPLZR")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.infos.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(100)
    )

### standard includes
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START60_V0::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = source

### validation-specific includes
process.load("SimTracker.TrackAssociation.TrackAssociatorByPosition_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("Validation.RecoTrack.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("DQMServices.Components.EDMtoMEConverter_cff")

process.load("Validation.Configuration.postValidation_cff")
process.load("RecoTracker.FinalTrackSelectors.selectHighPurity_cfi")

process.TrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
process.TrackAssociatorByHits.Purity_SimToReco = cms.double(0.75)
process.TrackAssociatorByHits.Cut_RecoToSim = cms.double(0.75)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.myAnalyzer = cms.EDAnalyzer("TrackNtuplizer",
                                    source=cms.string("generalTracks"),
                                    beamspot=cms.InputTag("offlineBeamSpot"),
                                    simSource=cms.InputTag("mergedtruth","MergedTrackTruth"),
                                    saveTrees=cms.bool(True),
                                    outfile=cms.string("trackAnalyzerOut.root")
                                    )
process.p = cms.Path(
#    process.selectHighPurity*
    process.myAnalyzer
)


