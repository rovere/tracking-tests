import FWCore.ParameterSet.Config as cms

forUpgrade = True

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import *
if forUpgrade:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup')

#process.GlobalTag.globaltag = 'START61_V11::All'

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
if forUpgrade:
    process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
else:
    process.load('Configuration.StandardSequences.Geometry_cff')
    process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/home/rovere/reco_tenmu_marcofix.root',
        'file:../../reco_tenmu_fixformarcofix.root'
    #'/store/relval/CMSSW_6_2_0_pre2-START61_V11/RelValSingleMuPt10/GEN-SIM-RECO/v1/00000/52A31EA2-7478-E211-9031-002590489C9E.root',
    # '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/FE48CBD7-437A-E211-AE8D-003048F1CA6E.root',
    # '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/F600059F-2E7A-E211-9399-003048F1186A.root',
    # '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/76BF1027-B17A-E211-92CC-003048F00AF8.root',
    # '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/6C0AE730-317A-E211-8416-003048F0E00A.root',
    # '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/22D46130-377A-E211-81F5-003048F0258A.root',
    # '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/0E566081-337A-E211-9E7F-003048CF97B2.root',
    # '/store/relval/CMSSW_6_2_0_pre2-PU_START61_V11/RelValTTbar/GEN-SIM-RECO/v1/00000/0A6F8ADE-2D7A-E211-BE3D-0025901D5E20.root'
    )
)

process.beamSpot = cms.Sequence(
    process.offlineBeamSpot
)
process.clustToHits = cms.Sequence(
    process.siPixelRecHits*process.siStripMatchedRecHits
)

process.tracking = cms.Sequence(
    process.trackingGlobalReco
)

#process.select = cms.EDFilter('TrackSelector',
#        src = cms.InputTag('generalTracks'),
#        cut = cms.string("abs(eta)<0.9")
#        #cut = cms.string("quality('highPurity') & (algo=9 ) & abs(eta)<0.9")
#        #cut = cms.string("quality('highPurity') & (algo=10) & abs(eta)<0.9")
#        #cut = cms.string("quality('highPurity') & (algo=9 ) & abs(eta)>0.9 & abs(eta)<1.6")
#        #cut = cms.string("quality('highPurity') & (algo=10) & abs(eta)>0.9 & abs(eta)<1.6")
#        #cut = cms.string("quality('highPurity') & (algo=9 ) & abs(eta)>1.6 & abs(eta)<2.5")
#        #cut = cms.string("quality('highPurity') & (algo=10) & abs(eta)>1.6 & abs(eta)<2.5")
#)

process.demo = cms.EDAnalyzer('Occupancy',
        tracks = cms.untracked.InputTag('generalTracks'),
        TTRHBuilder = cms.string('WithTrackAngle'),
        saveHists = cms.bool(True),
        printInfo = cms.bool(False)
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("occupancy_marcofix.root"),
        closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.clustToHits
    * process.demo
)

if forUpgrade:
    # Need these until pixel templates are used
    process.load("SLHCUpgradeSimulations.Geometry.recoFromSimDigis_cff")
    # PixelCPEGeneric #
    process.PixelCPEGenericESProducer.Upgrade = cms.bool(True)
    process.PixelCPEGenericESProducer.UseErrorsFromTemplates = cms.bool(False)
    process.PixelCPEGenericESProducer.LoadTemplatesFromDB = cms.bool(False)
    process.PixelCPEGenericESProducer.TruncatePixelCharge = cms.bool(False)
    process.PixelCPEGenericESProducer.IrradiationBiasCorrection = False
    process.PixelCPEGenericESProducer.DoCosmics = False
     # CPE for other steps
    process.siPixelRecHits.CPE = cms.string('PixelCPEGeneric')

    # Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
    # from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_phase2_BE5DPixel10D

    # Call to customisation function cust_phase2_BE5DPixel10D imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
     # process = cust_phase2_BE5DPixel10D(process)

     # Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
     # from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn

     # Call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
     # process = setCrossingFrameOn(process)
