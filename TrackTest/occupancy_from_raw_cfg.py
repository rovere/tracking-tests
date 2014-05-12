import FWCore.ParameterSet.Config as cms

forUpgrade = False

process = cms.Process("Occupancy")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/group/phys_tracking/samples_710pre7/DIGI/AVE_PU25_BX25/TTbar/DIGI_PU25_BX25_DIGI_L1_DIGI2RAW_HLT_PU_98_1_aEo.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('RECO nevts:-1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOEventContent.outputCommands,
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('RECO')
    )
)

# Additional output definition

# Condition statements
from Configuration.AlCa.GlobalTag import GlobalTag
if forUpgrade:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

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

process.occupancy = cms.EDAnalyzer('Occupancy',
                                   tracks = cms.untracked.InputTag('generalTracks'),
                                   TTRHBuilder = cms.string('WithTrackAngle'),
                                   saveHists = cms.bool(True),
                                   printInfo = cms.bool(False)
                                   )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("occupancy.root"),
        closeFileFast = cms.untracked.bool(True)
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)
process.occupancy_step = cms.Path(process.occupancy)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,
                                process.L1Reco_step,
                                process.reconstruction_step,
                                process.occupancy_step,
                                process.endjob_step,
                                process.RECOoutput_step
                                )
### UPGRADE CUSTOMIZATIONS
if forUpgrade:
    process.load('Configuration.Geometry.GeometryExtendedPhase2TkBE5DPixel10DReco_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

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
