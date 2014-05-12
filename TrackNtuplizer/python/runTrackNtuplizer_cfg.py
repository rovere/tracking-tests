import FWCore.ParameterSet.Config as cms

process = cms.Process('TRKNTPLZR')

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

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

### validation-specific includes
process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
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
                                    simSource=cms.InputTag("mix", "MergedTrackTruth"),
                                    saveTrees=cms.bool(True),
                                    outfile=cms.string("trackAnalyzerOut.root")
                                    )

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)

process.ntuplizer_step = cms.Path(process.simHitTPAssocProducer * process.myAnalyzer)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,
                                process.L1Reco_step,
                                process.reconstruction_step,
                                process.ntuplizer_step,
                                process.endjob_step,
                                process.RECOoutput_step)



