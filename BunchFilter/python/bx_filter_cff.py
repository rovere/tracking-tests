import FWCore.ParameterSet.Config as cms

bxFilter = cms.EDFilter("BunchFilter",
                        lowest_bx  = cms.int32(1),
                        highest_bx = cms.int32(1),
                        filter = cms.bool(True))

bxfilter = cms.Sequence(bxFilter)
