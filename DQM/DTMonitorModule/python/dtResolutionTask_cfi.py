import FWCore.ParameterSet.Config as cms

dtResolutionAnalysisMonitor = cms.EDAnalyzer("DTResolutionAnalysisTask",
                                             # labels of 4D hits
                                             recHits4DLabel = cms.string('dt4DSegments'),
                                             # interval of lumi block after which we reset the histos
                                             ResetCycle = cms.untracked.int32(10000),
                                             # top folder for the histograms in DQMStore
                                             topHistoFolder = cms.untracked.string("DT/02-Segments")
                                             )


