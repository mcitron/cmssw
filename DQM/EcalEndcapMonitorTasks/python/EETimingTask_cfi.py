import FWCore.ParameterSet.Config as cms

ecalEndcapTimingTask = cms.EDAnalyzer("EETimingTask",
    prefixME = cms.untracked.string('EcalEndcap'),
    enableCleanup = cms.untracked.bool(False),
    mergeRuns = cms.untracked.bool(False),    
    EcalRawDataCollection = cms.InputTag("ecalEBunpacker"),
    EcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE")
)

