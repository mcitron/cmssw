import FWCore.ParameterSet.Config as cms

process = cms.Process('TNPEFFICIENCYPOSTPROCESSOR')

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load("DQMOffline.Trigger.TnPEfficiencyPostProcessor_cff")
process.jpsiClient.SavePlotsInRootFileName = cms.untracked.string("/tmp/TnPEfficiencyJpsiPlots.root")
process.upsilonClient.SavePlotsInRootFileName = cms.untracked.string("/tmp/TnPEfficiencyUpsilonPlots.root")
process.zClient.SavePlotsInRootFileName = cms.untracked.string("/tmp/TnPEfficiencyZPlots.root")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    fileMode = cms.untracked.string('FULLMERGE')
)

process.source = cms.Source("PoolSource",
    processingMode = cms.untracked.string('RunsAndLumis'),
    fileNames = cms.untracked.vstring('file:/tmp/TnPEfficiency.root')
)

process.load('Configuration/StandardSequences/EDMtoMEAtJobEnd_cff')
process.dqmSaver.workflow = '/TnPEfficiency/Post/Processor'
process.dqmSaver.dirName = '/tmp/'
process.path = cms.Path(process.EDMtoME*process.tagAndProbeEfficiencyPostProcessor)

process.endpath = cms.EndPath(process.DQMSaver)

