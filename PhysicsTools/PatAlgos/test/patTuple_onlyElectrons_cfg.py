## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

## load tau sequences up to selectedPatElectrons
process.load("PhysicsTools.PatAlgos.producersLayer1.electronProducer_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi")

## make sure to keep the created objects
process.out.outputCommands = ['keep *_selectedPat*_*_*',]

## let it run
process.p = cms.Path(
    process.makePatElectrons *
    process.selectedPatElectrons
)



## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#                                         ##
#   process.source.fileNames = [          ##
#    '/store/relval/CMSSW_3_5_0_pre1/RelValTTbar/GEN-SIM-RECO/STARTUP3X_V14-v1/0006/14920B0A-0DE8-DE11-B138-002618943926.root'
#   ]                                     ##  (e.g. 'file:AOD.root')
#                                         ##
#   process.maxEvents.input = ...         ##  (e.g. -1 to run on all events)
#                                         ##
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#                                         ##
#   process.out.fileName = ...            ##  (e.g. 'myTuple.root')
#                                         ##
#   process.options.wantSummary = True    ##  (to suppress the long output at the end of the job)    
