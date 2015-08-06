import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

fastSim = True
if fastSim:
    modulename =  "famosSimHits"
    inputFile = 'file:../../Generation/gensim.root'
    inputFile = 'file:/user/kiesel/root-files/fastSimStudies/FlatPhotonGun_fast_noTracker.root'
    inputFile = 'file:../../SingleElectronValidation.root'
    inputFile = '/store/user/kiesel/SingleElectronGun_E30_fast/SingleElectronGun_E30_fast/ddef0922dd59c60724f16c5ac4acaff9/out_10_1_ie2.root'

    from fastFileList import filelist
    filelist = [ x for x in filelist if "_E30_" in x ]
    #inputFile = filelist

    outputFile = "fastsim_E30_muchStat.root"
    outputFile = "fastsim.root"
else:
    modulename =  "g4SimHits"
    inputFile = 'file:/user/kiesel/root-files/fastSimStudies/FlatPhotonGun_full_tracker2.root'
    inputFile = 'file:/user/kiesel/root-files/fastSimStudies/FlatPhotonGun_full_noTracker.root'

    from fullFileList import filelist
    filelist = [ x for x in filelist if "_E30_" in x ]
    #inputFile = filelist

    outputFile = "fullsim_E30_muchStat.root"
    outputFile = "fullsim.root"





process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        inputFile
    )
)

process.SimTreeProducer = cms.EDAnalyzer('SimTreeWriter',
    ProducerModule = cms.string( modulename )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string( outputFile )
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.p = cms.Path( process.SimTreeProducer )
