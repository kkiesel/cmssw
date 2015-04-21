import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")

options = VarParsing.VarParsing ('analysis')
options.register ('fast',
    '0',
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "process fastsim"
)
options.register ('half',
    '0',
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,
    "process fastsim"
)
options.maxEvents = -1

options.parseArguments()

isFast = int(options.fast) != 0

if isFast:
    options.outputFile = '3d_fast'
    from fastFileList import filelist
else:
    options.outputFile = '3d_full'
    from fullFileList import filelist

# only e=30:
#filelist = [ x for x in filelist if "_E30_" in x ]
#options.outputFile += "_E30"
### end only 30

half = int( options.half )
if half:
    firstHalf = filelist[:len(filelist)/2]
    secondHalf = filelist[len(filelist)/2:]
    if half == 1:
        filelist = firstHalf
        options.outputFile += "1"
    elif half == 2:
        filelist = secondHalf
        options.outputFile += "2"
    else:
        print "please use half=1 or half=2"
        exit

options.inputFiles = filelist

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( options.maxEvents ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( options.inputFiles )
)

process.load( "Analyzer.ECALScaleFactorCalculator.CfiFile_cfi" )
if not isFast:
    process.ecalScaleFactorCalculator.module = 'g4SimHits'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string( options.outputFile )
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.default.limit = 100000000


process.p = cms.Path( process.ecalScaleFactorCalculator )
