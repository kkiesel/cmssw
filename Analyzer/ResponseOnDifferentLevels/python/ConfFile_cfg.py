import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# user function

def getFileNamesFromDataset( datasetname, prefix="root://xrootd-cms.infn.it//" ):
    # Returns filenames found on dbs to the accroding dataset

    import subprocess
    out = subprocess.check_output( 'das_client.py --limit 0 --query="file dataset=%s"'%datasetname, shell=True )

    # Split the output after each linebreak, and attach the prefix
    # After the last linebreak, there is no file (obmitted)
    files = [ prefix + x for x in out.split("\n")[:-1] ]

    if not files:
        print "Warning: No files found for", datasetname
    return files


process = cms.Process("Analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:../../reproduceDQMplotsFastSim/Hgg/H130GGgluonfusion_13TeV_cfi_GEN_SIM_RECO_EI_HLT_VALIDATION.root'
        #'file:../../reproduceDQMplotsFastSim/Hgg/H130GGgluonfusion_13TeV_TuneCUETP8M1_cfi_GEN_SIM_RECO_EI_HLT_VALIDATION.root'
        #getFileNamesFromDataset( "/RelValH130GGgluonfusion_13/CMSSW_7_5_0_pre3-MCRUN2_74_V7-v2/GEN-SIM-DIGI-RAW-HLTDEBUG" )
        #getFileNamesFromDataset( "/RelValH130GGgluonfusion_13/CMSSW_7_5_0_pre3-MCRUN2_74_V7-v2/GEN-SIM-RECO" )
        getFileNamesFromDataset( "/RelValH130GGgluonfusion_13/CMSSW_7_5_0_pre3-MCRUN2_74_V7_FastSim-v1/GEN-SIM-DIGI-RECO" )
    )
)

process.ana = cms.EDAnalyzer('ResponseOnDifferentLevels',
    module = cms.untracked.string( 'famosSimHits' )
)

if not "FastSim" in process.source.fileNames[0]:
    process.ana.module = 'g4SimHits'


# set output name
process.TFileService = cms.Service("TFileService",
    fileName = cms.string( "test.root" )
)

outputReplacements = {
    "HLTDEBUG": "full_debug",
    "SIM-RECO": "full_reco",
    "reproduceDQMplotsFastSim": "mod",
    "_FastSim-": "fast"
}

for exp, name in outputReplacements.iteritems():
    if exp in process.source.fileNames[0]:
        process.TFileService.fileName = name + ".root"

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100000000


process.p = cms.Path( process.ana )
