import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("Analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

filenames_hgg_fast = [
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/1A44894A-5281-E411-99B4-0025905B85D6.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/30DC154C-5281-E411-BBBF-0025905B85AA.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/50DC506E-5281-E411-8425-0025905A6122.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/68F72B5B-5281-E411-A9C4-00259059642E.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/70C5914C-5281-E411-81A7-0025905B858C.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/8666BF49-5281-E411-8F62-0025905B8592.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/CC07FA68-5281-E411-B77C-0025905A48D8.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/DE5CC649-5281-E411-9C73-0025905B85A2.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RECO/MCRUN2_73_V7_FastSim-v1/00000/E074EC54-5281-E411-9090-0025905B8590.root"
]

filenames_hgg_full_debug = [
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_73_V7-v1/00000/2412742F-5A81-E411-97CE-00261894397F.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_73_V7-v1/00000/2ABF16C3-5781-E411-B1AB-0025905B8576.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_73_V7-v1/00000/50B250C0-5781-E411-945F-00261894388A.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_73_V7-v1/00000/B662322F-5A81-E411-8588-002618943810.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-DIGI-RAW-HLTDEBUG/MCRUN2_73_V7-v1/00000/DC17DAC7-5781-E411-8F78-0025905B85EE.root"
]
filenames_hgg_full_reco = [
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-RECO/MCRUN2_73_V7-v1/00000/9851956E-8281-E411-A6FC-00259059391E.root",
    "/store/relval/CMSSW_7_3_0/RelValH130GGgluonfusion_13/GEN-SIM-RECO/MCRUN2_73_V7-v1/00000/6E405963-8281-E411-BD96-0026189438B5.root"
]



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:RelValH130GGgluonfusion_13_CMSSW_7_3_0-MCRUN2_73_V7_FastSim-v1_GEN-SIM-DIGI-REC_1A44894A-5281-E411-99B4-0025905B85D6.root'
        'file:../../reproduceDQMplotsFastSim/Hgg/H130GGgluonfusion_13TeV_cfi_GEN_SIM_RECO_EI_HLT_VALIDATION.root'
        #filenames_hgg_fast
        #filenames_hgg_full_reco
        #filenames_hgg_full_debug
    )
)

process.ana = cms.EDAnalyzer('ResponseOnDifferentLevels',
    module = cms.untracked.string( 'famosSimHits' )
)

if not "FastSim" in process.source.fileNames[0]:
    process.ana.module = 'g4SimHits'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string( "test.root" )
)

if "FastSim" in process.source.fileNames[0] and "store" in process.source.fileNames[0]:
    process.TFileService.fileName = "fast.root"
if "HLTDEBUG" in process.source.fileNames[0]:
    process.TFileService.fileName = "full_debug.root"
if "SIM-RECO" in process.source.fileNames[0]:
    process.TFileService.fileName = "full_reco.root"


process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100000000


process.p = cms.Path( process.ana )
