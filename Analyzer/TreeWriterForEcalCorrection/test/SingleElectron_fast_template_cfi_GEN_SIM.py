# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleGammaE30Eta03Phi0_cfi --conditions auto:run2_mc --magField 38T_PostLS1 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 -s GEN,SIM --eventcontent FEVTDEBUG --datatier GEN-SIM --fileout gensim.root --fast -n 1000
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_NoPileUp_cff')
process.load('FastSimulation.Configuration.Geometries_MC_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 10000 )
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string(''),
    name = cms.untracked.string('Applications')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    # To save disk space, save only important collections
    outputCommands = cms.untracked.vstring(
            'drop *',
            'keep *_*_EcalHits*_SIM',
            'keep *_genParticles_*_SIM'
    ),
    fileName = cms.untracked.string('out.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(11), # electrons
        MinPhi = cms.double(0.0),
        MaxPhi = cms.double(0.0),
        MinEta = cms.double(0.0),
        MaxEta = cms.double(3.1),
        MinE = cms.double( $ENERGY ),
        MaxE = cms.double( $ENERGY )
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single electron'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)


# Disable tracker
process.famosSimHits.MaterialEffects.Bremsstrahlung = False
process.famosSimHits.MaterialEffects.NuclearInteraction = False
process.famosSimHits.MaterialEffects.PairProduction = False
process.famosSimHits.MaterialEffects.MuonBremsstrahlung = False
process.famosSimHits.MaterialEffects.MultipleScattering = False
process.famosSimHits.MaterialEffects.EnergyLoss = False

# Disable vertex smearing
process.VertexSmearing = cms.Sequence()

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.simulationWithFamos)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from FastSimulation.Configuration.MixingModule_Full2Fast
from FastSimulation.Configuration.MixingModule_Full2Fast import setVertexGeneratorPileUpProducer 

#call to customisation function setVertexGeneratorPileUpProducer imported from FastSimulation.Configuration.MixingModule_Full2Fast
process = setVertexGeneratorPileUpProducer(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions

