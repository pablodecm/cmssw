import FWCore.ParameterSet.Config as cms

# set default max events
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.setDefault('maxEvents', -1 )
options.parseArguments()

process = cms.Process("Validation")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

# Reasonable logging level
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Load DQM
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.Core.DQM_cfg")
process.load("DQMOffline.RecoB.bTagSequences_cff")

# Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents))

# Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

# MiniAOD collections 
primaryVertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
secondaryVertexCollection = cms.InputTag('slimmedSecondaryVertices')
jetCollection = cms.InputTag("slimmedJets")
genJetCollection  = cms.InputTag("slimmedGenJets")
candidateCollection = cms.InputTag("packedPFCandidates")
muonCollection = cms.InputTag("slimmedMuons")
electronCollection = cms.InputTag("slimmedElectrons")
genParticleCollection = cms.InputTag("prunedGenParticles")

# Load b-tagging modules
process.load("RecoBTag.Configuration.RecoBTag_cff")

# use pfBtagging sequence
process.btagSeq = cms.Sequence(process.pfBTagging)

# remove IVF modules which do not have to be re-run
process.btagSeq.remove(process.inclusiveCandidateVertexFinder)
process.btagSeq.remove(process.candidateVertexMerger)
process.btagSeq.remove(process.candidateVertexArbitrator)
process.btagSeq.remove(process.inclusiveCandidateSecondaryVertices)

# For MC-based pileup jet ID
process.ak4GenJetsForPUid = cms.EDFilter("GenJetSelector",
    src = genJetCollection, 
    cut = cms.string('pt > 8.'),
    filter = cms.bool(False)
)

# Gen matching for PAT jets
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
process.patJetGenJetMatch.src = jetCollection
process.patJetGenJetMatch.matched = cms.InputTag("ak4GenJetsForPUid")
process.patJetGenJetMatch.maxDeltaR = cms.double(0.25)
process.patJetGenJetMatch.resolveAmbiguities = cms.bool(True)

# Load the jet flavor(DQM uses the old parton-based flavor definition) 
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
process.AK4byRef.jets = jetCollection
process.flavourSeq = cms.Sequence(
    process.myPartons *
    process.AK4Flavour
)

# Load b-tag validation
process.load("Validation.RecoB.bTagAnalysis_cfi")
# Some common plot parameters 
flavPlots = "allbcl" 
ptRanges = cms.vdouble(50.0,80.0,120.0)
etaRanges = cms.vdouble(0.0,1.4,2.4)

# Specify taggers for which produce the validation plots
from DQMOffline.RecoB.bTagCommon_cff import *
# tagger configuration (from bTagCommon)
tagConfig = cms.VPSet(
    cms.PSet(
        bTagTrackIPAnalysisBlock,
        type = cms.string('CandIP'),
        label = cms.InputTag("pfImpactParameterTagInfos"),
        folder = cms.string("IPTag")
    ),
    cms.PSet(
        bTagCombinedSVAnalysisBlock,
        ipTagInfos = cms.InputTag("pfImpactParameterTagInfos"),
        type = cms.string('GenericMVA'),
        svTagInfos = cms.InputTag("pfSecondaryVertexTagInfos"),
        label = cms.InputTag("candidateCombinedSecondaryVertexComputer"),
        folder = cms.string("CSVTag")
    ),
    cms.PSet(
        bTagTrackCountingAnalysisBlock,
        label = cms.InputTag("pfTrackCountingHighEffBJetTags"),
        folder = cms.string("TCHE")
    ),
    cms.PSet(
        bTagTrackCountingAnalysisBlock,
        label = cms.InputTag("pfTrackCountingHighPurBJetTags"),
        folder = cms.string("TCHP")
    ),
    cms.PSet(
        bTagProbabilityAnalysisBlock,
        label = cms.InputTag("pfJetProbabilityBJetTags"),
        folder = cms.string("JP")
    ),
    cms.PSet(
        bTagBProbabilityAnalysisBlock,
        label = cms.InputTag("pfJetBProbabilityBJetTags"),
        folder = cms.string("JBP")
    ),
    cms.PSet(
        bTagSimpleSVAnalysisBlock,
        label = cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"),
        folder = cms.string("SSVHE")
    ),
    cms.PSet(
        bTagSimpleSVAnalysisBlock,
        label = cms.InputTag("pfSimpleSecondaryVertexHighPurBJetTags"),
        folder = cms.string("SSVHP")
    ),
    #cms.PSet(
        #bTagGenericAnalysisBlock,
        #label = cms.InputTag("combinedSecondaryVertexBJetTags"),
        #folder = cms.string("CSV_tkOnly")
    #),
    cms.PSet(
        bTagGenericAnalysisBlock,
        label = cms.InputTag("pfCombinedSecondaryVertexBJetTags"),
        folder = cms.string("CSV")
    ),
    cms.PSet(
        bTagGenericAnalysisBlock,
        label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
        folder = cms.string("CSVv2")
    ),
    cms.PSet(
        bTagSoftLeptonAnalysisBlock,
        label = cms.InputTag("softPFMuonBJetTags"),
        folder = cms.string("SMT")
    ),
    cms.PSet(
        bTagSoftLeptonAnalysisBlock,
        label = cms.InputTag("softPFElectronBJetTags"),
        folder = cms.string("SET")
    ),
)    

# Validation configuration
process.bTagValidation.jetMCSrc = 'AK4byValAlgo'
process.bTagValidation.applyPtHatWeight = False
process.bTagValidation.genJetsMatched = cms.InputTag("patJetGenJetMatch")
process.bTagValidation.softLeptonInfo = cms.InputTag("softPFElectronsTagInfos")
process.bTagValidation.doPUid = cms.bool(True)
process.bTagValidation.flavPlots = flavPlots
process.bTagValidation.ptRanges = ptRanges
process.bTagValidation.etaRanges = etaRanges
process.bTagValidation.tagConfig = tagConfig
# Harvesting parameters (same as validation)
process.bTagHarvestMC.flavPlots = flavPlots
process.bTagHarvestMC.ptRanges = ptRanges
process.bTagHarvestMC.etaRanges = etaRanges
process.bTagHarvestMC.tagConfig = tagConfig

# Setup saving DQM parameters 
process.dqmEnv.subSystemFolder = 'BTAG'
process.dqmSaver.producer = 'DQM'
process.dqmSaver.workflow = '/POG/BTAG/BJET'
process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd =cms.untracked.bool(True)
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)

# Adapt module configurations to MiniAOD input
process.pfImpactParameterTagInfos.jets = jetCollection
process.pfImpactParameterTagInfos.primaryVertex = primaryVertexCollection
process.pfImpactParameterTagInfos.candidates = candidateCollection
process.pfInclusiveSecondaryVertexFinderTagInfos.extSVCollection = secondaryVertexCollection
process.myPartons.src = genParticleCollection
process.softPFMuonsTagInfos.jets = jetCollection
process.softPFMuonsTagInfos.vertex = primaryVertexCollection # temporary solution
#process.softPFMuonsTagInfos.primaryVertex = primaryVertexCollection # this will work in the final 740 release
process.softPFMuonsTagInfos.muons = muonCollection
process.softPFElectronsTagInfos.jets = jetCollection
process.softPFElectronsTagInfos.primaryVertex = primaryVertexCollection
process.softPFElectronsTagInfos.electrons = electronCollection

# Run Path
process.dqm = cms.Path( process.btagSeq * process.ak4GenJetsForPUid * process.patJetGenJetMatch * process.flavourSeq * process.bTagValidation * process.bTagHarvestMC * process.dqmSaver)

# Set filenames
process.PoolSource.fileNames = [ 
    '/store/relval/CMSSW_7_4_0_pre8/RelValTTbar_13/MINIAODSIM/PU50ns_MCRUN2_74_V6-v1/00000/90A385BF-60BD-E411-85D8-0025905938A4.root',
    '/store/relval/CMSSW_7_4_0_pre8/RelValTTbar_13/MINIAODSIM/PU50ns_MCRUN2_74_V6-v1/00000/A8E0D9B1-60BD-E411-BFC1-0025905A6080.root'
    ] 


