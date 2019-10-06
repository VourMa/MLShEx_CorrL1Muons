import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras


process = cms.Process( "Analysis",eras.Run2_2018 ) # For 2018 datasets


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


#Standard
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#MuonPropagator
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2018_realistic_v10', '') # For 2018 datasets


process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.Ana = cms.EDAnalyzer('L1uGMTAnalyzer_Prop',
    saveTags = cms.bool( True ),
    CandTag = cms.InputTag( 'gmtStage2Digis','Muon' ),
    RecoTag = cms.InputTag( 'slimmedMuons' ),
    PVTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
    #muon track extrapolation to 1nd station
    muProp1st = cms.PSet(
        useTrack = cms.string("tracker"),  # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
        useState = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
        useSimpleGeometry = cms.bool(True),
        useStation2 = cms.bool(False),
    ),
    #muon track extrapolation to 2nd station
    muProp2nd = cms.PSet(
        useTrack = cms.string("tracker"),  # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
        useState = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
        useSimpleGeometry = cms.bool(True),
        useStation2 = cms.bool(True),
        fallbackToME1 = cms.bool(False),
    ),
)

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/data/Run2018C/ZeroBias/MINIAOD/17Sep2018-v1/80000/EAB199EB-692F-8042-BFC4-C1F42CF7FFE8.root', # MINIAOD 2018
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("outputL1uGMTAnalyzer_Prop.root")
)

process.p = cms.Path(process.Ana)
