import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras


#process = cms.Process( "TEST",eras.Run2_2017 ) # For 2017 datasets
process = cms.Process( "TEST",eras.Run2_2018 ) # For 2018 datasets


process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


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
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v12', '') # For 2017 datasets
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2018_realistic_v10', '') # For 2018 datasets


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.myFilter = hlt.hltHighLevel.clone(
    #HLTPaths =[ "HLT_ZeroBias_v5", "HLT_ZeroBias_v6" ],
    HLTPaths =[ "HLT_Mu7p5_Track2_Jpsi_v*", "HLT_Mu7p5_Track3p5_Jpsi_v*", "HLT_Mu7p5_Track7_Jpsi_v*" ],
    throw = False
    ) # Filter to choose events only from specific HLT paths

process.Ana = cms.EDAnalyzer('L1uGMTAnalyzer',
    saveTags = cms.bool( True ),
    CentralBxOnly = cms.bool( False ),
    SelectQualities = cms.vint32(  ),
    CandTag = cms.InputTag( 'gmtStage2Digis','Muon' ),
    RecoTag = cms.InputTag( 'slimmedMuons' ), # MiniAOD
    #RecoTag = cms.InputTag( 'muons' ), # AOD
    PVTag = cms.InputTag("offlineSlimmedPrimaryVertices"), # MiniAOD
    #PVTag = cms.InputTag("offlinePrimaryVertices"), # AOD
    triggerobjects = cms.InputTag('hltTriggerSummaryAOD','','HLT'),
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
# 'root://xrootd-cms.infn.it//store/data/Run2017B/ZeroBias/MINIAOD/17Nov2017-v1/50000/EAED5891-5CD4-E711-BCB4-A4BF0112E3E8.root', # MINIAOD 2017
 'root://xrootd-cms.infn.it//store/data/Run2018A/MuOnia/MINIAOD/17Sep2018-v1/60000/FEC487B5-D47E-034E-8BE5-92D87627A58B.root', # MINIAOD 2018
 #'root://xrootd-cms.infn.it//store/data/Run2017F/ZeroBias/AOD/17Nov2017-v1/60001/2CA7FFAB-41E2-E711-B9E2-0025905A60FE.root' # AOD 2017

    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("outputL1uGMTAnalyzer.root")
)

process.p = cms.Path(process.Ana)
#process.p = cms.Path(process.myFilter*process.Ana) # Enable HLT filter
