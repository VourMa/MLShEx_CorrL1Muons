import ROOT
import os
import sys

IDOpt = "tight"
debugOpt = False
onlytxtOpt = False

ROOT.gROOT.ProcessLine(".L %s/src/MLShEx_CorrL1Muons/Configuration/test/skimming/skimming_L1Reco.C+" % os.environ['CMSSW_BASE'])
#ROOT.gROOT.ProcessLine(".L %s/src/MLShEx_CorrL1Muons/Configuration/test/skimming/skimming_Prop.C+" % os.environ['CMSSW_BASE'])

inputFileName = "%s/src/MLShEx_CorrL1Muons/Configuration/python/outputL1uGMTAnalyzer_L1Reco.root" % os.environ['CMSSW_BASE']
#inputFileName = "%s/src/MLShEx_CorrL1Muons/Configuration/python/outputL1uGMTAnalyzer_Prop.root" % os.environ['CMSSW_BASE']
fileExists = os.path.isfile(inputFileName)
if fileExists:
     inputChain = ROOT.TChain("events")
     inputChain.Add(inputFileName)
     print "Got file with",inputChain.GetEntries(),"entries\n"
     outOpt = "skimmedL1uGMTAnalyzer_L1Reco.txt"
     #outOpt = "skimmedL1uGMTAnalyzer_Prop.txt"
     #Running
     L = ROOT.skimming(inputChain)
     L.Loop(IDOpt, outOpt, debugOpt, onlytxtOpt);
else: print "Input file doesn't exist..."
print ""
