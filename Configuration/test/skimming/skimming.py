import ROOT
import os
import sys

IDOpt = "tight"
debugOpt = False
onlytxtOpt = False

ROOT.gROOT.ProcessLine(".L %s/src/MLShEx_CorrL1Muons/Configuration/test/skimming/skimming.C+" % os.environ['CMSSW_BASE'])

inputFileName = "%s/src/MLShEx_CorrL1Muons/Configuration/python/outputL1uGMTAnalyzer.root" % os.environ['CMSSW_BASE']
fileExists = os.path.isfile(inputFileName)
if fileExists:
     inputChain = ROOT.TChain("events")
     inputChain.Add(inputFileName)
     print "Got file with",inputChain.GetEntries(),"entries\n"
     outOpt = "skimmedL1uGMTAnalyzer.txt"
     #Running
     L = ROOT.skimming(inputChain)
     L.Loop(IDOpt, outOpt, debugOpt, onlytxtOpt);
else: print "Input file doesn't exist..."
print ""
