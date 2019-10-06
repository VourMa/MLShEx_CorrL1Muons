import ROOT
import os
import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--dataset", dest="datasetOpt", default=None, help="Choose dataset")
parser.add_option("--year", dest="yearOpt", default=None, help="Choose year")
parser.add_option("--era", dest="eraOpt", default=None, help="Choose eras, separated by commas")
options, args = parser.parse_args()

CMSSW_BASE = os.environ['CMSSW_BASE']
inDirOpt = "/eos/cms/store/cmst3/user/evourlio/"
IDOpt = "tight"
TFBinMethodOpt = "Eta"
performOnOpt = "A"
particleOpt = "JPsi"
debugOpt = False

print "MVA Testing: Mass spectrum\n"
ROOT.gROOT.ProcessLine(".L %s/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/MassSpectrum.C+" % CMSSW_BASE)
inputFileName = inDirOpt+"L1uGMTAnalyzer_Trees/L1toRecoMatchPlots_"+options.datasetOpt+options.yearOpt+"_"+IDOpt+"_"+options.eraOpt +".root"
fileExists = os.path.isfile(inputFileName)
if fileExists:
    inputChain = ROOT.TChain("mytree")
    inputChain.Add(inputFileName)
    print "Got file with",inputChain.GetEntries(),"entries\n"

    outFileName = "./MassSpectrum_"+options.datasetOpt+"_"+options.eraOpt+"_"+particleOpt+"_"+TFBinMethodOpt+".root"
    outFile = ROOT.TFile(outFileName,"recreate")
    print "Created output file: "+outFileName

    L = ROOT.MassSpectrum(inputChain)
    L.Loop(outFile,particleOpt,performOnOpt,TFBinMethodOpt,debugOpt);
    outFile.Write();
    outFile.Close();
    print "Output file written and closed!\n"

else: print "Input file doesn't exist..."
