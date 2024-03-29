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
guysOpt = "A"
performOnOpt = "A"
debugOpt = False

print "MVA Testing: Resolution\n"
ROOT.gROOT.ProcessLine(".L %s/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/Resolutions.C+" % CMSSW_BASE)
inputFileName = inDirOpt+"L1uGMTAnalyzer_Trees/L1toRecoMatchPlots_"+options.datasetOpt+options.yearOpt+"_"+IDOpt+"_"+options.eraOpt +".root"
fileExists = os.path.isfile(inputFileName)
if fileExists:
    inputChain = ROOT.TChain("mytree")
    inputChain.Add(inputFileName)
    print "Got file with",inputChain.GetEntries(),"entries\n"

    outFileName = "./Resolutions_"+options.datasetOpt+"_"+options.eraOpt+"_"+TFBinMethodOpt+"_"+guysOpt+"_"+performOnOpt+".root"
    outFile = ROOT.TFile(outFileName,"recreate")
    print "Created output file: "+outFileName

    L = ROOT.Resolutions(inputChain)
    L.Loop(outFile,guysOpt,performOnOpt,TFBinMethodOpt,debugOpt);
    outFile.Write();
    outFile.Close();
    print "Output file written and closed!\n"

else: print "Input file doesn't exist..."
