import ROOT
import os
import sys

CMSSW_BASE = os.environ['CMSSW_BASE']
datasetOpt = "ZeroBias"
yearOpt = "2018"
IDOpt = "tight"
TFOpt = "B,O,E"
totalErasOpt = "ABC"
guysOpt = "A"
extraTextOpt = ""
TFBinMethodOpt = "Eta"

TFs = TFOpt.split(',')
guys = guysOpt.split(',')
for TF in TFs:
    for guy in guys:
        print ".x {CMSSW_BASE}/src/MLShEx_CorrL1Muons/Configuration/test/TMVATraining/TMVARegression.C(\"{datasetArg}\",\"{yearArg}\",\"{IDArg}\",\"{TFArg}\",\"{erasArg}\",\"{guyArg}\",\"{extraTextArg}\",\"{etaOrIndexArg}\")".format(CMSSW_BASE=CMSSW_BASE, datasetArg=datasetOpt, yearArg=yearOpt, IDArg=IDOpt, TFArg=TF, erasArg=totalErasOpt, guyArg=guy, extraTextArg=extraTextOpt, etaOrIndexArg=TFBinMethodOpt)
        ROOT.gROOT.ProcessLine(".x {CMSSW_BASE}/src/MLShEx_CorrL1Muons/Configuration/test/TMVATraining/TMVARegression.C(\"{datasetArg}\",\"{yearArg}\",\"{IDArg}\",\"{TFArg}\",\"{erasArg}\",\"{guyArg}\",\"{extraTextArg}\",\"{etaOrIndexArg}\")".format(CMSSW_BASE=CMSSW_BASE, datasetArg=datasetOpt, yearArg=yearOpt, IDArg=IDOpt, TFArg=TF, erasArg=totalErasOpt, guyArg=guy, extraTextArg=extraTextOpt, etaOrIndexArg=TFBinMethodOpt))

ROOT.TMVA.TMVARegGui("TMVARegression_TFB_EraABC_GuysA_Eta.root")
raw_input("Press Enter to continue...")
