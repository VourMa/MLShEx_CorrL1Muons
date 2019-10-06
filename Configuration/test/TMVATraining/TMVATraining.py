import ROOT
import os
import sys

CMSSW_BASE = os.environ['CMSSW_BASE']
datasetOpt = "ZeroBias"
yearOpt = "2018"
IDOpt = "tight"
TFOpt = "B,O,E"
totalErasOpt = "ABC"
useTotalErasOpt = True
guysOpt = "A"
extraTextOpt = ""
TFBinMethodOpt = "Eta"

TFs = TFOpt.split(',')
guys = guysOpt.split(',')
for TF in TFs:
    for guy in guys:
        print ".x {CMSSW_BASE}/src/MLShEx_CorrL1Muons/Configuration/test/TMVATraining/TMVARegression.C(\"{datasetArg}\",\"{yearArg}\",\"{IDArg}\",\"{TFArg}\",\"{erasArg}\",{useTotalErasArg},\"{guyArg}\",\"{extraTextArg}\",\"{etaOrIndexArg}\")".format(CMSSW_BASE=CMSSW_BASE, datasetArg=datasetOpt, yearArg=yearOpt, IDArg=IDOpt, TFArg=TF, erasArg=totalErasOpt, useTotalErasArg=int(useTotalErasOpt), guyArg=guy, extraTextArg=extraTextOpt, etaOrIndexArg=TFBinMethodOpt)
        ROOT.gROOT.ProcessLine(".x {CMSSW_BASE}/src/MLShEx_CorrL1Muons/Configuration/test/TMVATraining/TMVARegression.C(\"{datasetArg}\",\"{yearArg}\",\"{IDArg}\",\"{TFArg}\",\"{erasArg}\",{useTotalErasArg},\"{guyArg}\",\"{extraTextArg}\",\"{etaOrIndexArg}\")".format(CMSSW_BASE=CMSSW_BASE, datasetArg=datasetOpt, yearArg=yearOpt, IDArg=IDOpt, TFArg=TF, erasArg=totalErasOpt, useTotalErasArg=int(useTotalErasOpt), guyArg=guy, extraTextArg=extraTextOpt, etaOrIndexArg=TFBinMethodOpt))
