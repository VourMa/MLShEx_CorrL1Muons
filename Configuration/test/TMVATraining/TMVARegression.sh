#!/bin/bash
cd /afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/test/TMVATraining
root -l -b -q "TMVARegression.C(\"$1\",\"$2\",\"$3\")"
