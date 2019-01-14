cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd $CMSSW_BASE/src/L1uGMTAnalyzer/Configuration/test/TMVA_PhiExtrapolation_New
root -l -b -q TMVARegression_Impr_00_08.C
root -l -b -q TMVARegression_Impr_08_12.C
root -l -b -q TMVARegression_Impr_12_24.C