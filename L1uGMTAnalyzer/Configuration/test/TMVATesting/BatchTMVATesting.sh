cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd $CMSSW_BASE/src/L1uGMTAnalyzer/Configuration/test/TMVA_PhiExtrapolation_New
root -l -b -q TMVARegReader_Impr_Runner.C