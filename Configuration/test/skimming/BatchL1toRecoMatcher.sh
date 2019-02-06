cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd $CMSSW_BASE/src/L1uGMTAnalyzer/Configuration/test/skimming/
#root -l -b -q "L1toRecoMatcher_Runner.C(\"tight\",false)" #ZeroBias
root -l -b -q "L1toRecoMatcher_Runner.C(\"tight\",true)" #Charm