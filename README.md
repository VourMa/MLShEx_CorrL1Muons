# CMS Data Analysis School (DAS) 2019, Beijing - Machine Learning Short Exercise: Using NN to Correct L1 Muons

## Setup
* _cmsrel CMSSW_10_2_11_  
   If there is an architecture problem, set SCRAM_ARCH=slc7_amd64_gcc700.
* _cd CMSSW_10_2_11/_
* _cmsenv_
* _git cms-init_
* _cd src/_
* _git clone git@github.com:VourMa/MLShEx_CorrL1Muons.git MLShEx_CorrL1Muons --single-branch --branch exercise_
* _scram b -j 4_
* _cd MLShEx_CorrL1Muons/Configuration/_

## Exercise
* Understand [_plugins/L1uGMTAnalyzer.cc_](../exercise/Configuration/plugins/L1uGMTAnalyzer.cc).
* Understand [_python/ana.py_](../exercise/Configuration/python/ana.py).
* Run  
   _cd python_  
   _cmsRun ana.py_

* Inspect _outputL1uGMTAnalyzer.root_.
* Inspect [_crabConfig.py_](../exercise/Configuration/python/crabConfig.py).

* Understand [_../test/skimming/skimming.py_](../exercise/Configuration/test/skimming/skimming.py), [_../test/skimming/skimming.h_](../exercise/Configuration/test/skimming/skimming.h), [_../test/skimming/skimming.C_](../exercise/Configuration/test/skimming/skimming.C).
* Run  
   _cd ../test/skimming_  
   _python skimming.py_
* Inspect _skimmedL1uGMTAnalyzer.root_.

* Discussion on the TMVA.

* **Build _../../plugins/L1uGMTAnalyzer.cc_**:  
   After line:
   ```c++
   #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
   ```
   add lines:
   ```c++
   #include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
   #include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
   ```

   After line:
   ```c++
   edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;
   ```
   add lines:
   ```c++
   PropagateToMuon muPropagator1st_;
   PropagateToMuon muPropagator2nd_;
   ```

   After line:
   ```c++
   vector<bool> recomuon_isStandAloneMuon;
   ```
   add lines:
   ```c++
   vector<float> recomuon_etaAtSt1;
   vector<float> recomuon_phiAtSt1;
   vector<float> recomuon_etaAtSt2;
   vector<float> recomuon_phiAtSt2;
   ```

   Replace line:
   ```c++
   PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ) {
   ```
   with lines:
   ```c++
   PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ),
   muPropagator1st_( iConfig.getParameter<edm::ParameterSet>("muProp1st") ),
   muPropagator2nd_( iConfig.getParameter<edm::ParameterSet>("muProp2nd") ) {
   ```

   After line:
   ```c++
   recomuon_pt.clear(); recomuon_eta.clear(); recomuon_phi.clear(); recomuon_isTight.clear(); recomuon_isMedium.clear(); recomuon_isLoose.clear(); recomuon_isSoft.clear(); recomuon_isGlobalMuon.clear(); recomuon_isTrackerMuon.clear(); recomuon_isStandAloneMuon.clear();
   ```
   add line:
   ```c++
   recomuon_etaAtSt1.clear(); recomuon_phiAtSt1.clear(); recomuon_etaAtSt2.clear(); recomuon_phiAtSt2.clear();
   ```

   After line:
   ```c++
   iEvent.getByToken(recoToken_, recoMuons);
   ```
   add lines:
   ```c++
   muPropagator1st_.init(iSetup);
   muPropagator2nd_.init(iSetup);
   ```

   After line:
   ```c++
   recomuon_pt.push_back(muon->pt()); recomuon_eta.push_back(muon->eta()); recomuon_phi.push_back(muon->phi());
   ```
   add lines:
   ```c++
   TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*muon);
   if (stateAtMuSt1.isValid()) {
      recomuon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta()); recomuon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
   } else {
      recomuon_etaAtSt1.push_back(-9999); recomuon_phiAtSt1.push_back(-9999);
   }

   TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*muon);
   if (stateAtMuSt2.isValid()) {
      recomuon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta()); recomuon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
   } else {
      recomuon_etaAtSt2.push_back(-9999); recomuon_phiAtSt2.push_back(-9999);
   }
   ```

   After line:
   ```c++
   events->Branch("recomuon_phi",&recomuon_phi);
   ```
   add lines:
   ```c++
   events->Branch("recomuon_etaAtSt1",&recomuon_etaAtSt1);
   events->Branch("recomuon_phiAtSt1",&recomuon_phiAtSt1);
   events->Branch("recomuon_etaAtSt2",&recomuon_etaAtSt2);
   events->Branch("recomuon_phiAtSt2",&recomuon_phiAtSt2);
   ```
* Run  
   _cd $CMSSW_BASE/src_  
   _scram b -j 4_  
   _cd MLShEx_CorrL1Muons/Configuration/python_

* **Build _ana.py_**:  
   After line:
   ```c++
   PVTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
   ```
   add lines:
   ```c++
   #muon track extrapolation to 1nd station
   muProp1st = cms.PSet(
      useTrack = cms.string("tracker"),  # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
      useState = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
      useSimpleGeometry = cms.bool(True),
      useStation2 = cms.bool(False),
   ),
   #muon track extrapolation to 2nd station
   muProp2nd = cms.PSet(
      useTrack = cms.string("tracker"),  # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
      useState = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
      useSimpleGeometry = cms.bool(True),
      useStation2 = cms.bool(True),
      fallbackToME1 = cms.bool(False),
   ),
   ```

* **Build _../test/skimming/skimming.h_**:  
   After line:
   ```c++
   vector<float>   *recomuon_phi;
   ```
   add lines:
   ```c++
   vector<float>   *recomuon_etaAtSt1;
   vector<float>   *recomuon_phiAtSt1;
   vector<float>   *recomuon_etaAtSt2;
   vector<float>   *recomuon_phiAtSt2;
   ```

   After line:
   ```c++
   TBranch        *b_recomuon_phi;
   ```
   add lines:
   ```c++
   TBranch        *b_recomuon_etaAtSt1;   //!
   TBranch        *b_recomuon_phiAtSt1;   //!
   TBranch        *b_recomuon_etaAtSt2;   //!
   TBranch        *b_recomuon_phiAtSt2;   //!
   ```

   After line:
   ```c++
   recomuon_phi = 0;
   ```
   add lines:
   ```c++
   recomuon_etaAtSt1 = 0;
   recomuon_phiAtSt1 = 0;
   recomuon_etaAtSt2 = 0;
   recomuon_phiAtSt2 = 0;
   ```

   After line:
   ```c++
   fChain->SetBranchAddress("recomuon_phi", &recomuon_phi, &b_recomuon_phi);
   ```
   add lines:
   ```c++
   fChain->SetBranchAddress("recomuon_etaAtSt1", &recomuon_etaAtSt1, &b_recomuon_etaAtSt1);
   fChain->SetBranchAddress("recomuon_phiAtSt1", &recomuon_phiAtSt1, &b_recomuon_phiAtSt1);
   fChain->SetBranchAddress("recomuon_etaAtSt2", &recomuon_etaAtSt2, &b_recomuon_etaAtSt2);
   fChain->SetBranchAddress("recomuon_phiAtSt2", &recomuon_phiAtSt2, &b_recomuon_phiAtSt2);
   ```

* **Build _../test/skimming/skimming.C_**:  
   Replace line:
   ```c++
   float dr = DR( recomuon_eta->at(i),recomuon_phi->at(i),muon_eta->at(j),muon_phi->at(j) );
   ```
   with lines:
   ```c++
   float dr = -9999;
   if (recomuon_phiAtSt2->at(i) != -9999) dr = DR( recomuon_eta->at(i),recomuon_phiAtSt2->at(i),muon_eta->at(j),muon_phi->at(j) );
   else if (recomuon_phiAtSt1->at(i) != -9999) dr = DR( recomuon_eta->at(i),recomuon_phiAtSt1->at(i),muon_eta->at(j),muon_phi->at(j) );
   else NotMatchedAtSt = true;
   ```

* Understand [_../test/TMVATraining/TMVATraining.py_](../exercise/Configuration/test/TMVATraining/TMVATraining.py), [_../test/TMVATraining/TMVARegression.C_](../exercise/Configuration/test/TMVATraining/TMVARegression.C).
* Run  
   _cd ../test/TMVATraining_  
   _python TMVARegression.py_
* Discuss on the results.

* Understand [_../TMVATesting/Resolutions.py_](../exercise/Configuration/test/TMVATesting/Resolutions.py), [_../TMVATesting/Resolutions.h_](../exercise/Configuration/test/TMVATesting/Resolutions.h), [_../TMVATesting/Resolutions.C_](../exercise/Configuration/test/TMVATesting/Resolutions.C).
* Understand [_../TMVATesting/MassSpectrum.py_](../exercise/Configuration/test/TMVATesting/MassSpectrum.py), [_../TMVATesting/MassSpectrum.h_](../exercise/Configuration/test/TMVATesting/MassSpectrum.h), [_../TMVATesting/MassSpectrum.C_](../exercise/Configuration/test/TMVATesting/MassSpectrum.C).
* **Change _trainingDir_ in _../TMVATesting/Resolutions.h_ and _../TMVATesting/MassSpectrum.h_**.
* Run  
   _cd ../TMVATesting_  
   _python Resolutions.py --dataset ZeroBias --year 2017 --era BCEF_  
   _python Resolutions.py --dataset Charmonium --year 2017 --era B_  

   _python MassSpectrum.py --dataset ZeroBias --year 2017 --era BCEF_  
   _python MassSpectrum.py --dataset Charmonium --year 2017 --era B_  
* Understand [_../plotter.py_](../exercise/Configuration/test/plotter/plotter.py).
* Run  
   _cd ../plotter_  
   _python plotter.py --test Resolutions —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/Resolutions_ZeroBias_BCEF_Eta_A_A.root_  
   _python plotter.py --test Resolutions —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/Resolutions_Charmonium_B_Eta_A_A.root_  

   _python plotter.py --test MassSpectrum —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/MassSpectrum_ZeroBias_BCEF_JPsi_Eta.root_  
   _python plotter.py --test MassSpectrum —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/MassSpectrum_Charmonium_B_JPsi_Eta.root_  
