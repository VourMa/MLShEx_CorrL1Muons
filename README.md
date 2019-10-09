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
* Understand _plugins/L1uGMTAnalyzer.cc_.
* Understand _python/ana.py_.
* Run  
   _cd python_  
   _cmsRun ana.py_
* Inspect _outputL1uGMTAnalyzer.root_.

* **Build _../plugins/L1uGMTAnalyzer.cc_**:  
   After line:
   ```c++
   #include "DataFormats/L1Trigger/interface/Muon.h"
   ```
   add lines:
   ```c++
   #include "DataFormats/PatCandidates/interface/Muon.h”
   
   #include "DataFormats/VertexReco/interface/Vertex.h"
   ```

   After line:
   ```c++
   edm::EDGetTokenT<l1t::MuonBxCollection>      candToken_;
   ```
   add lines:
   ```c++
   edm::InputTag                                recoTag_;
   edm::EDGetTokenT< std::vector<pat::Muon> >   recoToken_;
   edm::InputTag                                PVTag_;
   edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;
   ```

   After line:
   ```c++
   vector<float> muon_hwQual;
   ```
   add lines:
   ```c++
   vector<float> recomuon_pt;
   vector<float> recomuon_eta;
   vector<float> recomuon_phi;
   vector<float> recomuon_charge;
   vector<bool> recomuon_isTight;
   vector<bool> recomuon_isMedium;
   vector<bool> recomuon_isLoose;
   vector<bool> recomuon_isSoft;
   vector<bool> recomuon_isGlobalMuon;
   vector<bool> recomuon_isTrackerMuon;
   vector<bool> recomuon_isStandAloneMuon;
   ```

   Replace line:
   ```c++
   candToken_( consumes<l1t::MuonBxCollection>(candTag_)) {
   ```
   with lines:
   ```c++
   candToken_( consumes<l1t::MuonBxCollection>(candTag_)),
   recoTag_( iConfig.getParameter<edm::InputTag>("RecoTag") ),
   recoToken_( consumes< std::vector<pat::Muon> >(recoTag_)),
   PVTag_( iConfig.getParameter<edm::InputTag>("PVTag") ),
   PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ) {
   ```

   After line:
   ```c++
   using namespace l1t;
   ```
   add line:
   ```c++
   using namespace reco;
   ```

   After line:
   ```c++
   muon_hwPt.clear(); muon_hwEta.clear(); muon_hwEtaAtVtx.clear(); muon_hwPhi.clear(); muon_hwPhiAtVtx.clear(); muon_hwCharge.clear(); muon_hwChargeValid.clear(); muon_hwQual.clear();
   ```
   add line:
   ```c++
   recomuon_pt.clear(); recomuon_eta.clear(); recomuon_phi.clear(); recomuon_isTight.clear(); recomuon_isMedium.clear(); recomuon_isLoose.clear(); recomuon_isSoft.clear(); recomuon_isGlobalMuon.clear(); recomuon_isTrackerMuon.clear(); recomuon_isStandAloneMuon.clear();
   ```

   Before line:
   ```c++
   events->Fill();
   ```
   add lines:
   ```c++
   edm::Handle<std::vector<Vertex> > theVertices;
   iEvent.getByToken(PV_token,theVertices) ;

   int nvertex = theVertices->size();
   Vertex::Point PV(0,0,0);
   reco::Vertex TheVertex;
   if( nvertex) {
      PV = theVertices->begin()->position();
      TheVertex = * (theVertices->begin());
   }
   else {
      recomuon_pt.push_back(-99); recomuon_eta.push_back(-99); recomuon_phi.push_back(-99); recomuon_isTight.push_back(-99); recomuon_isMedium.push_back(-99); recomuon_isLoose.push_back(-99); recomuon_isSoft.push_back(-99); recomuon_isGlobalMuon.push_back(-99); recomuon_isTrackerMuon.push_back(-99); recomuon_isStandAloneMuon.push_back(-99);
      events->Fill();
      return;
   }


   Handle< vector<pat::Muon> > recoMuons;
   iEvent.getByToken(recoToken_, recoMuons);

   for( vector<pat::Muon>::const_iterator muon = (*recoMuons).begin(); muon != (*recoMuons).end(); muon++ ) {
      recomuon_pt.push_back(muon->pt()); recomuon_eta.push_back(muon->eta()); recomuon_phi.push_back(muon->phi());

      recomuon_isTight.push_back( muon->isTightMuon(TheVertex) ); recomuon_isMedium.push_back( muon->isMediumMuon() ); recomuon_isLoose.push_back( muon->isLooseMuon() ); recomuon_isSoft.push_back( muon->isSoftMuon(TheVertex) ); recomuon_isGlobalMuon.push_back( muon->isGlobalMuon() ); recomuon_isTrackerMuon.push_back( muon->isTrackerMuon() ); recomuon_isStandAloneMuon.push_back( muon->isStandAloneMuon() );
   }
   ```

   After line:
   ```c++
   events->Branch("muon_hwQual",&muon_hwQual);
   ```
   add lines:
   ```c++
   events->Branch("recomuon_pt",&recomuon_pt);
   events->Branch("recomuon_eta",&recomuon_eta);
   events->Branch("recomuon_phi",&recomuon_phi);
   events->Branch("recomuon_isTightMuon",&recomuon_isTight);
   events->Branch("recomuon_isMediumMuon",&recomuon_isMedium);
   events->Branch("recomuon_isLooseMuon",&recomuon_isLoose);
   events->Branch("recomuon_isSoftMuon",&recomuon_isSoft);
   events->Branch("recomuon_isGlobalMuon",&recomuon_isGlobalMuon);
   events->Branch("recomuon_isTrackerMuon",&recomuon_isTrackerMuon);
   events->Branch("recomuon_isStandAloneMuon",&recomuon_isStandAloneMuon);
   ```
* Run  
   _cd $CMSSW_BASE/src_  
   _scram b -j 4_  
   _cd MLShEx_CorrL1Muons/Configuration/python_

* **Build _ana.py_**:  
   After line:
   ```c++
   CandTag = cms.InputTag( 'gmtStage2Digis','Muon' ),
   ```
   add lines:
   ```c++
   RecoTag = cms.InputTag( 'slimmedMuons' ),
   PVTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
   ```
* Run  
   _cmsRun ana.py_

* Inspect _outputL1uGMTAnalyzer.root_.
* Inspect _crabConfig.py_.

* Understand _../test/skimming.py_, _../test/skimming.h_, _../test/skimming.C_.
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

* Understand _../TMVATraining/TMVATraining.py_, _../TMVATraining/TMVARegression.C_.
* Run  
   _cd ../TMVATraining_  
   _python TMVARegression.py_

* Understand _../TMVATesting/Resolutions.py_, _../TMVATesting/Resolutions.h_, _../TMVATesting/Resolutions.C_.
* Understand _../TMVATesting/MassSpectrum.py_, _../TMVATesting/MassSpectrum.h_, _../TMVATesting/MassSpectrum.C_.
* **Change _trainingDir_ in _../TMVATesting/Resolutions.h_ and _../TMVATesting/MassSpectrum.h_**.
* Run  
   _cd ../TMVATesting_  
   _python Resolutions.py --dataset ZeroBias --year 2017 --era BCEF_  
   _python Resolutions.py --dataset Charmonium --year 2017 --era B_  

   _python MassSpectrum.py --dataset ZeroBias --year 2017 --era BCEF_  
   _python MassSpectrum.py --dataset Charmonium --year 2017 --era B_  
* Understand _plotter.py_.
* Run  
   _python plotter.py --test Resolutions —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/Resolutions_ZeroBias_BCEF_Eta_A_A.root_  
   _python plotter.py --test Resolutions —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/Resolutions_Charmonium_B_Eta_A_A.root_  

   _python plotter.py --test MassSpectrum —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/MassSpectrum_ZeroBias_BCEF_JPsi_Eta.root_  
   _python plotter.py --test MassSpectrum —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/MassSpectrum_Charmonium_B_JPsi_Eta.root_  
