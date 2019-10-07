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

Exercise
* Understand plugins/L1uGMTAnalyzer.cc.
* Understand python/ana.py.
* Run
cd python
cmsRun ana.py
* Inspect outputL1uGMTAnalyzer.root.

* Build ../plugins/L1uGMTAnalyzer.cc:
After line:
#include "DataFormats/L1Trigger/interface/Muon.h"
add lines:
#include "DataFormats/PatCandidates/interface/Muon.h”

#include "DataFormats/VertexReco/interface/Vertex.h"

After line:
edm::EDGetTokenT<l1t::MuonBxCollection>      candToken_;
add lines:
edm::InputTag                                recoTag_;
edm::EDGetTokenT< std::vector<pat::Muon> >   recoToken_;
edm::InputTag                                PVTag_;
edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;

After line:
 
vector<float> muon_hwQual;
add lines:
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

Replace line:
candToken_( consumes<l1t::MuonBxCollection>(candTag_)) {
with lines:
candToken_( consumes<l1t::MuonBxCollection>(candTag_)),
recoTag_( iConfig.getParameter<edm::InputTag>("RecoTag") ),
recoToken_( consumes< std::vector<pat::Muon> >(recoTag_)),
PVTag_( iConfig.getParameter<edm::InputTag>("PVTag") ),
PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ) {

After line:
using namespace l1t;
add line:
using namespace reco;

After line:
muon_hwPt.clear(); muon_hwEta.clear(); muon_hwEtaAtVtx.clear(); muon_hwPhi.clear(); muon_hwPhiAtVtx.clear(); muon_hwCharge.clear(); muon_hwChargeValid.clear(); muon_hwQual.clear();
add line:
recomuon_pt.clear(); recomuon_eta.clear(); recomuon_phi.clear(); recomuon_isTight.clear(); recomuon_isMedium.clear(); recomuon_isLoose.clear(); recomuon_isSoft.clear(); recomuon_isGlobalMuon.clear(); recomuon_isTrackerMuon.clear(); recomuon_isStandAloneMuon.clear();

Before line:
events->Fill();
add lines:
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

After line:
events->Branch("muon_hwQual",&muon_hwQual);
add lines:
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
* Run
cd $CMSSW_BASE/src
scram b -j 4
cd MLShEx_CorrL1Muons/Configuration/python

* Build ana.py:
After line:
CandTag = cms.InputTag( 'gmtStage2Digis','Muon' ),
add lines:
RecoTag = cms.InputTag( 'slimmedMuons' ),
PVTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
* Run
cmsRun ana.py

* Inspect outputL1uGMTAnalyzer.root.
* Inspect crabConfig.py.

* Understand ../test/skimming.py, ../test/skimming.h, ../test/skimming.C.
* Run
cd ../test/skimming
python skimming.py
* Inspect skimmedL1uGMTAnalyzer.root.

* Discussion on the TMVA.

* Build ../../plugins/L1uGMTAnalyzer.cc:
After line:
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
add lines:
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

After line:
edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;
add lines:
PropagateToMuon muPropagator1st_;
PropagateToMuon muPropagator2nd_;

After line:
vector<bool> recomuon_isStandAloneMuon;
add lines:
vector<float> recomuon_etaAtSt1;
vector<float> recomuon_phiAtSt1;
vector<float> recomuon_etaAtSt2;
vector<float> recomuon_phiAtSt2;

Replace line:
PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ) {
with lines:
PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ),
muPropagator1st_( iConfig.getParameter<edm::ParameterSet>("muProp1st") ),
muPropagator2nd_( iConfig.getParameter<edm::ParameterSet>("muProp2nd") ) {

After line:
recomuon_pt.clear(); recomuon_eta.clear(); recomuon_phi.clear(); recomuon_isTight.clear(); recomuon_isMedium.clear(); recomuon_isLoose.clear(); recomuon_isSoft.clear(); recomuon_isGlobalMuon.clear(); recomuon_isTrackerMuon.clear(); recomuon_isStandAloneMuon.clear();
add line:
recomuon_etaAtSt1.clear(); recomuon_phiAtSt1.clear(); recomuon_etaAtSt2.clear(); recomuon_phiAtSt2.clear();

After line:
iEvent.getByToken(recoToken_, recoMuons);
add lines:
muPropagator1st_.init(iSetup);
muPropagator2nd_.init(iSetup);

After line:
recomuon_pt.push_back(muon->pt()); recomuon_eta.push_back(muon->eta()); recomuon_phi.push_back(muon->phi());
add lines:
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

After line:
events->Branch("recomuon_phi",&recomuon_phi);
add lines:
events->Branch("recomuon_etaAtSt1",&recomuon_etaAtSt1);
events->Branch("recomuon_phiAtSt1",&recomuon_phiAtSt1);
events->Branch("recomuon_etaAtSt2",&recomuon_etaAtSt2);
events->Branch("recomuon_phiAtSt2",&recomuon_phiAtSt2);
* Run
cd $CMSSW_BASE/src
scram b -j 4
cd MLShEx_CorrL1Muons/Configuration/python

* Build ana.py:
After line:
PVTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
add lines:
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

* Build ../test/skimming/skimming.h:
After line:
vector<float>   *recomuon_phi;
add lines:
vector<float>   *recomuon_etaAtSt1;
vector<float>   *recomuon_phiAtSt1;
vector<float>   *recomuon_etaAtSt2;
vector<float>   *recomuon_phiAtSt2;

After line:
TBranch        *b_recomuon_phi;
add lines:
TBranch        *b_recomuon_etaAtSt1;   //!
TBranch        *b_recomuon_phiAtSt1;   //!
TBranch        *b_recomuon_etaAtSt2;   //!
TBranch        *b_recomuon_phiAtSt2;   //!

After line:
recomuon_phi = 0;
add lines:
recomuon_etaAtSt1 = 0;
recomuon_phiAtSt1 = 0;
recomuon_etaAtSt2 = 0;
recomuon_phiAtSt2 = 0;

After line:
fChain->SetBranchAddress("recomuon_phi", &recomuon_phi, &b_recomuon_phi);
add lines:
fChain->SetBranchAddress("recomuon_etaAtSt1", &recomuon_etaAtSt1, &b_recomuon_etaAtSt1);
fChain->SetBranchAddress("recomuon_phiAtSt1", &recomuon_phiAtSt1, &b_recomuon_phiAtSt1);
fChain->SetBranchAddress("recomuon_etaAtSt2", &recomuon_etaAtSt2, &b_recomuon_etaAtSt2);
fChain->SetBranchAddress("recomuon_phiAtSt2", &recomuon_phiAtSt2, &b_recomuon_phiAtSt2);

* Build ../test/skimming/skimming.C:
Replace line:
float dr = DR( recomuon_eta->at(i),recomuon_phi->at(i),muon_eta->at(j),muon_phi->at(j) );
with lines:
float dr = -9999;
if (recomuon_phiAtSt2->at(i) != -9999) dr = DR( recomuon_eta->at(i),recomuon_phiAtSt2->at(i),muon_eta->at(j),muon_phi->at(j) );
else if (recomuon_phiAtSt1->at(i) != -9999) dr = DR( recomuon_eta->at(i),recomuon_phiAtSt1->at(i),muon_eta->at(j),muon_phi->at(j) );
else NotMatchedAtSt = true;

* Understand ../TMVATraining/TMVATraining.py, ../TMVATraining/TMVARegression.C.
* Run
cd ../TMVATraining
python TMVARegression.py

* Understand ../TMVATesting/Resolutions.py, ../TMVATesting/Resolutions.h, ../TMVATesting/Resolutions.C.
* Understand ../TMVATesting/MassSpectrum.py, ../TMVATesting/MassSpectrum.h, ../TMVATesting/MassSpectrum.C.
* Change trainingDir in ../TMVATesting/Resolutions.h and ../TMVATesting/MassSpectrum.h.
* Run
cd ../TMVATesting
python Resolutions.py --dataset ZeroBias --year 2017 --era BCEF
python Resolutions.py --dataset Charmonium --year 2017 --era B

python MassSpectrum.py --dataset ZeroBias --year 2017 --era BCEF
python MassSpectrum.py --dataset Charmonium --year 2017 --era B
* Run
python plotter.py --test Resolutions —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/Resolutions_ZeroBias_BCEF_Eta_A_A.root
python plotter.py --test Resolutions —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/Resolutions_Charmonium_B_Eta_A_A.root

python plotter.py --test MassSpectrum —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/MassSpectrum_ZeroBias_BCEF_JPsi_Eta.root
python plotter.py --test MassSpectrum —dir $CMSSW_BASE/src/MLShEx_CorrL1Muons/Configuration/test/TMVATesting/MassSpectrum_Charmonium_B_JPsi_Eta.root
