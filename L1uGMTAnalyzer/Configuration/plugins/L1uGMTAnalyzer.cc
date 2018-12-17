// Class:      L1uGMTAnalyzer
//
// Original Author:  Manos Vourliotis
//     Inspired by:  L1uGMTAnalyzer by David Sperka
//         Created:  Mon, 17 Dec 2018 
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// track extrapolation
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"


#include "TTree.h"
//
// class declaration
//

class L1uGMTAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit L1uGMTAnalyzer(const edm::ParameterSet&);
      ~L1uGMTAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      bool isLooseMuon(const pat::Muon *mu);
      bool isMediumMuon(const pat::Muon *mu);
      bool isTightMuon(const pat::Muon *mu, reco::Vertex::Point PV);
      //bool RecoHLTMatching(const edm::Event& iEvent, double recoeta, double recophi, edm::InputTag IT_filter, double dRmatching);


      // ----------member data ---------------------------

      TTree *events;

      //input tag identifying the product containing muons
      edm::InputTag                                candTag_;
      edm::EDGetTokenT<l1t::MuonBxCollection>      candToken_;
      edm::InputTag                                recoTag_;
      edm::EDGetTokenT< std::vector<pat::Muon> >   recoToken_;
      edm::InputTag                                PVTag_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;
      //edm::InputTag                                trigobjectsTag_;
      //edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsToken_;


      /// Quality codes:
      /// to be updated with new L1 quality definitions
      int qualityBitMask_;    
      /// use central bx only muons
      bool centralBxOnly_;
      PropagateToMuon muPropagator1st_;
      PropagateToMuon muPropagator2nd_;

      //TFileService
      edm::Service<TFileService> fs;
      
      int counter = 0;

      ULong64_t Run, Event, Lumi, BX, Orbit;

      vector<float> muon_BX;
      vector<float> muon_hwPt;
      vector<float> muon_hwEta;
      vector<float> muon_hwPhi;
      vector<float> muon_hwCharge;
      vector<float> muon_hwChargeValid;
      vector<float> muon_hwQual;
   
      vector<float> muon_pt;
      vector<float> muon_eta;
      vector<float> muon_phi;
      vector<float> muon_charge;
      vector<float> muon_tfMuonIndex;
   
      vector<float> muon_hwEtaAtVtx;
      vector<float> muon_hwPhiAtVtx;
      vector<float> muon_EtaAtVtx;
      vector<float> muon_PhiAtVtx;
      
      vector<float> recomuon_pt;
      vector<float> recomuon_eta;
      vector<float> recomuon_phi;
      vector<float> recomuon_etaAtSt1;
      vector<float> recomuon_phiAtSt1;
      vector<float> recomuon_etaAtSt2;
      vector<float> recomuon_phiAtSt2;
      vector<float> recomuon_charge;
      vector<int> recomuon_quality;
      vector<bool> recomuon_isTight;
      vector<bool> recomuon_isMedium;
      vector<bool> recomuon_isLoose;
      vector<bool> recomuon_isSoft;
      vector<bool> recomuon_isGlobalMuon;
      vector<bool> recomuon_isTrackerMuon;
      vector<bool> recomuon_isStandAloneMuon;
      
      //bool L1Seed_SingleMu;
      
   
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1uGMTAnalyzer::L1uGMTAnalyzer(const edm::ParameterSet& iConfig)
 :
  candTag_( iConfig.getParameter<edm::InputTag>("CandTag") ),
  candToken_( consumes<l1t::MuonBxCollection>(candTag_)),
  recoTag_( iConfig.getParameter<edm::InputTag>("RecoTag") ),
  recoToken_( consumes< std::vector<pat::Muon> >(recoTag_)),
  PVTag_( iConfig.getParameter<edm::InputTag>("PVTag") ),
  PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ),
  //trigobjectsTag_( iConfig.getParameter<edm::InputTag> ("triggerobjects")),
  //trigobjectsToken_( consumes<trigger::TriggerEvent>(trigobjectsTag_)),
  centralBxOnly_( iConfig.getParameter<bool>("CentralBxOnly") ),
  muPropagator1st_( iConfig.getParameter<edm::ParameterSet>("muProp1st") ),
  muPropagator2nd_( iConfig.getParameter<edm::ParameterSet>("muProp2nd") )
{


    events = new TTree("events","events");

    //set the quality bit mask
    qualityBitMask_ = 0;
    vector<int> selectQualities = iConfig.getParameter<vector<int> >("SelectQualities");
    for(int selectQualitie : selectQualities){
//     if(selectQualities[i] > 7){  // FIXME: this will be updated once we have info from L1
//       throw edm::Exception(edm::errors::Configuration) << "QualityBits must be smaller than 8!";
//     }
        qualityBitMask_ |= 1<<selectQualitie;
    }

}


L1uGMTAnalyzer::~L1uGMTAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1uGMTAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace trigger;
   using namespace l1t;
   using namespace reco;
   
   //cout<< "Event Counter = " << counter << endl;
   counter++;

   Run = counter;
   Run = iEvent.id().run();
   Event = iEvent.id().event();
   Lumi = iEvent.id().luminosityBlock();
   BX = iEvent.bunchCrossing();
   Orbit = iEvent.orbitNumber();

   muon_BX.clear();
   muon_hwPt.clear();
   muon_hwEta.clear();
   muon_hwPhi.clear();
   muon_hwCharge.clear();
   muon_hwChargeValid.clear();
   muon_hwQual.clear();
   
   muon_pt.clear();
   muon_eta.clear();
   muon_phi.clear();
   muon_charge.clear();
   muon_tfMuonIndex.clear();
   
   muon_hwEtaAtVtx.clear();
   muon_hwPhiAtVtx.clear();
   muon_EtaAtVtx.clear();
   muon_PhiAtVtx.clear();
   
   recomuon_pt.clear();
   recomuon_eta.clear();
   recomuon_phi.clear();
   recomuon_etaAtSt1.clear();
   recomuon_phiAtSt1.clear();
   recomuon_etaAtSt2.clear();
   recomuon_phiAtSt2.clear();
   recomuon_quality.clear();
   recomuon_isTight.clear();
   recomuon_isMedium.clear();
   recomuon_isLoose.clear();
   recomuon_isSoft.clear();
   recomuon_isGlobalMuon.clear();
   recomuon_isTrackerMuon.clear();
   recomuon_isStandAloneMuon.clear();
   
   //L1Seed_SingleMu = 0;


   Handle<MuonBxCollection> allMuons;
   iEvent.getByToken(candToken_, allMuons);

   //cout << "L1 muons" << endl;
   for (int ibx = allMuons->getFirstBX(); ibx <= allMuons->getLastBX(); ++ibx) {
       if (centralBxOnly_ && (ibx != 0)) continue;
       for (auto it = allMuons->begin(ibx); it != allMuons->end(ibx); it++){

           //MuonRef muon(allMuons, distance(allMuons->begin(allMuons->getFirstBX()),it) );
           auto muon = it;

           muon_BX.push_back(ibx);
           muon_hwPt.push_back(muon->hwPt());
           muon_hwEta.push_back(muon->hwEta());
           muon_hwPhi.push_back(muon->hwPhi());
           muon_hwCharge.push_back(muon->hwCharge());
           muon_hwChargeValid.push_back(muon->hwChargeValid());
           muon_hwQual.push_back(muon->hwQual());

           muon_pt.push_back(muon->pt());
           muon_eta.push_back(muon->eta());
           muon_phi.push_back(muon->phi());
           muon_charge.push_back(muon->charge());
           muon_tfMuonIndex.push_back(muon->tfMuonIndex());
           
           //cout << "pt = " << muon_pt.back() << endl;
           //cout << "eta = " << muon_eta.back() << endl;
           //cout << "phi = " << muon_phi.back() << endl;
           //cout << endl;
           
           muon_hwEtaAtVtx.push_back(muon->hwEtaAtVtx());
           muon_hwPhiAtVtx.push_back(muon->hwPhiAtVtx());
           muon_EtaAtVtx.push_back(muon->etaAtVtx());
           muon_PhiAtVtx.push_back(muon->phiAtVtx());

       }
   }
   
   Handle< vector<pat::Muon> > recoMuons;
   iEvent.getByToken(recoToken_, recoMuons);
   edm::Handle<std::vector<Vertex> > theVertices;
   iEvent.getByToken(PV_token,theVertices) ;
   int nvertex = theVertices->size();
   Vertex::Point PV(0,0,0);
   /*const */reco::Vertex TheVertex;
   if( nvertex) {
	   PV = theVertices->begin()->position();
	   TheVertex = * (theVertices->begin());
   }
   else return;
   
   
   muPropagator1st_.init(iSetup);
   muPropagator2nd_.init(iSetup);
   
   //cout << "RECO muons" << endl;
   for( vector<pat::Muon>::const_iterator muon = (*recoMuons).begin(); muon != (*recoMuons).end(); muon++ ) {   
	   recomuon_pt.push_back(muon->pt());
	   recomuon_eta.push_back(muon->eta());
	   recomuon_phi.push_back(muon->phi());
	   
	   //cout<< "pt = " << recomuon_pt.back() << endl;
	   TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*muon);
	   if (stateAtMuSt1.isValid()) {
		   recomuon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta());
		   recomuon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
	   } else {
		   recomuon_etaAtSt1.push_back(-9999);
		   recomuon_phiAtSt1.push_back(-9999);
	   }
	   //cout<< "eta(St1) = " << recomuon_etaAtSt1.back() << endl;
	   //cout<< "phi(St1) = " << recomuon_phiAtSt1.back() << endl;
	   TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*muon);
	   if (stateAtMuSt2.isValid()) {
		   recomuon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta());
		   recomuon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
	   } else {
		   recomuon_etaAtSt2.push_back(-9999);
		   recomuon_phiAtSt2.push_back(-9999);
	   }
	   //cout<< "eta(St2) = " << recomuon_etaAtSt2.back() << endl;
	   //cout<< "phi(St2) = " << recomuon_phiAtSt2.back() << endl;
	   
	   //recomuon_isTight.push_back(isTightMuon(&*muon,PV)); //Custom
	   recomuon_isTight.push_back( muon->isTightMuon(TheVertex) ); //Official
	   //recomuon_isMedium.push_back(isMediumMuon(&*muon)); //Custom
	   recomuon_isMedium.push_back( muon->isMediumMuon() ); //Official
	   //recomuon_isLoose.push_back(isLooseMuon(&*muon)); //Custom
	   recomuon_isLoose.push_back( muon->isLooseMuon() ); //Official
	   recomuon_isSoft.push_back( muon->isSoftMuon(TheVertex) );
	   recomuon_quality.push_back(muon->bestTrack()->qualityMask());
	   recomuon_isGlobalMuon.push_back(muon->isGlobalMuon());
	   recomuon_isTrackerMuon.push_back(muon->isTrackerMuon());
	   recomuon_isStandAloneMuon.push_back(muon->isStandAloneMuon());
	   
	   //if ( RecoHLTMatching(iEvent,muon->eta(),muon->phi(),InputTag("L1_SingleMu7","","HLT"),0.4) ) L1Seed_SingleMu = 1;
   }
   
   events->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
L1uGMTAnalyzer::beginJob()
{

    events->Branch("counter",&counter);
    
    events->Branch("Run",&Run,"Run/l");
    events->Branch("Event",&Event,"Event/l");
    events->Branch("Lumi",&Lumi,"Lumi/l");
    events->Branch("BX",&BX,"BX/l");
    events->Branch("Orbit",&Orbit,"Orbit/l");

    events->Branch("muon_BX",&muon_BX);
    events->Branch("muon_hwPt",&muon_hwPt);
    events->Branch("muon_hwEta",&muon_hwEta);
    events->Branch("muon_hwPhi",&muon_hwPhi);
    events->Branch("muon_hwCharge",&muon_hwCharge);
    events->Branch("muon_hwChargeValid",&muon_hwChargeValid);
    events->Branch("muon_hwQual",&muon_hwQual);
  
    events->Branch("muon_pt",&muon_pt);
    events->Branch("muon_eta",&muon_eta);
    events->Branch("muon_phi",&muon_phi);
    events->Branch("muon_charge",&muon_charge);
    events->Branch("muon_tfMuonIndex",&muon_tfMuonIndex);
  
    events->Branch("muon_hwEtaAtVtx",&muon_hwEtaAtVtx);
    events->Branch("muon_hwPhiAtVtx",&muon_hwPhiAtVtx);
    events->Branch("muon_EtaAtVtx",&muon_EtaAtVtx);
    events->Branch("muon_PhiAtVtx",&muon_PhiAtVtx);
    
    events->Branch("recomuon_pt",&recomuon_pt);
    events->Branch("recomuon_eta",&recomuon_eta);
    events->Branch("recomuon_phi",&recomuon_phi);
    events->Branch("recomuon_etaAtSt1",&recomuon_etaAtSt1);
    events->Branch("recomuon_phiAtSt1",&recomuon_phiAtSt1);
    events->Branch("recomuon_etaAtSt2",&recomuon_etaAtSt2);
    events->Branch("recomuon_phiAtSt2",&recomuon_phiAtSt2);
    events->Branch("recomuon_quality",&recomuon_quality);
    events->Branch("recomuon_isTightMuon",&recomuon_isTight);
    events->Branch("recomuon_isMediumMuon",&recomuon_isMedium);
    events->Branch("recomuon_isLooseMuon",&recomuon_isLoose);
    events->Branch("recomuon_isSoftMuon",&recomuon_isSoft);
    events->Branch("recomuon_isGlobalMuon",&recomuon_isGlobalMuon);
    events->Branch("recomuon_isTrackerMuon",&recomuon_isTrackerMuon);
    events->Branch("recomuon_isStandAloneMuon",&recomuon_isStandAloneMuon);
    
    //events->Branch("L1Seed_SingleMu",&L1Seed_SingleMu);

}

// ------------ method called once each job just after ending the event loop  ------------
void
L1uGMTAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1uGMTAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


bool
L1uGMTAnalyzer::isLooseMuon(const pat::Muon *mu){
	bool flag = false ;
	if(mu->isPFMuon() && (mu->isGlobalMuon() || mu->isTrackerMuon())) flag = true;

	return flag;
}

bool
L1uGMTAnalyzer::isMediumMuon(const pat::Muon *mu){
	bool goodGlob = mu->isGlobalMuon()
		&& mu->globalTrack()->normalizedChi2() < 3
		&& mu->combinedQuality().chi2LocalPosition < 12
		&& mu->combinedQuality().trkKink < 20; 
	//reco::Muon& recomu = (reco::Muon&)mu;
	bool isMedium = isLooseMuon(mu)
		&& mu->innerTrack()->validFraction() > 0.49 
		&& mu->segmentCompatibility() > (goodGlob ? 0.303 : 0.451);

	return isMedium; 
}

bool
L1uGMTAnalyzer::isTightMuon(const pat::Muon *mu, reco::Vertex::Point PV){
	bool isTight = mu->isGlobalMuon() && mu->isPFMuon()
		&& mu->globalTrack()->normalizedChi2() < 10.
		&& mu->globalTrack()->hitPattern().numberOfValidMuonHits() > 0
		&& mu->numberOfMatchedStations() > 1
		&& fabs(mu->muonBestTrack()->dxy(PV)) < 0.2
		&& fabs(mu->muonBestTrack()->dz(PV)) < 0.5
		&& mu->innerTrack()->hitPattern().numberOfValidPixelHits() > 0
		&& mu->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
		&& mu->globalTrack()->normalizedChi2() < 1;

	return isTight;
}

//bool L1uGMTAnalyzer::RecoHLTMatching(const edm::Event& iEvent, double recoeta, double recophi, edm::InputTag IT_filter, double dRmatching) {
	//edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
	//iEvent.getByToken(trigobjectsToken_ ,triggerObjectsSummary);
	//trigger::TriggerObjectCollection selectedObjects;
	
	//if (triggerObjectsSummary.isValid()) {
		//size_t filterIndex = (*triggerObjectsSummary).filterIndex(IT_filter);
		//trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();
		//if (filterIndex < (*triggerObjectsSummary).sizeFilters()) {
			//const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
			//for (size_t j = 0; j < keys.size(); j++) {
				//trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
				//if(deltaR(recoeta,recophi,foundObject.eta(), foundObject.phi() ) < dRmatching )return true;
			//}
		//}
	//}
	//return false;
//}

//define this as a plug-in
DEFINE_FWK_MODULE(L1uGMTAnalyzer);
