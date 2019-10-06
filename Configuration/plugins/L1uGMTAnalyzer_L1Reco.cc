// Class:      L1uGMTAnalyzer_L1Reco
//
// Original Author:  Manos Vourliotis
//         Created:  Mon, 26 Aug 2019 
//


// System include files
#include <memory>

// User include files
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

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


// Class declaration
class L1uGMTAnalyzer_L1Reco : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit L1uGMTAnalyzer_L1Reco(const edm::ParameterSet&);
		~L1uGMTAnalyzer_L1Reco();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		
		
		// Class variables
		TTree *events;
		
		edm::InputTag                                candTag_;
		edm::EDGetTokenT<l1t::MuonBxCollection>      candToken_;
		edm::InputTag                                recoTag_;
		edm::EDGetTokenT< std::vector<pat::Muon> >   recoToken_;
		edm::InputTag                                PVTag_;
		edm::EDGetTokenT<std::vector<reco::Vertex> > PV_token;
		
		edm::Service<TFileService> fs;
		
		int counter = 0;
		ULong64_t Run, Event, Lumi, BX, Orbit;
		
		vector<float> muon_BX;
		
		vector<float> muon_pt;
		vector<float> muon_eta;
		vector<float> muon_etaAtVtx;
		vector<float> muon_phi;
		vector<float> muon_phiAtVtx;
		vector<float> muon_charge;
		vector<float> muon_tfMuonIndex;
		
		vector<float> muon_hwPt;
		vector<float> muon_hwEta;
		vector<float> muon_hwEtaAtVtx;
		vector<float> muon_hwPhi;
		vector<float> muon_hwPhiAtVtx;
		vector<float> muon_hwCharge;
		vector<float> muon_hwChargeValid;
		vector<float> muon_hwQual;

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
};

// Constructor
L1uGMTAnalyzer_L1Reco::L1uGMTAnalyzer_L1Reco(const edm::ParameterSet& iConfig) :
candTag_( iConfig.getParameter<edm::InputTag>("CandTag") ),
candToken_( consumes<l1t::MuonBxCollection>(candTag_)),
recoTag_( iConfig.getParameter<edm::InputTag>("RecoTag") ),
recoToken_( consumes< std::vector<pat::Muon> >(recoTag_)),
PVTag_( iConfig.getParameter<edm::InputTag>("PVTag") ),
PV_token( consumes<std::vector<reco::Vertex> > (PVTag_) ) {
    events = new TTree("events","events");
}

// Destructor
L1uGMTAnalyzer_L1Reco::~L1uGMTAnalyzer_L1Reco() {
}

// Method called for each event
void L1uGMTAnalyzer_L1Reco::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace std;
	using namespace edm;
	using namespace l1t;
	using namespace reco;
	
	counter++;
	
	Run = iEvent.id().run();
	Event = iEvent.id().event();
	Lumi = iEvent.id().luminosityBlock();
	BX = iEvent.bunchCrossing();
	Orbit = iEvent.orbitNumber();
	
	muon_BX.clear();

	muon_pt.clear(); muon_eta.clear(); muon_etaAtVtx.clear(); muon_phi.clear(); muon_phiAtVtx.clear(); muon_charge.clear(); muon_tfMuonIndex.clear();

	muon_hwPt.clear(); muon_hwEta.clear(); muon_hwEtaAtVtx.clear(); muon_hwPhi.clear(); muon_hwPhiAtVtx.clear(); muon_hwCharge.clear(); muon_hwChargeValid.clear(); muon_hwQual.clear();
	
	recomuon_pt.clear(); recomuon_eta.clear(); recomuon_phi.clear(); recomuon_isTight.clear(); recomuon_isMedium.clear(); recomuon_isLoose.clear(); recomuon_isSoft.clear(); recomuon_isGlobalMuon.clear(); recomuon_isTrackerMuon.clear(); recomuon_isStandAloneMuon.clear();


	Handle<MuonBxCollection> L1Muons;
	iEvent.getByToken(candToken_, L1Muons);
	
	for (int ibx = L1Muons->getFirstBX(); ibx <= L1Muons->getLastBX(); ++ibx) {
		for (auto it = L1Muons->begin(ibx); it != L1Muons->end(ibx); it++) {
			auto muon = it;
		
			muon_BX.push_back(ibx);
		
			muon_pt.push_back(muon->pt()); muon_eta.push_back(muon->eta()); muon_etaAtVtx.push_back(muon->etaAtVtx()); muon_phi.push_back(muon->phi()); muon_phiAtVtx.push_back(muon->phiAtVtx()); muon_charge.push_back(muon->charge()); muon_tfMuonIndex.push_back(muon->tfMuonIndex());

			muon_hwPt.push_back(muon->hwPt()); muon_hwEta.push_back(muon->hwEta()); muon_hwEtaAtVtx.push_back(muon->hwEtaAtVtx()); muon_hwPhi.push_back(muon->hwPhi()); muon_hwPhiAtVtx.push_back(muon->hwPhiAtVtx()); muon_hwCharge.push_back(muon->hwCharge()); muon_hwChargeValid.push_back(muon->hwChargeValid()); muon_hwQual.push_back(muon->hwQual());
		}
	}
	

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

	events->Fill();
}

// Method called once, just before starting the event loop 
void L1uGMTAnalyzer_L1Reco::beginJob() {
	events->Branch("counter",&counter);
	
	events->Branch("Run",&Run,"Run/l");
	events->Branch("Event",&Event,"Event/l");
	events->Branch("Lumi",&Lumi,"Lumi/l");
	events->Branch("BX",&BX,"BX/l");
	events->Branch("Orbit",&Orbit,"Orbit/l");
	
	events->Branch("muon_BX",&muon_BX);
	
	events->Branch("muon_pt",&muon_pt);
	events->Branch("muon_eta",&muon_eta);
	events->Branch("muon_etaAtVtx",&muon_etaAtVtx);
	events->Branch("muon_phi",&muon_phi);
	events->Branch("muon_phiAtVtx",&muon_phiAtVtx);
	events->Branch("muon_charge",&muon_charge);
	events->Branch("muon_tfMuonIndex",&muon_tfMuonIndex);
	
	events->Branch("muon_hwPt",&muon_hwPt);
	events->Branch("muon_hwEta",&muon_hwEta);
	events->Branch("muon_hwEtaAtVtx",&muon_hwEtaAtVtx);
	events->Branch("muon_hwPhi",&muon_hwPhi);
	events->Branch("muon_hwPhiAtVtx",&muon_hwPhiAtVtx);
	events->Branch("muon_hwCharge",&muon_hwCharge);
	events->Branch("muon_hwChargeValid",&muon_hwChargeValid);
	events->Branch("muon_hwQual",&muon_hwQual);
    
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
}

// Method called once, just after ending the event loop  ------------
void L1uGMTAnalyzer_L1Reco::endJob() {
}

// Method that fills 'descriptions' with the allowed parameters for the module
void
L1uGMTAnalyzer_L1Reco::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1uGMTAnalyzer_L1Reco);
