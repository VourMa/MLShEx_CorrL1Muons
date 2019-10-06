// Class:      L1uGMTAnalyzer_OnlyL1
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

#include "TTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


// Class declaration
class L1uGMTAnalyzer_OnlyL1 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit L1uGMTAnalyzer_OnlyL1(const edm::ParameterSet&);
		~L1uGMTAnalyzer_OnlyL1();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		
		
		// Class variables
		TTree *events;
		
		edm::InputTag                                candTag_;
		edm::EDGetTokenT<l1t::MuonBxCollection>      candToken_;
		
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
};

// Constructor
L1uGMTAnalyzer_OnlyL1::L1uGMTAnalyzer_OnlyL1(const edm::ParameterSet& iConfig) :
candTag_( iConfig.getParameter<edm::InputTag>("CandTag") ),
candToken_( consumes<l1t::MuonBxCollection>(candTag_)) {
    events = new TTree("events","events");
}

// Destructor
L1uGMTAnalyzer_OnlyL1::~L1uGMTAnalyzer_OnlyL1() {
}

// Method called for each event
void L1uGMTAnalyzer_OnlyL1::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace std;
	using namespace edm;
	using namespace l1t;
	
	counter++;
	
	Run = iEvent.id().run();
	Event = iEvent.id().event();
	Lumi = iEvent.id().luminosityBlock();
	BX = iEvent.bunchCrossing();
	Orbit = iEvent.orbitNumber();
	
	muon_BX.clear();

	muon_pt.clear(); muon_eta.clear(); muon_etaAtVtx.clear(); muon_phi.clear(); muon_phiAtVtx.clear(); muon_charge.clear(); muon_tfMuonIndex.clear();

	muon_hwPt.clear(); muon_hwEta.clear(); muon_hwEtaAtVtx.clear(); muon_hwPhi.clear(); muon_hwPhiAtVtx.clear(); muon_hwCharge.clear(); muon_hwChargeValid.clear(); muon_hwQual.clear();
	
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
	
	events->Fill();
}

// Method called once, just before starting the event loop 
void L1uGMTAnalyzer_OnlyL1::beginJob() {
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
}

// Method called once, just after ending the event loop  ------------
void L1uGMTAnalyzer_OnlyL1::endJob() {
}

// Method that fills 'descriptions' with the allowed parameters for the module
void
L1uGMTAnalyzer_OnlyL1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1uGMTAnalyzer_OnlyL1);
