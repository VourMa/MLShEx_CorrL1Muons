#define Resolutions_cxx
#include "Resolutions.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

bool YouWantNoCloseMuons = false; //with
double DR_CloseMuons = 0.4;


double DPhi(double phi1, double phi2) {
	double result = phi1 - phi2;
	if(result>3.14) result -= 6.28;
	if(result<=-3.14) result += 6.28;
	return result;
}

float DR(float eta1, float phi1, float eta2, float phi2) {
	float result = TMath::Sqrt((eta1-eta2)*(eta1-eta2)+DPhi(phi1,phi2)*DPhi(phi1,phi2));
	return result;
}


class Histos {
	private:
	TH1D * h_DphiL1Reco, * h_PhiRecoReg, * h_DphiRegReco, * h_DphiExtReco, * h_DphiRegReco_B, * h_DphiExtReco_B, * h_DphiRegReco_O, * h_DphiExtReco_O, * h_DphiRegReco_E, * h_DphiExtReco_E;
	
	public:
	void CreateHistos();
	void FillHistos(float DphiL1Reco, float PhiRecoReg, float L1muon_eta_, float L1PhiAtVtx, float RecoPhi);
};

void Histos::CreateHistos () {
	h_DphiL1Reco = new TH1D("h_DphiL1Reco", "#Delta(#phi_{L1},#phi_{reco})_{Reg}", 300, -1.5, 1.5);
	h_PhiRecoReg = new TH1D("h_PhiRecoReg", "#phi_{Reg}", 64, -3.2, 3.2);
	h_DphiRegReco = new TH1D("h_DphiRegReco", "#Delta(#phi_{Reg},#phi_{reco})", 80, -0.4, 0.4);
	h_DphiExtReco = new TH1D("h_DphiExtReco", "#Delta(#phi_{ext},#phi_{reco})", 80, -0.4, 0.4);
	h_DphiRegReco_B = new TH1D("h_DphiRegReco_B", "#Delta(#phi_{Reg},#phi_{reco}), Barrel", 80, -0.4, 0.4);
	h_DphiExtReco_B = new TH1D("h_DphiExtReco_B", "#Delta(#phi_{ext},#phi_{reco}), Barrel", 80, -0.4, 0.4);
	h_DphiRegReco_O = new TH1D("h_DphiRegReco_O", "#Delta(#phi_{Reg},#phi_{reco}), Overlap", 80, -0.4, 0.4);
	h_DphiExtReco_O = new TH1D("h_DphiExtReco_O", "#Delta(#phi_{ext},#phi_{reco}), Overlap", 80, -0.4, 0.4);
	h_DphiRegReco_E = new TH1D("h_DphiRegReco_E", "#Delta(#phi_{Reg},#phi_{reco}), Endcap", 80, -0.4, 0.4);
	h_DphiExtReco_E = new TH1D("h_DphiExtReco_E", "#Delta(#phi_{ext},#phi_{reco}), Endcap", 80, -0.4, 0.4);
}

void Histos::FillHistos(float DphiL1Reco, float PhiRecoReg, float L1muon_eta_, float L1PhiAtVtx, float RecoPhi) {
	h_DphiL1Reco->Fill( DphiL1Reco );
	h_PhiRecoReg->Fill( PhiRecoReg );
	h_DphiRegReco->Fill( DPhi( PhiRecoReg,RecoPhi ) );
	h_DphiExtReco->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	if( fabs(L1muon_eta_) < 0.8 ) {
		h_DphiRegReco_B->Fill( DPhi( PhiRecoReg,RecoPhi ) );
		h_DphiExtReco_B->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	}
	else if( fabs(L1muon_eta_) < 1.2 ) {
		h_DphiRegReco_O->Fill( DPhi( PhiRecoReg,RecoPhi ) );
		h_DphiExtReco_O->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	}
	else {
		h_DphiRegReco_E->Fill( DPhi( PhiRecoReg,RecoPhi ) );
		h_DphiExtReco_E->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	}
}


void Resolutions::correctThePhi(float L1muon_eta_, double & DphiL1Reco, double & PhiRecoReg) {
	if( fabs(L1muon_eta_) < 0.8 )  {
				DphiL1Reco = (readerBMTF->EvaluateRegression( "MLP" ))[0];
				PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
			}
			
			if( fabs(L1muon_eta_) > 0.8 && fabs(L1muon_eta_) < 1.2 ) {
				DphiL1Reco = (readerOMTF->EvaluateRegression( "MLP" ))[0];
				PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
			}
			
			if( fabs(L1muon_eta_) > 1.2 ) {
				DphiL1Reco = (readerEMTF->EvaluateRegression( "MLP" ))[0];
				PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
			}
}

void Resolutions::Loop(TFile * out, TString performOnWhichGuys, bool debug)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	//cout << "Number of entries:" << nentries << endl;
	if(debug) nentries = 10000;
	
	out->cd();
	
	//Definitions
	Histos hist;
	hist.CreateHistos();
	
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if( ientry % 1000000 == 0 ) cout << "I have processed " << ientry << " events!" << endl;
		
		if( YouWantNoCloseMuons ){
			bool YouGotCloseMuons = false;
			for(int i = 0; i<L1muon_ptCorr->size(); i++) {
				for(int j = i+1; j < L1muon_ptCorr->size(); j++) {
					if( DR(L1muon_eta->at(i),L1muon_phiAtVtx->at(i),L1muon_eta->at(j),L1muon_phiAtVtx->at(j)) < DR_CloseMuons ) YouGotCloseMuons = true;
				}
			}
			if( YouGotCloseMuons ) continue;
		}
		
		for(int i = 0; i<L1muon_ptCorr->size(); i++) {
			//Cuts
			if(recomuon_dr->at(i) < 0.0 || recomuon_dr->at(i) > 0.2) continue;
			if(performOnWhichGuys == "G") if(L1muon_pt->at(i) < 10 || L1muon_ptCorr->at(i) > 40) continue;
			if(performOnWhichGuys == "B") if(L1muon_pt->at(i) > 10 && L1muon_ptCorr->at(i) < 40) continue;
			
			L1muon_ptCorr_ = (float) L1muon_ptCorr->at(i); L1muon_eta_ = (float) L1muon_eta->at(i); L1muon_phi_ = (float) L1muon_phiAtVtx->at(i); L1muon_charge_ = (float) L1muon_charge->at(i);
			
			double DphiL1Reco, PhiRecoReg;
			correctThePhi(L1muon_eta_, DphiL1Reco, PhiRecoReg);
			
			hist.FillHistos(DphiL1Reco,PhiRecoReg,L1muon_eta_,L1muon_phiAtVtx->at(i),recomuon_phi->at(i));
		}
	}
	
	return;
}
