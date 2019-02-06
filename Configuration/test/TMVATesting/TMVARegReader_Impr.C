#define TMVARegReader_Impr_cxx
#include "TMVARegReader_Impr.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

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

bool YouWantNoCloseMuons = false; //with
double DR_CloseMuons = 0.4;

void TMVARegReader_Impr::Loop(TH1D * h_DphiL1Reco,TH1D * h_PhiRecoReg,TH1D * h_DphiRegReco, TH1D * h_DphiExtReco)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	//cout << "Number of entries:" << nentries << endl;
	//nentries = 100;
	
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
					if( DR(L1muon_eta->at(i),L1muon_phi->at(i),L1muon_eta->at(j),L1muon_phi->at(j)) < DR_CloseMuons ) YouGotCloseMuons = true;
				}
			}
			if( YouGotCloseMuons ) continue;
			
		}
		
		//cout << L1muon_pt->size() << endl;
		for(int i = 0; i<L1muon_ptCorr->size(); i++) {
			//Cuts
			//cout << recomuon_dr->at(i) << endl;
			if(recomuon_dr->at(i) < 0.0 || recomuon_dr->at(i) > 0.2) continue;
			//if(fabs(L1muon_eta->at(i)) < 1.2) continue; // TESTING
			//if( !( (L1muon_tfMuonIndex->at(i) >= 0 && L1muon_tfMuonIndex->at(i) <= 16) || (L1muon_tfMuonIndex->at(i) >= 90 && L1muon_tfMuonIndex->at(i) <= 107) ) ) continue; // TESTING
			
			
			L1muon_ptCorr_ = (float) L1muon_ptCorr->at(i); L1muon_eta_ = (float) L1muon_eta->at(i); L1muon_phi_ = (float) L1muon_phiAtVtx->at(i); L1muon_charge_ = (float) L1muon_charge->at(i);
			
			double DphiL1Reco, PhiRecoReg;
			
			//if( L1muon_tfMuonIndex->at(i) >= 36 && L1muon_tfMuonIndex->at(i) <= 70 )  {
			if( fabs(L1muon_eta_) < 0.8 )  {
				DphiL1Reco = (readerBMTF->EvaluateRegression( "MLP" ))[0];
				PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
			}
			//if( (L1muon_tfMuonIndex->at(i) >= 17 && L1muon_tfMuonIndex->at(i) <= 35) || (L1muon_tfMuonIndex->at(i) >= 71 && L1muon_tfMuonIndex->at(i) <= 89) ) {
			if( fabs(L1muon_eta_) > 0.8 && fabs(L1muon_eta_) < 1.2 ) {
				DphiL1Reco = (readerOMTF->EvaluateRegression( "MLP" ))[0];
				PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
			}
			//if( (L1muon_tfMuonIndex->at(i) >= 0 && L1muon_tfMuonIndex->at(i) <= 16) || (L1muon_tfMuonIndex->at(i) >= 90 && L1muon_tfMuonIndex->at(i) <= 107) ) {
			if( fabs(L1muon_eta_) > 1.2 ) {
				DphiL1Reco = (readerEMTF->EvaluateRegression( "MLP" ))[0];
				PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
			}
			
			h_DphiL1Reco->Fill( DphiL1Reco );
			h_PhiRecoReg->Fill( PhiRecoReg );
			h_DphiRegReco->Fill( DPhi( PhiRecoReg,recomuon_phi->at(i) ) );
			h_DphiExtReco->Fill( DPhi( L1muon_phiAtVtx->at(i),recomuon_phi->at(i) ) );
		}
	}
	
	return;
}