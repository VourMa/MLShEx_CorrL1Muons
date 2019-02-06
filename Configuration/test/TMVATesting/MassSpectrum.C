#define MassSpectrum_cxx
#include "MassSpectrum.h"
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

float Mll(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2){
		TLorentzVector mu1, mu2;
		mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.104);
		mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.104);
		return (mu1+mu2).M();
}

void MassSpectrum::correctThePhi(float L1muon_eta_, double & DphiL1Reco, double & PhiRecoReg) {
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


void MassSpectrum::Loop(TH1D * h_recomll,TH1D * h_L1mll,TH1D * h_L1mllCorr, TH1D * h_L1mll_JPsi, TH1D * h_L1mllCorr_JPsi)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	//cout << "Number of entries:" << nentries << endl;
//	nentries = 10000;
	
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
//		if( ientry % 1000000 == 0 ) cout << "I have processed " << ientry << " events!" << endl;
		
		if( recomuon_N < 2 ) continue;
		
		for(int i = 0; i<recomuon_N; i++) {
			if(recomuon_dr->at(i) < 0.0 || recomuon_dr->at(i) > 0.2) continue;
			
			L1muon_ptCorr_ = (float) L1muon_ptCorr->at(i); L1muon_eta_ = (float) L1muon_eta->at(i); L1muon_phi_ = (float) L1muon_phiAtVtx->at(i); L1muon_charge_ = (float) L1muon_charge->at(i);
			double DphiL1Reco_i, PhiRecoReg_i;
			correctThePhi(L1muon_eta_, DphiL1Reco_i, PhiRecoReg_i);
			
			for(int j = i+1; j<recomuon_N; j++) {
				L1muon_ptCorr_ = (float) L1muon_ptCorr->at(j); L1muon_eta_ = (float) L1muon_eta->at(j); L1muon_phi_ = (float) L1muon_phiAtVtx->at(j); L1muon_charge_ = (float) L1muon_charge->at(j);
				double DphiL1Reco_j, PhiRecoReg_j;
				correctThePhi(L1muon_eta_, DphiL1Reco_j, PhiRecoReg_j);
				
				float recomll = Mll( recomuon_pt->at(i), recomuon_eta->at(i), recomuon_phi->at(i), recomuon_pt->at(j), recomuon_eta->at(j), recomuon_phi->at(j) );
				h_recomll->Fill(recomll);
				float L1mll = Mll(L1muon_ptCorr->at(i), L1muon_eta->at(i), L1muon_phi->at(i), L1muon_ptCorr->at(j), L1muon_eta->at(j), L1muon_phi->at(j));
				h_L1mll->Fill(L1mll);
				float L1mllCorr = Mll(L1muon_ptCorr->at(i), L1muon_eta->at(i), PhiRecoReg_i, L1muon_ptCorr->at(j), L1muon_eta->at(j), PhiRecoReg_j);
				h_L1mllCorr->Fill(L1mllCorr);
				if( recomll > 3.046 && recomll < 3.146 ) {
					h_L1mll_JPsi->Fill(L1mll);
					h_L1mllCorr_JPsi->Fill(L1mllCorr);
				}
			}
		}
	}
	return;
}