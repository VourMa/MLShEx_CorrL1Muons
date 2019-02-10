#define MassSpectrum_cxx
#include "MassSpectrum.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

const float JPsi_M = 3.096;

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

TString MatchWithTF(float eta) {
	if( fabs(eta) < 0.8 ) return "B";
	else if( fabs(eta) < 1.2 ) return "O";
	else return "E";
}


class TFHistos {
	private:
	TH1D * h_L1mllCorr_JPsi_BB, * h_L1mllCorr_JPsi_BO, * h_L1mllCorr_JPsi_BE, * h_L1mllCorr_JPsi_OO, * h_L1mllCorr_JPsi_OE, * h_L1mllCorr_JPsi_EE, * h_L1DmllCorr_o_mllCorr_JPsi_BB, * h_L1DmllCorr_o_mllCorr_JPsi_BO, * h_L1DmllCorr_o_mllCorr_JPsi_BE, * h_L1DmllCorr_o_mllCorr_JPsi_OO, * h_L1DmllCorr_o_mllCorr_JPsi_OE, * h_L1DmllCorr_o_mllCorr_JPsi_EE;
	
	public:
	void CreateHistos();
	void FillHistos(TString TF1, TString TF2, float mass);
};

void TFHistos::CreateHistos () {
	h_L1mllCorr_JPsi_BB = new TH1D("h_L1mllCorr_JPsi_BB", "M_{ll}(L1(Corrected),J/Psi) Barrel-Barrel", 100, 0.0, 10.0);
	h_L1mllCorr_JPsi_BO = new TH1D("h_L1mllCorr_JPsi_BO", "M_{ll}(L1(Corrected),J/Psi) Barrel-Overlap", 100, 0.0, 10.0);
	h_L1mllCorr_JPsi_BE = new TH1D("h_L1mllCorr_JPsi_BE", "M_{ll}(L1(Corrected),J/Psi) Barrel-Endcap", 100, 0.0, 10.0);
	h_L1mllCorr_JPsi_OO = new TH1D("h_L1mllCorr_JPsi_OO", "M_{ll}(L1(Corrected),J/Psi) Overlap-Overlap", 100, 0.0, 10.0);
	h_L1mllCorr_JPsi_OE = new TH1D("h_L1mllCorr_JPsi_OE", "M_{ll}(L1(Corrected),J/Psi) Overlap-Endcap", 100, 0.0, 10.0);
	h_L1mllCorr_JPsi_EE = new TH1D("h_L1mllCorr_JPsi_EE", "M_{ll}(L1(Corrected),J/Psi) Endcap-Endcap", 100, 0.0, 10.0);
	
	h_L1DmllCorr_o_mllCorr_JPsi_BB = new TH1D("h_L1DmllCorr_o_mllCorr_JPsi_BB", "( M_{ll}(L1(Corrected),J/Psi) - M(J/#psi) ) / M(J/#psi) Barrel-Barrel", 60, -1.0, 3.0);
	h_L1DmllCorr_o_mllCorr_JPsi_BO = new TH1D("h_L1DmllCorr_o_mllCorr_JPsi_BO", "( M_{ll}(L1(Corrected),J/Psi) - M(J/#psi) ) / M(J/#psi) Barrel-Overlap", 60, -1.0, 3.0);
	h_L1DmllCorr_o_mllCorr_JPsi_BE = new TH1D("h_L1DmllCorr_o_mllCorr_JPsi_BE", "( M_{ll}(L1(Corrected),J/Psi) - M(J/#psi) ) / M(J/#psi) Barrel-Endcap", 60, -1.0, 3.0);
	h_L1DmllCorr_o_mllCorr_JPsi_OO = new TH1D("h_L1DmllCorr_o_mllCorr_JPsi_OO", "( M_{ll}(L1(Corrected),J/Psi) - M(J/#psi) ) / M(J/#psi) Overlap-Overlap", 60, -1.0, 3.0);
	h_L1DmllCorr_o_mllCorr_JPsi_OE = new TH1D("h_L1DmllCorr_o_mllCorr_JPsi_OE", "( M_{ll}(L1(Corrected),J/Psi) - M(J/#psi) ) / M(J/#psi) Overlap-Endcap", 60, -1.0, 3.0);
	h_L1DmllCorr_o_mllCorr_JPsi_EE = new TH1D("h_L1DmllCorr_o_mllCorr_JPsi_EE", "( M_{ll}(L1(Corrected),J/Psi) - M(J/#psi) ) / M(J/#psi) Endcap-Endcap", 60, -1.0, 3.0);
}

void TFHistos::FillHistos(TString TF1, TString TF2, float mass) {
	if( TF1.Contains("B") ) {
		if( TF2.Contains("B") ) { h_L1mllCorr_JPsi_BB->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_BB->Fill( (mass - JPsi_M) / JPsi_M ); }
		else if( TF2.Contains("O") ) { h_L1mllCorr_JPsi_BO->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_BO->Fill( (mass - JPsi_M) / JPsi_M ); }
		else if( TF2.Contains("E") ) { h_L1mllCorr_JPsi_BE->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_BE->Fill( (mass - JPsi_M) / JPsi_M ); }
	}
	else if( TF1.Contains("O") ) {
		if( TF2.Contains("B") ) { h_L1mllCorr_JPsi_BO->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_BO->Fill( (mass - JPsi_M) / JPsi_M ); }
		else if( TF2.Contains("O") ) { h_L1mllCorr_JPsi_OO->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_OO->Fill( (mass - JPsi_M) / JPsi_M ); }
		else if( TF2.Contains("E") ) { h_L1mllCorr_JPsi_OE->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_OE->Fill( (mass - JPsi_M) / JPsi_M ); }
	}
	else if( TF1.Contains("E") ) {
		if( TF2.Contains("B") ) { h_L1mllCorr_JPsi_BE->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_BE->Fill( (mass - JPsi_M) / JPsi_M ); }
		else if( TF2.Contains("O") ) { h_L1mllCorr_JPsi_OE->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_OE->Fill( (mass - JPsi_M) / JPsi_M ); }
		else if( TF2.Contains("E") ) { h_L1mllCorr_JPsi_EE->Fill(mass); h_L1DmllCorr_o_mllCorr_JPsi_EE->Fill( (mass - JPsi_M) / JPsi_M ); }
	}
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


void MassSpectrum::Loop(TFile * out, bool SplitTFs, bool debug)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
//	cout << "Number of entries:" << nentries << endl;
	if(debug) nentries = 10000;
	
	out->cd();
	
	//Definitions
	TH1D * h_recomll = new TH1D("h_recomll", "M_{ll}(reco)", 300, 0.0, 30.0);
	TH1D * h_L1mll = new TH1D("h_L1mll", "M_{ll}(L1)", 300, 0.0, 30.0);
	TH1D * h_L1mllCorr = new TH1D("h_L1mllCorr", "M_{ll}(L1(Corrected))", 300, 0.0, 30.0);
	
	TH1D * h_L1mll_JPsi = new TH1D("h_L1mll_JPsi", "M_{ll}(L1,J/Psi)", 100, 0.0, 10.0);
	TH1D * h_L1Dmll_o_mll_JPsi = new TH1D("h_L1Dmll_o_mll_JPsi", "( M_{ll}(L1,J/Psi) - M(J/#psi) ) / M(J/#psi)", 60, -1.0, 3.0);
	TH1D * h_L1mllCorr_JPsi = new TH1D("h_L1mllCorr_JPsi", "M_{ll}(L1(Corrected),J/Psi)", 100, 0.0, 10.0);
	TH1D * h_L1Dmll_o_mllCorr_JPsi = new TH1D("h_L1Dmll_o_mllCorr_JPsi", "( M_{ll}(L1(Corrected),J/Psi) - M(J/#psi) ) / M(J/#psi)", 60, -1.0, 3.0);
	
	//Only pT corrected plots
	TH1D * h_L1mllCorrOnlyPT = new TH1D("h_L1mllCorrOnlyPT", "M_{ll}(L1(Corrected only p_{T}))", 300, 0.0, 30.0);
	TH1D * h_L1mllCorrOnlyPT_JPsi = new TH1D("h_L1mllCorrOnlyPT_JPsi", "M_{ll}(L1(Corrected only p_{T}),J/Psi)", 100, 0.0, 10.0);
	
	TH1D * h_L1Dmll_o_mllCorrOnlyPT_JPsi = new TH1D("h_L1Dmll_o_mllCorrOnlyPT_JPsi", "( M_{ll}(L1(Corrected only p_{T}),J/Psi) - M(J/#psi) ) / M(J/#psi)", 60, -1.0, 3.0);
	
	//Only phi corrected plots
	TH1D * h_L1mllCorrOnlyPhi = new TH1D("h_L1mllCorrOnlyPhi", "M_{ll}(L1(Corrected only #phi))", 300, 0.0, 30.0);
	TH1D * h_L1mllCorrOnlyPhi_JPsi = new TH1D("h_L1mllCorrOnlyPhi_JPsi", "M_{ll}(L1(Corrected only #phi),J/Psi)", 100, 0.0, 10.0);
	
	TH1D * h_L1Dmll_o_mllCorrOnlyPhi_JPsi = new TH1D("h_L1Dmll_o_mllCorrOnlyPhi_JPsi", "( M_{ll}(L1(Corrected only #phi),J/Psi) - M(J/#psi) ) / M(J/#psi)", 60, -1.0, 3.0);
	
	TFHistos SplitTFHistos;
	if( SplitTFs == 1 ) SplitTFHistos.CreateHistos();
	TString TF1, TF2;
	
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if( ientry % 1000000 == 0 ) cout << "I have processed " << ientry << " events!" << endl;
		
		if( recomuon_N < 2 ) continue;
		
		for(int i = 0; i<recomuon_N; i++) {
			if(recomuon_dr->at(i) < 0.0 || recomuon_dr->at(i) > 0.2) continue;
			
			L1muon_ptCorr_ = (float) L1muon_ptCorr->at(i); L1muon_eta_ = (float) L1muon_eta->at(i); L1muon_phi_ = (float) L1muon_phiAtVtx->at(i); L1muon_charge_ = (float) L1muon_charge->at(i);
			double DphiL1Reco_i, PhiRecoReg_i;
			correctThePhi(L1muon_eta_, DphiL1Reco_i, PhiRecoReg_i);
			
			if( SplitTFs == 1 ) TF1 = MatchWithTF(L1muon_eta_);
			
			for(int j = i+1; j<recomuon_N; j++) {
				L1muon_ptCorr_ = (float) L1muon_ptCorr->at(j); L1muon_eta_ = (float) L1muon_eta->at(j); L1muon_phi_ = (float) L1muon_phiAtVtx->at(j); L1muon_charge_ = (float) L1muon_charge->at(j);
				double DphiL1Reco_j, PhiRecoReg_j;
				correctThePhi(L1muon_eta_, DphiL1Reco_j, PhiRecoReg_j);
				
				float recomll = Mll( recomuon_pt->at(i), recomuon_eta->at(i), recomuon_phi->at(i), recomuon_pt->at(j), recomuon_eta->at(j), recomuon_phi->at(j) );
				h_recomll->Fill(recomll);
				float L1mll = Mll(L1muon_ptCorr->at(i), L1muon_eta->at(i), L1muon_phi->at(i), L1muon_ptCorr->at(j), L1muon_eta->at(j), L1muon_phi->at(j));
				h_L1mll->Fill(L1mll);
				float L1mllCorrOnlyPT = Mll(L1muon_ptCorr->at(i), L1muon_eta->at(i), L1muon_phiAtVtx->at(i), L1muon_ptCorr->at(j), L1muon_eta->at(j), L1muon_phiAtVtx->at(i));
				h_L1mllCorrOnlyPT->Fill(L1mllCorrOnlyPT);
				float L1mllCorrOnlyPhi = Mll(L1muon_pt->at(i), L1muon_eta->at(i), PhiRecoReg_i, L1muon_pt->at(j), L1muon_eta->at(j), PhiRecoReg_j);
				h_L1mllCorrOnlyPhi->Fill(L1mllCorrOnlyPhi);
				float L1mllCorr = Mll(L1muon_ptCorr->at(i), L1muon_eta->at(i), PhiRecoReg_i, L1muon_ptCorr->at(j), L1muon_eta->at(j), PhiRecoReg_j);
				h_L1mllCorr->Fill(L1mllCorr);
				if( recomll > 3.046 && recomll < 3.146 ) {
					h_L1mll_JPsi->Fill(L1mll);
					h_L1Dmll_o_mll_JPsi->Fill( (L1mll - JPsi_M) / JPsi_M );
					h_L1mllCorrOnlyPT_JPsi->Fill(L1mllCorrOnlyPT);
					h_L1Dmll_o_mllCorrOnlyPT_JPsi->Fill( (L1mllCorrOnlyPT - JPsi_M) / JPsi_M );
					h_L1mllCorrOnlyPhi_JPsi->Fill(L1mllCorrOnlyPhi);
					h_L1Dmll_o_mllCorrOnlyPhi_JPsi->Fill( (L1mllCorrOnlyPhi - JPsi_M) / JPsi_M );
					h_L1mllCorr_JPsi->Fill(L1mllCorr);
					h_L1Dmll_o_mllCorr_JPsi->Fill( (L1mllCorr - JPsi_M) / JPsi_M );
					if( SplitTFs == 1 ) {
						TF2 = MatchWithTF(L1muon_eta_);
						SplitTFHistos.FillHistos(TF1, TF2, L1mllCorr);
					}
				}
			}
		}
	}
	
	return;
}
