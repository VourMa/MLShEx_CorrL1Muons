#define ScoutedData_cxx
#include "ScoutedData.h"
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

void getProf(TFile & infile, TString histoname_1, TProfile & prof_1, TString histoname_2, TProfile & prof_2, TString histoname_3, TProfile & prof_3, TString histoname_4, TProfile & prof_4) {
	TProfile * h_temp = (TProfile*)infile.Get(histoname_1);
	h_temp->SetName(histoname_1);
	prof_1 = *h_temp;

	h_temp = (TProfile*)infile.Get(histoname_2);
	h_temp->SetName(histoname_2);
	prof_2 = *h_temp;

	h_temp = (TProfile*)infile.Get(histoname_3);
	h_temp->SetName(histoname_3);
	prof_3 = *h_temp;

	h_temp = (TProfile*)infile.Get(histoname_4);
	h_temp->SetName(histoname_4);
	prof_4 = *h_temp;

	return;
}

void correctpT(int i, vector<float> * pt, vector<float> * eta, TProfile BMTF, TProfile OMTF, TProfile EMTF, vector<float> & pt_corr) {
	double bin, corrFactor;
	if( fabs( eta->at(i) ) <= 0.8 ) {
		bin = BMTF.FindBin( pt->at(i) );
		corrFactor = BMTF.GetBinContent( bin );
		if( corrFactor == 0.0 ) corrFactor = 2.0; //Correct for empty bins
		pt_corr.push_back( pt->at(i) / corrFactor );
	}
	else if( fabs( eta->at(i) ) <= 1.2 ) {
		bin = OMTF.FindBin( pt->at(i) );
		corrFactor = OMTF.GetBinContent( bin );
		if( corrFactor == 0.0 ) corrFactor = 2.0; //Correct for empty bins
		pt_corr.push_back( pt->at(i) / corrFactor );
	}
	else {
		bin = EMTF.FindBin( pt->at(i) );
		corrFactor = EMTF.GetBinContent( bin );
		if( corrFactor == 0.0 ) corrFactor = 3.0; //Correct for empty bins
		pt_corr.push_back( pt->at(i) / corrFactor );
	}
}


void ScoutedData::correctThePhi(float L1muon_eta_, double & DphiL1Reco, double & PhiRecoReg) {
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

void ScoutedData::Loop(TFile * out, bool debug)
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	cout << "Number of entries:" << fChain->GetEntries() << endl;
	if(debug) nentries = 10000000;
	
	//___Get the profile histograms___
	TFile profFile("/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer/CMSSW_10_1_9_patch1/src/L1uGMTAnalyzer/Configuration/test/plots_tight_test_profile.root","read");
	TProfile AllTF, BMTF, OMTF, EMTF;
	getProf(profFile, "pt_rat_vs_L1muon_pt_pfx", AllTF, "pt_rat_vs_L1muon_pt_BMTF_pfx", BMTF, "pt_rat_vs_L1muon_pt_OMTF_pfx", OMTF, "pt_rat_vs_L1muon_pt_EMTF_pfx", EMTF);
	profFile.Close();
	
	out->cd();
	
	//Definitions
	TH1D * h_L1mll = new TH1D("h_L1mll", "M_{ll}(L1)", 1500, 0.0, 150.0);
	TH1D * h_L1mllOS = new TH1D("h_L1mllOS", "M_{ll}(L1,OS)", 1500, 0.0, 150.0);
	TH1D * h_L1mllSS = new TH1D("h_L1mllSS", "M_{ll}(L1,SS)", 1500, 0.0, 150.0);
	TH1D * h_L1mll_B = new TH1D("h_L1mll_B", "M_{ll}(L1,BMTF)", 1500, 0.0, 150.0);
	TH1D * h_L1mllOS_B = new TH1D("h_L1mllOS_B", "M_{ll}(L1,OS,BMTF)", 1500, 0.0, 150.0);
	TH1D * h_L1mllSS_B = new TH1D("h_L1mllSS_B", "M_{ll}(L1,SS,BMTF)", 1500, 0.0, 150.0);
	
	TH1D * h_L1mllCorr = new TH1D("h_L1mllCorr", "M_{ll}(L1(Corrected))", 1500, 0.0, 150.0);
	TH1D * h_L1mllCorrOS = new TH1D("h_L1mllCorrOS", "M_{ll}(L1(Corrected),OS)", 1500, 0.0, 150.0);
	TH1D * h_L1mllCorrSS = new TH1D("h_L1mllCorrSS", "M_{ll}(L1(Corrected),SS)", 1500, 0.0, 150.0);
	TH1D * h_L1mllCorr_B = new TH1D("h_L1mllCorr_B", "M_{ll}(L1(Corrected),BMTF)", 1500, 0.0, 150.0);
	TH1D * h_L1mllCorrOS_B = new TH1D("h_L1mllCorrOS_B", "M_{ll}(L1(Corrected),OS,,BMTF)", 1500, 0.0, 150.0);
	TH1D * h_L1mllCorrSS_B = new TH1D("h_L1mllCorrSS_B", "M_{ll}(L1(Corrected),SS,BMTF)", 1500, 0.0, 150.0);
	
	//Only pT corrected plots
	TH1D * h_L1mllCorrOnlyPT = new TH1D("h_L1mllCorrOnlyPT", "M_{ll}(L1(Corrected only p_{T}))", 300, 0.0, 150.0);
	
	//Only phi corrected plots
	TH1D * h_L1mllCorrOnlyPhi = new TH1D("h_L1mllCorrOnlyPhi", "M_{ll}(L1(Corrected only #phi))", 300, 0.0, 150.0);
	
	
	vector<float> pt_corr;
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if( ientry % 1000000 == 0 ) cout << "I have processed " << ientry << " events!" << endl;
		
		if(pt->size() < 2) continue;
		
		pt_corr.clear();
		
		for(int i = 0; i<pt->size(); i++) {
			if(i==0) correctpT(i, pt, eta, BMTF, OMTF, EMTF, pt_corr);
			
			L1muon_ptCorr_ = (float) pt_corr.at(i); L1muon_pt_ = (float) pt->at(i); L1muon_eta_ = (float) eta->at(i); L1muon_phi_ = (float) phip->at(i); L1muon_charge_ = (float) charge->at(i);
			double DphiL1Reco_i, PhiRecoReg_i;
			correctThePhi(L1muon_eta_, DphiL1Reco_i, PhiRecoReg_i);
			
			for(int j = i+1; j<pt->size(); j++) {
				if(i==0) correctpT(j, pt, eta, BMTF, OMTF, EMTF, pt_corr);
//				if (pt->at(i) < 5 || pt->at(j) < 5) continue;
				
				L1muon_ptCorr_ = (float) pt_corr.at(j); L1muon_pt_ = (float) pt->at(j); L1muon_eta_ = (float) eta->at(j); L1muon_phi_ = (float) phip->at(j); L1muon_charge_ = (float) charge->at(j);
				double DphiL1Reco_j, PhiRecoReg_j;
				correctThePhi(L1muon_eta_, DphiL1Reco_j, PhiRecoReg_j);
				
				float L1mll = Mll(pt->at(i), etap->at(i), phip->at(i), pt->at(j), etap->at(j), phip->at(j));
				h_L1mll->Fill(L1mll);
				
				float L1mllCorrOnlyPT = Mll(pt_corr.at(i), etap->at(i), phip->at(i), pt_corr.at(j), etap->at(j), phip->at(j));
				h_L1mllCorrOnlyPT->Fill(L1mllCorrOnlyPT);
				
				float L1mllCorrOnlyPhi = Mll(pt->at(i), etap->at(i), PhiRecoReg_i, pt->at(j), etap->at(j), PhiRecoReg_j);
				h_L1mllCorrOnlyPhi->Fill(L1mllCorrOnlyPhi);
				
				float L1mllCorr = Mll(pt_corr.at(i), etap->at(i), PhiRecoReg_i, pt_corr.at(j), etap->at(j), PhiRecoReg_j);
				h_L1mllCorr->Fill(L1mllCorr);
				
				if( charge->at(i) * charge->at(j) < 0) {
					h_L1mllOS->Fill(L1mll);
					h_L1mllCorrOS->Fill(L1mllCorr);
				}
				else {
					h_L1mllSS->Fill(L1mll);
					h_L1mllCorrSS->Fill(L1mllCorr);
				}
				
				if( fabs(eta->at(i)) < 0.8 && fabs(eta->at(j)) < 0.8 ) {
					h_L1mll_B->Fill(L1mll);
					h_L1mllCorr_B->Fill(L1mllCorr);
					
					if( charge->at(i) * charge->at(j) < 0) {
						h_L1mllOS_B->Fill(L1mll);
						h_L1mllCorrOS_B->Fill(L1mllCorr);
					}
					else {
						h_L1mllSS_B->Fill(L1mll);
						h_L1mllCorrSS_B->Fill(L1mllCorr);
					}
				}
				//FIXME Other TF can be added and a more efficient way to fill histos should be found
			}
		}
	}
}
