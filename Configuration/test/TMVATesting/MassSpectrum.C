#define MassSpectrum_cxx
#include "MassSpectrum.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

const float JPsi_M = 3.096;
const float Upsilon_M = 9.460;

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

void findTFandGuys(TString etaOrIndex, float L1muon_pt_, float L1muon_eta_, float L1muon_index_, string &readerTFTemp, string &readerpTtemp) {
	//Label the muon by TF
	if(etaOrIndex == "Eta") {
		if( fabs(L1muon_eta_) < 0.8 ) readerTFTemp = "B";
		else if( fabs(L1muon_eta_) < 1.2 ) readerTFTemp = "O";
		else readerTFTemp = "E";
	}
	else if(etaOrIndex == "Index") {
		if( (L1muon_index_ >= 36 && L1muon_index_ <= 70) ) readerTFTemp = "B";
		else if( (L1muon_index_ >= 17 && L1muon_index_ <= 35) || (L1muon_index_ >= 71 && L1muon_index_ <= 89) ) readerTFTemp = "O";
		else if( (L1muon_index_ >= 0 && L1muon_index_ <= 16) || (L1muon_index_ >= 90 && L1muon_index_ <= 107) ) readerTFTemp = "E";
		else { cout << "Index out of bounds, skipping muon..." << endl; readerTFTemp = "Invalid"; }
	}
	else { cout << "Invalid TF binning, exiting..." << endl; exit(1); }
	
	//Label the muon by pT
	if( ( L1muon_pt_ <= 10 || L1muon_pt_ >= 40 ) ) readerpTtemp = "B";
	else readerpTtemp = "G";
}


class TFHistos {
	private:
	TH1D * histoPtrs[6] = {};
	TH1D * histoDPtrs[6] = {};
	
	void FillHisto(string TF1, string TF2, unordered_map<string, TH1D *> histosMap, float mass, bool difference, string particle);
	
	public:
	string TFCombinations[6] = {"BB", "BO", "BE", "OO", "OE", "EE"};
	unordered_map<string, TH1D *> histoMap;
	unordered_map<string, TH1D *> histoDMap;
	
	void CreateHistos(string particle);
	void FillHistos(string TF1, string TF2, float mass, string particle);
};

void TFHistos::CreateHistos (string particle) {
	for( int i = 0; i < 6; i++) {
		TString histoName = "h_L1mllCorr_"+particle+"_"+TFCombinations[i];
		TString histoTitle = "M_{ll}(L1(Corrected),"+particle+") "+TFCombinations[i];
		TString histoDName = "h_L1DmllCorr_o_mllCorr_"+particle+"_"+TFCombinations[i];
		TString histoDTitle = "( M_{ll}(L1(Corrected),"+particle+") - M("+particle+") ) / M("+particle+") "+TFCombinations[i];
		
		histoPtrs[i] = new TH1D(histoName, histoTitle, 150, 0.0, 15.0);
		histoDPtrs[i] = new TH1D(histoDName, histoDTitle, 80, -1.0, 3.0);
		histoMap.insert( pair<string,TH1D *> (TFCombinations[i],histoPtrs[i]) );
		histoDMap.insert( pair<string,TH1D *> (TFCombinations[i],histoDPtrs[i]) );
	}
}

void TFHistos::FillHisto(string TF1, string TF2, unordered_map<string, TH1D *> histosMap, float mass, bool difference, string particle) {
	float M;
	if(particle == "Upsilon") M = Upsilon_M;
	else M = JPsi_M;
	
	string tempKey = TF1; tempKey.append(TF2);
	unordered_map<string, TH1D *>::iterator histosit = histosMap.find( tempKey );
	if( histosit == histosMap.end() ) {
		tempKey = TF2; tempKey.append(TF1);
		histosit = histosMap.find( tempKey );
		
	}
	if(difference == 0) (histosit->second)->Fill(mass);
	else (histosit->second)->Fill( (mass - M) / M );
}

void TFHistos::FillHistos(string TF1, string TF2, float mass, string particle) {
	FillHisto( TF1, TF2, histoMap, mass, 0, particle);
	FillHisto( TF1, TF2, histoDMap, mass, 1, particle);
}


void MassSpectrum::correctThePhi(string readerTFTemp, string whichGuys, double & DphiL1Reco, double & PhiRecoReg) {
	string stringTemp = readerTFTemp; stringTemp.append("MTF_"); stringTemp.append(whichGuys);
	unordered_map<string,TMVA::Reader *>::iterator it = readerMap.find(stringTemp);
	TMVA::Reader * readerTemp = it->second;
	
	DphiL1Reco = (readerTemp->EvaluateRegression( "MLP" ))[0];
	PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
}

void MassSpectrum::Loop(TFile * out, string particle, TString whichGuys, TString etaOrIndex, bool debug)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
//	cout << "Number of entries:" << nentries << endl;
	if(debug) nentries = 100000;
	
	out->cd();
	
	//FIXME Definition of histograms could be moved to a class
	//Definitions
	TH1D * h_recomll = new TH1D("h_recomll", "M_{ll}(reco)", 300, 0.0, 30.0);
	TH1D * h_L1mll = new TH1D("h_L1mll", "M_{ll}(L1)", 300, 0.0, 30.0);
	TH1D * h_L1mllCorr = new TH1D("h_L1mllCorr", "M_{ll}(L1(Corrected))", 300, 0.0, 30.0);
	
	TH1D * h_L1mll_particle = new TH1D("h_L1mll_particle", "M_{ll}(L1,particle)", 150, 0.0, 15.0);
	TH1D * h_L1Dmll_o_mll_particle = new TH1D("h_L1Dmll_o_mll_particle", "( M_{ll}(L1,particle) - M(J/#psi) ) / M(J/#psi)", 60, -1.0, 3.0);
	TH1D * h_L1mllCorr_particle = new TH1D("h_L1mllCorr_particle", "M_{ll}(L1(Corrected),particle)", 150, 0.0, 15.0);
	TH1D * h_L1Dmll_o_mllCorr_particle = new TH1D("h_L1Dmll_o_mllCorr_particle", "( M_{ll}(L1(Corrected),particle) - M(particle) ) / M(particle)", 60, -1.0, 3.0);
	
	//Only pT corrected plots
	TH1D * h_L1mllCorrOnlyPT = new TH1D("h_L1mllCorrOnlyPT", "M_{ll}(L1(Corrected only p_{T}))", 300, 0.0, 30.0);
	TH1D * h_L1mllCorrOnlyPT_particle = new TH1D("h_L1mllCorrOnlyPT_particle", "M_{ll}(L1(Corrected only p_{T}),particle)", 150, 0.0, 15.0);
	
	TH1D * h_L1Dmll_o_mllCorrOnlyPT_particle = new TH1D("h_L1Dmll_o_mllCorrOnlyPT_particle", "( M_{ll}(L1(Corrected only p_{T}),particle) - M(particle) ) / M(particle)", 60, -1.0, 3.0);
	
	//Only phi corrected plots
	TH1D * h_L1mllCorrOnlyPhi = new TH1D("h_L1mllCorrOnlyPhi", "M_{ll}(L1(Corrected only #phi))", 300, 0.0, 30.0);
	TH1D * h_L1mllCorrOnlyPhi_particle = new TH1D("h_L1mllCorrOnlyPhi_particle", "M_{ll}(L1(Corrected only #phi),particle)", 150, 0.0, 15.0);
	
	TH1D * h_L1Dmll_o_mllCorrOnlyPhi_particle = new TH1D("h_L1Dmll_o_mllCorrOnlyPhi_particle", "( M_{ll}(L1(Corrected only #phi),particle) - M(particle) ) / M(particle)", 60, -1.0, 3.0);
	
	TFHistos SplitTFHistos; SplitTFHistos.CreateHistos(particle);
	float M;
	if( particle == "Upsilon" ) M = Upsilon_M;
	else M = JPsi_M;
	
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		if( ientry % 1000000 == 0 ) cout << "I have processed " << ientry << " events!" << endl;
		
		if( recomuon_N < 2 ) continue;
		
		for(int i = 0; i<recomuon_N; i++) {
			if(recomuon_dr->at(i) < 0.0 || recomuon_dr->at(i) > 0.2) continue;
			
			L1muon_ptCorr_ = (float) L1muon_ptCorr->at(i); L1muon_pt_ = (float) L1muon_pt->at(i); L1muon_eta_ = (float) L1muon_eta->at(i); L1muon_phi_ = (float) L1muon_phiAtVtx->at(i); L1muon_charge_ = (float) L1muon_charge->at(i); L1muon_index_ = (float) L1muon_tfMuonIndex->at(i);
			
			string TF1, TF2, G1, G2;
			findTFandGuys(etaOrIndex, L1muon_pt_, L1muon_eta_, L1muon_index_, TF1, G1);
			if( TF1 == "Invalid" ) continue;
			
			double DphiL1Reco_i, PhiRecoReg_i;
			if( whichGuys == "Combined" ) {
				if( TF1 == "B" ) correctThePhi(TF1, "A", DphiL1Reco_i, PhiRecoReg_i);
				else correctThePhi(TF1, "G", DphiL1Reco_i, PhiRecoReg_i);
				//Change "B" and "G" according to the best combination
			}
			else correctThePhi(TF1, string(whichGuys), DphiL1Reco_i, PhiRecoReg_i);
			
			
			for(int j = i+1; j<recomuon_N; j++) {
				if(recomuon_dr->at(j) < 0.0 || recomuon_dr->at(j) > 0.2) continue;
				
				
				L1muon_ptCorr_ = (float) L1muon_ptCorr->at(j); L1muon_pt_ = (float) L1muon_pt->at(j); L1muon_eta_ = (float) L1muon_eta->at(j); L1muon_phi_ = (float) L1muon_phiAtVtx->at(j); L1muon_charge_ = (float) L1muon_charge->at(j); L1muon_index_ = (float) L1muon_tfMuonIndex->at(j);
				
				findTFandGuys(etaOrIndex, L1muon_pt_, L1muon_eta_, L1muon_index_, TF2, G2);
				if( TF2 == "Invalid" ) continue;
				
				double DphiL1Reco_j, PhiRecoReg_j;
				if( whichGuys == "Combined" ) {
					if( TF2 == "B" ) correctThePhi(TF2, "A", DphiL1Reco_j, PhiRecoReg_j);
					else correctThePhi(TF2, "G", DphiL1Reco_j, PhiRecoReg_j);
					//Change "B" and "G" according to the best combination
				}
				else correctThePhi(TF2, string(whichGuys), DphiL1Reco_j, PhiRecoReg_j);
				
				float recomll = Mll( recomuon_pt->at(i), recomuon_eta->at(i), recomuon_phi->at(i), recomuon_pt->at(j), recomuon_eta->at(j), recomuon_phi->at(j) );
				h_recomll->Fill(recomll);
				
				float L1mll = Mll(L1muon_pt->at(i), L1muon_etaAtVtx->at(i), L1muon_phiAtVtx->at(i), L1muon_pt->at(j), L1muon_etaAtVtx->at(j), L1muon_phiAtVtx->at(j));
				h_L1mll->Fill(L1mll);
				
				float L1mllCorrOnlyPT = Mll(L1muon_ptCorr->at(i), L1muon_etaAtVtx->at(i), L1muon_phiAtVtx->at(i), L1muon_ptCorr->at(j), L1muon_etaAtVtx->at(j), L1muon_phiAtVtx->at(i));
				h_L1mllCorrOnlyPT->Fill(L1mllCorrOnlyPT);
				
				float L1mllCorrOnlyPhi = Mll(L1muon_pt->at(i), L1muon_etaAtVtx->at(i), PhiRecoReg_i, L1muon_pt->at(j), L1muon_etaAtVtx->at(j), PhiRecoReg_j);
				h_L1mllCorrOnlyPhi->Fill(L1mllCorrOnlyPhi);
				
				float L1mllCorr = Mll(L1muon_ptCorr->at(i), L1muon_etaAtVtx->at(i), PhiRecoReg_i, L1muon_ptCorr->at(j), L1muon_etaAtVtx->at(j), PhiRecoReg_j);
				h_L1mllCorr->Fill(L1mllCorr);
				
				if( recomll > (M-0.050) && recomll < (M+0.050) ) {
					h_L1mll_particle->Fill(L1mll);
					h_L1Dmll_o_mll_particle->Fill( (L1mll - M) / M );
					h_L1mllCorrOnlyPT_particle->Fill(L1mllCorrOnlyPT);
					h_L1Dmll_o_mllCorrOnlyPT_particle->Fill( (L1mllCorrOnlyPT - M) / M );
					h_L1mllCorrOnlyPhi_particle->Fill(L1mllCorrOnlyPhi);
					h_L1Dmll_o_mllCorrOnlyPhi_particle->Fill( (L1mllCorrOnlyPhi - M) / M );
					h_L1mllCorr_particle->Fill(L1mllCorr);
					h_L1Dmll_o_mllCorr_particle->Fill( (L1mllCorr - M) / M );
					
					SplitTFHistos.FillHistos(TF1, TF2, L1mllCorr, particle);
				}
			}
		}
	}
	
	return;
}
