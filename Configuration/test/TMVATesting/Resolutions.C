//#define Resolutions_cxx
#include "Resolutions.h"
#include <unordered_map>
#include <iostream>
#include <string>
#include <string>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"


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


class Histos {
	private:
	TH1D * h_DphiL1Reco, * h_PhiRecoReg, * h_DphiRegReco, * h_DphiExtReco, * h_DphiRegReco_B, * h_DphiExtReco_B, * h_DphiRegReco_O, * h_DphiExtReco_O, * h_DphiRegReco_E, * h_DphiExtReco_E; //Add extra plots here
	
	public:
	void CreateHistos();
	void FillHistos(float DphiL1Reco, float PhiRecoReg, string readerTFTemp, float L1PhiAtVtx, float RecoPhi);
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
	//Define extra plots here
}

void Histos::FillHistos(float DphiL1Reco, float PhiRecoReg, string readerTFTemp, float L1PhiAtVtx, float RecoPhi) {
	h_DphiL1Reco->Fill( DphiL1Reco );
	h_PhiRecoReg->Fill( PhiRecoReg );
	h_DphiRegReco->Fill( DPhi( PhiRecoReg,RecoPhi ) );
	h_DphiExtReco->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	if( readerTFTemp == "B" ) {
		h_DphiRegReco_B->Fill( DPhi( PhiRecoReg,RecoPhi ) );
		h_DphiExtReco_B->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	}
	else if( readerTFTemp == "O" ) {
		h_DphiRegReco_O->Fill( DPhi( PhiRecoReg,RecoPhi ) );
		h_DphiExtReco_O->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	}
	else if( readerTFTemp == "E" ) {
		h_DphiRegReco_E->Fill( DPhi( PhiRecoReg,RecoPhi ) );
		h_DphiExtReco_E->Fill( DPhi( L1PhiAtVtx,RecoPhi ) );
	}
	else { cout << "Index out of bounds, skipping muon in histos!" << endl; } //This should never be printed, the code should have skipped the muon before!
	
	//Fill extra plots here
}


void Resolutions::correctThePhi(string readerTFTemp, string whichGuys, double & DphiL1Reco, double & PhiRecoReg) {
	string stringTemp = readerTFTemp; stringTemp.append("MTF_"); stringTemp.append(whichGuys);
	unordered_map<string,TMVA::Reader *>::iterator it = readerMap.find(stringTemp);
	TMVA::Reader * readerTemp = it->second;
	
	DphiL1Reco = (readerTemp->EvaluateRegression( "MLP" ))[0];
	PhiRecoReg = DPhi( L1muon_phi_,DphiL1Reco );
}

void Resolutions::Loop(TFile * out, TString whichGuys, TString performOnWhichGuys, TString etaOrIndex, Long64_t maxEvents)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	//cout << "Number of entries:" << nentries << endl;
	nentries = maxEvents;
	
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
		
		for(unsigned int i = 0; i<L1muon_ptCorr->size(); i++) {
			//Cuts
			if(recomuon_dr->at(i) < 0.0 || recomuon_dr->at(i) > 0.2) continue;
			
			//The following lines define the subset of muons on which the training will be performed on
			if(performOnWhichGuys == "G") if(L1muon_pt->at(i) <= 10 || L1muon_pt->at(i) >= 40) continue;
			if(performOnWhichGuys == "B") if(L1muon_pt->at(i) > 10 && L1muon_pt->at(i) < 40) continue;
			
			L1muon_ptCorr_ = (float) L1muon_ptCorr->at(i); L1muon_pt_ = (float) L1muon_pt->at(i); L1muon_eta_ = (float) L1muon_eta->at(i); L1muon_phi_ = (float) L1muon_phiAtVtx->at(i); L1muon_charge_ = (float) L1muon_charge->at(i); L1muon_index_ = (float) L1muon_tfMuonIndex->at(i);
			
			//The following lines label the muon
			string readerTFTemp, readerpTtemp;
			findTFandGuys(etaOrIndex,L1muon_pt_,L1muon_eta_,L1muon_index_,readerTFTemp,readerpTtemp);
			if( readerTFTemp == "Invalid" ) continue;
			
			//The following lines apply the weights, based on whichGuys were used for the training
			double DphiL1Reco, PhiRecoReg;
			if( whichGuys == "Combined" ) {
				if( readerpTtemp == "B" ) correctThePhi(readerTFTemp, "A", DphiL1Reco, PhiRecoReg);
				else correctThePhi(readerTFTemp, "G", DphiL1Reco, PhiRecoReg);
				//Change "A", "B" and "G" according to the best combination
			}
			else correctThePhi(readerTFTemp, string(whichGuys), DphiL1Reco, PhiRecoReg);
			
			hist.FillHistos(DphiL1Reco,PhiRecoReg,readerTFTemp,L1muon_phiAtVtx->at(i),recomuon_phi->at(i));
		}
	}
	
	return;
}
