#define L1toRecoMatcher_cxx
#include "L1toRecoMatcher.h"
#include "TLorentzVector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <utility>


float DPhi(float phi1, float phi2) {
	float result = phi1 - phi2;
	if(result>3.14) result -= 6.28;
	if(result<=-3.14) result += 6.28;
	return result;
}

float SPhi(float phi1, float phi2) {
	float result = phi1 + phi2;
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



void L1toRecoMatcher::Loop(TString ID,TFile * out)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	//cout << "Total number of events: " << nentries << endl;
	//nentries = 1000000;


	//___Get the profile histograms___
	TFile profFile("plots_tight_test_profile.root","read");
	TProfile AllTF, BMTF, OMTF, EMTF;
	getProf(profFile, "pt_rat_vs_L1muon_pt_pfx", AllTF, "pt_rat_vs_L1muon_pt_BMTF_pfx", BMTF, "pt_rat_vs_L1muon_pt_OMTF_pfx", OMTF, "pt_rat_vs_L1muon_pt_EMTF_pfx", EMTF);
	profFile.Close();


	//___Definitions__//
	out->cd();
	TTree * mytree =new TTree("mytree","mytree");
	
	int event_ = 0, recomuon_N, L1muon_N;
	MatchedMuons recomuon, L1muon;
	
	mytree->Branch("event",&event_);
	mytree->Branch("recomuon_N",&recomuon_N);
	mytree->Branch("recomuon_pt",&recomuon.Pt);
	mytree->Branch("recomuon_eta",&recomuon.Eta);
	mytree->Branch("recomuon_phi",&recomuon.Phi);
	mytree->Branch("recomuon_dr",&recomuon.Dr);
	//mytree->Branch("recomuon_isTight",&recomuon.IsTight);
	//mytree->Branch("recomuon_isMedium",&recomuon.IsMedium);
	//mytree->Branch("recomuon_isLoose",&recomuon.IsLoose);
	//mytree->Branch("recomuon_isGlobal",&recomuon.IsGlobal);
	//mytree->Branch("recomuon_isStandAlone",&recomuon.IsStandAlone);
	mytree->Branch("L1muon_N",&L1muon_N);
	mytree->Branch("L1muon_pt",&L1muon.Pt);
	mytree->Branch("L1muon_ptCorr",&L1muon.PtCorr);
	mytree->Branch("L1muon_eta",&L1muon.Eta);
	mytree->Branch("L1muon_etaAtVtx",&L1muon.EtaAtVtx);
	mytree->Branch("L1muon_phi",&L1muon.Phi);
	mytree->Branch("L1muon_phiAtVtx",&L1muon.PhiAtVtx);
	mytree->Branch("L1muon_charge",&L1muon.Charge);
	mytree->Branch("L1muon_tfMuonIndex",&L1muon.TfMuonIndex);
	mytree->Branch("L1muon_hwQual",&L1muon.HwQual);



	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;


		int nReco = 0;
		std::map<float, pair<int, int> > dr_reco_L1;
		recomuon.Clear(); L1muon.Clear();


		for(int i = 0; i < recomuon_pt->size(); i++) {
			if(ID == "tight") if(recomuon_isTightMuon->at(i) == 0) continue;
			if(ID == "medium") if(recomuon_isMediumMuon->at(i) == 0) continue;
			if(ID == "loose") if(recomuon_isLooseMuon->at(i) == 0) continue;
			if(ID == "global") if(recomuon_isGlobalMuon->at(i) == 0) continue;
			if(ID == "standalone") if(recomuon_isStandAloneMuon->at(i) == 0) continue;
			if(ID == "tracker") if(recomuon_isTrackerMuon->at(i) == 0) continue;
			
			nReco++;
			
			bool NotMatchedAtSt = false;
			for(int j = 0; j < muon_pt->size(); j++) {
				float dr = -9999;
				if (recomuon_phiAtSt2->at(i) != -9999) dr = DR( recomuon_eta->at(i),recomuon_phiAtSt2->at(i),muon_eta->at(j),muon_phi->at(j) );
				else if (recomuon_phiAtSt1->at(i) != -9999) dr = DR( recomuon_eta->at(i),recomuon_phiAtSt1->at(i),muon_eta->at(j),muon_phi->at(j) );
				else NotMatchedAtSt = true;
				pair<int,int> reco_L1 (i,j);
				dr_reco_L1.insert(pair<float, pair<int, int> > (dr,reco_L1));
			}
			if(muon_pt->size() == 0 || NotMatchedAtSt) {
				pair<int,int> reco_L1 (i,(i+1)*(-1));
				dr_reco_L1.insert(pair<float, pair<int, int> > ((i+1)*(-1),reco_L1));
			}
		}
		
		for (map<float, pair<int, int> >::iterator it=dr_reco_L1.begin(); it!=dr_reco_L1.end(); it++) {
			for (map<float, pair<int, int> >::iterator it1=next(it); it1!=dr_reco_L1.end(); ) {
				
				if(it1->second.first == it->second.first || it1->second.second == it->second.second) {
					it1++;
					
					dr_reco_L1.erase(prev(it1));
				}
				else it1++;
			}
		}
		
		
		for (map<float, pair<int, int> >::iterator it=dr_reco_L1.begin(); it!=dr_reco_L1.end(); it++) {
			float dr_map = it->first;
			int recoIndex_map = it->second.first;
			int L1Index_map = it->second.second;
			
			recomuon.Pt.push_back( recomuon_pt->at(recoIndex_map) );
			recomuon.Eta.push_back( recomuon_eta->at(recoIndex_map) );
			recomuon.Phi.push_back( recomuon_phi->at(recoIndex_map) );
			recomuon.EtaAtVtx.push_back( recomuon_eta->at(recoIndex_map) );
			recomuon.PhiAtVtx.push_back( recomuon_phi->at(recoIndex_map) );
			recomuon.Dr.push_back( dr_map );
			if(L1Index_map >= 0) {
				L1muon.Pt.push_back( muon_pt->at(L1Index_map) );
				L1muon.Eta.push_back( muon_eta->at(L1Index_map) );
				
				double bin, corrFactor;
				if( fabs( L1muon.Eta.back() ) <= 0.8)
				{
					bin = BMTF.FindBin( L1muon.Pt.back() );
					corrFactor = BMTF.GetBinContent( bin );
					if( corrFactor == 0.0 ) corrFactor = 2.0; //Correct for empty bins
					L1muon.PtCorr.push_back( L1muon.Pt.back() / corrFactor );
				}
				
				if( fabs( L1muon.Eta.back() ) > 0.8 && fabs( L1muon.Eta.back() ) <= 1.2)
				{
					bin = OMTF.FindBin( L1muon.Pt.back() );
					corrFactor = OMTF.GetBinContent( bin );
					if( corrFactor == 0.0 ) corrFactor = 2.0; //Correct for empty bins
					L1muon.PtCorr.push_back( L1muon.Pt.back() / corrFactor );
				}
				
				if( fabs( L1muon.Eta.back() ) > 1.2)
				{
					bin = EMTF.FindBin( L1muon.Pt.back() );
					corrFactor = EMTF.GetBinContent( bin );
					if( corrFactor == 0.0 ) corrFactor = 3.0; //Correct for empty bins
					L1muon.PtCorr.push_back( L1muon.Pt.back() / corrFactor );
				}
				//cout << "pt = " << L1muon.Pt.back() << " eta = " << L1muon.Eta.back() << endl;
				//cout << "bin = " << bin << " corrFactor = " << corrFactor << endl;
				//cout << "ptCorr = " << L1muon.PtCorr.back() << endl;
				
				L1muon.Phi.push_back( muon_phi->at(L1Index_map) );
				L1muon.EtaAtVtx.push_back( muon_EtaAtVtx->at(L1Index_map) );
				L1muon.PhiAtVtx.push_back( muon_PhiAtVtx->at(L1Index_map) );
				L1muon.Charge.push_back( muon_charge->at(L1Index_map) );
				L1muon.TfMuonIndex.push_back( muon_tfMuonIndex->at(L1Index_map) );
				L1muon.HwQual.push_back( muon_hwQual->at(L1Index_map) );
			}
			else {
				L1muon.Pt.push_back( -100.0 );
				L1muon.PtCorr.push_back( -100.0 );
				L1muon.Eta.push_back( -100.0 );
				L1muon.Phi.push_back( -100.0 );
				L1muon.EtaAtVtx.push_back( -100.0 );
				L1muon.PhiAtVtx.push_back( -100.0 );
				L1muon.Charge.push_back( -100.0 );
				L1muon.TfMuonIndex.push_back( -100.0 );
				L1muon.HwQual.push_back( -100.0 );
			}
			recomuon.Charge.push_back( L1muon.Charge.back() );
			recomuon.HwQual.push_back( L1muon.HwQual.back() );
		}
		event_++;
		recomuon_N = recomuon.Size(); L1muon_N = L1muon.Size(); 
		mytree->Fill();
	}
}