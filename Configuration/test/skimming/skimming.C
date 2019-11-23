//#define skimming_cxx
#include <iomanip>
#include <iostream>
#include <fstream>
#include "skimming.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TFile.h"
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



void skimming::Loop(TString ID,TString out, bool debug, bool onlytxt)
{
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	//cout << "Total number of events: " << nentries << endl;
	if( debug == 1 ) nentries = 100000;
	
	//___Get the profile histograms___
	TFile profFile("/home/cmsdas/public/store/MLShortExercise/L1uGMTAnalyzer_Trees/plots_tight_test_profile.root","read");
	TProfile AllTF, BMTF, OMTF, EMTF;
	getProf(profFile, "pt_rat_vs_L1muon_pt_pfx", AllTF, "pt_rat_vs_L1muon_pt_BMTF_pfx", BMTF, "pt_rat_vs_L1muon_pt_OMTF_pfx", OMTF, "pt_rat_vs_L1muon_pt_EMTF_pfx", EMTF);
	profFile.Close();

	//Create txt output file
	std::ofstream outfile;
	outfile.open(out);	
	std::cout<<"Created output file: "<<out<<std::endl;

	//Write name of variables in txt file (only once)
	outfile << std::setw(12) << "reco_pt"   << std::setw(12) << "reco_eta"  << std::setw(12)  << "reco_phi"       << std::setw(15) << "reco_etaAtVtx" <<  std::setw(15) << "reco_phiAtVtx"
		<< std::setw(10) << "l1_pt"     << std::setw(10) << "l1_eta"    << std::setw(12)  << "l1_phi"         << std::setw(15) << "l1_etaAtVtx"   <<  std::setw(15) << "l1_phiAtVtx"
		<< std::setw(12) << "l1_ptCorr" << std::setw(12) << "l1_charge" << std::setw(15)  << "l1_tfMuonIndex" << std::setw(10) << "l1_hwQual" 
		<< std::setw(10) << "dr_map"    << std::setw(10) << "dphi" << std::endl;						
	
	//___Definitions__//	
	int event_, recomuon_N, L1muon_N;
	MatchedMuons recomuon, L1muon;
	
	TTree * mytree =new TTree("mytree","mytree");	
	mytree->Branch("event",             &event_);
	mytree->Branch("recomuon_N",        &recomuon_N);
	mytree->Branch("recomuon_pt",       &recomuon.Pt);
	mytree->Branch("recomuon_eta",      &recomuon.Eta);
	mytree->Branch("recomuon_etaAtVtx", &recomuon.EtaAtVtx);
	mytree->Branch("recomuon_phi",      &recomuon.Phi);
	mytree->Branch("recomuon_phiAtVtx", &recomuon.PhiAtVtx);
	mytree->Branch("recomuon_dr",       &recomuon.Dr);
	mytree->Branch("L1muon_N",          &L1muon_N);
	mytree->Branch("L1muon_pt",         &L1muon.Pt);
	mytree->Branch("L1muon_ptCorr",     &L1muon.PtCorr);
	mytree->Branch("L1muon_eta",        &L1muon.Eta);
	mytree->Branch("L1muon_etaAtVtx",   &L1muon.EtaAtVtx);
	mytree->Branch("L1muon_phi",        &L1muon.Phi);
	mytree->Branch("L1muon_phiAtVtx",   &L1muon.PhiAtVtx);
	mytree->Branch("L1muon_charge",     &L1muon.Charge);
	mytree->Branch("L1muon_tfMuonIndex",&L1muon.TfMuonIndex);
	mytree->Branch("L1muon_hwQual",     &L1muon.HwQual);

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);		
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		bool bprint = (debug) ? (ientry % (nentries/10) == 0) : (ientry % 1000000 == 0);
		if(bprint) cout << "I have processed " << ientry << " events!" << endl;

		int nReco = 0;
		std::map<float, pair<int, int> > dr_reco_L1;
		recomuon.Clear(); L1muon.Clear();

		for(unsigned int i = 0; i < recomuon_pt->size(); i++) {
			if(ID == "tight") if(recomuon_isTightMuon->at(i) == 0) continue;
			if(ID == "medium") if(recomuon_isMediumMuon->at(i) == 0) continue;
			if(ID == "loose") if(recomuon_isLooseMuon->at(i) == 0) continue;
			if(ID == "global") if(recomuon_isGlobalMuon->at(i) == 0) continue;
			if(ID == "standalone") if(recomuon_isStandAloneMuon->at(i) == 0) continue;
			if(ID == "tracker") if(recomuon_isTrackerMuon->at(i) == 0) continue;
			
			nReco++;
			
			bool NotMatchedAtSt = false;
			for(unsigned int j = 0; j < muon_pt->size(); j++) {
				float dr = DR( recomuon_eta->at(i),recomuon_phi->at(i),muon_eta->at(j),muon_phi->at(j) );
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
		  float dr_map      = it->first;
		  int recoIndex_map = it->second.first;
		  int L1Index_map   = it->second.second;
			
		  //Reco muon
		  float reco_pt       = recomuon_pt->at(recoIndex_map);
		  float reco_eta      = recomuon_eta->at(recoIndex_map);
		  float reco_phi      = recomuon_phi->at(recoIndex_map);
		  float reco_etaAtVtx = recomuon_eta->at(recoIndex_map);
		  float reco_phiAtVtx = recomuon_phi->at(recoIndex_map);
			
		  //L1 muon
		  float l1_pt, l1_eta, l1_ptCorr, l1_phi, l1_etaAtVtx, l1_phiAtVtx, l1_charge, l1_tfMuonIndex, l1_hwQual;
		  l1_pt = l1_eta = l1_ptCorr = l1_phi = l1_etaAtVtx = l1_phiAtVtx = l1_charge = l1_tfMuonIndex = l1_hwQual = -100.0;
		  
		  //DeltaPhi(l1muon, recomuon)
		  float dphi = -100.0;

		  if(L1Index_map >= 0) {
		    l1_pt          = muon_pt->at(L1Index_map);
		    l1_eta         = muon_eta->at(L1Index_map);
		    l1_phi         = muon_phi->at(L1Index_map);
		    l1_etaAtVtx    = muon_etaAtVtx->at(L1Index_map);
		    l1_phiAtVtx    = muon_phiAtVtx->at(L1Index_map);
		    l1_charge      = muon_charge->at(L1Index_map);
		    l1_tfMuonIndex = muon_tfMuonIndex->at(L1Index_map);
		    l1_hwQual      = muon_hwQual->at(L1Index_map);

		    dphi = DPhi(muon_phiAtVtx->at(L1Index_map), recomuon_phi->at(recoIndex_map));
		    
		    double bin, corrFactor = 2.0;
		    if( fabs( l1_eta ) <= 0.8)
		      {
			bin = BMTF.FindBin( l1_pt );
			corrFactor = BMTF.GetBinContent( bin );
			if( corrFactor == 0.0 ) corrFactor = 2.0; //Correct for empty bins
		      }
				
		    else if( fabs( l1_eta ) <= 1.2)
		      {
			bin = OMTF.FindBin( l1_pt );
			corrFactor = OMTF.GetBinContent( bin );
			if( corrFactor == 0.0 ) corrFactor = 2.0; //Correct for empty bins
		      }
				
		    else //( fabs( l1_eta ) > 1.2)
		      {
			bin = EMTF.FindBin( l1_pt );
			corrFactor = EMTF.GetBinContent( bin );
			if( corrFactor == 0.0 ) corrFactor = 3.0; //Correct for empty bins
		      }
			  
		    l1_ptCorr = l1_pt / corrFactor;
		  } //if(L1Index_map >= 0)
			
		  //Store variables in txt file
		  // reco_pt, reco_eta, reco_phi, reco_etaAtVtx, reco_phiAtVtx, l1_pt, l1_eta, l1_phi, l1_etaAtVtx, l1_phiAtVtx, l1_ptCorr, l1_charge, l1_tfMuonIndex, l1_hwQual, dr_map, dphi
		  outfile << std::setw(12) << reco_pt   << std::setw(12) << reco_eta  << std::setw(12)  << reco_phi       << std::setw(12) << reco_etaAtVtx <<  std::setw(15) << reco_phiAtVtx
			  << std::setw(12) << l1_pt     << std::setw(12) << l1_eta    << std::setw(12)  << l1_phi         << std::setw(12) << l1_etaAtVtx   <<  std::setw(15) << l1_phiAtVtx
			  << std::setw(12) << l1_ptCorr << std::setw(12) << l1_charge << std::setw(12)  << l1_tfMuonIndex << std::setw(12) << l1_hwQual 
			  << std::setw(14) << dr_map    << std::setw(12) << dphi << std::endl;						
		
		  if (!onlytxt){
		    //Store variables in vectors
		    recomuon.Pt.push_back( reco_pt );
		    recomuon.Eta.push_back( reco_eta );
		    recomuon.Phi.push_back( reco_phi );
		    recomuon.EtaAtVtx.push_back( reco_etaAtVtx );
		    recomuon.PhiAtVtx.push_back( reco_phiAtVtx );
		    recomuon.Charge.push_back( l1_charge );
		    recomuon.HwQual.push_back( l1_hwQual );
		    recomuon.Dr.push_back( dr_map );
		    L1muon.Pt.push_back( l1_pt );
		    L1muon.Eta.push_back( l1_eta );
		    L1muon.PtCorr.push_back( l1_ptCorr );
		    L1muon.Phi.push_back( l1_phi );
		    L1muon.EtaAtVtx.push_back( l1_etaAtVtx );
		    L1muon.PhiAtVtx.push_back( l1_phiAtVtx );
		    L1muon.Charge.push_back( l1_charge );
		    L1muon.TfMuonIndex.push_back( l1_tfMuonIndex );
		    L1muon.HwQual.push_back( l1_hwQual );
		  }
		}
		
		if (!onlytxt){
		  event_ ++;
		  recomuon_N = recomuon.Size(); L1muon_N = L1muon.Size();
		  mytree->Fill();
		}
	}
	outfile.close();
	
	//Write TTree in root file 
	if (!onlytxt){
	  TString outName = out;
	  outName.ReplaceAll(".txt",".root");
	  TFile *outFile = new TFile(outName,"RECREATE");
	  std::cout<<"Created output file: "<<outName<<std::endl;
	  outFile -> cd();
	  mytree  -> Write();
	  outFile -> Close();
	  std::cout<<"Output file written and closed!\n"<<std::endl;
	}
	
}

