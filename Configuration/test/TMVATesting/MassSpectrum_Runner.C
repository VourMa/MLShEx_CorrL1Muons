#include "MassSpectrum.h"
#include "MassSpectrum.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MassSpectrum_Runner()
{
	//___ Input - Output ___//
	TChain *Input = new TChain("mytree");
	//Input->Add("../FriendTree_ptCorr_test.root");
	Input->Add("../L1toRecoMatchPlots_tight_B.root");
	Input->Add("../L1toRecoMatchPlots_tight_C.root");
	Input->Add("../L1toRecoMatchPlots_tight_E.root");
	Input->Add("../L1toRecoMatchPlots_tight_F.root");
	
//	Input->Add("../L1toRecoMatchPlots_Charm_tight_B.root");
//	Input->Add("../L1toRecoMatchPlots_Charm_tight_C.root");
	
	TFile * out = new TFile("MassSpectrum_BCEF.root","recreate");
	
//	TFile * out = new TFile("MassSpectrum_Charm_BC.root","recreate");


	//Definitions	
	TH1D * h_recomll = new TH1D("h_recomll", "M_{ll}(reco)", 300, 0.0, 30.0);
	TH1D * h_L1mll = new TH1D("h_L1mll", "M_{ll}(L1)", 300, 0.0, 30.0);
	TH1D * h_L1mllCorr = new TH1D("h_L1mllCorr", "M_{ll}(L1(Corrected))", 300, 0.0, 30.0);
	TH1D * h_L1mll_JPsi = new TH1D("h_L1mll_JPsi", "M_{ll}(L1,J/Psi)", 100, 0.0, 10.0);
	TH1D * h_L1mllCorr_JPsi = new TH1D("h_L1mllCorr_JPsi", "M_{ll}(L1(Corrected),J/Psi)", 100, 0.0, 10.0);


	//Running
	MassSpectrum Read(Input);
	Read.Loop(h_recomll, h_L1mll, h_L1mllCorr, h_L1mll_JPsi, h_L1mllCorr_JPsi);
	
	out->Write();
	out->Close();
	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	
	return;
}