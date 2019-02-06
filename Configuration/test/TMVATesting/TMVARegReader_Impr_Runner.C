#include "TMVARegReader_Impr.h"
#include "TMVARegReader_Impr.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void TMVARegReader_Impr_Runner()
{
	//___ Input - Output ___//
	TChain *Input = new TChain("mytree");
	//Input->Add("../L1toRecoMatchPlots_tight_B.root");
	//Input->Add("../L1toRecoMatchPlots_tight_C.root");
	//Input->Add("../L1toRecoMatchPlots_tight_E.root");
	//Input->Add("../L1toRecoMatchPlots_tight_F.root");
	
	Input->Add("../L1toRecoMatchPlots_Charm_tight_B.root");
	Input->Add("../L1toRecoMatchPlots_Charm_tight_C.root");
	
	//TFile * out = new TFile("TMVATest_Impr_BCEF.root","recreate");
	
	TFile * out = new TFile("TMVATest_Impr_Charm_BC.root","recreate");


	//Definitions	
	TH1D * h_DphiL1Reco = new TH1D("h_DphiL1Reco", "#Delta(#phi_{L1},#phi_{reco})_{Reg}", 300, -1.5, 1.5);
	TH1D * h_PhiRecoReg = new TH1D("h_PhiRecoReg", "#phi_{Reg}", 64, -3.2, 3.2);
	TH1D * h_DphiRegReco = new TH1D("h_DphiRegReco", "#Delta(#phi_{Reg},#phi_{reco})", 80, -0.4, 0.4);
	TH1D * h_DphiExtReco = new TH1D("h_DphiExtReco", "#Delta(#phi_{ext},#phi_{reco})", 80, -0.4, 0.4);


	//Running
	TMVARegReader_Impr Read(Input);
	Read.Loop(h_DphiL1Reco,h_PhiRecoReg,h_DphiRegReco,h_DphiExtReco);
	
	out->Write();
	out->Close();
	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	
	return;
}