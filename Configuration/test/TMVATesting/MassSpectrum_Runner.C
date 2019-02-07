#include "MassSpectrum.h"
#include "MassSpectrum.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MassSpectrum_Runner(TString dataset, bool debug = 0)
{
	//___ Input - Output ___//
	TChain *Input = new TChain("mytree");
	TString outSuffix = dataset;
	TString NTupleDir = "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer/CMSSW_10_1_9_patch1/src/L1uGMTAnalyzer/Configuration/test/";
	
	if( dataset == "ZeroBias2017" ) {
		Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_B.root");
		Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_C.root");
		Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_E.root");
		Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_F.root");
		outSuffix += "_BCEF";
	}
	else if( dataset == "Charmonium2017" ) {
		Input->Add(NTupleDir + "L1toRecoMatchPlots_Charm_tight_B.root");
		Input->Add(NTupleDir + "L1toRecoMatchPlots_Charm_tight_C.root");
		outSuffix += "_BC";
	}
	else if( dataset == "ZeroBias2018" ) {
		Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_B.root");
		outSuffix += "_B";
	}
	else {
		cout << "Non valid dataset, exiting." << endl;
		return;
	}
	
	TFile * out = new TFile("MassSpectrum_"+outSuffix+".root","recreate");


	//Running
	MassSpectrum Read(Input);
	Read.Loop(out, debug);
	
	out->Write();
	out->Close();
	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	
	return;
}
