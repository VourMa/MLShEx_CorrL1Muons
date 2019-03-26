#include "ScoutedData.h"
#include "ScoutedData.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ScoutedData_Runner(TString extraText = "", bool debug = 0)
{
	//___ Input - Output ___//
	TChain *Input = new TChain("t");
	TString NTupleDir = "/eos/cms/store/group/cmst3/group/daql1scout/run2/rootfiles/run_324970/";
	Input->Add(NTupleDir + "scout_324970_*.root");
	
	if(Input == NULL) {
		cout << "Input is empty, exiting." << endl;
		return;
	}
	
	TFile * out = new TFile("ScoutedData"+extraText+".root","recreate");


	//Running
	ScoutedData Read(Input);
	Read.Loop(out, debug);
	
	out->Write();
	out->Close();
	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	
	return;
}
