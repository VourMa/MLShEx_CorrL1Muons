#include "MassSpectrum.h"
#include "MassSpectrum.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MassSpectrum_Runner(TString dataset, TString era, string particle, string whichGuys = "A", bool debug = 0)
{
	TString etaOrIndex = "Eta"; // or "Index"
	
	//___ Input - Output ___//
	TChain *Input = new TChain("mytree");
	TString NTupleDir = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/";
	
	if( dataset == "ZeroBias2017" ) {
		if( era.Contains("B") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_ZeroBias2017_tight_B.root");
		if( era.Contains("C") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_ZeroBias2017_tight_C.root");
		if( era.Contains("E") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_ZeroBias2017_tight_E.root");
		if( era.Contains("F") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_ZeroBias2017_tight_F.root");
	}
	else if( dataset == "MuOnia2017" ) {
		if( era.Contains("A") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_MuOnia2017_tight_A.root");
		if( era.Contains("B") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_MuOnia2017_tight_B.root");
		if( era.Contains("C") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_MuOnia2017_tight_C.root");
	}
	else if( dataset == "Charmonium2017" ) {
		if( era.Contains("B") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_Charmonium2017_tight_B.root");
		if( era.Contains("C") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_Charmonium2017_tight_C.root");
		if( era.Contains("D") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_Charmonium2017_tight_D.root");
	}
	else if( dataset == "ZeroBias2018" ) {
		Input->Add(NTupleDir + "L1toRecoMatchPlots_ZeroBias2018_tight_ABC.root");
		era = "ABC";
	}
	else {
		cout << "Non valid dataset, exiting." << endl;
		return;
	}
	
	if(Input == NULL) {
		cout << "Input is empty, exiting." << endl;
		return;
	}
	
	TString outSuffix = dataset + "_" + era + "_" + particle + "_" + etaOrIndex;
	if( debug == 1 ) outSuffix += "_debugging";
	TFile * out = new TFile("MassSpectrum_"+outSuffix+".root","recreate");


	//Running
	MassSpectrum Read(Input, etaOrIndex);
	Read.Loop(out, particle, whichGuys, etaOrIndex, debug);
	
	out->Write();
	out->Close();
	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	
	return;
}
