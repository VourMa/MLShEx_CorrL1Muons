#include "Resolutions.h"
#include "Resolutions.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Resolutions_Runner(TString dataset, TString era, TString whichGuys = "A", TString performOnWhichGuys = "A", bool debug = 0)
{
	TString etaOrIndex = "Eta"; // or "Index"
	
	//Input and output
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
	
	//Safguards
	if(Input == NULL) { cout << "Input is empty, exiting." << endl; return; }
	if( whichGuys != "A" && whichGuys != "B" && whichGuys != "G" && whichGuys != "Combined" ) 
	{ cout << "Invalid \"whichGuys\" parameter:\nMLP was trained on \"A\"(ll), \"B\"(ad), \"G\"(ood) or \"Combined\" subsets"<<endl; cout << "Exiting." << endl; return; }
	if( performOnWhichGuys != "A" && performOnWhichGuys != "B" && performOnWhichGuys != "G") 
	{ cout << "Invalid \"performOnWhichGuys\" parameter:\nMLP can be applied on \"A\"(ll), \"B\"(ad) or \"G\"(ood) subsets"<<endl; cout << "Exiting." << endl; return; }
	
	
	TString outSuffix = dataset + "_" + era + "_" + etaOrIndex + "_" + whichGuys + "_" + performOnWhichGuys;
	if( debug == 1) outSuffix += "_debugging";
	TFile * out = new TFile("Resolutions_"+outSuffix+".root","recreate");
	
	//Running
	Resolutions Read(Input,etaOrIndex);
	Read.Loop(out, whichGuys, performOnWhichGuys, etaOrIndex, debug);
	
	out->Write();
	out->Close();
	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	
	return;
}
