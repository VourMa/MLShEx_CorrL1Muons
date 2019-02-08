#include "MassSpectrum.h"
#include "MassSpectrum.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MassSpectrum_Runner(TString dataset, TString era, bool SplitTFs, bool debug = 0)
{
	//___ Input - Output ___//
	TChain *Input = new TChain("mytree");
	TString NTupleDir = "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer/CMSSW_10_1_9_patch1/src/L1uGMTAnalyzer/Configuration/test/";
	
	if( dataset == "ZeroBias2017" ) {
		if( era.Contains("B") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_B.root");
		if( era.Contains("C") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_C.root");
		if( era.Contains("E") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_E.root");
		if( era.Contains("F") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_tight_F.root");
	}
	else if( dataset == "Charmonium2017" ) {
		if( era.Contains("B") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_Charm_tight_B.root");
		if( era.Contains("C") ) Input->Add(NTupleDir + "L1toRecoMatchPlots_Charm_tight_C.root");
		if( era.Contains("D") ) Input->Add("/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/test/skimming/L1toRecoMatchPlots_Charmonium2017_tight_D.root");
	}
	else if( dataset == "ZeroBias2018" ) {
		Input->Add(NTupleDir + "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/test/skimming/L1toRecoMatchPlots_ZeroBias2018_tight_ABC.root");
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
	
	TString outSuffix = dataset + "_" + era;
	if( SplitTFs == 1 ) outSuffix += "_TFSplit";
	TFile * out = new TFile("MassSpectrum_"+outSuffix+".root","recreate");


	//Running
	MassSpectrum Read(Input);
	Read.Loop(out, SplitTFs, debug);
	
	out->Write();
	out->Close();
	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	
	return;
}
