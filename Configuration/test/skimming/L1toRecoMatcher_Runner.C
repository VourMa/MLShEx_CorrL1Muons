#include "L1toRecoMatcher.h"
#include "L1toRecoMatcher.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TChain.h"
#include "TFile.h"
#include <TGraphAsymmErrors.h>



void L1toRecoMatcher_Runner(TString dataset, TString era, TString ID="tight", bool debug = 0)
{
	cout << "Dataset = " << dataset << endl;
	cout << "ID = " << ID << endl;
	cout << "Eras = " << era << endl;


	//___ Inputs ___//
	TChain * InputChain = new TChain("events");
	
	if( dataset == "ZeroBias2017" ) {
		if( era.Contains("B") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017B_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017B_2018_12_17/181217_134608/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("C") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017C_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017C_2018_12_17/181218_000708/0000/outputL1uGMTAnalyzer.root");
//		if( era.Contains("D") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017D_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017D_2018_12_17/181225_184134/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("E") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017E_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017E_2018_12_17/181225_184329/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("F") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017F_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017F_2018_12_17/181225_184458/0000/outputL1uGMTAnalyzer.root");
	}
	else if( dataset == "Charmonium2017" ) {
		if( era.Contains("B") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Charmonium_Run2017B_2018_12_17/Charmonium/crab_L1uGMTAnalyzer_Charmonium_Run2017B_2018_12_17/181217_163536/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("C") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Charmonium_Run2017C_2018_12_17/Charmonium/crab_L1uGMTAnalyzer_Charmonium_Run2017C_2018_12_17/181225_201704/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("D") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Charmonium_Run2017D_2018_12_17/Charmonium/crab_L1uGMTAnalyzer_Charmonium_Run2017D_2018_12_17/181225_201812/0000/outputL1uGMTAnalyzer.root");
	}
	else if( dataset == "MuOnia2017" ) {
		if( era.Contains("A") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_MuOnia_Run2018A_2019_02_26/MuOnia/crab_L1uGMTAnalyzer_MuOnia_Run2018A_2019_02_26/190226_141001/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("B") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_MuOnia_Run2018B_2019_02_26/MuOnia/crab_L1uGMTAnalyzer_MuOnia_Run2018B_2019_02_26/190227_081713/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("C") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_MuOnia_Run2018C_2019_02_26/MuOnia/crab_L1uGMTAnalyzer_MuOnia_Run2018C_2019_02_26/190227_081856/0000/outputL1uGMTAnalyzer.root");
	}
	else if( dataset == "ZeroBias2018" ) {
		if( era.Contains("A") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2018A_2019_02_06/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2018A_2019_02_06/190206_231114/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("B") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2018B_2019_02_06/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2018B_2019_02_06/190206_231152/0000/outputL1uGMTAnalyzer.root");
		if( era.Contains("C") ) InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2018C_2019_02_06/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2018C_2019_02_06/190206_175901/0000/outputL1uGMTAnalyzer.root");
	}
	else {
		cout << "Non valid dataset, exiting." << endl;
		return;
	}
	
	if(InputChain == NULL) {
		cout << "InputChain is empty, exiting." << endl;
		return;
	}
	
	TFile * out  = new TFile("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/L1toRecoMatchPlots_"+dataset+"_"+ID+"_"+era+".root","recreate");
	
	//___Running___//
	L1toRecoMatcher L;
	L.Init(InputChain);
	L.Loop(ID,out,debug);


	//Closing procedures
	out->Write();
	out->Close();


	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	return;
}
