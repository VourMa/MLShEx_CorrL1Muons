#include "L1toRecoMatcher.h"
#include "L1toRecoMatcher.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TChain.h"
#include "TFile.h"
#include <TGraphAsymmErrors.h>



void L1toRecoMatcher_Runner(TString ID="tight",bool Charm = false)
{
	cout << "Charm = " << Charm << endl;
	cout << "ID = " << ID << endl;


	//___ Inputs ___//
	TChain * InputChain = new TChain("events");
	TFile * out = nullptr;
	
	if(Charm) {
		//___Inclusive
		InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Charmonium_Run2017B_2018_12_17/Charmonium/crab_L1uGMTAnalyzer_Charmonium_Run2017B_2018_12_17/181217_163536/0000/outputL1uGMTAnalyzer.root");
		//InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Charmonium_Run2017C_2018_12_17/Charmonium/crab_L1uGMTAnalyzer_Charmonium_Run2017C_2018_12_17/181225_201704/0000/outputL1uGMTAnalyzer.root");
		//InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Charmonium_Run2017D_2018_12_17/Charmonium/crab_L1uGMTAnalyzer_Charmonium_Run2017D_2018_12_17/181225_201812/0000/outputL1uGMTAnalyzer.root");
		
		out  = new TFile("L1toRecoMatchPlots_Charm_"+ID+"_B.root","recreate");
	}
	if(!Charm) {
		//___2017
		//InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017B_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017B_2018_12_17/181217_134608/0000/outputL1uGMTAnalyzer.root");
		//InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017C_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017C_2018_12_17/181218_000708/0000/outputL1uGMTAnalyzer.root");
		//InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017D_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017D_2018_12_17/181225_184134/0000/outputL1uGMTAnalyzer.root");
		//InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017E_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017E_2018_12_17/181225_184329/0000/outputL1uGMTAnalyzer.root");
		InputChain->Add("/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_ZeroBias_Run2017F_2018_12_17/ZeroBias/crab_L1uGMTAnalyzer_ZeroBias_Run2017F_2018_12_17/181225_184458/0000/outputL1uGMTAnalyzer.root");
		
		out  = new TFile("L1toRecoMatchPlots_"+ID+"_F.root","recreate");
	}


	//___Running___//
	L1toRecoMatcher L;
	L.Init(InputChain);
	L.Loop(ID,out);


	//Closing procedures
	out->Write();
	out->Close();


	cout << "...Finished! File " << out->GetName() << " created!" << endl;
	return;
}