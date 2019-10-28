#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

#include "deltaPhi.h"


void TFCut(TCut &mycut, TString TF, TString etaOrIndex = "Eta") {
	if( etaOrIndex == "Index" ) {
		if(TF == "B") mycut += " L1muon_tfMuonIndex >= 36 && L1muon_tfMuonIndex <= 70 ";
		else if(TF == "O") mycut += " (L1muon_tfMuonIndex >= 17 && L1muon_tfMuonIndex <= 35) || (L1muon_tfMuonIndex >= 71 && L1muon_tfMuonIndex <= 89) ";
		else if(TF == "E") mycut += " (L1muon_tfMuonIndex >= 0 && L1muon_tfMuonIndex <= 16) || (L1muon_tfMuonIndex >= 90 && L1muon_tfMuonIndex <= 107) ";
		else { cout << "Invalid TF, exiting." << endl; exit(1); }
	}
	else if( etaOrIndex == "Eta" ) {
		if(TF == "B") mycut += " fabs(L1muon_eta) < 0.8 ";
		else if(TF == "O") mycut += " fabs(L1muon_eta) >= 0.8 && fabs(L1muon_eta) < 1.2 ";
		else if(TF == "E") mycut += " fabs(L1muon_eta) >= 1.2 ";
		else { cout << "Invalid TF, exiting." << endl; exit(1); }
	}
	else { cout << "Invalid TF binning, neither \"Eta\" nor \"Index\". Exiting..." << endl; exit(1); }
	
	return;
}

void GuysCut(TCut &mycut, TString guys) {
	if(guys == "G") mycut += " L1muon_pt >= 10 && L1muon_pt <= 40 ";
	else if(guys == "B") mycut += " L1muon_pt <= 10 || L1muon_pt >= 40 ";
	else if(guys == "L") mycut += " L1muon_pt <= 10 ";
	else if(guys == "H") mycut += " L1muon_pt >= 40 ";
	else cout << "Neither G(ood) nor B(ad) guys selected ==> Inclusive training." << endl;
	
	return;
}


TFile * InitializeFile(TString dataset, TString year, TString ID, TString fileEras) {
	TFile *input(0);
	TString NTupleDir = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/";
	TString fname = NTupleDir+"L1toRecoMatchPlots_"+dataset+year+"_"+ID+"_"+fileEras+".root";
	if (!gSystem->AccessPathName( fname )) {
		input = TFile::Open( fname );
	}
	if (!input) cout << "ERROR: could not open data file " << fileEras << endl;
	return input;
}


using namespace TMVA;

void TMVARegression( TString dataset, TString year, TString ID, TString TF, TString fileEras, TString guys, TString extraText, TString etaOrIndex)
{
	Tools::Instance();
	
	
	cout << endl;
	cout << "==> Start TMVARegression" << endl;
	
	
	TString CMSSW_BASE = getenv("CMSSW_BASE");
	TString dataloaderName( "TMVARegression_TF"+TF+"_Era"+fileEras+"_Guys"+guys+"_"+etaOrIndex+extraText );
	TString outfileName( CMSSW_BASE+"/src/MLShEx_CorrL1Muons/Configuration/test/TMVATraining/"+dataloaderName );
	TFile* outputFile = TFile::Open( outfileName+".root", "RECREATE" );
	
	TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
	
	TMVA::DataLoader *dataloader = new TMVA::DataLoader(dataloaderName);
	
	dataloader->AddVariable( "L1muon_ptCorr", "p_{T,corrected}(L1 #mu)", "GeV", 'F' );
	dataloader->AddVariable( "L1muon_eta", "#eta(L1 #mu)", "", 'F' );
	dataloader->AddVariable( "L1muon_phiAtVtx", "#phiAtVtx(L1 #mu)", "", 'F' );
	dataloader->AddVariable( "L1muon_charge", "charge(L1 #mu)", "", 'I' );
	
	dataloader->AddTarget( "deltaPhi(L1muon_phiAtVtx,recomuon_phi)" );
	
	TFile * inputFile = InitializeFile(dataset,year,ID,fileEras);
	TTree * inputTree;
	if( inputFile != NULL ) {
		inputTree =  (TTree*)inputFile->Get("mytree");
		dataloader->AddRegressionTree( inputTree );
	}
	else return;
	
	
	TCut mycut = " recomuon_dr >= 0.0 && recomuon_dr < 0.2 ";
	TFCut(mycut, TF, etaOrIndex);
	GuysCut(mycut, guys);
	cout << mycut << endl;
	
	
	TString preparationString = "nTrain_Regression=25000:nTest_Regression=25000:SplitMode=Random:NormMode=NumEvents:!V";
	dataloader->PrepareTrainingAndTestTree( mycut, preparationString );
	
	
	factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=300:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.6:SamplingEpoch=0.8:ConvergenceImprove=1e-4:ConvergenceTests=15" );
	
	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();
	
	
	outputFile->Close();
	
	cout << "==> Wrote root file: " << outputFile->GetName() << endl;
	cout << "==> TMVARegression is done!" << endl;
	
	delete factory;
	delete dataloader;
	
	return;
}
