#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

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


TString NTupleDir = "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer/CMSSW_10_1_9_patch1/src/L1uGMTAnalyzer/Configuration/test/";


float deltaPhi(float phi1, float phi2) {
	float result=phi1-phi2;
	if (result>3.14) result-=6.28;
	if (result<=-3.14) result+=6.28;
	return result;
}

void TFCut(TCut &mycut, TString TF) {
	if(TF == "B") mycut += " (L1muon_tfMuonIndex >= 36 && L1muon_tfMuonIndex <= 70) ";
	else if(TF == "O") mycut += " ( (L1muon_tfMuonIndex >= 17 && L1muon_tfMuonIndex <= 35) || (L1muon_tfMuonIndex >= 71 && L1muon_tfMuonIndex <= 89) ) ";
	else if(TF == "E") mycut += " ( (L1muon_tfMuonIndex >= 0 && L1muon_tfMuonIndex <= 16) || (L1muon_tfMuonIndex >= 90 && L1muon_tfMuonIndex <= 107) ) ";
	else { cout << "Invalid TF, exiting." << endl; exit(1); }
	
	return;
}

void GuysCut(TCut &mycut, TString guys) {
	if(guys == "G") mycut += " ( L1muon_pt >= 10 && L1muon_pt <= 40 ) ";
	else if(guys == "B") mycut += " ( L1muon_pt <= 10 || L1muon_pt >= 40 ) ";
	else cout << "Neither G(ood) nor B(ad) guys selected ==> Inclusive training." << endl;
	
	return;
}


class LoadFiles {
	private:
	static const int NEras = 6;
	TString eras[NEras] = {"A","B","C","D","E","F"};
	vector<TFile *> InitializedFiles;
	vector<TTree *> regTree; 
	
	TFile * InitializeFile(TString era);
	void DataLoadFiles(vector<TFile *> Files, TMVA::DataLoader * dataloader);
	
	public:
	void FilesInReg(TString fileEras, TMVA::DataLoader * dataloader);
};

void LoadFiles::FilesInReg(TString fileEras, TMVA::DataLoader * dataloader) {
	for(int i = 0; i < NEras; i++) {
		if( fileEras.Contains( eras[i] ) ) {
			if( InitializeFile( eras[i] ) != NULL ) InitializedFiles.push_back( InitializeFile( eras[i] ) );
		}
	}
	DataLoadFiles(InitializedFiles, dataloader);
}

TFile * LoadFiles::InitializeFile(TString era) {
	TFile *input(0);
	TString fname = NTupleDir+"L1toRecoMatchPlots_tight_"+era+".root";
	if (!gSystem->AccessPathName( fname )) {
		input = TFile::Open( fname );
		//cout << "--- TMVARegression           : Using input file: " << input->GetName() << endl;
	}
	if (!input) {
		cout << "ERROR: could not open data file " << era << endl;
		input = NULL;
	}
	
	return input;
}

void LoadFiles::DataLoadFiles(vector<TFile *> Files, TMVA::DataLoader * dataloader) {
	regTree.clear();
	for(int i = 0; i < Files.size(); i++) {
			regTree.push_back( (TTree*)Files.at(i)->Get("mytree") );
			dataloader->AddRegressionTree( regTree.at(i) );
	}
	return;
}


using namespace TMVA;

void TMVARegression( TString TF, TString fileEras, TString guys )
{
	TMVA::Tools::Instance();
	
	
	cout << endl;
	cout << "==> Start TMVARegression" << endl;
	
	
	TString outfileName( "TMVARegression_TF"+TF+"_Era"+fileEras+"_Guys"+guys );
	TFile* outputFile = TFile::Open( outfileName+".root", "RECREATE" );
	
	TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
	
	
	TMVA::DataLoader *dataloader=new TMVA::DataLoader(outfileName);
	
	dataloader->AddVariable( "L1muon_ptCorr", "p_{T,corrected}(L1 #mu)", "GeV", 'F' );
	dataloader->AddSpectator( "L1muon_pt", "p_{T}(L1 #mu)", "GeV", 'F' );
	dataloader->AddVariable( "L1muon_eta", "#eta(L1 #mu)", "", 'F' );
	dataloader->AddVariable( "L1muon_phiAtVtx", "#phiAtVtx(L1 #mu)", "", 'F' );
	dataloader->AddVariable( "L1muon_charge", "charge(L1 #mu)", "", 'I' );
	
	dataloader->AddTarget( "deltaPhi(L1muon_phiAtVtx,recomuon_phi)" );
	
	
	LoadFiles TMVAInputFiles; TMVAInputFiles.FilesInReg(fileEras,dataloader);
	
	
	TCut mycut = " recomuon_dr >= 0.0 && recomuon_dr < 0.2 ";
	TFCut(mycut, TF);
	GuysCut(mycut, guys);
	cout << mycut << endl;
	
	
	TString preparationString;
	if(TF == "E") preparationString = "nTrain_Regression=250000:nTest_Regression=250000:SplitMode=Random:NormMode=NumEvents:!V";
	else preparationString = "SplitMode=Random:NormMode=NumEvents:!V";
	dataloader->PrepareTrainingAndTestTree( mycut, preparationString );
	
	
	factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15" );
	//factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "");
	
	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();
	
	
	outputFile->Close();
	
	cout << "==> Wrote root file: " << outputFile->GetName() << ".root" << endl;
	cout << "==> TMVARegression is done!" << endl;
	
	delete factory;
	delete dataloader;
	
	if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName+".root" );
}
