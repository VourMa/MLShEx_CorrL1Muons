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


float deltaPhi(float phi1, float phi2){
	float result=phi1-phi2;
	if (result>3.14) result-=6.28;
	if (result<=-3.14) result+=6.28;
	return result;
}


using namespace TMVA;

void TMVARegression_Impr_08_12( TString myMethodList = "" )
{
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // Neural Network
   Use["MLP"] = 1;
   Use["DNN_CPU"] = 0;

   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegression" << std::endl;

   // ---------------------------------------------------------------

   // Create a new root output file
   TString outfileName( "TMVAReg_Impr_08_12.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
   "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );


   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset_Impr_08_12");
   
   dataloader->AddVariable( "L1muon_ptCorr", "p_{T}(L1 #mu)", "GeV", 'F' );
   dataloader->AddVariable( "L1muon_eta", "#eta(L1 #mu)", "", 'F' );
   dataloader->AddVariable( "L1muon_phiAtVtx", "#phiAtVtx(L1 #mu)", "", 'F' );
   dataloader->AddVariable( "L1muon_charge", "charge(L1 #mu)", "", 'I' );

   dataloader->AddTarget( "deltaPhi(L1muon_phiAtVtx,recomuon_phi)" );


   TFile *input_B(0);
   TString fname_B = "../L1toRecoMatchPlots_tight_B.root";
   if (!gSystem->AccessPathName( fname_B )) {
      input_B = TFile::Open( fname_B ); // check if file in local directory exists
   }
   if (!input_B) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input_B->GetName() << std::endl;
   
   TFile *input_C(0);
   TString fname_C = "../L1toRecoMatchPlots_tight_C.root";
   if (!gSystem->AccessPathName( fname_C )) {
      input_C = TFile::Open( fname_C ); // check if file in local directory exists
   }
   if (!input_C) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input_C->GetName() << std::endl;
   
   TFile *input_E(0);
   TString fname_E = "../L1toRecoMatchPlots_tight_E.root";
   if (!gSystem->AccessPathName( fname_E )) {
      input_E = TFile::Open( fname_E ); // check if file in local directory exists
   }
   if (!input_E) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input_E->GetName() << std::endl;
   
   TFile *input_F(0);
   TString fname_F = "../L1toRecoMatchPlots_tight_F.root";
   if (!gSystem->AccessPathName( fname_F )) {
      input_F = TFile::Open( fname_F ); // check if file in local directory exists
   }
   if (!input_F) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegression           : Using input file: " << input_F->GetName() << std::endl;


   //Double_t regWeight  = 1.0;
   TTree *regTree_B = (TTree*)input_B->Get("mytree");
   TTree *regTree_C = (TTree*)input_C->Get("mytree");
   TTree *regTree_E = (TTree*)input_E->Get("mytree");
   TTree *regTree_F = (TTree*)input_F->Get("mytree");
   dataloader->AddRegressionTree( regTree_B/*, regWeight*/ );
   dataloader->AddRegressionTree( regTree_C/*, regWeight*/ );
   dataloader->AddRegressionTree( regTree_E/*, regWeight*/ );
   dataloader->AddRegressionTree( regTree_F/*, regWeight*/ );

   // This would set individual event weights (the variables defined in the
   // expression need to exist in the original TTree)
   //dataloader->SetWeightExpression( "var1", "Regression" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = " recomuon_dr >= 0.0 && recomuon_dr < 0.2 ";
   //mycut += " ( (L1muon_tfMuonIndex >= 17 && L1muon_tfMuonIndex <= 35) || (L1muon_tfMuonIndex >= 71 && L1muon_tfMuonIndex <= 89) )";
   mycut += " ( fabs( L1muon_eta ) > 0.8 && fabs( L1muon_eta ) <= 1.2 ) ";
   

   dataloader->PrepareTrainingAndTestTree( mycut,"SplitMode=Random:NormMode=NumEvents:!V" );
   

   // Book MVA methods
   
   // Neural network (MLP)
   if (Use["MLP"])
      factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15" );
      //factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP", "");

   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;

   delete factory;
   delete dataloader;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
}