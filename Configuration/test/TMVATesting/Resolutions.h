//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 20 12:28:44 2018 by ROOT version 6.10/09
// from TTree mytree/mytree
// found on file: L1uGMTPlots_tight_tree.root
//////////////////////////////////////////////////////////

#ifndef Resolutions_h
#define Resolutions_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class Resolutions {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<double>  *recomuon_phi;
   vector<double>  *recomuon_dr;
   vector<double>  *L1muon_ptCorr;
   vector<double>  *L1muon_pt;
   vector<double>  *L1muon_eta;
   vector<double>  *L1muon_phi;
   vector<double>  *L1muon_phiAtVtx;
   vector<double>  *L1muon_charge;
   //vector<double>  *L1muon_tfMuonIndex;

   // List of branches
   TBranch        *b_recomuon_phi;   //!
   TBranch        *b_recomuon_dr;   //!
   TBranch        *b_L1muon_ptCorr;   //!
   TBranch        *b_L1muon_pt;   //!
   TBranch        *b_L1muon_eta;   //!
   TBranch        *b_L1muon_phi;   //!
   TBranch        *b_L1muon_phiAtVtx;   //!
   TBranch        *b_L1muon_charge;   //!
   //TBranch        *b_L1muon_tfMuonIndex;   //!
   
   TMVA::Reader* readerBMTF = new TMVA::Reader();
   TMVA::Reader* readerOMTF = new TMVA::Reader();
   TMVA::Reader* readerEMTF = new TMVA::Reader();
   float L1muon_ptCorr_,L1muon_pt_, L1muon_eta_,L1muon_phi_, L1muon_charge_;

   Resolutions(TTree *tree=0, TString etaOrIndex = "eta", TString whichGuys = "A");
   virtual ~Resolutions();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   void correctThePhi(float L1muon_eta_, double & DphiL1Reco, double & PhiRecoReg);
   virtual void     Loop(TFile * out, TString performOnWhichGuys, bool debug);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Resolutions_cxx
Resolutions::Resolutions(TTree *tree, TString etaOrIndex, TString whichGuys) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   //if (tree == 0) {
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("L1uGMTPlots_tight_tree.root");
      //if (!f || !f->IsOpen()) {
         //f = new TFile("L1uGMTPlots_tight_tree.root");
      //}
      //f->GetObject("mytree",tree);

   //}
   Init(tree);
   
   readerBMTF->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
   if(etaOrIndex == "index") readerBMTF->AddSpectator("L1muon_pt", &L1muon_pt_);
   readerBMTF->AddVariable("L1muon_eta", &L1muon_eta_);
   readerBMTF->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
   readerBMTF->AddVariable("L1muon_charge", &L1muon_charge_);
   
   readerOMTF->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
   if(etaOrIndex == "index") readerOMTF->AddSpectator("L1muon_pt", &L1muon_pt_);
   readerOMTF->AddVariable("L1muon_eta", &L1muon_eta_);
   readerOMTF->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
   readerOMTF->AddVariable("L1muon_charge", &L1muon_charge_);
   
   readerEMTF->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
   if(etaOrIndex == "index") readerEMTF->AddSpectator("L1muon_pt", &L1muon_pt_);
   readerEMTF->AddVariable("L1muon_eta", &L1muon_eta_);
   readerEMTF->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
   readerEMTF->AddVariable("L1muon_charge", &L1muon_charge_);
   
   if(etaOrIndex == "eta") {
	   readerBMTF->BookMVA("MLP", "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer/CMSSW_10_1_9_patch1/src/L1uGMTAnalyzer/Configuration/test/TMVA_PhiExtrapolation_New/dataset_Impr_00_08/weights/TMVARegression_MLP.weights.xml");
	   readerOMTF->BookMVA("MLP", "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer/CMSSW_10_1_9_patch1/src/L1uGMTAnalyzer/Configuration/test/TMVA_PhiExtrapolation_New/dataset_Impr_08_12/weights/TMVARegression_MLP.weights.xml");
	   readerEMTF->BookMVA("MLP", "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer/CMSSW_10_1_9_patch1/src/L1uGMTAnalyzer/Configuration/test/TMVA_PhiExtrapolation_New/dataset_Impr_12_24/weights/TMVARegression_MLP.weights.xml");
   }
   else if(etaOrIndex == "index") {
	   readerBMTF->BookMVA("MLP", "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/test/TMVATraining/TMVARegression_TFB_EraBCEF_Guys"+whichGuys+"/weights/TMVARegression_MLP.weights.xml");
	   readerOMTF->BookMVA("MLP", "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/test/TMVATraining/TMVARegression_TFO_EraBCEF_Guys"+whichGuys+"/weights/TMVARegression_MLP.weights.xml");
	   readerEMTF->BookMVA("MLP", "/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/test/TMVATraining/TMVARegression_TFE_EraBCEF_Guys"+whichGuys+"/weights/TMVARegression_MLP.weights.xml");
   }
}

Resolutions::~Resolutions()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Resolutions::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Resolutions::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Resolutions::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   recomuon_phi = 0;
   recomuon_dr = 0;
   L1muon_ptCorr = 0;
   L1muon_pt = 0;
   L1muon_eta = 0;
   L1muon_phi = 0;
   L1muon_phiAtVtx = 0;
   L1muon_charge = 0;
   //L1muon_tfMuonIndex = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("recomuon_phi", &recomuon_phi, &b_recomuon_phi);
   fChain->SetBranchAddress("recomuon_dr", &recomuon_dr, &b_recomuon_dr);
   fChain->SetBranchAddress("L1muon_ptCorr", &L1muon_ptCorr, &b_L1muon_ptCorr);
   fChain->SetBranchAddress("L1muon_pt", &L1muon_pt, &b_L1muon_pt);
   fChain->SetBranchAddress("L1muon_eta", &L1muon_eta, &b_L1muon_eta);
   fChain->SetBranchAddress("L1muon_phi", &L1muon_phi, &b_L1muon_phi);
   fChain->SetBranchAddress("L1muon_phiAtVtx", &L1muon_phiAtVtx, &b_L1muon_phiAtVtx);
   fChain->SetBranchAddress("L1muon_charge", &L1muon_charge, &b_L1muon_charge);
   //fChain->SetBranchAddress("L1muon_tfMuonIndex", &L1muon_tfMuonIndex, &b_L1muon_tfMuonIndex);
   Notify();
}

Bool_t Resolutions::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Resolutions::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Resolutions::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Resolutions_cxx
