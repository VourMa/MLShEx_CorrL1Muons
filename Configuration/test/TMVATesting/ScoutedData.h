//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 24 22:38:04 2019 by ROOT version 6.12/07
// from TTree t/Analysis
// found on file: /eos/cms/store/group/cmst3/group/daql1scout/run2/rootfiles/run_324970/scout_324970_000000.root
//////////////////////////////////////////////////////////

#ifndef ScoutedData_h
#define ScoutedData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class ScoutedData {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          orbit;
   UInt_t          bx;
   vector<float>   *pt;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<float>   *etap;
   vector<float>   *phip;
   vector<int>     *charge;
   vector<int>     *index;
   vector<int>     *qual;

   // List of branches
   TBranch        *b_orbit;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_etap;   //!
   TBranch        *b_phip;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_index;   //!
   TBranch        *b_qual;   //!
   
   TString trainingDir = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/TMVATrainingFiles/";
   TMVA::Reader* readerBMTF = new TMVA::Reader();
   TMVA::Reader* readerOMTF = new TMVA::Reader();
   TMVA::Reader* readerEMTF = new TMVA::Reader();
   float L1muon_ptCorr_, L1muon_pt_, L1muon_eta_,L1muon_phi_, L1muon_charge_;

   ScoutedData(TTree *tree=0, TString whichGuys = "A");
   virtual ~ScoutedData();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   void correctThePhi(float L1muon_eta_, double & DphiL1Reco, double & PhiRecoReg);
   virtual void     Loop(TFile * out, bool debug);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ScoutedData_cxx
ScoutedData::ScoutedData(TTree *tree, TString whichGuys = "A") : fChain(0) 
{
   Init(tree);
   
   readerBMTF->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
   readerBMTF->AddVariable("L1muon_eta", &L1muon_eta_);
   readerBMTF->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
   readerBMTF->AddVariable("L1muon_charge", &L1muon_charge_);
   
   readerOMTF->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
   readerOMTF->AddVariable("L1muon_eta", &L1muon_eta_);
   readerOMTF->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
   readerOMTF->AddVariable("L1muon_charge", &L1muon_charge_);
   
   readerEMTF->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
   readerEMTF->AddVariable("L1muon_eta", &L1muon_eta_);
   readerEMTF->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
   readerEMTF->AddVariable("L1muon_charge", &L1muon_charge_);
   
   readerBMTF->BookMVA("MLP", trainingDir+"TMVARegression_TFB_EraBCEF_Guys"+whichGuys+"_Eta/weights/TMVARegression_MLP.weights.xml");
   readerOMTF->BookMVA("MLP", trainingDir+"TMVARegression_TFO_EraBCEF_Guys"+whichGuys+"_Eta/weights/TMVARegression_MLP.weights.xml");
   readerEMTF->BookMVA("MLP", trainingDir+"TMVARegression_TFE_EraBCEF_Guys"+whichGuys+"_Eta/weights/TMVARegression_MLP.weights.xml");
}

ScoutedData::~ScoutedData()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ScoutedData::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ScoutedData::LoadTree(Long64_t entry)
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

void ScoutedData::Init(TTree *tree)
{
   // Set object pointer
   pt = 0;
   eta = 0;
   phi = 0;
   etap = 0;
   phip = 0;
   charge = 0;
   index = 0;
   qual = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("etap", &etap, &b_etap);
   fChain->SetBranchAddress("phip", &phip, &b_phip);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("index", &index, &b_index);
   fChain->SetBranchAddress("qual", &qual, &b_qual);
   Notify();
}

Bool_t ScoutedData::Notify()
{
   return kTRUE;
}

void ScoutedData::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ScoutedData::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef ScoutedData_cxx
