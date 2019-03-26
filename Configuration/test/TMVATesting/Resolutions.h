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
   vector<double>  *L1muon_tfMuonIndex;

   // List of branches
   TBranch        *b_recomuon_phi;   //!
   TBranch        *b_recomuon_dr;   //!
   TBranch        *b_L1muon_ptCorr;   //!
   TBranch        *b_L1muon_pt;   //!
   TBranch        *b_L1muon_eta;   //!
   TBranch        *b_L1muon_phi;   //!
   TBranch        *b_L1muon_phiAtVtx;   //!
   TBranch        *b_L1muon_charge;   //!
   TBranch        *b_L1muon_tfMuonIndex;   //!
   
   string readerTFpT[9] = {"BMTF_A","BMTF_B","BMTF_G","OMTF_A","OMTF_B","OMTF_G","EMTF_A","EMTF_B","EMTF_G"};
   TString readerTF[3] = {"B","O","E"};
   TString readerpT[3] = {"A","B","G"};
   unordered_map<string,TMVA::Reader *> readerMap;
   TMVA::Reader* readerBMTF_A = new TMVA::Reader(); TMVA::Reader* readerBMTF_B = new TMVA::Reader(); TMVA::Reader* readerBMTF_G = new TMVA::Reader();
   TMVA::Reader* readerOMTF_A = new TMVA::Reader(); TMVA::Reader* readerOMTF_B = new TMVA::Reader(); TMVA::Reader* readerOMTF_G = new TMVA::Reader();
   TMVA::Reader* readerEMTF_A = new TMVA::Reader(); TMVA::Reader* readerEMTF_B = new TMVA::Reader(); TMVA::Reader* readerEMTF_G = new TMVA::Reader();
   
   TString trainingDir = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/TMVATrainingFiles/";
   float L1muon_ptCorr_, L1muon_pt_, L1muon_eta_,L1muon_phi_, L1muon_charge_, L1muon_index_;

   Resolutions(TTree *tree=0, TString etaOrIndex = "Eta");
   virtual ~Resolutions();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   void correctThePhi(string readerTFTemp, string readerpTtemp,/*TString etaOrIndex, float L1muon_eta_, float L1muon_index_,*/ double & DphiL1Reco, double & PhiRecoReg);
   virtual void     Loop(TFile * out, TString whichGuys, TString performOnWhichGuys, TString etaOrIndex, bool debug);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Resolutions_cxx
Resolutions::Resolutions(TTree *tree, TString etaOrIndex) : fChain(0) 
{
   Init(tree);
   
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[0], readerBMTF_A) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[1], readerBMTF_B) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[2], readerBMTF_G) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[3], readerOMTF_A) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[4], readerOMTF_B) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[5], readerOMTF_G) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[6], readerEMTF_A) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[7], readerEMTF_B) );
   readerMap.insert(pair<string ,TMVA::Reader *> (readerTFpT[8], readerEMTF_G) );
   
   for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
         unordered_map<string,TMVA::Reader *>::iterator it = readerMap.find(readerTFpT[3*i + j]);
         TMVA::Reader * readerTemp = it->second;
         
         readerTemp->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
         readerTemp->AddVariable("L1muon_eta", &L1muon_eta_);
         readerTemp->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
         readerTemp->AddVariable("L1muon_charge", &L1muon_charge_);
         
         readerTemp->BookMVA("MLP", trainingDir+"TMVARegression_TF"+readerTF[i]+"_EraBCEF_Guys"+readerpT[j]+"_"+etaOrIndex+"/weights/TMVARegression_MLP.weights.xml");
      }
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
   // Set object pointer
   recomuon_phi = 0;
   recomuon_dr = 0;
   L1muon_ptCorr = 0;
   L1muon_pt = 0;
   L1muon_eta = 0;
   L1muon_phi = 0;
   L1muon_phiAtVtx = 0;
   L1muon_charge = 0;
   L1muon_tfMuonIndex = 0;
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
   fChain->SetBranchAddress("L1muon_tfMuonIndex", &L1muon_tfMuonIndex, &b_L1muon_tfMuonIndex);
   Notify();
}

Bool_t Resolutions::Notify()
{
   return kTRUE;
}

void Resolutions::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Resolutions::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef Resolutions_cxx
