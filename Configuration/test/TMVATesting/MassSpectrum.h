//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 20 12:28:44 2018 by ROOT version 6.10/09
// from TTree mytree/mytree
// found on file: L1uGMTPlots_tight_tree.root
//////////////////////////////////////////////////////////

//#ifndef MassSpectrum_h
//#define MassSpectrum_h

#include <unordered_map>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"
#include "TMVA/PyMethodBase.h"

// Header file for the classes stored in the TTree if any.
#include "vector"

class MassSpectrum {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   int  recomuon_N;
   vector<double>  *recomuon_pt;
   vector<double>  *recomuon_eta;
   vector<double>  *recomuon_phi;
   vector<double>  *recomuon_dr;
   vector<double>  *L1muon_pt;
   vector<double>  *L1muon_ptCorr;
   vector<double>  *L1muon_eta;
   vector<double>  *L1muon_etaAtVtx;
   vector<double>  *L1muon_phi;
   vector<double>  *L1muon_phiAtVtx;
   vector<double>  *L1muon_charge;
   vector<double>  *L1muon_tfMuonIndex;

   // List of branches
   TBranch        *b_recomuon_N;   //!
   TBranch        *b_recomuon_pt;   //!
   TBranch        *b_recomuon_eta;   //!
   TBranch        *b_recomuon_phi;   //!
   TBranch        *b_recomuon_dr;   //!
   TBranch        *b_L1muon_pt;   //!
   TBranch        *b_L1muon_ptCorr;   //!
   TBranch        *b_L1muon_eta;   //!
   TBranch        *b_L1muon_etaAtVtx;   //!
   TBranch        *b_L1muon_phi;   //!
   TBranch        *b_L1muon_phiAtVtx;   //!
   TBranch        *b_L1muon_charge;   //!
   TBranch        *b_L1muon_tfMuonIndex;   //!
   
   string readerTF[3] = {"B","O","E"};
   TString readerpT[1] = {"A"};
   unordered_map<string,TMVA::Reader *> readerMap;
   TMVA::Reader* reader[3] = {};
   
<<<<<<< HEAD
   TString trainingDir = "/home/cmsdas/public/store/MLShortExercise/TMVATrainingFiles/"; // Change accordingly
=======
   TString trainingDir = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/TMVATrainingFiles/"; // Change accordingly
   TString era = "BCEF"; //Change accordingly
>>>>>>> 5aaa8d3... Modifications to main scripts to accommodate the bonus exercise
   float L1muon_ptCorr_, L1muon_pt_, L1muon_eta_,L1muon_phi_, L1muon_charge_, L1muon_index_;

   MassSpectrum(TTree *tree=0, TString etaOrIndex = "Eta");
   virtual ~MassSpectrum();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   void correctThePhi(string readerTFTemp, string whichGuys, double & DphiL1Reco, double & PhiRecoReg);
   virtual void     Loop(TFile * out, string particle, TString whichGuys, TString etaOrIndex, Long64_t maxEvents);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef MassSpectrum_cxx
MassSpectrum::MassSpectrum(TTree *tree, TString etaOrIndex) : fChain(0)
{
   Init(tree);
   
   for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 1; j++) {
         string key = readerTF[i]; key.append("MTF_"); key.append(readerpT[j]);
         reader[3*j+i] = new TMVA::Reader();
         
         reader[3*j+i]->AddVariable("L1muon_ptCorr", &L1muon_ptCorr_);
         reader[3*j+i]->AddVariable("L1muon_eta", &L1muon_eta_);
         reader[3*j+i]->AddVariable("L1muon_phiAtVtx", &L1muon_phi_);
         reader[3*j+i]->AddVariable("L1muon_charge", &L1muon_charge_);
         
         reader[3*j+i]->BookMVA("MLP", trainingDir+"TMVARegression_TF"+readerTF[i]+"_Era"+era+"_Guys"+readerpT[j]+"_"+etaOrIndex+"/weights/TMVARegression_MLP.weights.xml");
         
         readerMap.insert(pair<string,TMVA::Reader *> (key, reader[3*j+i]) );
      }
   }
}

MassSpectrum::~MassSpectrum()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MassSpectrum::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MassSpectrum::LoadTree(Long64_t entry)
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

void MassSpectrum::Init(TTree *tree)
{
   // Set object pointer
   recomuon_pt = 0;
   recomuon_eta = 0;
   recomuon_phi = 0;
   recomuon_dr = 0;
   L1muon_pt = 0;
   L1muon_ptCorr = 0;
   L1muon_eta = 0;
   L1muon_etaAtVtx = 0;
   L1muon_phi = 0;
   L1muon_phiAtVtx = 0;
   L1muon_charge = 0;
   L1muon_tfMuonIndex = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("recomuon_N", &recomuon_N, &b_recomuon_N);
   fChain->SetBranchAddress("recomuon_pt", &recomuon_pt, &b_recomuon_pt);
   fChain->SetBranchAddress("recomuon_eta", &recomuon_eta, &b_recomuon_eta);
   fChain->SetBranchAddress("recomuon_phi", &recomuon_phi, &b_recomuon_phi);
   fChain->SetBranchAddress("recomuon_dr", &recomuon_dr, &b_recomuon_dr);
   fChain->SetBranchAddress("L1muon_pt", &L1muon_pt, &b_L1muon_pt);
   fChain->SetBranchAddress("L1muon_ptCorr", &L1muon_ptCorr, &b_L1muon_ptCorr);
   fChain->SetBranchAddress("L1muon_eta", &L1muon_eta, &b_L1muon_eta);
   fChain->SetBranchAddress("L1muon_etaAtVtx", &L1muon_etaAtVtx, &b_L1muon_etaAtVtx);
   fChain->SetBranchAddress("L1muon_phi", &L1muon_phi, &b_L1muon_phi);
   fChain->SetBranchAddress("L1muon_phiAtVtx", &L1muon_phiAtVtx, &b_L1muon_phiAtVtx);
   fChain->SetBranchAddress("L1muon_charge", &L1muon_charge, &b_L1muon_charge);
   fChain->SetBranchAddress("L1muon_tfMuonIndex", &L1muon_tfMuonIndex, &b_L1muon_tfMuonIndex);
   Notify();
}

Bool_t MassSpectrum::Notify()
{
   return kTRUE;
}

void MassSpectrum::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MassSpectrum::Cut(Long64_t entry)
{
   return 1;
}
//#endif // #ifdef MassSpectrum_cxx
