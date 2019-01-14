//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  5 15:22:33 2018 by ROOT version 6.10/09
// from TTree events/events
// found on file: /eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_test/ZeroBias/crab_L1uGMTAnalyzer_test/181005_104201/0000/outputL1uGMTAnalyzer.root
//////////////////////////////////////////////////////////

#ifndef L1toRecoMatcher_h
#define L1toRecoMatcher_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class MatchedMuons {
  public:
    vector<double> Pt, PtCorr, Eta, Phi, EtaAtVtx, PhiAtVtx, TfMuonIndex, EtaAtSt1, PhiAtSt1, EtaAtSt2, PhiAtSt2, Charge, HwQual, Dr/*, IsTight, IsMedium, IsLoose, IsGlobal, IsStandAlone*/;
    MatchedMuons ();
    int Size();
    void Clear();
};

MatchedMuons::MatchedMuons () {
  Pt.clear(); PtCorr.clear();
  Eta.clear(); EtaAtVtx.clear(); EtaAtSt1.clear(); EtaAtSt2.clear();
  Phi.clear(); PhiAtVtx.clear(); PhiAtSt1.clear(); PhiAtSt2.clear();
  TfMuonIndex.clear();
  Charge.clear();
  HwQual.clear();
  Dr.clear();
  //IsTight.clear(); IsMedium.clear(); IsLoose.clear(); IsGlobal.clear(); IsStandAlone.clear();
}

int MatchedMuons::Size () {
	return Pt.size();
}

void MatchedMuons::Clear () {
  Pt.clear(); PtCorr.clear();
  Eta.clear(); EtaAtVtx.clear(); EtaAtSt1.clear(); EtaAtSt2.clear();
  Phi.clear(); PhiAtVtx.clear(); PhiAtSt1.clear(); PhiAtSt2.clear();
  TfMuonIndex.clear();
  Charge.clear();
  HwQual.clear();
  Dr.clear();
  //IsTight.clear(); IsMedium.clear(); IsLoose.clear(); IsGlobal.clear(); IsStandAlone.clear();
}

class L1toRecoMatcher {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           counter;
   ULong64_t       Run;
   ULong64_t       Event;
   ULong64_t       Lumi;
   ULong64_t       BX;
   ULong64_t       Orbit;
   vector<float>   *muon_BX;
   vector<float>   *muon_hwPt;
   vector<float>   *muon_hwEta;
   vector<float>   *muon_hwPhi;
   vector<float>   *muon_hwCharge;
   vector<float>   *muon_hwChargeValid;
   vector<float>   *muon_hwQual;
   vector<float>   *muon_pt;
   vector<float>   *muon_eta;
   vector<float>   *muon_phi;
   vector<float>   *muon_charge;
   vector<float>   *muon_tfMuonIndex;
   vector<float>   *muon_hwEtaAtVtx;
   vector<float>   *muon_hwPhiAtVtx;
   vector<float>   *muon_EtaAtVtx;
   vector<float>   *muon_PhiAtVtx;
   vector<float>   *recomuon_pt;
   vector<float>   *recomuon_eta;
   vector<float>   *recomuon_phi;
   vector<float>   *recomuon_etaAtSt1;
   vector<float>   *recomuon_phiAtSt1;
   vector<float>   *recomuon_etaAtSt2;
   vector<float>   *recomuon_phiAtSt2;
   vector<int>     *recomuon_quality;
   vector<bool>    *recomuon_isTightMuon;
   vector<bool>    *recomuon_isMediumMuon;
   vector<bool>    *recomuon_isLooseMuon;
   vector<bool>    *recomuon_isGlobalMuon;
   vector<bool>    *recomuon_isTrackerMuon;
   vector<bool>    *recomuon_isStandAloneMuon;

   // List of branches
   TBranch        *b_counter;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_BX;   //!
   TBranch        *b_Orbit;   //!
   TBranch        *b_muon_BX;   //!
   TBranch        *b_muon_hwPt;   //!
   TBranch        *b_muon_hwEta;   //!
   TBranch        *b_muon_hwPhi;   //!
   TBranch        *b_muon_hwCharge;   //!
   TBranch        *b_muon_hwChargeValid;   //!
   TBranch        *b_muon_hwQual;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_tfMuonIndex;   //!
   TBranch        *b_muon_hwEtaAtVtx;   //!
   TBranch        *b_muon_hwPhiAtVtx;   //!
   TBranch        *b_muon_EtaAtVtx;   //!
   TBranch        *b_muon_PhiAtVtx;   //!
   TBranch        *b_recomuon_pt;   //!
   TBranch        *b_recomuon_eta;   //!
   TBranch        *b_recomuon_phi;   //!
   TBranch        *b_recomuon_etaAtSt1;   //!
   TBranch        *b_recomuon_phiAtSt1;   //!
   TBranch        *b_recomuon_etaAtSt2;   //!
   TBranch        *b_recomuon_phiAtSt2;   //!
   TBranch        *b_recomuon_quality;   //!
   TBranch        *b_recomuon_isTightMuon;   //!
   TBranch        *b_recomuon_isMediumMuon;   //!
   TBranch        *b_recomuon_isLooseMuon;   //!
   TBranch        *b_recomuon_isGlobalMuon;   //!
   TBranch        *b_recomuon_isTrackerMuon;   //!
   TBranch        *b_recomuon_isStandAloneMuon;   //!

   L1toRecoMatcher(TTree *tree=0);
   virtual ~L1toRecoMatcher();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString ID,TFile * out);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef L1toRecoMatcher_cxx
L1toRecoMatcher::L1toRecoMatcher(TTree *tree) : fChain(0) 
{
   Init(tree);
}

L1toRecoMatcher::~L1toRecoMatcher()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t L1toRecoMatcher::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t L1toRecoMatcher::LoadTree(Long64_t entry)
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

void L1toRecoMatcher::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   muon_BX = 0;
   muon_hwPt = 0;
   muon_hwEta = 0;
   muon_hwPhi = 0;
   muon_hwCharge = 0;
   muon_hwChargeValid = 0;
   muon_hwQual = 0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_charge = 0;
   muon_tfMuonIndex = 0;
   muon_hwEtaAtVtx = 0;
   muon_hwPhiAtVtx = 0;
   muon_EtaAtVtx = 0;
   muon_PhiAtVtx = 0;
   recomuon_pt = 0;
   recomuon_eta = 0;
   recomuon_phi = 0;
   recomuon_etaAtSt1 = 0;
   recomuon_phiAtSt1 = 0;
   recomuon_etaAtSt2 = 0;
   recomuon_phiAtSt2 = 0;
   recomuon_quality = 0;
   recomuon_isTightMuon = 0;
   recomuon_isMediumMuon = 0;
   recomuon_isLooseMuon = 0;
   recomuon_isGlobalMuon = 0;
   recomuon_isTrackerMuon = 0;
   recomuon_isStandAloneMuon = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("counter", &counter, &b_counter);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("BX", &BX, &b_BX);
   fChain->SetBranchAddress("Orbit", &Orbit, &b_Orbit);
   fChain->SetBranchAddress("muon_BX", &muon_BX, &b_muon_BX);
   fChain->SetBranchAddress("muon_hwPt", &muon_hwPt, &b_muon_hwPt);
   fChain->SetBranchAddress("muon_hwEta", &muon_hwEta, &b_muon_hwEta);
   fChain->SetBranchAddress("muon_hwPhi", &muon_hwPhi, &b_muon_hwPhi);
   fChain->SetBranchAddress("muon_hwCharge", &muon_hwCharge, &b_muon_hwCharge);
   fChain->SetBranchAddress("muon_hwChargeValid", &muon_hwChargeValid, &b_muon_hwChargeValid);
   fChain->SetBranchAddress("muon_hwQual", &muon_hwQual, &b_muon_hwQual);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_tfMuonIndex", &muon_tfMuonIndex, &b_muon_tfMuonIndex);
   fChain->SetBranchAddress("muon_hwEtaAtVtx", &muon_hwEtaAtVtx, &b_muon_hwEtaAtVtx);
   fChain->SetBranchAddress("muon_hwPhiAtVtx", &muon_hwPhiAtVtx, &b_muon_hwPhiAtVtx);
   fChain->SetBranchAddress("muon_EtaAtVtx", &muon_EtaAtVtx, &b_muon_EtaAtVtx);
   fChain->SetBranchAddress("muon_PhiAtVtx", &muon_PhiAtVtx, &b_muon_PhiAtVtx);
   fChain->SetBranchAddress("recomuon_pt", &recomuon_pt, &b_recomuon_pt);
   fChain->SetBranchAddress("recomuon_eta", &recomuon_eta, &b_recomuon_eta);
   fChain->SetBranchAddress("recomuon_phi", &recomuon_phi, &b_recomuon_phi);
   fChain->SetBranchAddress("recomuon_etaAtSt1", &recomuon_etaAtSt1, &b_recomuon_etaAtSt1);
   fChain->SetBranchAddress("recomuon_phiAtSt1", &recomuon_phiAtSt1, &b_recomuon_phiAtSt1);
   fChain->SetBranchAddress("recomuon_etaAtSt2", &recomuon_etaAtSt2, &b_recomuon_etaAtSt2);
   fChain->SetBranchAddress("recomuon_phiAtSt2", &recomuon_phiAtSt2, &b_recomuon_phiAtSt2);
   fChain->SetBranchAddress("recomuon_quality", &recomuon_quality, &b_recomuon_quality);
   fChain->SetBranchAddress("recomuon_isTightMuon", &recomuon_isTightMuon, &b_recomuon_isTightMuon);
   fChain->SetBranchAddress("recomuon_isMediumMuon", &recomuon_isMediumMuon, &b_recomuon_isMediumMuon);
   fChain->SetBranchAddress("recomuon_isLooseMuon", &recomuon_isLooseMuon, &b_recomuon_isLooseMuon);
   fChain->SetBranchAddress("recomuon_isGlobalMuon", &recomuon_isGlobalMuon, &b_recomuon_isGlobalMuon);
   fChain->SetBranchAddress("recomuon_isTrackerMuon", &recomuon_isTrackerMuon, &b_recomuon_isTrackerMuon);
   fChain->SetBranchAddress("recomuon_isStandAloneMuon", &recomuon_isStandAloneMuon, &b_recomuon_isStandAloneMuon);
   Notify();
}

Bool_t L1toRecoMatcher::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void L1toRecoMatcher::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t L1toRecoMatcher::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef L1toRecoMatcher_cxx