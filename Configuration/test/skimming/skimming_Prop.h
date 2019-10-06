//#ifndef skimming_h
//#define skimming_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "vector"

class MatchedMuons {
  public:
    vector<double> Pt, PtCorr, Eta, Phi, EtaAtVtx, PhiAtVtx, TfMuonIndex, EtaAtSt1, PhiAtSt1, EtaAtSt2, PhiAtSt2, Charge, HwQual, Dr;
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
}

class skimming {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   vector<float>   *muon_etaAtVtx;
   vector<float>   *muon_phiAtVtx;
   vector<float>   *recomuon_pt;
   vector<float>   *recomuon_eta;
   vector<float>   *recomuon_phi;
   vector<float>   *recomuon_etaAtSt1;
   vector<float>   *recomuon_phiAtSt1;
   vector<float>   *recomuon_etaAtSt2;
   vector<float>   *recomuon_phiAtSt2;
   vector<bool>    *recomuon_isTightMuon;
   vector<bool>    *recomuon_isMediumMuon;
   vector<bool>    *recomuon_isLooseMuon;
   vector<bool>    *recomuon_isGlobalMuon;
   vector<bool>    *recomuon_isTrackerMuon;
   vector<bool>    *recomuon_isStandAloneMuon;

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
   TBranch        *b_muon_etaAtVtx;   //!
   TBranch        *b_muon_phiAtVtx;   //!
   TBranch        *b_recomuon_pt;   //!
   TBranch        *b_recomuon_eta;   //!
   TBranch        *b_recomuon_phi;   //!
   TBranch        *b_recomuon_etaAtSt1;   //!
   TBranch        *b_recomuon_phiAtSt1;   //!
   TBranch        *b_recomuon_etaAtSt2;   //!
   TBranch        *b_recomuon_phiAtSt2;   //!
   TBranch        *b_recomuon_isTightMuon;   //!
   TBranch        *b_recomuon_isMediumMuon;   //!
   TBranch        *b_recomuon_isLooseMuon;   //!
   TBranch        *b_recomuon_isGlobalMuon;   //!
   TBranch        *b_recomuon_isTrackerMuon;   //!
   TBranch        *b_recomuon_isStandAloneMuon;   //!

   skimming(TTree *tree=0);
   virtual ~skimming();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString ID,TString out,bool debug, bool onlytxt);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef skimming_cxx
skimming::skimming(TTree *tree) : fChain(0) 
{
   Init(tree);
}

skimming::~skimming()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skimming::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skimming::LoadTree(Long64_t entry)
{
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skimming::Init(TTree *tree)
{
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
   muon_etaAtVtx = 0;
   muon_phiAtVtx = 0;
   recomuon_pt = 0;
   recomuon_eta = 0;
   recomuon_phi = 0;
   recomuon_etaAtSt1 = 0;
   recomuon_phiAtSt1 = 0;
   recomuon_etaAtSt2 = 0;
   recomuon_phiAtSt2 = 0;
   recomuon_isTightMuon = 0;
   recomuon_isMediumMuon = 0;
   recomuon_isLooseMuon = 0;
   recomuon_isGlobalMuon = 0;
   recomuon_isTrackerMuon = 0;
   recomuon_isStandAloneMuon = 0;

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
   fChain->SetBranchAddress("muon_etaAtVtx", &muon_etaAtVtx, &b_muon_etaAtVtx);
   fChain->SetBranchAddress("muon_phiAtVtx", &muon_phiAtVtx, &b_muon_phiAtVtx);
   fChain->SetBranchAddress("recomuon_pt", &recomuon_pt, &b_recomuon_pt);
   fChain->SetBranchAddress("recomuon_eta", &recomuon_eta, &b_recomuon_eta);
   fChain->SetBranchAddress("recomuon_phi", &recomuon_phi, &b_recomuon_phi);
   fChain->SetBranchAddress("recomuon_etaAtSt1", &recomuon_etaAtSt1, &b_recomuon_etaAtSt1);
   fChain->SetBranchAddress("recomuon_phiAtSt1", &recomuon_phiAtSt1, &b_recomuon_phiAtSt1);
   fChain->SetBranchAddress("recomuon_etaAtSt2", &recomuon_etaAtSt2, &b_recomuon_etaAtSt2);
   fChain->SetBranchAddress("recomuon_phiAtSt2", &recomuon_phiAtSt2, &b_recomuon_phiAtSt2);
   fChain->SetBranchAddress("recomuon_isTightMuon", &recomuon_isTightMuon, &b_recomuon_isTightMuon);
   fChain->SetBranchAddress("recomuon_isMediumMuon", &recomuon_isMediumMuon, &b_recomuon_isMediumMuon);
   fChain->SetBranchAddress("recomuon_isLooseMuon", &recomuon_isLooseMuon, &b_recomuon_isLooseMuon);
   fChain->SetBranchAddress("recomuon_isGlobalMuon", &recomuon_isGlobalMuon, &b_recomuon_isGlobalMuon);
   fChain->SetBranchAddress("recomuon_isTrackerMuon", &recomuon_isTrackerMuon, &b_recomuon_isTrackerMuon);
   fChain->SetBranchAddress("recomuon_isStandAloneMuon", &recomuon_isStandAloneMuon, &b_recomuon_isStandAloneMuon);
   Notify();
}

Bool_t skimming::Notify()
{
   return kTRUE;
}

void skimming::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skimming::Cut(Long64_t entry)
{
   return 1;
}
//#endif // #ifdef skimming_cxx
