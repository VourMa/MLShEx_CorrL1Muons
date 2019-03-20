#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <sstream>


float DPhi(float phi1, float phi2) {
	float result = phi1 - phi2;
	if(result>3.14) result -= 6.28;
	if(result<=-3.14) result += 6.28;
	return result;
}

float DR(float eta1, float phi1, float eta2, float phi2){
	float result= TMath::Sqrt((eta1-eta2)*(eta1-eta2)+DPhi(phi1,phi2)*DPhi(phi1,phi2));
	return result;
}

float Mll(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2){
	TLorentzVector mu1, mu2;
	mu1.SetPtEtaPhiM(pt1,eta1,phi1,0.104);
	mu2.SetPtEtaPhiM(pt2,eta2,phi2,0.104);
	return (mu1+mu2).M();
}


TH1D* GetHisto(TChain *chain, TString histo, TString histo_name, TString cuts, int nbins=0, float lowedge=0.0, float highedge=0.0) {
	TString nbins_str = std::to_string(nbins);
	TString lowedge_str = std::to_string(lowedge);
	TString highedge_str = std::to_string(highedge);
	
	if(nbins==0) {
		chain->Draw(histo+">>"+histo_name,cuts,"goff");
		TH1D *res = (TH1D*)gDirectory->Get(histo_name);
		return res;
	}
	if(nbins>0) {
		chain->Draw(histo+">>"+histo_name+"("+nbins_str+","+lowedge_str+","+highedge_str+")",cuts,"goff");
		TH1D *res = (TH1D*)gDirectory->Get(histo_name);
		return res;
	}
	
	return NULL;
}

TH1D* GetHisto(TChain *chain, TString histo, TString histo_name, TString cuts, int nbins, float bins[]) {
	TH1D * res = new TH1D(histo_name,"",nbins,bins);
	chain->Draw(histo+">>"+histo_name,cuts,"goff");
	return res;
}

void CustomPlot(TH1D *histo, TString name, EColor color, TString xtitle) {
	histo->SetName(name);
	histo->SetLineColor(color);
	histo->SetLineWidth(3);
	histo->GetXaxis()->SetTitle(xtitle);
	histo->GetXaxis()->SetTitleSize(0.048);
	histo->GetXaxis()->SetTitleOffset(0.9);
	histo->GetYaxis()->SetTitle("Normalized events");
	histo->GetYaxis()->SetTitleSize(0.048);
	histo->GetYaxis()->SetTitleOffset(0.9);
	//Create underflow bin
	histo->SetBinContent(1,histo->GetBinContent(0) + histo->GetBinContent(1) );
	//Create overflow bin
	histo->SetBinContent(histo->GetNbinsX(),histo->GetBinContent(histo->GetNbinsX()) + histo->GetBinContent(histo->GetNbinsX() + 1) );
}

void DrawTogether(TH1D *histo, TString name, bool setlogy, TLegend *leg) {
	TCanvas *canvas = new TCanvas("canvas_"+name,"",800,600);
	
	histo->DrawNormalized();
	
	if(setlogy) canvas->SetLogy();
	
	//No legends in this implementation of the plotter but the snippets below hint at the implementation, if need be
	if(leg != nullptr) {
		leg->AddEntry(histo,"NeutrinoGun@PU200","l");
		leg->Draw();
	}
	
	canvas->Write(name);
	
	//Uncomment below to create .png files directly
//	TImage *img = TImage::Create();
//	img->FromPad(canvas);
//	img->WriteImage(name+".png");
//	delete img; //FIXME Verify that it works, never tested!
	
	if(leg != nullptr) leg->Clear();
}

void GetVariablePlots(TChain *Input, TString name, TString histo_name, TString xtitle, TString cuts, bool setlogy, TLegend *leg, int nbins=0, float lowedge=0.0, float highedge=0.0) {
	TH1D *h_signal_temp_PU0;
	if(nbins==0) h_signal_temp_PU0 = GetHisto(Input,name,histo_name,cuts);
	if(nbins>0) h_signal_temp_PU0 = GetHisto(Input,name,histo_name,cuts,nbins,lowedge,highedge);
	CustomPlot(h_signal_temp_PU0,histo_name,kBlue,xtitle);
	
	DrawTogether(h_signal_temp_PU0,histo_name,setlogy,leg);
	
	cout << "Canvas " << histo_name << " created!" << endl;
}

void GetVariablePlots(TChain *Input, TString name, TString histo_name, TString xtitle, TString cuts, bool setlogy, TLegend *leg, int nbins,float bins[]) {
	TH1D *h_signal_temp_PU0;
	h_signal_temp_PU0 = GetHisto(Input,name,histo_name,cuts,nbins,bins);
	CustomPlot(h_signal_temp_PU0,histo_name,kBlue,xtitle);
	
	DrawTogether(h_signal_temp_PU0,histo_name,setlogy,leg);
	
	cout << "Canvas " << histo_name << " created!" << endl;
}


TH2D* GetHisto2D(TChain *chain, TString histo,TString histo_name, TString cuts, TString bins) { //(xbins,xlow,xhigh,ybins,ylow,yhigh)
	chain->Draw(histo+">>"+histo_name+bins,cuts,"COLZ");
	TH2D *res = (TH2D*)gDirectory->Get(histo_name);
	return res;
}

TH2D* GetHisto2D(TChain *chain, TString histo,TString histo_name, TString cuts, int nbinsx, double binsx[], int nbinsy, double lowy, double highy) { 
	TH2D * res = new TH2D(histo_name,"",nbinsx,binsx,nbinsy,lowy,highy);
	chain->Draw(histo+">>"+histo_name,cuts,"COLZ");
	return res;
}

void CustomPlot2D(TH2D *histo, TString name, TString xtitle, TString ytitle) {
	histo->SetName(name);
	histo->GetXaxis()->SetTitle(xtitle);
	histo->GetXaxis()->SetTitleSize(0.048);
	histo->GetXaxis()->SetTitleOffset(0.9);
	histo->GetYaxis()->SetTitle(ytitle);
	histo->GetYaxis()->SetTitleSize(0.048);
	histo->GetYaxis()->SetTitleOffset(0.9);
}

void GetVariablePlots2D(TChain *Input, TString name, TString histo_name, int nbinsx, double binsx[], int nbinsy, float lowy, float highy, TString xtitle, TString ytitle, TString cuts) {
	TH2D *h_signal_temp_PU0;
	h_signal_temp_PU0 = GetHisto2D(Input,name,histo_name,cuts,nbinsx,binsx,nbinsy,lowy,highy);
	CustomPlot2D(h_signal_temp_PU0,histo_name,xtitle,ytitle);
	
	TCanvas *canvas = new TCanvas("canvas_"+histo_name,"",800,600);
	h_signal_temp_PU0->Draw("COLZ");
	
	canvas->Write(histo_name);
	
	cout << "Canvas " << histo_name << " created!" << endl;
}

void PlotEfficiency(TH1D * num,TH1D * den,TString name,TString xtitle,double low,double high) {
	TGraphAsymmErrors * eff = new TGraphAsymmErrors(num,den,"cp");
	eff->SetName(name+"_eff");
	
	TCanvas *c_eff = new TCanvas("c_eff","",200,10,700,500);
	
	eff->Draw("AP");
	eff->GetXaxis()->SetTitle(xtitle);
	eff->GetYaxis()->SetTitle("Efficiency");
	eff->GetXaxis()->SetLimits(low,high);
	eff->GetYaxis()->SetRangeUser(0.0,1.1);
	eff->Write(name+"_eff");
	
	cout << "Canvas " << name << "_eff created!" << endl;
	
	return;
}



//----------------------------Main----------------------------
void Plotter(TString fileName)
{
	
TChain *Input = new TChain("mytree");
Input->Add(fileName+".root");


TLegend * leg = nullptr; //TLegend *leg = new TLegend(0.6,0.75,0.9,0.9);
TFile * HistoFile  = new TFile(fileName+"_profilePlots.root","recreate");
TString cuts = "recomuon_dr > 0.0 && recomuon_dr < 0.2";


cout<<"File \""+fileName+"\": "<<Input->GetEntries()<<" events"<<endl;
gStyle->SetOptStat(0);

const int pt_size = 21;
float pt_bins[pt_size] = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0,12.0,14.0,16.0,18.0,20.0,30.0,50.0};

//GetVariablePlots(Input,"recomuon_dr","recomuon_dr","#DeltaR(L1,reco)","",false,leg,100, 0.0, 5.0);
//GetVariablePlots(Input,"L1muon_hwQual","L1muon_hwQual","L1 Quality",cuts,false,leg,15, 0.0, 15.0);

GetVariablePlots(Input,"recomuon_pt","recomuon_pt_den","p_{T}(reco) [GeV]","",false,leg,pt_size-1,pt_bins);
TH1D * recomuon_pt_den = (TH1D*)gDirectory->Get("recomuon_pt_den");
GetVariablePlots(Input,"recomuon_pt","recomuon_pt_num","p_{T}(reco) [GeV]",cuts,false,leg,pt_size-1,pt_bins);
TH1D * recomuon_pt_num = (TH1D*)gDirectory->Get("recomuon_pt_num");
PlotEfficiency(recomuon_pt_num,recomuon_pt_den,"recomuon_pt","p_{T}(reco)",0.0,50.0);

GetVariablePlots(Input,"L1muon_pt","L1muon_pt_den","p_{T}(L1) [GeV]","",false,leg,pt_size-1,pt_bins);
GetVariablePlots(Input,"L1muon_pt","L1muon_pt_num","p_{T}(L1) [GeV]",cuts,false,leg,pt_size-1,pt_bins);
//GetVariablePlots(Input,"(L1muon_pt - recomuon_pt) / recomuon_pt","pt_res","(p_{T}(L1) - p_{T}(reco)) / p_{T}(reco)",cuts,false,leg,600,-1.0,5.0);
//GetVariablePlots(Input,"L1muon_pt / recomuon_pt","pt_rat","p_{T}(L1) / p_{T}(reco)",cuts,false,leg,500,0.0,5.0);
//GetVariablePlots(Input,"(L1muon_ptCorr - recomuon_pt) / recomuon_pt","pt_res","(p_{T}(L1) - p_{T}(reco)) / p_{T}(reco)",cuts,false,leg,600,-1.0,5.0);
//GetVariablePlots(Input,"L1muon_ptCorr / recomuon_pt","pt_rat","p_{T}(L1) / p_{T}(reco)",cuts,false,leg,500,0.0,5.0);

//GetVariablePlots(Input,"recomuon_eta","recomuon_eta_den","#eta(reco)","",false,leg,50,-2.5,2.5);
//TH1D * recomuon_eta_den = (TH1D*)gDirectory->Get("recomuon_eta_den");
//GetVariablePlots(Input,"recomuon_eta","recomuon_eta_num","#eta(reco)",cuts,false,leg,50,-2.5,2.5);
//TH1D * recomuon_eta_num = (TH1D*)gDirectory->Get("recomuon_eta_num");
//PlotEfficiency(recomuon_eta_num,recomuon_eta_den,"recomuon_eta","#eta(reco)",-2.5,2.5);


//GetVariablePlots(Input,"L1muon_etaAtVtx - recomuon_eta","eta_dif","#eta(L1) - #eta(reco)",cuts,false,leg,80,-0.4,0.4);
/*GetVariablePlots(Input,"(L1muon_etaAtVtx - recomuon_eta) / recomuon_eta","eta_res","(#eta(L1) - #eta(reco)) / #eta(reco)",cuts,false,leg,200,-1.0,1.0);
GetVariablePlots(Input,"L1muon_etaAtVtx / recomuon_eta","eta_rat","#eta(L1) / #eta(reco)",cuts,false,leg,200,0.0,2.0);*/

//GetVariablePlots(Input,"recomuon_phi","recomuon_phi_den","#phi(reco)","",false,leg,64,-3.2,3.2);
//TH1D * recomuon_phi_den = (TH1D*)gDirectory->Get("recomuon_phi_den");
//GetVariablePlots(Input,"recomuon_phi","recomuon_phi_num","#phi(reco)",cuts,false,leg,64,-3.2,3.2);
//TH1D * recomuon_phi_num = (TH1D*)gDirectory->Get("recomuon_phi_num");
//PlotEfficiency(recomuon_phi_num,recomuon_phi_den,"recomuon_phi","#phi(reco)",-3.2,3.2);


//GetVariablePlots(Input,"DPhi(L1muon_phiAtVtx,recomuon_phi)","phi_dif","#phi(L1) - #phi(reco)",cuts,false,leg,80,-0.4,0.4);
/*GetVariablePlots(Input,"DPhi(L1muon_phiAtVtx,recomuon_phi) / recomuon_phi","phi_res","(#phi(L1) - #phi(reco)) / #phi(reco)",cuts,false,leg,600,-3.0,3.0);
GetVariablePlots(Input,"L1muon_phiAtVtx / recomuon_phi","phi_rat","#phi(L1) / #phi(reco)",cuts,false,leg,600,0.0,6.0);*/

//GetVariablePlots(Input,"Mll(recomuon_pt[0], recomuon_eta[0], recomuon_phi[0], recomuon_pt[1], recomuon_eta[1], recomuon_phi[1])","reco_mll","M_{ll}(reco)",cuts,false,leg,200,0.0,100.0);
//GetVariablePlots(Input,"Mll(L1muon_pt[0], L1muon_etaAtVtx[0], L1muon_phiAtVtx[0], L1muon_pt[1], L1muon_etaAtVtx[1], L1muon_phiAtVtx[1])","L1muon_mll","M_{ll}(L1)",cuts,false,leg,200,0.0,100.0);

//GetVariablePlots2D(Input, "L1muon_pt:recomuon_pt","L1muon_pt_vs_recomuon_pt","(50,0.0,50.0,50,0.0,50.0)","p_{T}(reco) [GeV]","p_{T}(L1) [GeV]",cuts);
//GetVariablePlots2D(Input, "recomuon_phi:recomuon_eta","recomuon_phi_vs_recomuon_eta_den","(20,-2.5,2.5,32,-3.2,3.2)","#eta(reco)","#phi(reco)","");
//GetVariablePlots2D(Input, "recomuon_phi:recomuon_eta","recomuon_phi_vs_recomuon_eta_num","(20,-2.5,2.5,32,-3.2,3.2)","#eta(reco)","#phi(reco)",cuts);
//GetVariablePlots2D(Input, "DPhi(L1muon_phiAtVtx,recomuon_phi):(L1muon_pt / recomuon_pt)","phi_dif_vs_pt_rat","(500,0.0,5.0,200,-1.0,1.0)","p_{T}(L1) / p_{T}(reco)","#phi(L1) - #phi(reco)",cuts);
//GetVariablePlots2D(Input, "(L1muon_etaAtVtx - recomuon_eta):(L1muon_pt / recomuon_pt)","eta_dif_vs_pt_rat","(500,0.0,5.0,80,-0.4,0.4)","p_{T}(L1) / p_{T}(reco)","#eta(L1) - #eta(reco)",cuts);
//GetVariablePlots2D(Input, "DR(L1muon_etaAtVtx, L1muon_phiAtVtx, recomuon_eta, recomuon_phi):(L1muon_pt / recomuon_pt)","dr_vs_pt_rat","(500,0.0,5.0,40,0.0,0.4)","p_{T}(L1) / p_{T}(reco)","#DeltaR(L1,reco)",cuts);


////___ratio_pt in 2D to get correction factos from profile plots___
const int ptCorr_size = 57;
double ptCorr_bins[ptCorr_size] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,16.5,17.0,17.5,18.0,18.5,19.0,19.5,20.0,22.0,24.0,26.0,28.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,80.0,100.0,120.0,140.0,160.0};

GetVariablePlots2D(Input, "(L1muon_pt / recomuon_pt):L1muon_pt","pt_rat_vs_L1muon_pt",ptCorr_size-1,ptCorr_bins,50,0.0,5.0,"p_{T}(L1) [GeV]","p_{T}(L1) / p_{T}(reco)",cuts);
TH2D * pt_rat_vs_L1muon_pt = (TH2D*)gDirectory->Get("pt_rat_vs_L1muon_pt");
pt_rat_vs_L1muon_pt->ProfileX();
//___BMTF___
GetVariablePlots2D(Input, "(L1muon_pt / recomuon_pt):L1muon_pt","pt_rat_vs_L1muon_pt_BMTF",ptCorr_size-1,ptCorr_bins,50,0.0,5.0,"p_{T}(L1) [GeV]","p_{T}(L1) / p_{T}(reco)",cuts+" && fabs(L1muon_eta) < 0.8");
TH2D * pt_rat_vs_L1muon_pt_BMTF = (TH2D*)gDirectory->Get("pt_rat_vs_L1muon_pt_BMTF");
pt_rat_vs_L1muon_pt_BMTF->ProfileX();
//___OMTF___
GetVariablePlots2D(Input, "(L1muon_pt / recomuon_pt):L1muon_pt","pt_rat_vs_L1muon_pt_OMTF",ptCorr_size-1,ptCorr_bins,50,0.0,5.0,"p_{T}(L1) [GeV]","p_{T}(L1) / p_{T}(reco)",cuts+" && fabs(L1muon_eta) > 0.8 && fabs(L1muon_eta) < 1.2");
TH2D * pt_rat_vs_L1muon_pt_OMTF = (TH2D*)gDirectory->Get("pt_rat_vs_L1muon_pt_OMTF");
pt_rat_vs_L1muon_pt_OMTF->ProfileX();
//___EMTF___
GetVariablePlots2D(Input, "(L1muon_pt / recomuon_pt):L1muon_pt","pt_rat_vs_L1muon_pt_EMTF",ptCorr_size-1,ptCorr_bins,50,0.0,5.0,"p_{T}(L1) [GeV]","p_{T}(L1) / p_{T}(reco)",cuts+" && fabs(L1muon_eta) > 1.2");
TH2D * pt_rat_vs_L1muon_pt_EMTF = (TH2D*)gDirectory->Get("pt_rat_vs_L1muon_pt_EMTF");
pt_rat_vs_L1muon_pt_EMTF->ProfileX();

const int L1pt_size = 51;
double L1pt_bins[L1pt_size] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0};
GetVariablePlots2D(Input, "recomuon_pt:L1muon_pt","recomuon_pt_vs_L1muon_pt_fixedBinning",L1pt_size-1,L1pt_bins,50,0.0,50.0,"p_{T}(L1) [GeV]","p_{T}(reco) [GeV]",cuts);
GetVariablePlots2D(Input, "(L1muon_pt / recomuon_pt):L1muon_pt","pt_rat_vs_L1muon_pt_fixedBinning",L1pt_size-1,L1pt_bins,50,0.0,5.0,"p_{T}(L1) [GeV]","p_{T}(L1) / p_{T}(reco)",cuts);

//pT bump study
//GetVariablePlots(Input,"recomuon_dr","recomuon_dr_ratio","#DeltaR(L1,reco)",cuts+" && L1muon_pt / recomuon_pt < 0.6 && L1muon_pt / recomuon_pt > 0.4",false,leg,100, 0.0, 5.0);
//GetVariablePlots(Input,"recomuon_pt","recomuon_pt_num_ratio","p_{T}(reco) [GeV]",cuts+" && L1muon_pt / recomuon_pt < 0.6 && L1muon_pt / recomuon_pt > 0.4",false,leg,pt_size-1,pt_bins);
//GetVariablePlots(Input,"L1muon_pt","L1muon_pt_num_ratio","p_{T}(L1) [GeV]",cuts+" && L1muon_pt / recomuon_pt < 0.6 && L1muon_pt / recomuon_pt > 0.4",false,leg,pt_size-1,pt_bins);


HistoFile->Write();
cout << HistoFile->GetName() << " created!"<<endl;

return;

}
