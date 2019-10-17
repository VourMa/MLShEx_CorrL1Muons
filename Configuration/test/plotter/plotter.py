#!/usr/bin/env python
import sys
import os
import ROOT
import argparse

# Useful functions
def canvasToImage(canvas):
    img = ROOT.TImage.Create()
    img.FromPad(canvas)
    img.WriteImage(str( canvas.GetName() )+".png")

color = { "Reco" : 1, "L1" : 2, "L1Corr" : 4 }
def histoCosmetics(histo,kind,test,xAxisTitle):
    histo.SetLineColor(color[kind]); histo.SetLineWidth(2)
    histo.GetXaxis().SetTitle(xAxisTitle); histo.GetXaxis().SetTitleSize(0.05); histo.GetXaxis().SetTitleOffset(0.90)
    histo.GetYaxis().SetTitle(""); histo.GetYaxis().SetTitleSize(0.05); histo.GetYaxis().SetTitleOffset(0.90) 
    return histo.GetBinContent(histo.GetMaximumBin())
    

# Argument/Flag definition
parser = argparse.ArgumentParser()
parser.add_argument("--test", dest="test", default=None, required=True, help="Choose test to plot: 'Resolutions' or 'MassSpectrum'")
parser.add_argument("--dir", dest="dir", default=None, required=True, help="Provide file directory")
parser.add_argument("--sample", dest="sample", default=None, required=True, help="Provide name of the sample")
args = parser.parse_args()

# Argument/Flag checking
if args.test not in ["Resolutions", "MassSpectrum"]: raise RuntimeError("The argument 'test' should be either 'Resolutions' or 'MassSpectrum'")
fileExists = os.path.isfile(args.dir)
if not fileExists: raise RuntimeError("File cannot be found")


# I/O and general ROOT settings
inFile = ROOT.TFile( args.dir, "read")
outFile = ROOT.TFile( args.test+"_"+args.sample+".root", "recreate")
ROOT.gROOT.SetBatch(1); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)


# Get, draw and superimpose plots
if args.test == "Resolutions":
    TFs = {"" : "Inclusive", "B" : "Barrel", "O" : "Overlap", "E" : "Endcap"}
    for TF in TFs.keys():
        name = "DeltaPhi_Resolution_"+args.sample+"_"+TFs[TF]+""
        canvas = ROOT.TCanvas(name,name)

        histoExt = inFile.Get("h_DphiExtReco") if TF == "" else inFile.Get("h_DphiExtReco_"+TF)
        histoExt.SetName("L1"); histoExt.SetTitle("L1"); histoExt.Draw("same")
        maxExt = histoCosmetics(histoExt,"L1",args.test,"#Delta#phi(L1-Reco)")

        histoReg = inFile.Get("h_DphiRegReco") if TF == "" else inFile.Get("h_DphiRegReco_"+TF)
        histoReg.SetName("L1Corr"); histoReg.SetTitle("L1Corr"); histoReg.Draw("same")
        maxReg = histoCosmetics(histoReg,"L1Corr",args.test,"#Delta#phi(L1-Reco)")

        histoExt.SetMaximum( 1.1*max(maxExt,maxReg) )
        canvas.BuildLegend(0.7,0.7,0.9,0.9)
        outFile.cd(); canvas.Write(); canvasToImage(canvas)

elif args.test == "MassSpectrum":
    canvasSpectrum = ROOT.TCanvas("Mass_Spectrum_"+args.sample,"Mass_Spectrum_"+args.sample)
    histoReco = inFile.Get("h_recomll"); histoReco.SetName("Reco"); histoReco.SetTitle("Reco"); histoReco.Draw("same"); maxReco = histoCosmetics(histoReco,"Reco",args.test,"M_{ll} [GeV]")
    histoL1 = inFile.Get("h_L1mll"); histoL1.SetName("L1"); histoL1.SetTitle("L1"); histoL1.Draw("same"); maxL1 = histoCosmetics(histoL1,"L1",args.test,"M_{ll} [GeV]")
    histoL1Corr = inFile.Get("h_L1mllCorr"); histoL1Corr.SetName("L1Corr"); histoL1Corr.SetTitle("L1Corr"); histoL1Corr.Draw("same"); maxL1Corr = histoCosmetics(histoL1Corr,"L1Corr",args.test,"M_{ll} [GeV]")
    histoReco.SetMaximum( 1.1*max(maxReco,maxL1,maxL1Corr) ); histoReco.GetXaxis().SetRangeUser(0.0,10.0); canvasSpectrum.BuildLegend(0.7,0.7,0.9,0.9); outFile.cd(); canvasSpectrum.Write(); canvasToImage(canvasSpectrum)

    canvasResolution = ROOT.TCanvas("Mass_Resolution_"+args.sample,"Mass_Resolution_"+args.sample); canvasResolution.cd()
    histoL1 = inFile.Get("h_L1Dmll_o_mll_particle"); histoL1.SetName("L1"); histoL1.SetTitle("L1"); histoL1.Draw("same"); maxL1 = histoCosmetics(histoL1,"L1",args.test,"(M_{reco} - M_{L1}) / M_{reco} (J/#psi)")
    histoL1Corr = inFile.Get("h_L1Dmll_o_mllCorr_particle"); histoL1Corr.SetName("L1Corr"); histoL1Corr.SetTitle("L1Corr"); histoL1Corr.Draw("same"); maxL1Corr = histoCosmetics(histoL1Corr,"L1Corr",args.test,"(M_{reco} - M_{L1}) / M_{reco} (J/#psi)")
    histoL1.SetMaximum( 1.1*max(maxL1,maxL1Corr) ); canvasResolution.BuildLegend(0.7,0.7,0.9,0.9); outFile.cd(); canvasResolution.Write(); canvasToImage(canvasResolution)


# Closing procedures
outFile.Write()
outFile.Close()
inFile.Close()
