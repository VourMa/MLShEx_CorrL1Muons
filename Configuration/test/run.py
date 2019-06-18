import ROOT
import os
import sys

helpText = "\n Choose what to run: --skim or --TMVATrain or --MVATest or --scout (not implemented yet)\n \
\n \
Mandatory arguments for all but --scout: --dataset, --year, --era\n \
Optional arguments for all but --scout: --useTotalEras, --ID, --extraText, --debug\n \
\n \
Mandatory arguments for --TMVATrain: --TF\n \
Optional arguments for --TMVATrain:  --guys\n \
\n \
Optional arguments for --MVATest mass: --performOn, --particle\n \
\n \
Optional arguments for --MVATest resolution: --guys, --performOn"

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(helpText)

    parser.add_option("--skim", dest="skimOpt", default=False, action="store_true", help="Do the skimming")
    parser.add_option("--TMVATrain", dest="TMVATrainOpt", default=False, action="store_true", help="Do the TMVA training")
    parser.add_option("--MVATest", dest="MVATestOpt", default=None, help="Do the MVA testing. Type \"mass\" for computation of mass spectrums, \"resolution\" for computation of resolutions or \"scout\" for application on the scouted data")

    parser.add_option("--dataset", dest="datasetOpt", default=None, help="Choose dataset")
    parser.add_option("--year", dest="yearOpt", default=None, help="Choose year")
    parser.add_option("--era", dest="eraOpt", default=None, help="Choose eras, separated by commas")
    parser.add_option("--TF", dest="TFOpt", default=None, help="Choose TF for the training, separated by commas. B for BMTF, O for OMTF, E for EMTF")
    parser.add_option("--useTotalEras", dest="useTotalErasOpt", default=False, action="store_true", help="Create/Use files with multiple eras [Default: %default]")
    parser.add_option("--ID", dest="IDOpt", default="tight", help="Choose muon ID [Default: %default]")
    parser.add_option("--extraText", dest="extraTextOpt", default=None, help="Text to be added at the end of the output file name [Default: %default]")
    parser.add_option("--debug", dest="debugOpt", default=False, action="store_true", help="Enable debugging [Default: %default]")

    parser.add_option("--guys", dest="guysOpt", default="A", help="Choose muon category for the training, separated by commas. A for all, B for pt <= 10 or pt >= 40, G for pt >= 10 and pt <= 40 [Default: %default]")
    parser.add_option("--performOn", dest="performOnOpt", default="A", help="Choose muon category to apply the training, separated by commas. A for all, B for pt <= 10 or pt >= 40, G for pt >= 10 and pt <= 40 [Default: %default]")
    parser.add_option("--TFBinMethod", dest="TFBinMethodOpt", default="Eta", help="Choose how to separate the different Track Finders (TF), Eta or Index [Default: %default]")
    parser.add_option("--particle", dest="particleOpt", default="JPsi", help="Choose particle for mass comparison, JPsi or Upsilon [Default: %default]")

    options, args = parser.parse_args()

    def interpretEras(eras,totalEras):
	    eras = options.eraOpt.split(',')
	    for era in eras:
		    totalEras += era
	    return eras, totalEras

    eras = []
    totalEras = ""
    if options.extraTextOpt: extraText = "_"+options.extraTextOpt
    else: extraText = ""

    # Skimming
    if options.skimOpt:
	    eras = options.eraOpt.split(','); totalEras = ""
	    for era in eras:
		    totalEras += era
	    ROOT.gROOT.ProcessLine(".L %s/src/L1uGMTAnalyzer/Configuration/test/skimming/L1toRecoMatcher.C+" % os.environ['CMSSW_BASE'])

	    inputFileNames = []

	    for iera, era in enumerate(eras):
		    print "Getting file of era "+era
		    inputFileNames.append('/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_NTuples/'+options.datasetOpt+'/Run'+options.yearOpt+era+'/outputL1uGMTAnalyzer.root')
		    fileExists = os.path.isfile(inputFileNames[iera])
		    if fileExists:
			    if not options.useTotalErasOpt:
				    inputChain = ROOT.TChain("events")
				    inputChain.Add(inputFileNames[iera])
				    print "with",inputChain.GetEntries(),"entries\n"
				    outFileName = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/L1toRecoMatchPlots_"+options.datasetOpt+"_"+options.IDOpt+"_"+era+extraText+".root"
				    outFile = ROOT.TFile(outFileName,"recreate")
				    print "Created output file: "+outFileName
				    #Running
				    L = ROOT.L1toRecoMatcher(inputChain)
				    L.Loop(options.IDOpt,outFile,options.debugOpt);
				    #Closing procedures
				    outFile.Write();
				    outFile.Close();
				    print "Output file written and closed!\n"
		    else: print "Input file doesn't exist..."
		    print ""

	    if options.useTotalErasOpt:
		    inputChain = ROOT.TChain("events")
		    for inputFileName in inputFileNames:
			    fileExists = os.path.isfile(inputFileName)
			    if fileExists: inputChain.Add(inputFileName)
		    print "with",inputChain.GetEntries(),"entries in total\n"
		    outFileName = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/L1toRecoMatchPlots_"+options.datasetOpt+"_"+options.IDOpt+"_"+totalEras+extraText+".root"
		    outFile = ROOT.TFile(outFileName,"recreate");
		    print "Created output file: "+outFileName
		    #Running
		    L = ROOT.L1toRecoMatcher(inputChain)
		    L.Loop(options.IDOpt,outFile,options.debugOpt);
		    #Closing procedures
		    outFile.Write();
		    outFile.Close();
		    print "Output file written and closed!"


    # TMVA Training
    elif options.TMVATrainOpt:
	    eras = options.eraOpt.split(','); totalEras = ""
	    for era in eras:
		    totalEras += era
	    TFs = options.TFOpt.split(',')
	    guys = options.guysOpt.split(',')
	    for TF in TFs:
		    for guy in guys:
			    #print ".x {CMSSW_BASE}/src/L1uGMTAnalyzer/Configuration/test/TMVATraining/TMVARegression.C(\"{datasetArg}\",\"{yearArg}\",\"{IDArg}\",\"{TFArg}\",\"{erasArg}\",{useTotalErasArg},\"{guyArg}\",\"{extraTextArg}\",\"{etaOrIndexArg}\")".format(CMSSW_BASE=os.environ['CMSSW_BASE'], datasetArg=options.datasetOpt, yearArg=options.yearOpt, IDArg=options.IDOpt, TFArg=TF, erasArg=totalEras, useTotalErasArg=int(options.useTotalErasOpt), guyArg=guy, extraTextArg=options.extraTextOpt, etaOrIndexArg=options.TFBinMethodOpt)
			    ROOT.gROOT.ProcessLine(".x {CMSSW_BASE}/src/L1uGMTAnalyzer/Configuration/test/TMVATraining/TMVARegression.C(\"{datasetArg}\",\"{yearArg}\",\"{IDArg}\",\"{TFArg}\",\"{erasArg}\",{useTotalErasArg},\"{guyArg}\",\"{extraTextArg}\",\"{etaOrIndexArg}\")".format(CMSSW_BASE=os.environ['CMSSW_BASE'], datasetArg=options.datasetOpt, yearArg=options.yearOpt, IDArg=options.IDOpt, TFArg=TF, erasArg=totalEras, useTotalErasArg=int(options.useTotalErasOpt), guyArg=guy, extraTextArg=options.extraTextOpt, etaOrIndexArg=options.TFBinMethodOpt))


    # MVA Testing
    elif options.MVATestOpt:
	    if options.MVATestOpt == "scout":
		    print "MVA Testing: Scouted Data mass spectrum\n"
		    #ROOT.gROOT.ProcessLine(".L %s/src/L1uGMTAnalyzer/Configuration/test/TMVATesting/ScoutedData.C+" % os.environ['CMSSW_BASE'])

	    elif options.MVATestOpt == "mass" or options.MVATestOpt == "resolution":
		    eras, totalEras = interpretEras(eras,totalEras)

		    if options.MVATestOpt == "mass":
			    print "MVA Testing: Mass spectrum\n"
			    ROOT.gROOT.ProcessLine(".L %s/src/L1uGMTAnalyzer/Configuration/test/TMVATesting/MassSpectrum.C+" % os.environ['CMSSW_BASE'])
		    elif options.MVATestOpt == "resolution":
			    print "MVA Testing: Resolution\n"
			    ROOT.gROOT.ProcessLine(".L %s/src/L1uGMTAnalyzer/Configuration/test/TMVATesting/Resolutions.C+" % os.environ['CMSSW_BASE'])

		    inputFileNames = []

		    for iera, era in enumerate(eras):
			    print "Getting file of era "+era
			    inputFileNames.append('/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/L1toRecoMatchPlots_'+options.datasetOpt+options.yearOpt+'_'+options.IDOpt+'_'+era+'.root')
			    fileExists = os.path.isfile(inputFileNames[iera])
			    if fileExists:
				    if not options.useTotalErasOpt:
					    inputChain = ROOT.TChain("mytree")
					    inputChain.Add(inputFileNames[iera])
					    print "with",inputChain.GetEntries(),"entries\n"

					    if options.MVATestOpt == "mass":
						    outFileName = "./TMVATesting/MassSpectrum_"+options.datasetOpt+"_"+era+"_"+options.particleOpt+"_"+options.TFBinMethodOpt+extraText+".root"
						    outFile = ROOT.TFile(outFileName,"recreate")
						    print "Created output file: "+outFileName

						    #Running
						    L = ROOT.MassSpectrum(inputChain,options.TFBinMethodOpt)
						    L.Loop(outFile,options.particleOpt,options.performOnOpt,options.TFBinMethodOpt,options.debugOpt);

					    elif options.MVATestOpt == "resolution":
						    outFileName = "./TMVATesting/Resolutions_"+options.datasetOpt+"_"+era+"_"+options.TFBinMethodOpt+"_"+options.guysOpt+"_"+options.performOnOpt+extraText+".root"
						    outFile = ROOT.TFile(outFileName,"recreate")
						    print "Created output file: "+outFileName

						    #Running
						    L = ROOT.Resolutions(inputChain)
						    L.Loop(outFile,options.guysOpt,options.performOnOpt,options.TFBinMethodOpt,options.debugOpt);

					    #Closing procedures
					    outFile.Write();
					    outFile.Close();
					    print "Output file written and closed!\n"
			    else: print "Input file doesn't exist..."
			    print ""

		    if options.useTotalErasOpt:
			    inputChain = ROOT.TChain("mytree")
			    for inputFileName in inputFileNames:
				    fileExists = os.path.isfile(inputFileName)
				    if fileExists: inputChain.Add(inputFileName)
			    print "with",inputChain.GetEntries(),"entries in total\n"

			    if options.MVATestOpt == "mass":
				    outFileName = "./TMVATesting/MassSpectrum_"+options.datasetOpt+"_"+totalEras+"_"+options.particleOpt+"_"+options.TFBinMethodOpt+extraText+".root"
				    outFile = ROOT.TFile(outFileName,"recreate")
				    print "Created output file: "+outFileName

				    #Running
				    L = ROOT.MassSpectrum(inputChain,options.TFBinMethodOpt)
				    L.Loop(outFile,options.particleOpt,options.performOnOpt,options.TFBinMethodOpt,options.debugOpt);

			    elif options.MVATestOpt == "resolution":
				    outFileName = "./TMVATesting/Resolutions_"+options.datasetOpt+"_"+totalEras+"_"+options.TFBinMethodOpt+"_"+options.guysOpt+"_"+options.performOnOpt+extraText+".root"
				    outFile = ROOT.TFile(outFileName,"recreate")
				    print "Created output file: "+outFileName

				    #Running
				    L = ROOT.Resolutions(inputChain)
				    L.Loop(outFile,options.guysOpt,options.performOnOpt,options.TFBinMethodOpt,options.debugOpt);

			    #Closing procedures
			    outFile.Write();
			    outFile.Close();
			    print "Output file written and closed!"

	    else:
		    sys.exit("Invalid MVA Test")

    else:
	    sys.exit("No script function")
