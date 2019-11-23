import os
import argparse
from ROOT import TMVA, TFile, TTree, TCut, TChain, gROOT

from keras.models import Sequential
from keras.layers.core import Dense,Reshape,Flatten
from keras.layers.convolutional import Conv1D
from keras.optimizers import Adam,SGD,RMSprop

TMVA.PyMethodBase.PyInitialize()


def TFCut(mycutOld, TF, etaOrIndex = "Eta"):
    mycut = mycutOld
    if etaOrIndex == "Index":
        if TF == "B": mycut = mycut + TCut(" L1muon_tfMuonIndex >= 36 && L1muon_tfMuonIndex <= 70 ")
        elif TF == "O": mycut = mycut + TCut(" (L1muon_tfMuonIndex >= 17 && L1muon_tfMuonIndex <= 35) || (L1muon_tfMuonIndex >= 71 && L1muon_tfMuonIndex <= 89) ")
        elif TF == "E": mycut = mycut + TCut(" (L1muon_tfMuonIndex >= 0 && L1muon_tfMuonIndex <= 16) || (L1muon_tfMuonIndex >= 90 && L1muon_tfMuonIndex <= 107) ")
        else: raise RuntimeError("Invalid TF, exiting.")
    elif etaOrIndex == "Eta":
        if TF == "B": mycut = mycut + TCut(" fabs(L1muon_eta) < 0.8 ")
        elif TF == "O": mycut = mycut + TCut(" fabs(L1muon_eta) >= 0.8 && fabs(L1muon_eta) < 1.2 ")
        elif TF == "E": mycut = mycut + TCut(" fabs(L1muon_eta) >= 1.2 ")
	else: raise RuntimeError("Invalid TF, exiting.")	
    else: raise RuntimeError("Invalid TF binning, neither 'Eta' nor 'Index'. Exiting...")
    return mycut

def GuysCut(mycutOld, guys):
    mycut = mycutOld
    if guys == "G": mycut = mycut + TCut(" L1muon_pt >= 10 && L1muon_pt <= 40 ")
    elif guys == "B": mycut = mycut + TCut(" L1muon_pt <= 10 || L1muon_pt >= 40 ")
    elif guys == "L": mycut = mycut + TCut(" L1muon_pt <= 10 ")
    elif guys == "H": mycut = mycut + TCut(" L1muon_pt >= 40 ")
    else: print "Neither G(ood) nor B(ad) guys selected ==> Inclusive training."
    return mycut


def InitializeFile(dataset, year, ID, fileEras):
    inputFile = None
    NTupleDir = "/eos/cms/store/cmst3/user/evourlio/L1uGMTAnalyzer_Trees/"
    fname = NTupleDir+"L1toRecoMatchPlots_"+dataset+year+"_"+ID+"_"+fileEras+".root"
    if os.path.isfile(fname): inputFile = TFile.Open( fname );
    else: raise RuntimeError("Could not open data file %s") % fileEras
    return inputFile;


parser = argparse.ArgumentParser()
parser.add_argument("--dataset", default="ZeroBias")
parser.add_argument("--year", default="2017")
parser.add_argument("--ID", default="tight")
parser.add_argument("--TF", default=None, required=True)
parser.add_argument("--fileEras", default="B")
parser.add_argument("--guys", default="A")
parser.add_argument("--extraText", default="_Bonus")
parser.add_argument("--etaOrIndex", default="Eta")
parser.add_argument("--model", default="1LMedium")
args = parser.parse_args()


print ""
print "==> Start TMVARegression"

gROOT.ProcessLine(".L deltaPhi.h+")

CMSSW_BASE = os.environ["CMSSW_BASE"]
dataloaderName = "TMVARegression_TF"+args.TF+"_Era"+args.fileEras+"_Guys"+args.guys+"_"+args.model+"_"+args.etaOrIndex+args.extraText
outfileName = CMSSW_BASE+"/src/MLShEx_CorrL1Muons/Configuration/test/TMVATraining/"+dataloaderName
outputFile = TFile.Open( outfileName+".root", "RECREATE" )

factory = TMVA.Factory( "TMVARegression", outputFile,"!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" )

dataloader = TMVA.DataLoader(dataloaderName)

dataloader.AddVariable( "L1muon_ptCorr", "p_{T,corrected}(L1 #mu)", "GeV", "F" )
dataloader.AddVariable( "L1muon_eta", "#eta(L1 #mu)", "", "F" )
dataloader.AddVariable( "L1muon_phiAtVtx", "#phiAtVtx(L1 #mu)", "", "F" )
dataloader.AddVariable( "L1muon_charge", "charge(L1 #mu)", "", "I" )

dataloader.AddTarget( "deltaPhi(L1muon_phiAtVtx,recomuon_phi)" )

inputFile = InitializeFile(args.dataset,args.year,args.ID,args.fileEras)
inputTree = inputFile.Get("mytree")
dataloader.AddRegressionTree( inputTree )


mycut = TCut(" recomuon_dr >= 0.0 && recomuon_dr < 0.2 ")
mycut = TFCut(mycut, args.TF, args.etaOrIndex)
mycut = GuysCut(mycut, args.guys)
print mycut


preparationString = "nTrain_Regression=25000:nTest_Regression=25000:SplitMode=Random:NormMode=NumEvents:!V"
dataloader.PrepareTrainingAndTestTree( mycut, preparationString )

# Define model
model = Sequential()
if args.model == "1LMedium":
    model.add(Dense(16, input_dim=4,activation="relu"))
elif args.model == "3LWide":
    model.add(Dense(16,input_dim=4,activation="relu"))
    model.add(Dense(32, activation="relu"))
    model.add(Dense(12, activation="relu"))
else:
    raise RuntimeError("Unknown model")
model.add(Dense(1))

# Set loss and optimizer
model.compile(loss="mse", optimizer="rmsprop", metrics=["mae"])

# Store model to file
model.save(args.model+".h5")
model.summary()


factory.BookMethod( dataloader,  TMVA.Types.kMLP, "MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=300:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.6:SamplingEpoch=0.8:ConvergenceImprove=1e-4:ConvergenceTests=15" )
factory.BookMethod(dataloader, TMVA.Types.kPyKeras,"PyKeras","!V:VarTransform=G:FilenameModel="+args.model+".h5:NumEpochs=20:BatchSize=32")

factory.TrainAllMethods();
factory.TestAllMethods();
factory.EvaluateAllMethods();


outputFile.Close();

print "==> Wrote root file: %s" % outfileName
print "==> TMVARegression is done!"
