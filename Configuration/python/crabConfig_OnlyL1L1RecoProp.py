from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

folder = "L1uGMTAnalyzer_ZeroBias_Run2018B" #Change accordingly
config.General.requestName =  folder
config.General.workArea = folder
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ana.py' #Change accordingly
config.JobType.outputFiles = ['outputL1uGMTAnalyzer.root']

config.Data.inputDataset ='/ZeroBias/Run2018B-17Sep2018-v1/MINIAOD' #Change accordingly

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 35
config.Data.lumiMask = '/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/python/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt' # Lumi mask for 2018
config.Data.outLFNDirBase = '/store/group/cmst3/user/evourlio/'+folder # Change accordingly
config.Data.publication = False

config.Site.storageSite = 'T2_CH_CERN'
