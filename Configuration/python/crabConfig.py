from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#folder = "L1uGMTAnalyzer_ZeroBias_Run2017B_2018_12_17"
#folder = "L1uGMTAnalyzer_Charmonium_Run2017B_2018_12_17"
#folder = "L1uGMTAnalyzer_ZeroBias_Run2018B_2019_02_06"
folder = "L1uGMTAnalyzer_MuOnia_Run2018C_2019_02_26"
config.General.requestName =  folder
config.General.workArea = folder
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ana.py'
config.JobType.outputFiles = ['outputL1uGMTAnalyzer.root']

#config.Data.inputDataset ='/ZeroBias/Run2017B-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset ='/ZeroBias/Run2017C-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset ='/ZeroBias/Run2017D-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset ='/ZeroBias/Run2017E-17Nov2017-v1/MINIAOD'
#config.Data.inputDataset ='/ZeroBias/Run2017F-17Nov2017-v1/MINIAOD'

#config.Data.inputDataset ='/ZeroBias/Run2018A-17Sep2018-v1/MINIAOD'
#config.Data.inputDataset ='/ZeroBias/Run2018B-17Sep2018-v1/MINIAOD'
#config.Data.inputDataset ='/ZeroBias/Run2018C-17Sep2018-v1/MINIAOD'

#config.Data.inputDataset ='/Charmonium/Run2017B-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset ='/Charmonium/Run2017C-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset ='/Charmonium/Run2017D-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset ='/Charmonium/Run2017E-31Mar2018-v1/MINIAOD'
#config.Data.inputDataset ='/Charmonium/Run2017F-31Mar2018-v1/MINIAOD'

#config.Data.inputDataset ='/MuOnia/Run2018A-17Sep2018-v1/MINIAOD'
#config.Data.inputDataset ='/MuOnia/Run2018B-17Sep2018-v1/MINIAOD'
config.Data.inputDataset ='/MuOnia/Run2018C-17Sep2018-v1/MINIAOD'

config.Data.inputDBS = 'global'
#config.Data.userInputFiles=[]
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 35
#config.Data.totalUnits= 20
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt' # Lumi mask for 2017
config.Data.lumiMask = '/afs/cern.ch/work/e/evourlio/private/L1uGMTAnalyzer_v2/CMSSW_10_2_11/src/L1uGMTAnalyzer/Configuration/python/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt' # Lumi mask for 2018
config.Data.outLFNDirBase = '/store/group/cmst3/user/evourlio/'+folder # Change accordingly
config.Data.publication = False

config.Site.storageSite = 'T2_CH_CERN'
