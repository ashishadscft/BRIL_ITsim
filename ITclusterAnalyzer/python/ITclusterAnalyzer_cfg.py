# imports
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


# create a new CMS process
process = cms.Process("ITclusterAnalyzer")

# set up the options
options = VarParsing.VarParsing('analysis')
#set up the defaults
# options.inputFiles = 'file:/afs/cern.ch/user/g/gauzinge/ITsim/myPU_35sample/step3.root'
options.inputFiles = 'file:/afs/cern.ch/work/g/gauzinge/public/PU100/step3_100.0.root'
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.1.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.2.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.3.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.4.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.5.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.6.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.7.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.8.root',
# 'file:/afs/cern.ch/work/g/gauzinge/public/condorout/step3_100.9.root'
options.tag='0'
options.outputFile='summary.root'
options.maxEvents = -1 #all events

#get and parse command line arguments
options.parseArguments()

#not sure these are needed
# Import all the necessary files
# process.load('Configuration.StandardSequences.Services_cff')
# process.load('Configuration.EventContent.EventContent_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_cff')
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# load the geomtry that i modified
process.load('Configuration.Geometry.GeometryExtended2023D21Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('ITclusterAnalyzer')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# the input file
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(options.inputFiles)
                            )

# the config of my analyzer
process.BRIL_IT_Analysis = cms.EDAnalyzer('ITclusterAnalyzer',
                                         clusters=cms.InputTag("siPixelClusters"),
                                         maxBin=cms.untracked.uint32(5000),
                                         docoincidence=cms.untracked.bool(True),
                                         dx_cut=cms.double(.5),
                                         dy_cut=cms.double(.5),
                                         dz_cut=cms.double(0.9)
                                         )

# the TFIleService that produces the output root files
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile)
                                   )


process.p = cms.Path(process.BRIL_IT_Analysis)