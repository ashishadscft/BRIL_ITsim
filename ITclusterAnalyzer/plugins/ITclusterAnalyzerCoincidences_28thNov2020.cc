// -*- C++ -*-
//
// Package:    BRIL_ITsim/ITclusterAnalyzer
// Class:      ITclusterAnalyzerCoincidences
//
/**\class  ITclusterAnalyzerCoincidences.cc BRIL_ITsim/ITclusterAnalyzer/plugins/ITclusterAnalyzerCoincidences.cc
   Description: [one line class summary]
   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Ashish Sehrawat
//         Created:  Tue, 28 Nov 2020 13:24:06 GMT
//
//

// system include files
#include <algorithm>
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "DataFormats/Phase2TrackerDigi/interface/Phase2TrackerDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TStyle.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class ITclusterAnalyzerCoincidences : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ITclusterAnalyzerCoincidences(const edm::ParameterSet&);
  ~ITclusterAnalyzerCoincidences();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //bool findCoincidence(DetId, Global3DPoint, bool);
  const SiPixelCluster* findCoincidence2x(DetId, Global3DPoint, unsigned int&, edmNew::DetSet<SiPixelCluster>::const_iterator, unsigned int);
  
  const SiPixelCluster* findCoincidenceInR2x(DetId, Global3DPoint, edmNew::DetSet<SiPixelCluster>::const_iterator, unsigned int&);

  uint32_t get2xRModuleMap(unsigned int, unsigned int, unsigned int, int, unsigned int&, unsigned int&, unsigned int);


  edm::DetSetVector<PixelDigiSimLink>::const_iterator findSimLinkDetSet(unsigned int thedetid);
  std::set<unsigned int> getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator, edmNew::DetSet<SiPixelCluster>::const_iterator, bool print);
  bool areSameSimTrackId(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>&);
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> m_tokenClusters;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> m_tokenSimLinks;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> m_tokenDigis;
  
  // the pointers to geometry, topology and clusters
  // these are members so all functions can access them without passing as argument
  const TrackerTopology* tTopo = NULL;
  const TrackerGeometry* tkGeom = NULL;
  const edmNew::DetSetVector<SiPixelCluster>* clusters = NULL;
  const edm::DetSetVector<PixelDigiSimLink>* simlinks = NULL;
  const edm::DetSetVector<PixelDigi>* digis = NULL;  //defining pointer to digis - COB 26.02.19
  
  //max bins of Counting histogram
  uint32_t m_maxBin;
  //flag for checking coincidences
  bool m_docoincidence;
  
  //array of TH2F for clusters per disk per ring
  TH2F* m_diskHistosCluster[2][4];
  
  //tracker maps for clusters
  TH2F* m_trackerLayoutClustersZR;
  TH2F* m_trackerLayoutClustersYX;
  
  //array of TH2F for 2xcoinc per disk per ring
  //first all coincidences
  TH2F* m_diskHistos2x[2][4];
  TH2F* m_diskHistos2x_InR[2][4];

  //and the real ones
  TH2F* m_diskHistos2xreal[2][4];
  TH2F* m_diskHistos2xreal_InR[2][4];
  
  //tracker maps for 2xcoinc
  TH2F* m_trackerLayout2xZR;
  TH2F* m_trackerLayout2xYX;
  
  TH2F* m_trackerLayout2xZR_InR;
  TH2F* m_trackerLayout2xYX_InR;

  
  //simple residual histograms for the cuts

  TH1F* m_dX[2][4][5];
  TH1F* m_dY[2][4][5];
  TH1F* m_dR[2][4][5];
  TH1F* m_dr[2][4][5];
  TH1F* m_deltaphi[2][4][5];
  
  TH1F* m_dX_sametrack[2][4][5];
  TH1F* m_dY_sametrack[2][4][5];
  TH1F* m_dR_sametrack[2][4][5];
  TH1F* m_dr_sametrack[2][4][5];
  TH1F* m_deltaphi_sametrack[2][4][5];
  
  TH1F* m_dX_notsametrack[2][4][5];
  TH1F* m_dY_notsametrack[2][4][5];
  TH1F* m_dR_notsametrack[2][4][5];
  TH1F* m_dr_notsametrack[2][4][5];
  TH1F* m_deltaphi_notsametrack[2][4][5];
  
  TH1F* m_deltaphi_fabs[2][4][5];
  TH1F* m_deltaphi_fabs_sametrack[2][4][5];
  TH1F* m_deltaphi_fabs_notsametrack[2][4][5];




  TH1F* m_dX_InR[2][4][5];
  TH1F* m_dY_InR[2][4][5];
  TH1F* m_dR_InR[2][4][5];
  TH1F* m_dr_InR[2][4][5];
  TH1F* m_deltaphi_InR[2][4][5];

  TH1F* m_dX_sametrack_InR[2][4][5];
  TH1F* m_dY_sametrack_InR[2][4][5];
  TH1F* m_dR_sametrack_InR[2][4][5];
  TH1F* m_dr_sametrack_InR[2][4][5];
  TH1F* m_deltaphi_sametrack_InR[2][4][5];

  TH1F* m_dX_notsametrack_InR[2][4][5];
  TH1F* m_dY_notsametrack_InR[2][4][5];
  TH1F* m_dR_notsametrack_InR[2][4][5];
  TH1F* m_dr_notsametrack_InR[2][4][5];
  TH1F* m_deltaphi_notsametrack_InR[2][4][5];

  TH1F* m_deltaphi_fabs_InR[2][4][5];
  TH1F* m_deltaphi_fabs_sametrack_InR[2][4][5];
  TH1F* m_deltaphi_fabs_notsametrack_InR[2][4][5];


  
  //the number of clusters per module
  TH1F* m_nClusters;
  
  //cuts for the coincidence
  double m_dx;
  double m_dy;
  double m_dz;
  
  const float C1 = 2;
  const float C2 = 2;

  
  float m_dr_cuts[2][4][5] = {{
                               {0.00633554, 0.00709231, 0.00757662, 0.008659, 0.00872148},           //Disk -1 
                               {0.00594206, 0.00663493, 0.00718173, 0.00766371, 0.00797059},           //Disk -2
                               {0.00594206, 0.00663493, 0.00660603, 0.00713035, 0.00679445},           //Disk -3 
                               {0.00530235, 0.00587692, 0.00616465, 0.00661972, 0.00660599}},          //Disk -4
			       {{0.00635769, 0.0071538, 0.00780034, 0.00863764, 0.00883121},           //Disk 1
			       {0.00587583, 0.00656126, 0.00700576, 0.00773345, 0.00786885},           //Disk 2 
			       {0.00557387, 0.00617934, 0.00664346, 0.0070085, 0.00727192},           //Disk 3
			       {0.0053299, 0.00586655, 0.0061499, 0.00667791, 0.00678366}}};         //Disk 4
  
  float m_dphi_cuts[2][4][5] = {{
                                {0.000856727,0.00128078,0.00151297,0.00151297,0.0016561},         //Disk -1
                                {0.000765842,0.00110698,0.00131104,0.00143555,0.00150265},          //Disk -2
                                {0.000684645,0.000975385,0.00115899,0.00126858,0.00131598},          //Disk -3         
                                {0.000598839,0.000862728,0.00104327,0.00112621,0.00115609}},          //Disk -4
                              
                                {{0.000877912,0.00125898,0.0015378,0.00164165,0.00173015},          //Disk 1
				{0.00076339,0.00110755,0.00132174,0.00144755,0.00150547},          //Disk 2
				{0.000679543,0.000972429,0.00116735,0.00125955,0.00132885},          //Disk 3 
				{0.000609859,0.000881345,0.00101584,0.00112554,0.00118477}}};        //Disk 4
  
  
  
  float m_dr_cuts_offset[2][4][5] = {{
                                     {0.0342041, 0.0511416, 0.0676943, 0.0834355, 0.0984932},          //Disk -1                             
				     {0.0298276, 0.0445398, 0.0589828, 0.0727454, 0.0853587},          //Disk -2                             
                                     {0.0260326, 0.0386608, 0.05122, 0.0632478, 0.0746467},            //Disk -3                         
				     {0.0227089, 0.0336467, 0.044476, 0.0551065, 0.06484}},            //Disk -4                             
				     {{0.0341565, 0.0511018, 0.0676485, 0.0835835, 0.0984028},          //Disk 1 
				     {0.0298515, 0.0444016, 0.0589681, 0.0727789, 0.0854843},          //Disk 2
				     {0.0260239, 0.0387303, 0.0511824, 0.0633755, 0.0742994},          //Disk 3 
				     {0.022729, 0.0336624, 0.0445526, 0.0550466, 0.0644863}}};          //Disk 4 		


  float m_dphi_cuts_offset[2][4][5] = {{
                                         {0,0,0,0,0},
                                         {0,0,0,0,0},
                                         {0,0,0,0,0},
					 {0,0,0,0,0}},
				         {{0,0,0,0,0},
					 {0,0,0,0,0},
					 {0,0,0,0,0},
					 {0,0,0,0,0}}};


  float m_dr_cuts_InR[2][4][5] = {{
                                     {1,1,1,1,1},
                                     {1,1,1,1,1},
                                     {1,1,1,1,1},
                                    {1,1,1,1,1}},
				    {{1,1,1,1,1},
				    {1,1,1,1,1},
				    {1,1,1,1,1},
				    {1,1,1,1,1}}};


  float m_dphi_cuts_InR[2][4][5] = {{
                                     {1,1,1,1,1},
                                     {1,1,1,1,1},
                                     {1,1,1,1,1},
                                     {1,1,1,1,1}},
				    {{1,1,1,1,1},
				     {1,1,1,1,1},
				     {1,1,1,1,1},
				     {1,1,1,1,1}}};
  




  float m_dr_cuts_offset_InR[2][4][5];
  float m_dphi_cuts_offset_InR[2][4][5];

  
  //event counter
  uint32_t m_nevents;
  
};


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ITclusterAnalyzerCoincidences::ITclusterAnalyzerCoincidences(const edm::ParameterSet& iConfig)
  : //m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>> ("clusters"))
  m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
  , m_tokenSimLinks(consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("simlinks")))
  , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"))) //adding digis variable - COB 26.02.19
  , m_maxBin(iConfig.getUntrackedParameter<uint32_t>("maxBin"))
  , m_docoincidence(iConfig.getUntrackedParameter<bool>("docoincidence"))
  , m_dx(iConfig.getParameter<double>("dx_cut"))
  , m_dy(iConfig.getParameter<double>("dy_cut"))
  , m_dz(iConfig.getParameter<double>("dz_cut"))
  , C1(iConfig.getParameter<double>("C1_cut"))
  , C2(iConfig.getParameter<double>("C2_cut")) {
  
  
  //now do what ever initialization is needed
  m_nevents = 0;
  
}
ITclusterAnalyzerCoincidences::~ITclusterAnalyzerCoincidences() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void ITclusterAnalyzerCoincidences::beginJob() {
  
  edm::Service<TFileService> fs;
  
  fs->file().cd("/");
  TFileDirectory td = fs->mkdir("TEPX");
  
  fs->file().cd("/");
  td = fs->mkdir("TEPX/perModule");
  m_nClusters = td.make<TH1F>("Number of Clusters per module per event", "# of Clusters;# of Clusters; Occurence", 500, 0, 500);
  
  
  fs->file().cd("/");
  td = fs->mkdir("TEPX/Clusters");

  
  //now lets create the histograms
  for (unsigned int k = 0; k < 2; k++) {
    for (unsigned int i = 0; i < 4; i++) {
      
      unsigned int disk = i + 1;
      unsigned int side = k + 1;
      
      std::stringstream histoname;
      histoname << "Numberofclusters_side" << side <<"_Disk" << disk;
      
      //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
      m_diskHistosCluster[k][i] = td.make<TH2F>(histoname.str().c_str()," ", 5, .5, 5.5, m_maxBin, 0, m_maxBin);
      
    }  
  }
  
  m_trackerLayoutClustersZR = td.make<TH2F>("RVsZ", "RVsZposition_Clusters", 6000, -300.0, 300.0, 600, 0.0, 30.0);
  m_trackerLayoutClustersYX = td.make<TH2F>("XVsY", "XVsYposition_Clusters", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
  
  
  
  if (m_docoincidence) {
    fs->file().cd("/");
    td = fs->mkdir("TEPX/2xCoincidences");
    //now lets create the histograms
    
    for (unsigned int k = 0; k < 2; k++) {
      for (unsigned int i = 0; i < 4; i++) {
	
	int disk = i + 1;
	int side = k + 1;
	
	std::stringstream histoname;
	histoname << "Numberof2xCoincidences_Inphi_side" << side << "_Disk" << disk;
	
	//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
	m_diskHistos2x[k][i] = td.make<TH2F>(histoname.str().c_str()," ",  5, .5, 5.5, m_maxBin, 0, m_maxBin);
	histoname.str("");
	
	histoname << "Numberofreal2xCoincidences_Inphi_side" << side << "_Disk" << disk;
	
	//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
	m_diskHistos2xreal[k][i] = td.make<TH2F>(histoname.str().c_str()," ", 5, .5, 5.5, m_maxBin, 0, m_maxBin);
	histoname.str("");


	histoname << "Numberof2xCoincidences_InR_side" << side << "_Disk" << disk;  
   
        //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh                       

        m_diskHistos2x_InR[k][i] = td.make<TH2F>(histoname.str().c_str()," ",  5, .5, 5.5, m_maxBin, 0, m_maxBin);
        histoname.str("");

        histoname << "Numberofreal2xCoincidences_InR_side" << side << "_Disk" << disk;

        //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh                                                                                
 
        m_diskHistos2xreal_InR[k][i] = td.make<TH2F>(histoname.str().c_str()," ", 5, .5, 5.5, m_maxBin, 0, m_maxBin);
        histoname.str("");

	
      }
    }
  
  
    m_trackerLayout2xZR = td.make<TH2F>("2xCoincidencesInphi_RVsZ", "RvsZpositionInphi_2xCoincidences", 6000, -300.0, 300.0, 600, 0.0, 30.0);
    m_trackerLayout2xYX = td.make<TH2F>("2xCoincidencesInphi_XVsY", "XvsYpositionInphi_2xCoincidences", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    

    m_trackerLayout2xZR_InR = td.make<TH2F>("2xCoincidencesInR_RVsZ", "RvsZpositionInR_2xCoincidences", 6000, -300.0, 300.0, 600, 0.0, 30.0);
    m_trackerLayout2xYX_InR = td.make<TH2F>("2xCoincidencesInR_XVsY", "XvsYpositionInR_2xCoincidences", 1000, -50.0, 50.0, 1000, -50.0, 50.0);


    
    for (unsigned int k = 0; k < 2; k++) {
      for (unsigned int i = 0; i < 4; i++) {
	for (unsigned int j = 0; j < 5; j++) {
	  
	  unsigned int disk = i + 1;
	  unsigned int ring = j + 1;
	  unsigned int side = k + 1;
	  
	  std::stringstream histoname;
	  histoname << "m_dX_side" << side <<"_Disk" << disk  << "_Ring" << ring;
	  m_dX[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dY_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dY[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dR_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dr_side" << side << "_Disk" << disk <<"_Ring" << ring;
	  m_dr[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dphi_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_deltaphi[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.04, 0.04);
	  histoname.str("");
	  
	  
	  histoname << "m_dphi_abs_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_deltaphi_fabs[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.04);
	  histoname.str("");
	  
	  
	  histoname << "m_dX_sametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dX_sametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dY_sametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dY_sametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  histoname << "m_dR_sametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dR_sametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, 0, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dr_sametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dr_sametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dphi_sametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_deltaphi_sametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.04, 0.04);
	  histoname.str("");
	  
	  
	  histoname << "m_dphi_abs_sametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_deltaphi_fabs_sametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.04);
	  histoname.str("");
	  
	  
	  histoname << "m_dX_notsametrack_side" << side <<"_Disk" << disk << "_Ring" << ring;
	  m_dX_notsametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dY_notsametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dY_notsametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dR_notsametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dR_notsametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, 0, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dr_notsametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_dr_notsametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
	  histoname.str("");
	  
	  
	  histoname << "m_dphi_notsametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
	  m_deltaphi_notsametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.04, 0.04);
	  histoname.str("");
	  
	  
	  histoname << "m_dphi_abs_notsametrack_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_deltaphi_fabs_notsametrack[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.04);
          histoname.str("");
	  
	  




	  //Cut histograms for 2xinR	
	  histoname << "m_dX_InR_side" << side <<"_Disk" << disk  << "_Ring" << ring;
          m_dX_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
          histoname.str("");


          histoname << "m_dY_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dY_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.15, 0.15);
          histoname.str("");


          histoname << "m_dR_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dR_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.15);
          histoname.str("");

	  histoname << "m_dr_InR_side" << side << "_Disk" << disk <<"_Ring" << ring;
          m_dr_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
          histoname.str("");


          histoname << "m_dphi_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_deltaphi_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.1, 0.1);
          histoname.str("");


          histoname << "m_dphi_abs_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_deltaphi_fabs_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.1);
          histoname.str("");


          histoname << "m_dX_sametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dX_sametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
          histoname.str("");


          histoname << "m_dY_sametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dY_sametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.15, 0.15);
          histoname.str("");

          histoname << "m_dR_sametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dR_sametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, 0, 0.15);
          histoname.str("");


	  histoname << "m_dr_sametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dr_sametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
          histoname.str("");


          histoname << "m_dphi_sametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_deltaphi_sametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.1, 0.1);
          histoname.str("");


          histoname << "m_dphi_abs_sametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_deltaphi_fabs_sametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.04);
          histoname.str("");


          histoname << "m_dX_notsametrack_InR_side" << side <<"_Disk" << disk << "_Ring" << ring;
          m_dX_notsametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.15, 0.15);
          histoname.str("");


          histoname << "m_dY_notsametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dY_notsametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
          histoname.str("");


	  histoname << "m_dR_notsametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dR_notsametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, 0, 0.15);
          histoname.str("");


          histoname << "m_dr_notsametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_dr_notsametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, -0.15, 0.15);
          histoname.str("");


          histoname << "m_dphi_notsametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_deltaphi_notsametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(), "", 200, -0.1, 0.1);
          histoname.str("");


          histoname << "m_dphi_abs_notsametrack_InR_side" << side << "_Disk" << disk << "_Ring" << ring;
          m_deltaphi_fabs_notsametrack_InR[k][i][j] = td.make<TH1F>(histoname.str().c_str(),"", 200, 0, 0.1);
          histoname.str("");

	  
	}
      }
    }
  }
}




// ------------ method called for each event  ------------
void ITclusterAnalyzerCoincidences::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //get the digis - COB 26.02.19
  edm::Handle<edm::DetSetVector<PixelDigi>> tdigis;
  iEvent.getByToken(m_tokenDigis, tdigis);
  
  //get the clusters
  edm::Handle<edmNew::DetSetVector<SiPixelCluster>> tclusters;
  iEvent.getByToken(m_tokenClusters, tclusters);
  
  //get the simlinks
  edm::Handle<edm::DetSetVector<PixelDigiSimLink>> tsimlinks;
  iEvent.getByToken(m_tokenSimLinks, tsimlinks);
  
  // Get the geometry
  edm::ESHandle<TrackerGeometry> tgeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get("idealForDigi", tgeomHandle);
  
  // Get the topology
  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  
  //get the pointers to geometry, topology and clusters
  tTopo = tTopoHandle.product();
  //const TrackerGeometry* tkGeom = &(*geomHandle);
  tkGeom = tgeomHandle.product();
  clusters = tclusters.product();
  simlinks = tsimlinks.product();
  digis = tdigis.product();  //pointer to digis - COB 26.02.19
  
  //a 2D counter array to count the number of clusters per disk and per ring
  unsigned int cluCounter[2][4][5];
  memset(cluCounter, 0, sizeof(cluCounter));
  
  //counter for 2x coincidences
  unsigned int x2Counter[2][4][5];
  memset(x2Counter, 0, sizeof(x2Counter));
  
  unsigned int x2Counterreal[2][4][5];
  memset(x2Counterreal, 0, sizeof(x2Counterreal));
  
  unsigned int x2CounterInR[2][4][5];
  memset(x2CounterInR, 0, sizeof(x2CounterInR));

  unsigned int x2CounterrealInR[2][4][5];
  memset(x2CounterrealInR, 0, sizeof(x2CounterrealInR));


  
  //-------------------------------------------------------------
  
  //loop the modules in the cluster collection
  for (typename edmNew::DetSetVector<SiPixelCluster>::const_iterator DSVit = clusters->begin(); DSVit != clusters->end(); DSVit++) {
    
    //get the detid
    unsigned int rawid(DSVit->detId());
    DetId detId(rawid);
    
    //figure out the module type using the detID
    TrackerGeometry::ModuleType mType = tkGeom->getDetectorType(detId);
    if (mType != TrackerGeometry::ModuleType::Ph2PXF && detId.subdetId() != PixelSubdetector::PixelEndcap)
      continue;
    
    
    //find out which layer, side and ring
    unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
    unsigned int layer = (tTopo->pxfDisk(detId)); //values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
    unsigned int ring = (tTopo->pxfBlade(detId));
    unsigned int module = (tTopo->pxfModule(detId));    
    
    
    if (layer > 8) { // TEPX modules
      
      //the index in my histogram map
      int disk = 1;
      
      
      if (side == 1 || side == 2) {
	disk = layer - 8;
      }
      
      // Get the geomdet
      const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
      if (!geomDetUnit)
	continue;
      
      MeasurementPoint origin(0, 0);
      Local3DPoint localPosClu_origin = geomDetUnit->topology().localPosition(origin);
      Global3DPoint globalPosClu_origin = geomDetUnit->surface().toGlobal(localPosClu_origin);
      float phiangle = globalPosClu_origin.phi().value();
      float rvalue = globalPosClu_origin.perp();
      
      //std::cout << "      " << side << "      " << disk << "      " << ring << "      " << module << "      " << globalPosClu_origin.z() << std::endl;
      
      int disk_layer_type = -1;
      
      if(side == 1 && 1<= ring && ring <= 4 && module % 2 != 0) disk_layer_type = 1;
      if(side == 1 && 1<= ring && ring <= 4 && module % 2 == 0) disk_layer_type = 2;
      if(side == 1 && ring == 5 && module % 2 != 0) disk_layer_type = 2;
      if(side == 1 && ring == 5 && module % 2 == 0) disk_layer_type = 1;
      
      if(side == 2 && 1<= ring && ring <= 4 && module % 2 == 0) disk_layer_type = 1;
      if(side == 2 && 1<= ring && ring <= 4 && module % 2 != 0) disk_layer_type = 2;
      if(side == 2 && ring == 5 && module % 2 != 0) disk_layer_type = 1;
      if(side == 2 && ring == 5 && module % 2 == 0) disk_layer_type = 2;
      
      
      //std::cout << disk_layer_type << std::endl;
      
      //std::cout << " " << side << " " << globalPosClu_origin.z() << " " << disk << " " << ring << " " << module << " " << rawid << " " << ((rawid >> 2) & 0xFF) << " " << rvalue << " " << phiangle << " " << std::endl;
      
      //std::cout << " rawid " << rawid << " ((rawid >> 2) & 0xFF) " << ((rawid >> 2) & 0xFF) << " module " << module << " side " << side << 
      //  " disk " << disk << " ring " << ring << " r " << rvalue << " phi " << phiangle << " z " << globalPosClu_origin.z() << std::endl;
      
      //std::cout << " rawid " << rawid << " ((rawid >> 2) & 0xFF) " << ((rawid >> 2) & 0xFF) << " module " << module << " side " << side << " disk " << disk << " ring " << ring << " r " << rvalue << " phi " << phiangle << " z " << globalPosClu1.z() << std::endl;    
      
      unsigned int nClu = 0;
      
      //fill the number of clusters for this module
      m_nClusters->Fill(DSVit->size());
      
      
      //now loop the clusters for each detector
      for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit1 = DSVit->begin(); cluit1 != DSVit->end(); cluit1++) {
	
	
	//increment the counters
	nClu++;
	cluCounter[side-1][disk-1][ring-1]++;
	
	// determine the position
	MeasurementPoint mpClu1(cluit1->x(), cluit1->y());
	Local3DPoint localPosClu1 = geomDetUnit->topology().localPosition(mpClu1);
	Global3DPoint globalPosClu1 = geomDetUnit->surface().toGlobal(localPosClu1);
	
	
	m_trackerLayoutClustersZR->Fill(globalPosClu1.z(), globalPosClu1.perp());
	m_trackerLayoutClustersYX->Fill(globalPosClu1.x(), globalPosClu1.y());
	
	
	
	
	if (m_docoincidence) {
	  if (disk_layer_type == 1) {
	    
	    //MeasurementPoint mpClu1(cluit1->x(), cluit1->y());
	    //Local3DPoint localPosClu1 = geomDetUnit->topology().localPosition(mpClu1);
	    //Global3DPoint globalPosClu1 = geomDetUnit->surface().toGlobal(localPosClu1);
	    
	    unsigned int coincidenceId;
	    
	    const SiPixelCluster* found2xcoincidencecluster = this->findCoincidence2x(detId, globalPosClu1, coincidenceId, cluit1, 1);
	    
	    if (!found2xcoincidencecluster) 
	      found2xcoincidencecluster = this->findCoincidence2x(detId, globalPosClu1, coincidenceId, cluit1, 2);
	    
	    if (found2xcoincidencecluster) {
	      
	      x2Counter[side-1][disk-1][ring-1]++;
	      
	      
	      m_trackerLayout2xZR->Fill(globalPosClu1.z(), globalPosClu1.perp());
	      m_trackerLayout2xYX->Fill(globalPosClu1.x(), globalPosClu1.y());
	      
	      
	    }
	  }
	  
	  
	  
	  unsigned int coincidenceIdInR;
	  
	  const SiPixelCluster* found2xcoincidenceclusterInR = this->findCoincidenceInR2x(detId, globalPosClu1, cluit1, coincidenceIdInR);
	  
	  if (found2xcoincidenceclusterInR) {
	    
	    x2CounterInR[side-1][disk-1][ring-1]++;
	    
	    
	    m_trackerLayout2xZR_InR->Fill(globalPosClu1.z(), globalPosClu1.perp());
	    m_trackerLayout2xYX_InR->Fill(globalPosClu1.x(), globalPosClu1.y());
	    
	    
	  }
	}
      }
    }
  }
  
  
  //-----------------------------------------     
  // End of cluster loop
  //end of module loop
  
  
  //ok, now I know the number of clusters/hits per ring per disk and should fill the histogram once for this event
  for (unsigned int k = 0; k < 2; k++) {
    //loop over sides
    for (unsigned int i = 0; i < 4; i++) {  // TEPX
      //loop the disks
      for (unsigned int j = 0; j < 5; j++) {
	//and the rings
	m_diskHistosCluster[k][i]->Fill(j + 1, cluCounter[k][i][j]);
	if (m_docoincidence) {
	  m_diskHistos2x[k][i]->Fill(j + 1, x2Counter[k][i][j]);
	  m_diskHistos2xreal[k][i]->Fill(j + 1, x2Counterreal[k][i][j]);
	  
	  m_diskHistos2x_InR[k][i]->Fill(j + 1, x2CounterInR[k][i][j]);
          m_diskHistos2xreal_InR[k][i]->Fill(j + 1, x2CounterrealInR[k][i][j]);
	  
	  
	}
      }
    }
  }
  
  m_nevents++;
  
}


// ------------ method called once each job just after ending the event loop  ------------
void ITclusterAnalyzerCoincidences::endJob() {
  
  std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ITclusterAnalyzerCoincidences::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
  
  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


const SiPixelCluster* ITclusterAnalyzerCoincidences::findCoincidence2x(DetId thedetid, Global3DPoint globalPosClu1, unsigned int& foundDetId, edmNew::DetSet<SiPixelCluster>::const_iterator cluit1, unsigned int neighbor) {
  
  const SiPixelCluster* found2xcoincidencecluster = NULL; 
  uint32_t rawid = thedetid.rawId();
  uint32_t newmodule = - 1;  
  
  //now I have the raw ID and can mess with the bits
  //the side, layer and ring are the same and I just have to increment or decrement the module number
  
  unsigned int themodule = (tTopo->pxfModule(thedetid));
  unsigned int thering = (tTopo->pxfBlade(thedetid));
  unsigned int thelayer = (tTopo->pxfDisk(thedetid));
  unsigned int theside = (tTopo->pxfSide(thedetid));  
  unsigned int thedisk = thelayer - 8;
  
  
  if (neighbor == 1) { 
    
    newmodule = themodule + 1;
    
    if (thering == 1 && themodule == 20)
      newmodule = 1;
    else if (thering == 2 && themodule == 28)
      newmodule = 1;
    else if (thering == 3 && themodule == 36)
      newmodule = 1;
    else if (thering == 4 && themodule == 44)
      newmodule = 1;
    else if (thering == 5 && themodule == 48)
      newmodule = 1;
    
  }
  
  
  if (neighbor == 2) {
    
    newmodule = themodule - 1;
    
    if (thering == 1 && themodule == 1)
      newmodule = 20;
    if (thering == 2 && themodule == 1)
      newmodule = 28;
    if (thering == 3 && themodule == 1)
      newmodule = 36;
    if (thering == 4 && themodule == 1)
      newmodule = 44;
    if (thering == 5 && themodule == 1)
      newmodule = 48;
    
  }
  
  
  uint32_t newid = (rawid & 0xFFFFFC03) | ((newmodule & 0xFF) << 2);
  
  DetId id(newid);
  
  edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(id);
  
  if (theit == clusters->end()) {
    return found2xcoincidencecluster;
  }
  
  const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(id));
  
  
  //at the end of the day, need to find the closest coincidence hit, so store the minimum 2D distance in a temporary variable and a vector for all values
  double R_min = 1000.;
  
  
  for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit2 = theit->begin(); cluit2 != theit->end(); cluit2++) {
    
    
    // determine the position
    MeasurementPoint mpClu2(cluit2->x(), cluit2->y());
    Local3DPoint localPosClu2 = geomDetUnit->topology().localPosition(mpClu2);
    Global3DPoint globalPosClu2 = geomDetUnit->surface().toGlobal(localPosClu2);    
    
    double r1 = sqrt(pow(globalPosClu1.x(), 2) + pow(globalPosClu1.y(), 2));
    double r2 = sqrt(pow(globalPosClu2.x(), 2) + pow(globalPosClu2.y(), 2));
    
    double phi1 = TMath::ATan2(globalPosClu1.y(), globalPosClu1.x());
    double phi2 = TMath::ATan2(globalPosClu2.y(), globalPosClu2.x());
    
    double dr = r2-r1;
    double dphi = phi2-phi1;
    
    if(fabs(dphi-m_dphi_cuts_offset[theside-1][thedisk-1][thering-1]) < C1*m_dphi_cuts[theside-1][thedisk-1][thering-1] && 
       fabs(dr-m_dr_cuts_offset[theside-1][thedisk-1][thering-1]) < C2*m_dr_cuts[theside-1][thedisk-1][thering-1] && 
       fabs(globalPosClu1.z() - globalPosClu2.z()) < m_dz) {
      
      double dX = - globalPosClu1.x() + globalPosClu2.x();
      double dY = - globalPosClu1.y() + globalPosClu2.y();
      double dR = sqrt(pow(dr-m_dr_cuts_offset[theside-1][thedisk-1][thering-1], 2) + pow((r2+r1)/2,2) * pow(dphi,2));
      
      if (dR < R_min) {
	
	R_min = dR;
	foundDetId = newid;
	found2xcoincidencecluster = cluit2; 
	
      }
      
      m_dX[theside-1][thedisk-1][thering-1] -> Fill(dX);
      m_dY[theside-1][thedisk-1][thering-1] -> Fill(dY);
      m_dR[theside-1][thedisk-1][thering-1] -> Fill(dR);
      m_dr[theside-1][thedisk-1][thering-1] -> Fill(dr);
      m_deltaphi[theside-1][thedisk-1][thering-1] -> Fill(phi2-phi1);
      m_deltaphi_fabs[theside-1][thedisk-1][thering-1] -> Fill(fabs(phi2-phi1));
      
      
      edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = findSimLinkDetSet(rawid);                                       
      std::set<unsigned int> simTrackId = this->getSimTrackId(simLinkDSViter, cluit1, false);                                              
      
      //now get the simlink detset based on the coincidence hit detid
      simLinkDSViter = findSimLinkDetSet(newid);                                                                                           
      std::set<unsigned int> coincidencesimTrackId = this->getSimTrackId(simLinkDSViter, cluit2, false);                
      std::set<unsigned int> intersection;                                                                                                 
      bool areSame = areSameSimTrackId(simTrackId, coincidencesimTrackId, intersection);   
      
      
      if(areSame) {
	
	
	m_dX_sametrack[theside-1][thedisk-1][thering-1] -> Fill(dX);
	m_dY_sametrack[theside-1][thedisk-1][thering-1] -> Fill(dY);
	m_dR_sametrack[theside-1][thedisk-1][thering-1] -> Fill(dR);
	m_dr_sametrack[theside-1][thedisk-1][thering-1] -> Fill(dr);
	m_deltaphi_sametrack[theside-1][thedisk-1][thering-1] -> Fill(phi2-phi1);
	m_deltaphi_fabs_sametrack[theside-1][thedisk-1][thering-1] -> Fill(fabs(phi2-phi1));
	
	
      }
      
      
      else if(!areSame) {
	
	
	m_dX_notsametrack[theside-1][thedisk-1][thering-1] -> Fill(dX);
	m_dY_notsametrack[theside-1][thedisk-1][thering-1] -> Fill(dY);
	m_dR_notsametrack[theside-1][thedisk-1][thering-1] -> Fill(dR);
	m_dr_notsametrack[theside-1][thedisk-1][thering-1] -> Fill(dr);
	m_deltaphi_notsametrack[theside-1][thedisk-1][thering-1] -> Fill(phi2-phi1);
	m_deltaphi_fabs_notsametrack[theside-1][thedisk-1][thering-1] -> Fill(fabs(phi2-phi1));
	
      }    
    }
  }  
  
  
  return found2xcoincidencecluster;
  
}





//2x Coincidence in R function
const SiPixelCluster* ITclusterAnalyzerCoincidences::findCoincidenceInR2x(DetId thedetid, Global3DPoint globalPosClu1, edmNew::DetSet<SiPixelCluster>::const_iterator cluit1, unsigned int& foundDetId) {
  
  
  const SiPixelCluster* found2xcoincidenceclusterInR = NULL;
  
  unsigned int thering = (tTopo->pxfBlade(thedetid));
  unsigned int thelayer = (tTopo->pxfDisk(thedetid));
  unsigned int theside = (tTopo->pxfSide(thedetid));
  unsigned int themodule = (tTopo->pxfModule(thedetid));  
  unsigned int thedisk = thelayer - 8;
  
  unsigned int mapRing;
  unsigned int mapModule; 
  int type = -1;
  
 
  uint32_t ovModid = get2xRModuleMap(thedisk, thering, themodule, type, mapRing, mapModule, theside); 
  
  const GeomDetUnit* geomDetUnit_tmp(tkGeom->idToDetUnit(ovModid));
  
  
  //check if there are clusters in this new module
  edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(ovModid);
  

  //if (theit == clusters->end()) {
  //ovModid = (ovModid & 0xFFFFFC03) | (((mapModule) & 0xFF) << 2);
  //continue;
  // }
  
  
  double R_min = 1000.;
  
  for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit2 = theit->begin(); cluit2 != theit->end(); cluit2++) {
    
    
    // determine the position
    MeasurementPoint mpClu2(cluit2->x(), cluit2->y());
    Local3DPoint localPosClu2 = geomDetUnit_tmp->topology().localPosition(mpClu2);
    Global3DPoint globalPosClu2 = geomDetUnit_tmp->surface().toGlobal(localPosClu2);    
    
    double r1 = sqrt(pow(globalPosClu1.x(), 2) + pow(globalPosClu1.y(), 2));
    double r2 = sqrt(pow(globalPosClu2.x(), 2) + pow(globalPosClu2.y(), 2));
    
    double phi1 = TMath::ATan2(globalPosClu1.y(), globalPosClu1.x());
    double phi2 = TMath::ATan2(globalPosClu2.y(), globalPosClu2.x());
    
    double dr = r2-r1;
    double dphi = phi2-phi1;
    
    
    if (fabs(dphi) < C1 && fabs(dr) < C2 && fabs(globalPosClu1.z() - globalPosClu2.z()) < m_dz) {
      
    
      //if(fabs(dphi) < C1*m_dphi_cuts_InR[theside-1][thedisk-1][thering-1] && 
      //fabs(dr) < C2*m_dr_cuts_InR[theside-1][thedisk-1][thering-1] && 
      //fabs(globalPosClu1.z() - globalPosClu2.z()) < m_dz) {
      
      
      double dX = - globalPosClu1.x() + globalPosClu2.x();
      double dY = - globalPosClu1.y() + globalPosClu2.y();
      double dR = sqrt(pow(dr-m_dr_cuts_offset[theside-1][thedisk-1][thering-1], 2) + pow((r2+r1)/2,2) * pow(dphi,2));
      
	if (dR < R_min) {
	  
	  R_min = dR;
	  foundDetId = ovModid; 
	  found2xcoincidenceclusterInR = cluit2; 
	  
	}
        
	m_dX_InR[theside-1][thedisk-1][thering-1] -> Fill(dX);
	m_dY_InR[theside-1][thedisk-1][thering-1] -> Fill(dY);
	m_dR_InR[theside-1][thedisk-1][thering-1] -> Fill(dR);
	m_dr_InR[theside-1][thedisk-1][thering-1] -> Fill(dr);
	m_deltaphi_InR[theside-1][thedisk-1][thering-1] -> Fill(phi2-phi1);
	m_deltaphi_fabs_InR[theside-1][thedisk-1][thering-1] -> Fill(fabs(phi2-phi1));
	
	
	edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = findSimLinkDetSet(thedetid);                                    
	std::set<unsigned int> simTrackId = this->getSimTrackId(simLinkDSViter, cluit1, false);                                              
	
	//now get the simlink detset based on the coincidence hit detid
	simLinkDSViter = findSimLinkDetSet(ovModid);                                                                                         
	std::set<unsigned int> coincidencesimTrackId = this->getSimTrackId(simLinkDSViter, cluit2, false);                
	std::set<unsigned int> intersection;                                                                                                 
	bool areSame = areSameSimTrackId(simTrackId, coincidencesimTrackId, intersection);   
	
	
	if(areSame) {
	  
	  
	  m_dX_sametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dX);
	  m_dY_sametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dY);
	  m_dR_sametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dR);
	  m_dr_sametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dr);
	  m_deltaphi_sametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(phi2-phi1);
	  m_deltaphi_fabs_sametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(fabs(phi2-phi1));
	  
	  
       	}
	
	
	else if(!areSame) {
	  
	  
	  m_dX_notsametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dX);
	  m_dY_notsametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dY);
	  m_dR_notsametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dR);
	  m_dr_notsametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(dr);
	  m_deltaphi_notsametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(phi2-phi1);
	  m_deltaphi_fabs_notsametrack_InR[theside-1][thedisk-1][thering-1] -> Fill(fabs(phi2-phi1));
	  
	}    
    }
  }
  
  
  
  
  
  return found2xcoincidenceclusterInR;
  
}




edm::DetSetVector<PixelDigiSimLink>::const_iterator ITclusterAnalyzerCoincidences::findSimLinkDetSet(unsigned int thedetid) {
  ////basic template
  edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDS = simlinks->find(thedetid);
  return simLinkDS;
}

std::set<unsigned int> ITclusterAnalyzerCoincidences::getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter, edmNew::DetSet<SiPixelCluster>::const_iterator cluster, bool print) {
  int size = cluster->size();
  std::set<unsigned int> simTrackIds;
  
  for (int i = 0; i < size; i++) {
    
    SiPixelCluster::Pixel pix = cluster->pixel(i);
    unsigned int clusterChannel = PixelDigi::pixelToChannel(pix.x, pix.y);
    
    if (simLinkDSViter != simlinks->end()) {
      for (edm::DetSet<PixelDigiSimLink>::const_iterator it = simLinkDSViter->data.begin(); it != simLinkDSViter->data.end(); it++) {
	if (clusterChannel == it->channel()) {
	  simTrackIds.insert(it->SimTrackId());
	  if (print)
	    std::cout << "Channel: " << clusterChannel << " SimTrack ID: " << it->SimTrackId() << std::endl;
	}
      }
    }
  }
  
  
  return simTrackIds;
}


bool ITclusterAnalyzerCoincidences::areSameSimTrackId(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>& intersection) {
  //method to check if the sim Track id is present in both sets
  //std::set<unsigned int> intersection;
  std::set_intersection(first.begin(), first.end(), second.begin(), second.end(), std::inserter(intersection, intersection.begin()));
  if (!intersection.size()) {
    //std::cout << "WARNING, these clusters have not been caused by the same SimTrackID" << std::endl;
    return false;
  } else if (intersection.size() == 1) {
    return true;
  } else {
    //std::cout << "WARNING: both clusters caused by multiple tracks!" << std::endl;
    return true;
  }
}





uint32_t ITclusterAnalyzerCoincidences::get2xRModuleMap(unsigned int disk, unsigned int ring, unsigned int mod, int type, unsigned int& mapRing, unsigned int& mapModule, unsigned int side) {
  
  
  //do this using rawid and return rawid of module which contain 2xCoincidenceinR cluster
  
  if (side == 2) {
    if (disk == 1) {
      
      if(ring == 1 && mod == 354685956) { mapRing = 2; mapModule = 354690060;}
      if(ring == 1 && mod == 354685964) { mapRing = 2; mapModule = 354690068;}
      if(ring == 1 && mod == 354685972) { mapRing = 2; mapModule = 354690076;}
      if(ring == 1 && mod == 354685980) { mapRing = 2; mapModule = 354690084;}
      if(ring == 1 && mod == 354685988) { mapRing = 2; mapModule = 354690100;}
      if(ring == 1 && mod == 354685996) { mapRing = 2; mapModule = 354690108;}
      if(ring == 1 && mod == 354686004) { mapRing = 2; mapModule = 354690124;}
      if(ring == 1 && mod == 354686012) { mapRing = 2; mapModule = 354690132;}
      if(ring == 1 && mod == 354686020) { mapRing = 2; mapModule = 354690140;}
      if(ring == 1 && mod == 354686028) { mapRing = 2; mapModule = 354690156;}
      
      if(ring == 1 && mod == 354685960) { mapRing = 2; mapModule = 354690060;}
      if(ring == 1 && mod == 354685968) { mapRing = 2; mapModule = 354690068;}
      if(ring == 1 && mod == 354685976) { mapRing = 2; mapModule = 354690084;}
      if(ring == 1 && mod == 354685984) { mapRing = 2; mapModule = 354690092;}
      if(ring == 1 && mod == 354685992) { mapRing = 2; mapModule = 354690100;}
      if(ring == 1 && mod == 354686000) { mapRing = 2; mapModule = 354690116;}
      if(ring == 1 && mod == 354686008) { mapRing = 2; mapModule = 354690124;}
      if(ring == 1 && mod == 354686016) { mapRing = 2; mapModule = 354690140;}
      if(ring == 1 && mod == 354686024) { mapRing = 2; mapModule = 354690148;}
      if(ring == 1 && mod == 354686032) { mapRing = 2; mapModule = 354690156;}
      
      if(ring == 1 && mod == 354685960) { mapRing = 2; mapModule = 354690056;}
      if(ring == 1 && mod == 354685968) { mapRing = 2; mapModule = 354690072;}
      if(ring == 1 && mod == 354685976) { mapRing = 2; mapModule = 354690080;}
      if(ring == 1 && mod == 354685984) { mapRing = 2; mapModule = 354690096;}
      if(ring == 1 && mod == 354685992) { mapRing = 2; mapModule = 354690104;}
      if(ring == 1 && mod == 354686000) { mapRing = 2; mapModule = 354690112;}
      if(ring == 1 && mod == 354686008) { mapRing = 2; mapModule = 354690128;}
      if(ring == 1 && mod == 354686016) { mapRing = 2; mapModule = 354690136;}
      if(ring == 1 && mod == 354686024) { mapRing = 2; mapModule = 354690144;}
      if(ring == 1 && mod == 354686032) { mapRing = 2; mapModule = 354690160;}
      
      
      
      if(ring == 3 && mod == 354694148) { mapRing = 2; mapModule = 354690052;}
      if(ring == 3 && mod == 354694156) { mapRing = 2; mapModule = 354690060;}
      if(ring == 3 && mod == 354694164) { mapRing = 2; mapModule = 354690068;}
      if(ring == 3 && mod == 354694172) { mapRing = 2; mapModule = 354690068;}
      if(ring == 3 && mod == 354694180) { mapRing = 2; mapModule = 354690076;}
      if(ring == 3 && mod == 354694188) { mapRing = 2; mapModule = 354690084;}
      if(ring == 3 && mod == 354694196) { mapRing = 2; mapModule = 354690092;}
      if(ring == 3 && mod == 354694204) { mapRing = 2; mapModule = 354690092;}
      if(ring == 3 && mod == 354694212) { mapRing = 2; mapModule = 354690100;}
      if(ring == 3 && mod == 354694220) { mapRing = 2; mapModule = 354690108;}
      if(ring == 3 && mod == 354694228) { mapRing = 2; mapModule = 354690116;}
      if(ring == 3 && mod == 354694236) { mapRing = 2; mapModule = 354690124;}
      if(ring == 3 && mod == 354694244) { mapRing = 2; mapModule = 354690124;}
      if(ring == 3 && mod == 354694252) { mapRing = 2; mapModule = 354690132;}
      if(ring == 3 && mod == 354694260) { mapRing = 2; mapModule = 354690140;}
      if(ring == 3 && mod == 354694268) { mapRing = 2; mapModule = 354690148;}
      if(ring == 3 && mod == 354694276) { mapRing = 2; mapModule = 354690148;}
      if(ring == 3 && mod == 354694284) { mapRing = 2; mapModule = 354690156;}
      
      if(ring == 2 && mod == 354690056) { mapRing = 3; mapModule = 354694156;}
      if(ring == 2 && mod == 354690064) { mapRing = 3; mapModule = 354694164;}
      if(ring == 2 && mod == 354690072) { mapRing = 3; mapModule = 354694172;}
      if(ring == 2 && mod == 354690080) { mapRing = 3; mapModule = 354694188;}
      if(ring == 2 && mod == 354690088) { mapRing = 3; mapModule = 354694196;}
      if(ring == 2 && mod == 354690096) { mapRing = 3; mapModule = 354694204;}
      if(ring == 2 && mod == 354690104) { mapRing = 3; mapModule = 354694212;}
      if(ring == 2 && mod == 354690112) { mapRing = 3; mapModule = 354694228;}
      if(ring == 2 && mod == 354690120) { mapRing = 3; mapModule = 354694236;}
      if(ring == 2 && mod == 354690128) { mapRing = 3; mapModule = 354694244;}
      if(ring == 2 && mod == 354690136) { mapRing = 3; mapModule = 354694252;}
      if(ring == 2 && mod == 354690144) { mapRing = 3; mapModule = 354694268;}
      if(ring == 2 && mod == 354690152) { mapRing = 3; mapModule = 354694276;}
      if(ring == 2 && mod == 354690160) { mapRing = 3; mapModule = 354694284;}
      
      if(ring == 3 && mod == 354694152) { mapRing = 2; mapModule = 354690056;}
      if(ring == 3 && mod == 354694160) { mapRing = 2; mapModule = 354690064;}
      if(ring == 3 && mod == 354694168) { mapRing = 2; mapModule = 354690064;}
      if(ring == 3 && mod == 354694176) { mapRing = 2; mapModule = 354690072;}
      if(ring == 3 && mod == 354694184) { mapRing = 2; mapModule = 354690080;}
      if(ring == 3 && mod == 354694192) { mapRing = 2; mapModule = 354690088;}
      if(ring == 3 && mod == 354694200) { mapRing = 2; mapModule = 354690096;}
      if(ring == 3 && mod == 354694208) { mapRing = 2; mapModule = 354690096;}
      if(ring == 3 && mod == 354694216) { mapRing = 2; mapModule = 354690104;}
      if(ring == 3 && mod == 354694224) { mapRing = 2; mapModule = 354690112;}
      if(ring == 3 && mod == 354694232) { mapRing = 2; mapModule = 354690112;}
      if(ring == 3 && mod == 354694240) { mapRing = 2; mapModule = 354690120;}
      if(ring == 3 && mod == 354694248) { mapRing = 2; mapModule = 354690128;}
      if(ring == 3 && mod == 354694256) { mapRing = 2; mapModule = 354690136;}
      if(ring == 3 && mod == 354694264) { mapRing = 2; mapModule = 354690144;}
      if(ring == 3 && mod == 354694272) { mapRing = 2; mapModule = 354690152;}    
      if(ring == 3 && mod == 354694280) { mapRing = 2; mapModule = 354690152;}
      if(ring == 3 && mod == 354694288) { mapRing = 2; mapModule = 354690160;}
      
      
      
      if(ring == 3 && mod == 354694148) { mapRing = 4; mapModule = 354698244;}
      if(ring == 3 && mod == 354694156) { mapRing = 4; mapModule = 354698252;}
      if(ring == 3 && mod == 354694164) { mapRing = 4; mapModule = 354698268;}
      if(ring == 3 && mod == 354694172) { mapRing = 4; mapModule = 354698276;}
      if(ring == 3 && mod == 354694180) { mapRing = 4; mapModule = 354698284;}
      if(ring == 3 && mod == 354694188) { mapRing = 4; mapModule = 354698292;}
      if(ring == 3 && mod == 354694196) { mapRing = 4; mapModule = 354698300;}
      if(ring == 3 && mod == 354694204) { mapRing = 4; mapModule = 354698316;}
      if(ring == 3 && mod == 354694212) { mapRing = 4; mapModule = 354698324;}
      if(ring == 3 && mod == 354694220) { mapRing = 4; mapModule = 354698332;}
      if(ring == 3 && mod == 354694228) { mapRing = 4; mapModule = 354698340;}
      if(ring == 3 && mod == 354694236) { mapRing = 4; mapModule = 354698356;}
      if(ring == 3 && mod == 354694244) { mapRing = 4; mapModule = 354698364;}
      if(ring == 3 && mod == 354694252) { mapRing = 4; mapModule = 354698372;}
      if(ring == 3 && mod == 354694260) { mapRing = 4; mapModule = 354698380;}
      if(ring == 3 && mod == 354694268) { mapRing = 4; mapModule = 354698388;}
      if(ring == 3 && mod == 354694276) { mapRing = 4; mapModule = 354698404;}
      if(ring == 3 && mod == 354694284) { mapRing = 4; mapModule = 354698412;}
      
      if(ring == 3 && mod == 354694152) { mapRing = 4; mapModule = 354698252;}
      if(ring == 3 && mod == 354694160) { mapRing = 4; mapModule = 354698260;}
      if(ring == 3 && mod == 354694168) { mapRing = 4; mapModule = 354698268;}
      if(ring == 3 && mod == 354694176) { mapRing = 4; mapModule = 354698276;}
      if(ring == 3 && mod == 354694184) { mapRing = 4; mapModule = 354698284;}
      if(ring == 3 && mod == 354694192) { mapRing = 4; mapModule = 354698300;}
      if(ring == 3 && mod == 354694200) { mapRing = 4; mapModule = 354698308;}
      if(ring == 3 && mod == 354694208) { mapRing = 4; mapModule = 354698316;}
      if(ring == 3 && mod == 354694216) { mapRing = 4; mapModule = 354698324;}
      if(ring == 3 && mod == 354694224) { mapRing = 4; mapModule = 354698340;}
      if(ring == 3 && mod == 354694232) { mapRing = 4; mapModule = 354698348;}
      if(ring == 3 && mod == 354694240) { mapRing = 4; mapModule = 354698356;}
      if(ring == 3 && mod == 354694248) { mapRing = 4; mapModule = 354698364;}
      if(ring == 3 && mod == 354694256) { mapRing = 4; mapModule = 354698372;}
      if(ring == 3 && mod == 354694264) { mapRing = 4; mapModule = 354698388;}
      if(ring == 3 && mod == 354694272) { mapRing = 4; mapModule = 354698396;}
      if(ring == 3 && mod == 354694280) { mapRing = 4; mapModule = 354698404;}
      if(ring == 3 && mod == 354694288) { mapRing = 4; mapModule = 354698412;}
      
      if(ring == 3 && mod == 354694152) { mapRing = 4; mapModule = 354698248;}
      if(ring == 3 && mod == 354694160) { mapRing = 4; mapModule = 354698256;}
      if(ring == 3 && mod == 354694168) { mapRing = 4; mapModule = 354698272;}
      if(ring == 3 && mod == 354694176) { mapRing = 4; mapModule = 354698280;}
      if(ring == 3 && mod == 354694184) { mapRing = 4; mapModule = 354698288;}
      if(ring == 3 && mod == 354694192) { mapRing = 4; mapModule = 354698296;}
      if(ring == 3 && mod == 354694200) { mapRing = 4; mapModule = 354698304;}
      if(ring == 3 && mod == 354694208) { mapRing = 4; mapModule = 354698320;}
      if(ring == 3 && mod == 354694216) { mapRing = 4; mapModule = 354698328;}
      if(ring == 3 && mod == 354694224) { mapRing = 4; mapModule = 354698336;}
      if(ring == 3 && mod == 354694232) { mapRing = 4; mapModule = 354698344;}
      if(ring == 3 && mod == 354694240) { mapRing = 4; mapModule = 354698360;}
      if(ring == 3 && mod == 354694248) { mapRing = 4; mapModule = 354698368;}
      if(ring == 3 && mod == 354694256) { mapRing = 4; mapModule = 354698376;}
      if(ring == 3 && mod == 354694264) { mapRing = 4; mapModule = 354698384;}
      if(ring == 3 && mod == 354694272) { mapRing = 4; mapModule = 354698392;}
      if(ring == 3 && mod == 354694280) { mapRing = 4; mapModule = 354698408;}
      if(ring == 3 && mod == 354694288) { mapRing = 4; mapModule = 354698416;}
      

      if(ring == 5 && mod == 354702340) { mapRing = 4; mapModule = 354698416;}
      if(ring == 5 && mod == 354702348) { mapRing = 4; mapModule = 354698248;}
      if(ring == 5 && mod == 354702356) { mapRing = 4; mapModule = 354698256;}
      if(ring == 5 && mod == 354702364) { mapRing = 4; mapModule = 354698264;}
      if(ring == 5 && mod == 354702372) { mapRing = 4; mapModule = 354698272;}
      if(ring == 5 && mod == 354702380) { mapRing = 4; mapModule = 354698280;}
      if(ring == 5 && mod == 354702388) { mapRing = 4; mapModule = 354698288;}
      if(ring == 5 && mod == 354702396) { mapRing = 4; mapModule = 354698296;}
      if(ring == 5 && mod == 354702404) { mapRing = 4; mapModule = 354698304;}
      if(ring == 5 && mod == 354702412) { mapRing = 4; mapModule = 354698312;}
      if(ring == 5 && mod == 354702420) { mapRing = 4; mapModule = 354698320;}
      if(ring == 5 && mod == 354702428) { mapRing = 4; mapModule = 354698328;}
      if(ring == 5 && mod == 354702436) { mapRing = 4; mapModule = 354698328;}
      if(ring == 5 && mod == 354702444) { mapRing = 4; mapModule = 354698336;}
      if(ring == 5 && mod == 354702452) { mapRing = 4; mapModule = 354698344;}
      if(ring == 5 && mod == 354702460) { mapRing = 4; mapModule = 354698352;}
      if(ring == 5 && mod == 354702468) { mapRing = 4; mapModule = 354698360;}
      if(ring == 5 && mod == 354702476) { mapRing = 4; mapModule = 354698368;}
      if(ring == 5 && mod == 354702484) { mapRing = 4; mapModule = 354698376;}
      if(ring == 5 && mod == 354702492) { mapRing = 4; mapModule = 354698384;}
      if(ring == 5 && mod == 354702500) { mapRing = 4; mapModule = 354698392;}
      if(ring == 5 && mod == 354702508) { mapRing = 4; mapModule = 354698400;}
      if(ring == 5 && mod == 354702516) { mapRing = 4; mapModule = 354698408;}
      if(ring == 5 && mod == 354702524) { mapRing = 4; mapModule = 354698416;}
      

      if(ring == 4 && mod == 354698248) { mapRing = 5; mapModule = 354702344;}
      if(ring == 4 && mod == 354698256) { mapRing = 5; mapModule = 354702352;}
      if(ring == 4 && mod == 354698264) { mapRing = 5; mapModule = 354702360;}
      if(ring == 4 && mod == 354698272) { mapRing = 5; mapModule = 354702368;}
      if(ring == 4 && mod == 354698280) { mapRing = 5; mapModule = 354702376;}
      if(ring == 4 && mod == 354698288) { mapRing = 5; mapModule = 354702384;}
      if(ring == 4 && mod == 354698296) { mapRing = 5; mapModule = 354702392;}
      if(ring == 4 && mod == 354698304) { mapRing = 5; mapModule = 354702408;}
      if(ring == 4 && mod == 354698312) { mapRing = 5; mapModule = 354702416;}
      if(ring == 4 && mod == 354698320) { mapRing = 5; mapModule = 354702424;}
      if(ring == 4 && mod == 354698328) { mapRing = 5; mapModule = 354702432;}
      if(ring == 4 && mod == 354698336) { mapRing = 5; mapModule = 354702440;}
      if(ring == 4 && mod == 354698344) { mapRing = 5; mapModule = 354702448;}
      if(ring == 4 && mod == 354698352) { mapRing = 5; mapModule = 354702456;}
      if(ring == 4 && mod == 354698360) { mapRing = 5; mapModule = 354702464;}
      if(ring == 4 && mod == 354698368) { mapRing = 5; mapModule = 354702472;}
      if(ring == 4 && mod == 354698376) { mapRing = 5; mapModule = 354702488;}
      if(ring == 4 && mod == 354698384) { mapRing = 5; mapModule = 354702496;}
      if(ring == 4 && mod == 354698392) { mapRing = 5; mapModule = 354702504;}
      if(ring == 4 && mod == 354698400) { mapRing = 5; mapModule = 354702512;}
      if(ring == 4 && mod == 354698408) { mapRing = 5; mapModule = 354702520;}
      if(ring == 4 && mod == 354698416) { mapRing = 5; mapModule = 354702528;}
      

      if(ring == 5 && mod == 354702344) { mapRing = 4; mapModule = 354698244;}
      if(ring == 5 && mod == 354702352) { mapRing = 4; mapModule = 354698252;}
      if(ring == 5 && mod == 354702360) { mapRing = 4; mapModule = 354698260;}
      if(ring == 5 && mod == 354702368) { mapRing = 4; mapModule = 354698268;}
      if(ring == 5 && mod == 354702376) { mapRing = 4; mapModule = 354698276;}
      if(ring == 5 && mod == 354702384) { mapRing = 4; mapModule = 354698284;}
      if(ring == 5 && mod == 354702392) { mapRing = 4; mapModule = 354698292;}
      if(ring == 5 && mod == 354702400) { mapRing = 4; mapModule = 354698300;}
      if(ring == 5 && mod == 354702408) { mapRing = 4; mapModule = 354698308;}
      if(ring == 5 && mod == 354702416) { mapRing = 4; mapModule = 354698316;}
      if(ring == 5 && mod == 354702424) { mapRing = 4; mapModule = 354698324;}
      if(ring == 5 && mod == 354702432) { mapRing = 4; mapModule = 354698332;}
      if(ring == 5 && mod == 354702440) { mapRing = 4; mapModule = 354698332;}
      if(ring == 5 && mod == 354702448) { mapRing = 4; mapModule = 354698340;}
      if(ring == 5 && mod == 354702456) { mapRing = 4; mapModule = 354698348;}
      if(ring == 5 && mod == 354702464) { mapRing = 4; mapModule = 354698356;}
      if(ring == 5 && mod == 354702472) { mapRing = 4; mapModule = 354698364;}
      if(ring == 5 && mod == 354702480) { mapRing = 4; mapModule = 354698372;}
      if(ring == 5 && mod == 354702488) { mapRing = 4; mapModule = 354698380;}
      if(ring == 5 && mod == 354702496) { mapRing = 4; mapModule = 354698388;}
      if(ring == 5 && mod == 354702504) { mapRing = 4; mapModule = 354698396;}
      if(ring == 5 && mod == 354702512) { mapRing = 4; mapModule = 354698404;}
      if(ring == 5 && mod == 354702520) { mapRing = 4; mapModule = 354698412;}
      if(ring == 5 && mod == 354702528) { mapRing = 4; mapModule = 354698412;}
      
      
    }
    
    
    
    
    if(disk == 2) {
      
      
      
      if(ring == 1 && mod == 354948100) { mapRing = 2; mapModule = 354952204;}
      if(ring == 1 && mod == 354948108) { mapRing = 2; mapModule = 354952212;}
      if(ring == 1 && mod == 354948116) { mapRing = 2; mapModule = 354952220;}
      if(ring == 1 && mod == 354948124) { mapRing = 2; mapModule = 354952228;}
      if(ring == 1 && mod == 354948132) { mapRing = 2; mapModule = 354952244;}
      if(ring == 1 && mod == 354948140) { mapRing = 2; mapModule = 354952252;}
      if(ring == 1 && mod == 354948148) { mapRing = 2; mapModule = 354952268;}
      if(ring == 1 && mod == 354948156) { mapRing = 2; mapModule = 354952276;}
      if(ring == 1 && mod == 354948164) { mapRing = 2; mapModule = 354952284;}
      if(ring == 1 && mod == 354948172) { mapRing = 2; mapModule = 354952300;}
      
      
      
      if(ring == 1 && mod == 354948104) { mapRing = 2; mapModule = 354952204;}
      if(ring == 1 && mod == 354948112) { mapRing = 2; mapModule = 354952212;}
      if(ring == 1 && mod == 354948120) { mapRing = 2; mapModule = 354952228;}
      if(ring == 1 && mod == 354948128) { mapRing = 2; mapModule = 354952236;}
      if(ring == 1 && mod == 354948136) { mapRing = 2; mapModule = 354952244;}
      if(ring == 1 && mod == 354948144) { mapRing = 2; mapModule = 354952260;}
      if(ring == 1 && mod == 354948152) { mapRing = 2; mapModule = 354952268;}
      if(ring == 1 && mod == 354948160) { mapRing = 2; mapModule = 354952284;}
      if(ring == 1 && mod == 354948168) { mapRing = 2; mapModule = 354952292;}
      if(ring == 1 && mod == 354948176) { mapRing = 2; mapModule = 354952300;}
      
      
      
      if(ring == 1 && mod == 354948104) { mapRing = 2; mapModule = 354952200;}
      if(ring == 1 && mod == 354948112) { mapRing = 2; mapModule = 354952216;}
      if(ring == 1 && mod == 354948120) { mapRing = 2; mapModule = 354952224;}
      if(ring == 1 && mod == 354948128) { mapRing = 2; mapModule = 354952240;}
      if(ring == 1 && mod == 354948136) { mapRing = 2; mapModule = 354952248;}
      if(ring == 1 && mod == 354948144) { mapRing = 2; mapModule = 354952256;}
      if(ring == 1 && mod == 354948152) { mapRing = 2; mapModule = 354952272;}
      if(ring == 1 && mod == 354948160) { mapRing = 2; mapModule = 354952280;}
      if(ring == 1 && mod == 354948168) { mapRing = 2; mapModule = 354952288;}
      if(ring == 1 && mod == 354948176) { mapRing = 2; mapModule = 354952304;}
      
      
      
      if(ring == 3 && mod == 354956292) { mapRing = 2; mapModule = 354952196;}
      if(ring == 3 && mod == 354956300) { mapRing = 2; mapModule = 354952204;}
      if(ring == 3 && mod == 354956308) { mapRing = 2; mapModule = 354952212;}
      if(ring == 3 && mod == 354956316) { mapRing = 2; mapModule = 354952212;}
      if(ring == 3 && mod == 354956324) { mapRing = 2; mapModule = 354952220;}
      if(ring == 3 && mod == 354956332) { mapRing = 2; mapModule = 354952228;}
      if(ring == 3 && mod == 354956340) { mapRing = 2; mapModule = 354952236;}
      if(ring == 3 && mod == 354956348) { mapRing = 2; mapModule = 354952236;}
      if(ring == 3 && mod == 354956356) { mapRing = 2; mapModule = 354952244;}
      if(ring == 3 && mod == 354956364) { mapRing = 2; mapModule = 354952252;}
      if(ring == 3 && mod == 354956372) { mapRing = 2; mapModule = 354952260;}
      if(ring == 3 && mod == 354956380) { mapRing = 2; mapModule = 354952268;}
      if(ring == 3 && mod == 354956388) { mapRing = 2; mapModule = 354952268;}
      if(ring == 3 && mod == 354956396) { mapRing = 2; mapModule = 354952276;}
      if(ring == 3 && mod == 354956404) { mapRing = 2; mapModule = 354952284;}
      if(ring == 3 && mod == 354956412) { mapRing = 2; mapModule = 354952292;}
      if(ring == 3 && mod == 354956420) { mapRing = 2; mapModule = 354952292;}
      if(ring == 3 && mod == 354956428) { mapRing = 2; mapModule = 354952300;}
      
      
      
      if(ring == 2 && mod == 354952200) { mapRing = 3; mapModule = 354956300;}
      if(ring == 2 && mod == 354952208) { mapRing = 3; mapModule = 354956308;}
      if(ring == 2 && mod == 354952216) { mapRing = 3; mapModule = 354956316;}
      if(ring == 2 && mod == 354952224) { mapRing = 3; mapModule = 354956332;}
      if(ring == 2 && mod == 354952232) { mapRing = 3; mapModule = 354956340;}
      if(ring == 2 && mod == 354952240) { mapRing = 3; mapModule = 354956348;}
      if(ring == 2 && mod == 354952248) { mapRing = 3; mapModule = 354956356;}
      if(ring == 2 && mod == 354952256) { mapRing = 3; mapModule = 354956372;}
      if(ring == 2 && mod == 354952264) { mapRing = 3; mapModule = 354956380;}
      if(ring == 2 && mod == 354952272) { mapRing = 3; mapModule = 354956388;}
      if(ring == 2 && mod == 354952280) { mapRing = 3; mapModule = 354956396;}
      if(ring == 2 && mod == 354952288) { mapRing = 3; mapModule = 354956412;}
      if(ring == 2 && mod == 354952296) { mapRing = 3; mapModule = 354956420;}
      if(ring == 2 && mod == 354952304) { mapRing = 3; mapModule = 354956428;}
      
      
      
      if(ring == 3 && mod == 354956296) { mapRing = 2; mapModule = 354952200;}
      if(ring == 3 && mod == 354956304) { mapRing = 2; mapModule = 354952208;}
      if(ring == 3 && mod == 354956312) { mapRing = 2; mapModule = 354952208;}
      if(ring == 3 && mod == 354956320) { mapRing = 2; mapModule = 354952216;}
      if(ring == 3 && mod == 354956328) { mapRing = 2; mapModule = 354952224;}
      if(ring == 3 && mod == 354956336) { mapRing = 2; mapModule = 354952232;}
      if(ring == 3 && mod == 354956344) { mapRing = 2; mapModule = 354952240;}
      if(ring == 3 && mod == 354956352) { mapRing = 2; mapModule = 354952240;}
      if(ring == 3 && mod == 354956360) { mapRing = 2; mapModule = 354952248;}
      if(ring == 3 && mod == 354956368) { mapRing = 2; mapModule = 354952256;}
      if(ring == 3 && mod == 354956376) { mapRing = 2; mapModule = 354952256;}
      if(ring == 3 && mod == 354956384) { mapRing = 2; mapModule = 354952264;}
      if(ring == 3 && mod == 354956392) { mapRing = 2; mapModule = 354952272;}
      if(ring == 3 && mod == 354956400) { mapRing = 2; mapModule = 354952280;}
      if(ring == 3 && mod == 354956408) { mapRing = 2; mapModule = 354952288;}
      if(ring == 3 && mod == 354956416) { mapRing = 2; mapModule = 354952296;}    
      if(ring == 3 && mod == 354956424) { mapRing = 2; mapModule = 354952296;}
      if(ring == 3 && mod == 354956432) { mapRing = 2; mapModule = 354952304;}
      
      
      if(ring == 3 && mod == 354956292) { mapRing = 4; mapModule = 354960388;}
      if(ring == 3 && mod == 354956300) { mapRing = 4; mapModule = 354960396;}
      if(ring == 3 && mod == 354956308) { mapRing = 4; mapModule = 354960412;}
      if(ring == 3 && mod == 354956316) { mapRing = 4; mapModule = 354960420;}
      if(ring == 3 && mod == 354956324) { mapRing = 4; mapModule = 354960428;}
      if(ring == 3 && mod == 354956332) { mapRing = 4; mapModule = 354960436;}
      if(ring == 3 && mod == 354956340) { mapRing = 4; mapModule = 354960444;}
      if(ring == 3 && mod == 354956348) { mapRing = 4; mapModule = 354960460;}
      if(ring == 3 && mod == 354956356) { mapRing = 4; mapModule = 354960468;}
      if(ring == 3 && mod == 354956364) { mapRing = 4; mapModule = 354960476;}
      if(ring == 3 && mod == 354956372) { mapRing = 4; mapModule = 354960484;}
      if(ring == 3 && mod == 354956380) { mapRing = 4; mapModule = 354960500;}
      if(ring == 3 && mod == 354956388) { mapRing = 4; mapModule = 354960508;}
      if(ring == 3 && mod == 354956396) { mapRing = 4; mapModule = 354960516;}
      if(ring == 3 && mod == 354956404) { mapRing = 4; mapModule = 354960524;}
      if(ring == 3 && mod == 354956412) { mapRing = 4; mapModule = 354960532;}
      if(ring == 3 && mod == 354956420) { mapRing = 4; mapModule = 354960548;}
      if(ring == 3 && mod == 354956428) { mapRing = 4; mapModule = 354960556;}
      
      
      
      if(ring == 3 && mod == 354956296) { mapRing = 4; mapModule = 354960396;}
      if(ring == 3 && mod == 354956304) { mapRing = 4; mapModule = 354960404;}
      if(ring == 3 && mod == 354956312) { mapRing = 4; mapModule = 354960412;}
      if(ring == 3 && mod == 354956320) { mapRing = 4; mapModule = 354960420;}
      if(ring == 3 && mod == 354956328) { mapRing = 4; mapModule = 354960428;}
      if(ring == 3 && mod == 354956336) { mapRing = 4; mapModule = 354960444;}
      if(ring == 3 && mod == 354956344) { mapRing = 4; mapModule = 354960452;}
      if(ring == 3 && mod == 354956352) { mapRing = 4; mapModule = 354960460;}
      if(ring == 3 && mod == 354956360) { mapRing = 4; mapModule = 354960468;}
      if(ring == 3 && mod == 354956368) { mapRing = 4; mapModule = 354960484;}
      if(ring == 3 && mod == 354956376) { mapRing = 4; mapModule = 354960492;}
      if(ring == 3 && mod == 354956384) { mapRing = 4; mapModule = 354960500;}
      if(ring == 3 && mod == 354956392) { mapRing = 4; mapModule = 354960508;}
      if(ring == 3 && mod == 354956400) { mapRing = 4; mapModule = 354960516;}
      if(ring == 3 && mod == 354956408) { mapRing = 4; mapModule = 354960532;}
      if(ring == 3 && mod == 354956416) { mapRing = 4; mapModule = 354960540;}
      if(ring == 3 && mod == 354956424) { mapRing = 4; mapModule = 354960548;}
      if(ring == 3 && mod == 354956432) { mapRing = 4; mapModule = 354960556;}
      
      
      
      if(ring == 3 && mod == 354956296) { mapRing = 4; mapModule = 354960392;}
      if(ring == 3 && mod == 354956304) { mapRing = 4; mapModule = 354960400;}
      if(ring == 3 && mod == 354956312) { mapRing = 4; mapModule = 354960416;}
      if(ring == 3 && mod == 354956320) { mapRing = 4; mapModule = 354960424;}
      if(ring == 3 && mod == 354956328) { mapRing = 4; mapModule = 354960432;}
      if(ring == 3 && mod == 354956336) { mapRing = 4; mapModule = 354960440;}
      if(ring == 3 && mod == 354956344) { mapRing = 4; mapModule = 354960448;}
      if(ring == 3 && mod == 354956352) { mapRing = 4; mapModule = 354960464;}
      if(ring == 3 && mod == 354956360) { mapRing = 4; mapModule = 354960472;}
      if(ring == 3 && mod == 354956368) { mapRing = 4; mapModule = 354960480;}
      if(ring == 3 && mod == 354956376) { mapRing = 4; mapModule = 354960488;}
      if(ring == 3 && mod == 354956384) { mapRing = 4; mapModule = 354960504;}
      if(ring == 3 && mod == 354956392) { mapRing = 4; mapModule = 354960512;}
      if(ring == 3 && mod == 354956400) { mapRing = 4; mapModule = 354960520;}
      if(ring == 3 && mod == 354956408) { mapRing = 4; mapModule = 354960528;}
      if(ring == 3 && mod == 354956416) { mapRing = 4; mapModule = 354960536;}
      if(ring == 3 && mod == 354956424) { mapRing = 4; mapModule = 354960552;}
      if(ring == 3 && mod == 354956432) { mapRing = 4; mapModule = 354960560;}
      
      
      
      if(ring == 5 && mod == 354964484) { mapRing = 4; mapModule = 354960560;}
      if(ring == 5 && mod == 354964492) { mapRing = 4; mapModule = 354960392;}
      if(ring == 5 && mod == 354964500) { mapRing = 4; mapModule = 354960400;}
      if(ring == 5 && mod == 354964508) { mapRing = 4; mapModule = 354960408;}
      if(ring == 5 && mod == 354964516) { mapRing = 4; mapModule = 354960416;}
      if(ring == 5 && mod == 354964524) { mapRing = 4; mapModule = 354960424;}
      if(ring == 5 && mod == 354964532) { mapRing = 4; mapModule = 354960432;}
      if(ring == 5 && mod == 354964540) { mapRing = 4; mapModule = 354960440;}
      if(ring == 5 && mod == 354964548) { mapRing = 4; mapModule = 354960448;}
      if(ring == 5 && mod == 354964556) { mapRing = 4; mapModule = 354960456;}
      if(ring == 5 && mod == 354964564) { mapRing = 4; mapModule = 354960464;}
      if(ring == 5 && mod == 354964572) { mapRing = 4; mapModule = 354960472;}
      if(ring == 5 && mod == 354964580) { mapRing = 4; mapModule = 354960472;}
      if(ring == 5 && mod == 354964588) { mapRing = 4; mapModule = 354960480;}
      if(ring == 5 && mod == 354964596) { mapRing = 4; mapModule = 354960488;}
      if(ring == 5 && mod == 354964604) { mapRing = 4; mapModule = 354960496;}
      if(ring == 5 && mod == 354964612) { mapRing = 4; mapModule = 354960504;}
      if(ring == 5 && mod == 354964620) { mapRing = 4; mapModule = 354960512;}
      if(ring == 5 && mod == 354964628) { mapRing = 4; mapModule = 354960520;}
      if(ring == 5 && mod == 354964636) { mapRing = 4; mapModule = 354960528;}
      if(ring == 5 && mod == 354964644) { mapRing = 4; mapModule = 354960536;}
      if(ring == 5 && mod == 354964652) { mapRing = 4; mapModule = 354960544;}
      if(ring == 5 && mod == 354964660) { mapRing = 4; mapModule = 354960552;}
      if(ring == 5 && mod == 354964668) { mapRing = 4; mapModule = 354960560;}
      
      
      
      
      if(ring == 4 && mod == 354960392) { mapRing = 5; mapModule = 354964488;}
      if(ring == 4 && mod == 354960400) { mapRing = 5; mapModule = 354964496;}
      if(ring == 4 && mod == 354960408) { mapRing = 5; mapModule = 354964504;}
      if(ring == 4 && mod == 354960416) { mapRing = 5; mapModule = 354964512;}
      if(ring == 4 && mod == 354960424) { mapRing = 5; mapModule = 354964520;}
      if(ring == 4 && mod == 354960432) { mapRing = 5; mapModule = 354964528;}
      if(ring == 4 && mod == 354960440) { mapRing = 5; mapModule = 354964536;}
      if(ring == 4 && mod == 354960448) { mapRing = 5; mapModule = 354964552;}
      if(ring == 4 && mod == 354960456) { mapRing = 5; mapModule = 354964560;}
      if(ring == 4 && mod == 354960464) { mapRing = 5; mapModule = 354964568;}
      if(ring == 4 && mod == 354960472) { mapRing = 5; mapModule = 354964576;}
      if(ring == 4 && mod == 354960480) { mapRing = 5; mapModule = 354964584;}
      if(ring == 4 && mod == 354960488) { mapRing = 5; mapModule = 354964592;}
      if(ring == 4 && mod == 354960496) { mapRing = 5; mapModule = 354964600;}
      if(ring == 4 && mod == 354960504) { mapRing = 5; mapModule = 354964608;}
      if(ring == 4 && mod == 354960512) { mapRing = 5; mapModule = 354964616;}
      if(ring == 4 && mod == 354960520) { mapRing = 5; mapModule = 354964632;}
      if(ring == 4 && mod == 354960528) { mapRing = 5; mapModule = 354964640;}
      if(ring == 4 && mod == 354960536) { mapRing = 5; mapModule = 354964648;}
      if(ring == 4 && mod == 354960544) { mapRing = 5; mapModule = 354964656;}
      if(ring == 4 && mod == 354960552) { mapRing = 5; mapModule = 354964664;}
      if(ring == 4 && mod == 354960560) { mapRing = 5; mapModule = 354964672;}
      
      
      
      if(ring == 5 && mod == 354964488) { mapRing = 4; mapModule = 354960388;}
      if(ring == 5 && mod == 354964496) { mapRing = 4; mapModule = 354960396;}
      if(ring == 5 && mod == 354964504) { mapRing = 4; mapModule = 354960404;}
      if(ring == 5 && mod == 354964512) { mapRing = 4; mapModule = 354960412;}
      if(ring == 5 && mod == 354964520) { mapRing = 4; mapModule = 354960420;}
      if(ring == 5 && mod == 354964528) { mapRing = 4; mapModule = 354960428;}
      if(ring == 5 && mod == 354964536) { mapRing = 4; mapModule = 354960436;}
      if(ring == 5 && mod == 354964544) { mapRing = 4; mapModule = 354960444;}
      if(ring == 5 && mod == 354964552) { mapRing = 4; mapModule = 354960452;}
      if(ring == 5 && mod == 354964560) { mapRing = 4; mapModule = 354960460;}
      if(ring == 5 && mod == 354964568) { mapRing = 4; mapModule = 354960468;}
      if(ring == 5 && mod == 354964576) { mapRing = 4; mapModule = 354960476;}
      if(ring == 5 && mod == 354964584) { mapRing = 4; mapModule = 354960476;}
      if(ring == 5 && mod == 354964592) { mapRing = 4; mapModule = 354960484;}
      if(ring == 5 && mod == 354964600) { mapRing = 4; mapModule = 354960492;}
      if(ring == 5 && mod == 354964608) { mapRing = 4; mapModule = 354960500;}
      if(ring == 5 && mod == 354964616) { mapRing = 4; mapModule = 354960508;}
      if(ring == 5 && mod == 354964624) { mapRing = 4; mapModule = 354960516;}
      if(ring == 5 && mod == 354964632) { mapRing = 4; mapModule = 354960524;}
      if(ring == 5 && mod == 354964640) { mapRing = 4; mapModule = 354960532;}
      if(ring == 5 && mod == 354964648) { mapRing = 4; mapModule = 354960540;}
      if(ring == 5 && mod == 354964656) { mapRing = 4; mapModule = 354960548;}
      if(ring == 5 && mod == 354964664) { mapRing = 4; mapModule = 354960556;}
      if(ring == 5 && mod == 354964672) { mapRing = 4; mapModule = 354960556;}
      
      
      
    }
    
    
    if(disk == 3) {
      
      
      
      if(ring == 1 && mod == 355210244) { mapRing = 2; mapModule = 355214348;}
      if(ring == 1 && mod == 355210252) { mapRing = 2; mapModule = 355214356;}
      if(ring == 1 && mod == 355210260) { mapRing = 2; mapModule = 355214364;}
      if(ring == 1 && mod == 355210268) { mapRing = 2; mapModule = 355214372;}
      if(ring == 1 && mod == 355210276) { mapRing = 2; mapModule = 355214388;}
      if(ring == 1 && mod == 355210284) { mapRing = 2; mapModule = 355214396;}
      if(ring == 1 && mod == 355210292) { mapRing = 2; mapModule = 355214412;}
      if(ring == 1 && mod == 355210300) { mapRing = 2; mapModule = 355214420;}
      if(ring == 1 && mod == 355210308) { mapRing = 2; mapModule = 355214428;}
      if(ring == 1 && mod == 355210316) { mapRing = 2; mapModule = 355214444;}
      
      
      
      if(ring == 1 && mod == 355210248) { mapRing = 2; mapModule = 355214348;}
      if(ring == 1 && mod == 355210256) { mapRing = 2; mapModule = 355214356;}
      if(ring == 1 && mod == 355210264) { mapRing = 2; mapModule = 355214372;}
      if(ring == 1 && mod == 355210272) { mapRing = 2; mapModule = 355214380;}
      if(ring == 1 && mod == 355210280) { mapRing = 2; mapModule = 355214388;}
      if(ring == 1 && mod == 355210288) { mapRing = 2; mapModule = 355214404;}
      if(ring == 1 && mod == 355210296) { mapRing = 2; mapModule = 355214412;}
      if(ring == 1 && mod == 355210304) { mapRing = 2; mapModule = 355214428;}
      if(ring == 1 && mod == 355210312) { mapRing = 2; mapModule = 355214436;}
      if(ring == 1 && mod == 355210320) { mapRing = 2; mapModule = 355214444;}
      
      
      
      if(ring == 1 && mod == 355210248) { mapRing = 2; mapModule = 355214344;}
      if(ring == 1 && mod == 355210256) { mapRing = 2; mapModule = 355214360;}
      if(ring == 1 && mod == 355210264) { mapRing = 2; mapModule = 355214368;}
      if(ring == 1 && mod == 355210272) { mapRing = 2; mapModule = 355214384;}
      if(ring == 1 && mod == 355210280) { mapRing = 2; mapModule = 355214392;}
      if(ring == 1 && mod == 355210288) { mapRing = 2; mapModule = 355214400;}
      if(ring == 1 && mod == 355210296) { mapRing = 2; mapModule = 355214416;}
      if(ring == 1 && mod == 355210304) { mapRing = 2; mapModule = 355214424;}
      if(ring == 1 && mod == 355210312) { mapRing = 2; mapModule = 355214432;}
      if(ring == 1 && mod == 355210320) { mapRing = 2; mapModule = 355214448;}
      

      
      if(ring == 3 && mod == 355218436) { mapRing = 2; mapModule = 355214340;}
      if(ring == 3 && mod == 355218444) { mapRing = 2; mapModule = 355214348;}
      if(ring == 3 && mod == 355218452) { mapRing = 2; mapModule = 355214356;}
      if(ring == 3 && mod == 355218460) { mapRing = 2; mapModule = 355214356;}
      if(ring == 3 && mod == 355218468) { mapRing = 2; mapModule = 355214364;}
      if(ring == 3 && mod == 355218476) { mapRing = 2; mapModule = 355214372;}
      if(ring == 3 && mod == 355218484) { mapRing = 2; mapModule = 355214380;}
      if(ring == 3 && mod == 355218492) { mapRing = 2; mapModule = 355214380;}
      if(ring == 3 && mod == 355218500) { mapRing = 2; mapModule = 355214388;}
      if(ring == 3 && mod == 355218508) { mapRing = 2; mapModule = 355214396;}
      if(ring == 3 && mod == 355218516) { mapRing = 2; mapModule = 355214404;}
      if(ring == 3 && mod == 355218524) { mapRing = 2; mapModule = 355214412;}
      if(ring == 3 && mod == 355218532) { mapRing = 2; mapModule = 355214412;}
      if(ring == 3 && mod == 355218540) { mapRing = 2; mapModule = 355214420;}
      if(ring == 3 && mod == 355218548) { mapRing = 2; mapModule = 355214428;}
      if(ring == 3 && mod == 355218556) { mapRing = 2; mapModule = 355214436;}
      if(ring == 3 && mod == 355218564) { mapRing = 2; mapModule = 355214436;}
      if(ring == 3 && mod == 355218572) { mapRing = 2; mapModule = 355214444;}
      
      
      
      if(ring == 2 && mod == 355214344) { mapRing = 3; mapModule = 355218444;}
      if(ring == 2 && mod == 355214352) { mapRing = 3; mapModule = 355218452;}
      if(ring == 2 && mod == 355214360) { mapRing = 3; mapModule = 355218460;}
      if(ring == 2 && mod == 355214368) { mapRing = 3; mapModule = 355218476;}
      if(ring == 2 && mod == 355214376) { mapRing = 3; mapModule = 355218484;}
      if(ring == 2 && mod == 355214384) { mapRing = 3; mapModule = 355218492;}
      if(ring == 2 && mod == 355214392) { mapRing = 3; mapModule = 355218500;}
      if(ring == 2 && mod == 355214400) { mapRing = 3; mapModule = 355218516;}
      if(ring == 2 && mod == 355214408) { mapRing = 3; mapModule = 355218524;}
      if(ring == 2 && mod == 355214416) { mapRing = 3; mapModule = 355218532;}
      if(ring == 2 && mod == 355214424) { mapRing = 3; mapModule = 355218540;}
      if(ring == 2 && mod == 355214432) { mapRing = 3; mapModule = 355218556;}
      if(ring == 2 && mod == 355214440) { mapRing = 3; mapModule = 355218564;}
      if(ring == 2 && mod == 355214448) { mapRing = 3; mapModule = 355218572;}
      
      
      
      if(ring == 3 && mod == 355218440) { mapRing = 2; mapModule = 355214344;}
      if(ring == 3 && mod == 355218448) { mapRing = 2; mapModule = 355214352;}
      if(ring == 3 && mod == 355218456) { mapRing = 2; mapModule = 355214352;}
      if(ring == 3 && mod == 355218464) { mapRing = 2; mapModule = 355214360;}
      if(ring == 3 && mod == 355218472) { mapRing = 2; mapModule = 355214368;}
      if(ring == 3 && mod == 355218480) { mapRing = 2; mapModule = 355214376;}
      if(ring == 3 && mod == 355218488) { mapRing = 2; mapModule = 355214384;}
      if(ring == 3 && mod == 355218496) { mapRing = 2; mapModule = 355214384;}
      if(ring == 3 && mod == 355218504) { mapRing = 2; mapModule = 355214392;}
      if(ring == 3 && mod == 355218512) { mapRing = 2; mapModule = 355214400;}
      if(ring == 3 && mod == 355218520) { mapRing = 2; mapModule = 355214400;}
      if(ring == 3 && mod == 355218528) { mapRing = 2; mapModule = 355214408;}
      if(ring == 3 && mod == 355218536) { mapRing = 2; mapModule = 355214416;}
      if(ring == 3 && mod == 355218544) { mapRing = 2; mapModule = 355214424;}
      if(ring == 3 && mod == 355218552) { mapRing = 2; mapModule = 355214432;}
      if(ring == 3 && mod == 355218560) { mapRing = 2; mapModule = 355214440;}    
      if(ring == 3 && mod == 355218568) { mapRing = 2; mapModule = 355214440;}
      if(ring == 3 && mod == 355218576) { mapRing = 2; mapModule = 355214448;}

      
      if(ring == 3 && mod == 355218436) { mapRing = 4; mapModule = 355222532;}
      if(ring == 3 && mod == 355218444) { mapRing = 4; mapModule = 355222540;}
      if(ring == 3 && mod == 355218452) { mapRing = 4; mapModule = 355222556;}
      if(ring == 3 && mod == 355218460) { mapRing = 4; mapModule = 355222564;}
      if(ring == 3 && mod == 355218468) { mapRing = 4; mapModule = 355222572;}
      if(ring == 3 && mod == 355218476) { mapRing = 4; mapModule = 355222580;}
      if(ring == 3 && mod == 355218484) { mapRing = 4; mapModule = 355222588;}
      if(ring == 3 && mod == 355218492) { mapRing = 4; mapModule = 355222604;}
      if(ring == 3 && mod == 355218500) { mapRing = 4; mapModule = 355222612;}
      if(ring == 3 && mod == 355218508) { mapRing = 4; mapModule = 355222620;}
      if(ring == 3 && mod == 355218516) { mapRing = 4; mapModule = 355222628;}
      if(ring == 3 && mod == 355218524) { mapRing = 4; mapModule = 355222644;}
      if(ring == 3 && mod == 355218532) { mapRing = 4; mapModule = 355222652;}
      if(ring == 3 && mod == 355218540) { mapRing = 4; mapModule = 355222660;}
      if(ring == 3 && mod == 355218548) { mapRing = 4; mapModule = 355222668;}
      if(ring == 3 && mod == 355218556) { mapRing = 4; mapModule = 355222676;}
      if(ring == 3 && mod == 355218564) { mapRing = 4; mapModule = 355222692;}
      if(ring == 3 && mod == 355218572) { mapRing = 4; mapModule = 355222700;}
      

      
      if(ring == 3 && mod == 355218440) { mapRing = 4; mapModule = 355222540;}
      if(ring == 3 && mod == 355218448) { mapRing = 4; mapModule = 355222548;}
      if(ring == 3 && mod == 355218456) { mapRing = 4; mapModule = 355222556;}
      if(ring == 3 && mod == 355218464) { mapRing = 4; mapModule = 355222564;}
      if(ring == 3 && mod == 355218472) { mapRing = 4; mapModule = 355222572;}
      if(ring == 3 && mod == 355218480) { mapRing = 4; mapModule = 355222588;}
      if(ring == 3 && mod == 355218488) { mapRing = 4; mapModule = 355222596;}
      if(ring == 3 && mod == 355218496) { mapRing = 4; mapModule = 355222604;}
      if(ring == 3 && mod == 355218504) { mapRing = 4; mapModule = 355222612;}
      if(ring == 3 && mod == 355218512) { mapRing = 4; mapModule = 355222628;}
      if(ring == 3 && mod == 355218520) { mapRing = 4; mapModule = 355222636;}
      if(ring == 3 && mod == 355218528) { mapRing = 4; mapModule = 355222644;}
      if(ring == 3 && mod == 355218536) { mapRing = 4; mapModule = 355222652;}
      if(ring == 3 && mod == 355218544) { mapRing = 4; mapModule = 355222660;}
      if(ring == 3 && mod == 355218552) { mapRing = 4; mapModule = 355222676;}
      if(ring == 3 && mod == 355218560) { mapRing = 4; mapModule = 355222684;}
      if(ring == 3 && mod == 355218568) { mapRing = 4; mapModule = 355222692;}
      if(ring == 3 && mod == 355218576) { mapRing = 4; mapModule = 355222700;}

      
      
      if(ring == 3 && mod == 355218440) { mapRing = 4; mapModule = 355222536;}
      if(ring == 3 && mod == 355218448) { mapRing = 4; mapModule = 355222544;}
      if(ring == 3 && mod == 355218456) { mapRing = 4; mapModule = 355222560;}
      if(ring == 3 && mod == 355218464) { mapRing = 4; mapModule = 355222568;}
      if(ring == 3 && mod == 355218472) { mapRing = 4; mapModule = 355222576;}
      if(ring == 3 && mod == 355218480) { mapRing = 4; mapModule = 355222584;}
      if(ring == 3 && mod == 355218488) { mapRing = 4; mapModule = 355222592;}
      if(ring == 3 && mod == 355218496) { mapRing = 4; mapModule = 355222608;}
      if(ring == 3 && mod == 355218504) { mapRing = 4; mapModule = 355222616;}
      if(ring == 3 && mod == 355218512) { mapRing = 4; mapModule = 355222624;}
      if(ring == 3 && mod == 355218520) { mapRing = 4; mapModule = 355222632;}
      if(ring == 3 && mod == 355218528) { mapRing = 4; mapModule = 355222648;}
      if(ring == 3 && mod == 355218536) { mapRing = 4; mapModule = 355222656;}
      if(ring == 3 && mod == 355218544) { mapRing = 4; mapModule = 355222664;}
      if(ring == 3 && mod == 355218552) { mapRing = 4; mapModule = 355222672;}
      if(ring == 3 && mod == 355218560) { mapRing = 4; mapModule = 355222680;}
      if(ring == 3 && mod == 355218568) { mapRing = 4; mapModule = 355222696;}
      if(ring == 3 && mod == 355218576) { mapRing = 4; mapModule = 355222704;}
      
      
      
      if(ring == 5 && mod == 355226628) { mapRing = 4; mapModule = 355222704;}
      if(ring == 5 && mod == 355226636) { mapRing = 4; mapModule = 355222536;}
      if(ring == 5 && mod == 355226644) { mapRing = 4; mapModule = 355222544;}
      if(ring == 5 && mod == 355226652) { mapRing = 4; mapModule = 355222552;}
      if(ring == 5 && mod == 355226660) { mapRing = 4; mapModule = 355222560;}
      if(ring == 5 && mod == 355226668) { mapRing = 4; mapModule = 355222568;}
      if(ring == 5 && mod == 355226676) { mapRing = 4; mapModule = 355222576;}
      if(ring == 5 && mod == 355226684) { mapRing = 4; mapModule = 355222584;}
      if(ring == 5 && mod == 355226692) { mapRing = 4; mapModule = 355222592;}
      if(ring == 5 && mod == 355226700) { mapRing = 4; mapModule = 355222600;}
      if(ring == 5 && mod == 355226708) { mapRing = 4; mapModule = 355222608;}
      if(ring == 5 && mod == 355226716) { mapRing = 4; mapModule = 355222616;}
      if(ring == 5 && mod == 355226724) { mapRing = 4; mapModule = 355222616;}
      if(ring == 5 && mod == 355226732) { mapRing = 4; mapModule = 355222624;}
      if(ring == 5 && mod == 355226740) { mapRing = 4; mapModule = 355222632;}
      if(ring == 5 && mod == 355226748) { mapRing = 4; mapModule = 355222640;}
      if(ring == 5 && mod == 355226756) { mapRing = 4; mapModule = 355222648;}
      if(ring == 5 && mod == 355226764) { mapRing = 4; mapModule = 355222656;}
      if(ring == 5 && mod == 355226772) { mapRing = 4; mapModule = 355222664;}
      if(ring == 5 && mod == 355226780) { mapRing = 4; mapModule = 355222672;}
      if(ring == 5 && mod == 355226788) { mapRing = 4; mapModule = 355222680;}
      if(ring == 5 && mod == 355226796) { mapRing = 4; mapModule = 355222688;}
      if(ring == 5 && mod == 355226804) { mapRing = 4; mapModule = 355222696;}
      if(ring == 5 && mod == 355226812) { mapRing = 4; mapModule = 355222704;}
      
      
      

      if(ring == 4 && mod == 355222536) { mapRing = 5; mapModule = 355226632;}
      if(ring == 4 && mod == 355222544) { mapRing = 5; mapModule = 355226640;}
      if(ring == 4 && mod == 355222552) { mapRing = 5; mapModule = 355226648;}
      if(ring == 4 && mod == 355222560) { mapRing = 5; mapModule = 355226656;}
      if(ring == 4 && mod == 355222568) { mapRing = 5; mapModule = 355226664;}
      if(ring == 4 && mod == 355222576) { mapRing = 5; mapModule = 355226672;}
      if(ring == 4 && mod == 355222584) { mapRing = 5; mapModule = 355226680;}
      if(ring == 4 && mod == 355222592) { mapRing = 5; mapModule = 355226696;}
      if(ring == 4 && mod == 355222600) { mapRing = 5; mapModule = 355226704;}
      if(ring == 4 && mod == 355222608) { mapRing = 5; mapModule = 355226712;}
      if(ring == 4 && mod == 355222616) { mapRing = 5; mapModule = 355226720;}
      if(ring == 4 && mod == 355222624) { mapRing = 5; mapModule = 355226728;}
      if(ring == 4 && mod == 355222632) { mapRing = 5; mapModule = 355226736;}
      if(ring == 4 && mod == 355222640) { mapRing = 5; mapModule = 355226744;}
      if(ring == 4 && mod == 355222648) { mapRing = 5; mapModule = 355226752;}
      if(ring == 4 && mod == 355222656) { mapRing = 5; mapModule = 355226760;}
      if(ring == 4 && mod == 355222664) { mapRing = 5; mapModule = 355226776;}
      if(ring == 4 && mod == 355222672) { mapRing = 5; mapModule = 355226784;}
      if(ring == 4 && mod == 355222680) { mapRing = 5; mapModule = 355226792;}
      if(ring == 4 && mod == 355222688) { mapRing = 5; mapModule = 355226800;}
      if(ring == 4 && mod == 355222696) { mapRing = 5; mapModule = 355226808;}
      if(ring == 4 && mod == 355222704) { mapRing = 5; mapModule = 355226816;}
      

      
      if(ring == 5 && mod == 355226632) { mapRing = 4; mapModule = 355222532;}
      if(ring == 5 && mod == 355226640) { mapRing = 4; mapModule = 355222540;}
      if(ring == 5 && mod == 355226648) { mapRing = 4; mapModule = 355222548;}
      if(ring == 5 && mod == 355226656) { mapRing = 4; mapModule = 355222556;}
      if(ring == 5 && mod == 355226664) { mapRing = 4; mapModule = 355222564;}
      if(ring == 5 && mod == 355226672) { mapRing = 4; mapModule = 355222572;}
      if(ring == 5 && mod == 355226680) { mapRing = 4; mapModule = 355222580;}
      if(ring == 5 && mod == 355226688) { mapRing = 4; mapModule = 355222588;}
      if(ring == 5 && mod == 355226696) { mapRing = 4; mapModule = 355222596;}
      if(ring == 5 && mod == 355226704) { mapRing = 4; mapModule = 355222604;}
      if(ring == 5 && mod == 355226712) { mapRing = 4; mapModule = 355222612;}
      if(ring == 5 && mod == 355226720) { mapRing = 4; mapModule = 355222620;}
      if(ring == 5 && mod == 355226728) { mapRing = 4; mapModule = 355222620;}
      if(ring == 5 && mod == 355226736) { mapRing = 4; mapModule = 355222628;}
      if(ring == 5 && mod == 355226744) { mapRing = 4; mapModule = 355222636;}
      if(ring == 5 && mod == 355226752) { mapRing = 4; mapModule = 355222644;}
      if(ring == 5 && mod == 355226760) { mapRing = 4; mapModule = 355222652;}
      if(ring == 5 && mod == 355226768) { mapRing = 4; mapModule = 355222660;}
      if(ring == 5 && mod == 355226776) { mapRing = 4; mapModule = 355222668;}
      if(ring == 5 && mod == 355226784) { mapRing = 4; mapModule = 355222676;}
      if(ring == 5 && mod == 355226792) { mapRing = 4; mapModule = 355222684;}
      if(ring == 5 && mod == 355226800) { mapRing = 4; mapModule = 355222692;}
      if(ring == 5 && mod == 355226808) { mapRing = 4; mapModule = 355222700;}
      if(ring == 5 && mod == 355226816) { mapRing = 4; mapModule = 355222700;}
      
      
    }
    
    
    if (disk == 4) {
      
      
      if(ring == 1 && mod == 355472388) { mapRing = 2; mapModule = 355476492;}
      if(ring == 1 && mod == 355472396) { mapRing = 2; mapModule = 355476500;}
      if(ring == 1 && mod == 355472404) { mapRing = 2; mapModule = 355476508;}
      if(ring == 1 && mod == 355472412) { mapRing = 2; mapModule = 355476516;}
      if(ring == 1 && mod == 355472420) { mapRing = 2; mapModule = 355476532;}
      if(ring == 1 && mod == 355472428) { mapRing = 2; mapModule = 355476540;}
      if(ring == 1 && mod == 355472436) { mapRing = 2; mapModule = 355476556;}
      if(ring == 1 && mod == 355472444) { mapRing = 2; mapModule = 355476564;}
      if(ring == 1 && mod == 355472452) { mapRing = 2; mapModule = 355476572;}
      if(ring == 1 && mod == 355472460) { mapRing = 2; mapModule = 355476588;}
      
      
      if(ring == 1 && mod == 355472392) { mapRing = 2; mapModule = 355476492;}
      if(ring == 1 && mod == 355472400) { mapRing = 2; mapModule = 355476500;}
      if(ring == 1 && mod == 355472408) { mapRing = 2; mapModule = 355476516;}
      if(ring == 1 && mod == 355472416) { mapRing = 2; mapModule = 355476524;}
      if(ring == 1 && mod == 355472424) { mapRing = 2; mapModule = 355476532;}
      if(ring == 1 && mod == 355472432) { mapRing = 2; mapModule = 355476548;}
      if(ring == 1 && mod == 355472440) { mapRing = 2; mapModule = 355476556;}
      if(ring == 1 && mod == 355472448) { mapRing = 2; mapModule = 355476572;}
      if(ring == 1 && mod == 355472456) { mapRing = 2; mapModule = 355476580;}
      if(ring == 1 && mod == 355472464) { mapRing = 2; mapModule = 355476588;}
      
      
      if(ring == 1 && mod == 355472392) { mapRing = 2; mapModule = 355476488;}
      if(ring == 1 && mod == 355472400) { mapRing = 2; mapModule = 355476504;}
      if(ring == 1 && mod == 355472408) { mapRing = 2; mapModule = 355476512;}
      if(ring == 1 && mod == 355472416) { mapRing = 2; mapModule = 355476528;}
      if(ring == 1 && mod == 355472424) { mapRing = 2; mapModule = 355476536;}
      if(ring == 1 && mod == 355472432) { mapRing = 2; mapModule = 355476544;}
      if(ring == 1 && mod == 355472440) { mapRing = 2; mapModule = 355476560;}
      if(ring == 1 && mod == 355472448) { mapRing = 2; mapModule = 355476568;}
      if(ring == 1 && mod == 355472456) { mapRing = 2; mapModule = 355476576;}
      if(ring == 1 && mod == 355472464) { mapRing = 2; mapModule = 355476592;}
      
      
      if(ring == 3 && mod == 355480580) { mapRing = 2; mapModule = 355476484;}
      if(ring == 3 && mod == 355480588) { mapRing = 2; mapModule = 355476492;}
      if(ring == 3 && mod == 355480596) { mapRing = 2; mapModule = 355476500;}
      if(ring == 3 && mod == 355480604) { mapRing = 2; mapModule = 355476500;}
      if(ring == 3 && mod == 355480612) { mapRing = 2; mapModule = 355476508;}
      if(ring == 3 && mod == 355480620) { mapRing = 2; mapModule = 355476516;}
      if(ring == 3 && mod == 355480628) { mapRing = 2; mapModule = 355476524;}
      if(ring == 3 && mod == 355480636) { mapRing = 2; mapModule = 355476524;}
      if(ring == 3 && mod == 355480644) { mapRing = 2; mapModule = 355476532;}
      if(ring == 3 && mod == 355480652) { mapRing = 2; mapModule = 355476540;}
      if(ring == 3 && mod == 355480660) { mapRing = 2; mapModule = 355476548;}
      if(ring == 3 && mod == 355480668) { mapRing = 2; mapModule = 355476556;}
      if(ring == 3 && mod == 355480676) { mapRing = 2; mapModule = 355476556;}
      if(ring == 3 && mod == 355480684) { mapRing = 2; mapModule = 355476564;}
      if(ring == 3 && mod == 355480692) { mapRing = 2; mapModule = 355476572;}
      if(ring == 3 && mod == 355480700) { mapRing = 2; mapModule = 355476580;}
      if(ring == 3 && mod == 355480708) { mapRing = 2; mapModule = 355476580;}
      if(ring == 3 && mod == 355480716) { mapRing = 2; mapModule = 355476588;}
      
      
      if(ring == 2 && mod == 355476488) { mapRing = 3; mapModule = 355480588;}
      if(ring == 2 && mod == 355476496) { mapRing = 3; mapModule = 355480596;}
      if(ring == 2 && mod == 355476504) { mapRing = 3; mapModule = 355480604;}
      if(ring == 2 && mod == 355476512) { mapRing = 3; mapModule = 355480620;}
      if(ring == 2 && mod == 355476520) { mapRing = 3; mapModule = 355480628;}
      if(ring == 2 && mod == 355476528) { mapRing = 3; mapModule = 355480636;}
      if(ring == 2 && mod == 355476536) { mapRing = 3; mapModule = 355480644;}
      if(ring == 2 && mod == 355476544) { mapRing = 3; mapModule = 355480660;}
      if(ring == 2 && mod == 355476552) { mapRing = 3; mapModule = 355480668;}
      if(ring == 2 && mod == 355476560) { mapRing = 3; mapModule = 355480676;}
      if(ring == 2 && mod == 355476568) { mapRing = 3; mapModule = 355480684;}
      if(ring == 2 && mod == 355476576) { mapRing = 3; mapModule = 355480700;}
      if(ring == 2 && mod == 355476584) { mapRing = 3; mapModule = 355480708;}
      if(ring == 2 && mod == 355476592) { mapRing = 3; mapModule = 355480716;}
      
      
      if(ring == 3 && mod == 355480584) { mapRing = 2; mapModule = 355476488;}
      if(ring == 3 && mod == 355480592) { mapRing = 2; mapModule = 355476496;}
      if(ring == 3 && mod == 355480600) { mapRing = 2; mapModule = 355476496;}
      if(ring == 3 && mod == 355480608) { mapRing = 2; mapModule = 355476504;}
      if(ring == 3 && mod == 355480616) { mapRing = 2; mapModule = 355476512;}
      if(ring == 3 && mod == 355480624) { mapRing = 2; mapModule = 355476520;}
      if(ring == 3 && mod == 355480632) { mapRing = 2; mapModule = 355476528;}
      if(ring == 3 && mod == 355480640) { mapRing = 2; mapModule = 355476528;}
      if(ring == 3 && mod == 355480648) { mapRing = 2; mapModule = 355476536;}
      if(ring == 3 && mod == 355480656) { mapRing = 2; mapModule = 355476544;}
      if(ring == 3 && mod == 355480664) { mapRing = 2; mapModule = 355476544;}
      if(ring == 3 && mod == 355480672) { mapRing = 2; mapModule = 355476552;}
      if(ring == 3 && mod == 355480680) { mapRing = 2; mapModule = 355476560;}
      if(ring == 3 && mod == 355480688) { mapRing = 2; mapModule = 355476568;}
      if(ring == 3 && mod == 355480696) { mapRing = 2; mapModule = 355476576;}
      if(ring == 3 && mod == 355480704) { mapRing = 2; mapModule = 355476584;}    
      if(ring == 3 && mod == 355480712) { mapRing = 2; mapModule = 355476584;}
      if(ring == 3 && mod == 355480720) { mapRing = 2; mapModule = 355476592;}
      
      
      if(ring == 3 && mod == 355480580) { mapRing = 4; mapModule = 355484676;}
      if(ring == 3 && mod == 355480588) { mapRing = 4; mapModule = 355484684;}
      if(ring == 3 && mod == 355480596) { mapRing = 4; mapModule = 355484700;}
      if(ring == 3 && mod == 355480604) { mapRing = 4; mapModule = 355484708;}
      if(ring == 3 && mod == 355480612) { mapRing = 4; mapModule = 355484716;}
      if(ring == 3 && mod == 355480620) { mapRing = 4; mapModule = 355484724;}
      if(ring == 3 && mod == 355480628) { mapRing = 4; mapModule = 355484732;}
      if(ring == 3 && mod == 355480636) { mapRing = 4; mapModule = 355484748;}
      if(ring == 3 && mod == 355480644) { mapRing = 4; mapModule = 355484756;}
      if(ring == 3 && mod == 355480652) { mapRing = 4; mapModule = 355484764;}
      if(ring == 3 && mod == 355480660) { mapRing = 4; mapModule = 355484772;}
      if(ring == 3 && mod == 355480668) { mapRing = 4; mapModule = 355484788;}
      if(ring == 3 && mod == 355480676) { mapRing = 4; mapModule = 355484796;}
      if(ring == 3 && mod == 355480684) { mapRing = 4; mapModule = 355484804;}
      if(ring == 3 && mod == 355480692) { mapRing = 4; mapModule = 355484812;}
      if(ring == 3 && mod == 355480700) { mapRing = 4; mapModule = 355484820;}
      if(ring == 3 && mod == 355480708) { mapRing = 4; mapModule = 355484836;}
      if(ring == 3 && mod == 355480716) { mapRing = 4; mapModule = 355484844;}
      
      
      
      if(ring == 3 && mod == 355480584) { mapRing = 4; mapModule = 355484684;}
      if(ring == 3 && mod == 355480592) { mapRing = 4; mapModule = 355484692;}
      if(ring == 3 && mod == 355480600) { mapRing = 4; mapModule = 355484700;}
      if(ring == 3 && mod == 355480608) { mapRing = 4; mapModule = 355484708;}
      if(ring == 3 && mod == 355480616) { mapRing = 4; mapModule = 355484716;}
      if(ring == 3 && mod == 355480624) { mapRing = 4; mapModule = 355484732;}
      if(ring == 3 && mod == 355480632) { mapRing = 4; mapModule = 355484740;}
      if(ring == 3 && mod == 355480640) { mapRing = 4; mapModule = 355484748;}
      if(ring == 3 && mod == 355480648) { mapRing = 4; mapModule = 355484756;}
      if(ring == 3 && mod == 355480656) { mapRing = 4; mapModule = 355484772;}
      if(ring == 3 && mod == 355480664) { mapRing = 4; mapModule = 355484780;}
      if(ring == 3 && mod == 355480672) { mapRing = 4; mapModule = 355484788;}
      if(ring == 3 && mod == 355480680) { mapRing = 4; mapModule = 355484796;}
      if(ring == 3 && mod == 355480688) { mapRing = 4; mapModule = 355484804;}
      if(ring == 3 && mod == 355480696) { mapRing = 4; mapModule = 355484820;}
      if(ring == 3 && mod == 355480704) { mapRing = 4; mapModule = 355484828;}
      if(ring == 3 && mod == 355480712) { mapRing = 4; mapModule = 355484836;}
      if(ring == 3 && mod == 355480720) { mapRing = 4; mapModule = 355484844;}
      
      
      
      if(ring == 3 && mod == 355480584) { mapRing = 4; mapModule = 355484680;}
      if(ring == 3 && mod == 355480592) { mapRing = 4; mapModule = 355484688;}
      if(ring == 3 && mod == 355480600) { mapRing = 4; mapModule = 355484704;}
      if(ring == 3 && mod == 355480608) { mapRing = 4; mapModule = 355484712;}
      if(ring == 3 && mod == 355480616) { mapRing = 4; mapModule = 355484720;}
      if(ring == 3 && mod == 355480624) { mapRing = 4; mapModule = 355484728;}
      if(ring == 3 && mod == 355480632) { mapRing = 4; mapModule = 355484736;}
      if(ring == 3 && mod == 355480640) { mapRing = 4; mapModule = 355484752;}
      if(ring == 3 && mod == 355480648) { mapRing = 4; mapModule = 355484760;}
      if(ring == 3 && mod == 355480656) { mapRing = 4; mapModule = 355484768;}
      if(ring == 3 && mod == 355480664) { mapRing = 4; mapModule = 355484776;}
      if(ring == 3 && mod == 355480672) { mapRing = 4; mapModule = 355484792;}
      if(ring == 3 && mod == 355480680) { mapRing = 4; mapModule = 355484800;}
      if(ring == 3 && mod == 355480688) { mapRing = 4; mapModule = 355484808;}
      if(ring == 3 && mod == 355480696) { mapRing = 4; mapModule = 355484816;}
      if(ring == 3 && mod == 355480704) { mapRing = 4; mapModule = 355484824;}
      if(ring == 3 && mod == 355480712) { mapRing = 4; mapModule = 355484840;}
      if(ring == 3 && mod == 355480720) { mapRing = 4; mapModule = 355484848;}
      
      
      if(ring == 5 && mod == 355488772) { mapRing = 4; mapModule = 355484848;}
      if(ring == 5 && mod == 355488780) { mapRing = 4; mapModule = 355484680;}
      if(ring == 5 && mod == 355488788) { mapRing = 4; mapModule = 355484688;}
      if(ring == 5 && mod == 355488796) { mapRing = 4; mapModule = 355484696;}
      if(ring == 5 && mod == 355488804) { mapRing = 4; mapModule = 355484704;}
      if(ring == 5 && mod == 355488812) { mapRing = 4; mapModule = 355484712;}
      if(ring == 5 && mod == 355488820) { mapRing = 4; mapModule = 355484720;}
      if(ring == 5 && mod == 355488828) { mapRing = 4; mapModule = 355484728;}
      if(ring == 5 && mod == 355488836) { mapRing = 4; mapModule = 355484736;}
      if(ring == 5 && mod == 355488844) { mapRing = 4; mapModule = 355484744;}
      if(ring == 5 && mod == 355488852) { mapRing = 4; mapModule = 355484752;}
      if(ring == 5 && mod == 355488860) { mapRing = 4; mapModule = 355484760;}
      if(ring == 5 && mod == 355488868) { mapRing = 4; mapModule = 355484760;}
      if(ring == 5 && mod == 355488876) { mapRing = 4; mapModule = 355484768;}
      if(ring == 5 && mod == 355488884) { mapRing = 4; mapModule = 355484776;}
      if(ring == 5 && mod == 355488892) { mapRing = 4; mapModule = 355484784;}
      if(ring == 5 && mod == 355488900) { mapRing = 4; mapModule = 355484792;}
      if(ring == 5 && mod == 355488908) { mapRing = 4; mapModule = 355484800;}
      if(ring == 5 && mod == 355488916) { mapRing = 4; mapModule = 355484808;}
      if(ring == 5 && mod == 355488924) { mapRing = 4; mapModule = 355484816;}
      if(ring == 5 && mod == 355488932) { mapRing = 4; mapModule = 355484824;}
      if(ring == 5 && mod == 355488940) { mapRing = 4; mapModule = 355484832;}
      if(ring == 5 && mod == 355488948) { mapRing = 4; mapModule = 355484840;}
      if(ring == 5 && mod == 355488956) { mapRing = 4; mapModule = 355484848;}
      
      
      
      
      if(ring == 4 && mod == 355484680) { mapRing = 5; mapModule = 355488776;}
      if(ring == 4 && mod == 355484688) { mapRing = 5; mapModule = 355488784;}
      if(ring == 4 && mod == 355484696) { mapRing = 5; mapModule = 355488792;}
      if(ring == 4 && mod == 355484704) { mapRing = 5; mapModule = 355488800;}
      if(ring == 4 && mod == 355484712) { mapRing = 5; mapModule = 355488808;}
      if(ring == 4 && mod == 355484720) { mapRing = 5; mapModule = 355488816;}
      if(ring == 4 && mod == 355484728) { mapRing = 5; mapModule = 355488824;}
      if(ring == 4 && mod == 355484736) { mapRing = 5; mapModule = 355488840;}
      if(ring == 4 && mod == 355484744) { mapRing = 5; mapModule = 355488848;}
      if(ring == 4 && mod == 355484752) { mapRing = 5; mapModule = 355488856;}
      if(ring == 4 && mod == 355484760) { mapRing = 5; mapModule = 355488864;}
      if(ring == 4 && mod == 355484768) { mapRing = 5; mapModule = 355488872;}
      if(ring == 4 && mod == 355484776) { mapRing = 5; mapModule = 355488880;}
      if(ring == 4 && mod == 355484784) { mapRing = 5; mapModule = 355488888;}
      if(ring == 4 && mod == 355484792) { mapRing = 5; mapModule = 355488896;}
      if(ring == 4 && mod == 355484800) { mapRing = 5; mapModule = 355488904;}
      if(ring == 4 && mod == 355484808) { mapRing = 5; mapModule = 355488920;}
      if(ring == 4 && mod == 355484816) { mapRing = 5; mapModule = 355488928;}
      if(ring == 4 && mod == 355484824) { mapRing = 5; mapModule = 355488936;}
      if(ring == 4 && mod == 355484832) { mapRing = 5; mapModule = 355488944;}
      if(ring == 4 && mod == 355484840) { mapRing = 5; mapModule = 355488952;}
      if(ring == 4 && mod == 355484848) { mapRing = 5; mapModule = 355488960;}
      
      
      
      if(ring == 5 && mod == 355488776) { mapRing = 4; mapModule = 355484676;}
      if(ring == 5 && mod == 355488784) { mapRing = 4; mapModule = 355484684;}
      if(ring == 5 && mod == 355488792) { mapRing = 4; mapModule = 355484692;}
      if(ring == 5 && mod == 355488800) { mapRing = 4; mapModule = 355484700;}
      if(ring == 5 && mod == 355488808) { mapRing = 4; mapModule = 355484708;}
      if(ring == 5 && mod == 355488816) { mapRing = 4; mapModule = 355484716;}
      if(ring == 5 && mod == 355488824) { mapRing = 4; mapModule = 355484724;}
      if(ring == 5 && mod == 355488832) { mapRing = 4; mapModule = 355484732;}
      if(ring == 5 && mod == 355488840) { mapRing = 4; mapModule = 355484740;}
      if(ring == 5 && mod == 355488848) { mapRing = 4; mapModule = 355484748;}
      if(ring == 5 && mod == 355488856) { mapRing = 4; mapModule = 355484756;}
      if(ring == 5 && mod == 355488864) { mapRing = 4; mapModule = 355484764;}
      if(ring == 5 && mod == 355488872) { mapRing = 4; mapModule = 355484764;}
      if(ring == 5 && mod == 355488880) { mapRing = 4; mapModule = 355484772;}
      if(ring == 5 && mod == 355488888) { mapRing = 4; mapModule = 355484780;}
      if(ring == 5 && mod == 355488896) { mapRing = 4; mapModule = 355484788;}
      if(ring == 5 && mod == 355488904) { mapRing = 4; mapModule = 355484796;}
      if(ring == 5 && mod == 355488912) { mapRing = 4; mapModule = 355484804;}
      if(ring == 5 && mod == 355488920) { mapRing = 4; mapModule = 355484812;}
      if(ring == 5 && mod == 355488928) { mapRing = 4; mapModule = 355484820;}
      if(ring == 5 && mod == 355488936) { mapRing = 4; mapModule = 355484828;}
      if(ring == 5 && mod == 355488944) { mapRing = 4; mapModule = 355484836;}
      if(ring == 5 && mod == 355488952) { mapRing = 4; mapModule = 355484844;}
      if(ring == 5 && mod == 355488960) { mapRing = 4; mapModule = 355484844;}
     
     
    }
  
  }
  
  
  //uint32_t tmpid = getModuleID(side, disk, mapRing); 
  //uint32_t longid = (tmpid & 0xFFFFFC03) | (((mapModule) & 0xFF) << 2);
  //std::cout<< longid << std::endl;
  //return longid;
    
  return mapModule;
}


DEFINE_FWK_MODULE(ITclusterAnalyzerCoincidences);
