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
//         Created:  Fri, 23 Oct 2020 19:00:06 GMT
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

  //and the real ones
  TH2F* m_diskHistos2xreal[2][4];
  
  //tracker maps for 2xcoinc
  TH2F* m_trackerLayout2xZR;
  TH2F* m_trackerLayout2xYX;
  
  
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
	histoname << "Numberof2xCoincidences_side" << side << "_Disk" << disk;
	
	//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
	m_diskHistos2x[k][i] = td.make<TH2F>(histoname.str().c_str()," ",  5, .5, 5.5, m_maxBin, 0, m_maxBin);
	histoname.str("");
	
	histoname << "Numberofreal2xCoincidences_side" << side << "_Disk" << disk;
	
	//name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
	m_diskHistos2xreal[k][i] = td.make<TH2F>(histoname.str().c_str()," ", 5, .5, 5.5, m_maxBin, 0, m_maxBin);
	histoname.str("");
	
      }
    }
  
  
    m_trackerLayout2xZR = td.make<TH2F>("2xCoincidences_RVsZ", "RvsZposition_2xCoincidences", 6000, -300.0, 300.0, 600, 0.0, 30.0);
    m_trackerLayout2xYX = td.make<TH2F>("2xCoincidences_XVsY", "XvsYposition_2xCoincidences", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    
    
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
      
      if(side == 1 && 1<= ring && ring <= 4 && module % 2 != 0) disk_layer_type = 2;
      if(side == 1 && 1<= ring && ring <= 4 && module % 2 == 0) disk_layer_type = 1;
      if(side == 1 && ring == 5 && module % 2 != 0) disk_layer_type = 1;
      if(side == 1 && ring == 5 && module % 2 == 0) disk_layer_type = 2;
      
      if(side == 2 && 1<= ring && ring <= 4 && module % 2 == 0) disk_layer_type = 2;
      if(side == 2 && 1<= ring && ring <= 4 && module % 2 != 0) disk_layer_type = 1;
      if(side == 2 && ring == 5 && module % 2 != 0) disk_layer_type = 2;
      if(side == 2 && ring == 5 && module % 2 == 0) disk_layer_type = 1;
      
      
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
	if (disk_layer_type == 1) {   
	  
	  //increment the counters
	  nClu++;
	  cluCounter[side-1][disk-1][ring-1]++;
	  
	  // determine the position
	  MeasurementPoint mpClu1(cluit1->x(), cluit1->y());
	  Local3DPoint localPosClu1 = geomDetUnit->topology().localPosition(mpClu1);
	  Global3DPoint globalPosClu1 = geomDetUnit->surface().toGlobal(localPosClu1);
	  
	  m_trackerLayoutClustersZR->Fill(globalPosClu1.z(), globalPosClu1.perp());
	  m_trackerLayoutClustersYX->Fill(globalPosClu1.x(), globalPosClu1.y());
	  
	}
	
	if (m_docoincidence) {
	  if (disk_layer_type == 2) {
	    
	    MeasurementPoint mpClu1(cluit1->x(), cluit1->y());
	    Local3DPoint localPosClu1 = geomDetUnit->topology().localPosition(mpClu1);
	    Global3DPoint globalPosClu1 = geomDetUnit->surface().toGlobal(localPosClu1);
	    
	    unsigned int coincidenceId;
	    
	    const SiPixelCluster* found2xcoincidencecluster = this->findCoincidence2x(detId, globalPosClu1, coincidenceId, cluit1, 1);
	    
	    if (!found2xcoincidencecluster) 
	      found2xcoincidencecluster = this->findCoincidence2x(detId, globalPosClu1, coincidenceId, cluit1, 2);
	    
	    if (found2xcoincidencecluster) {
	      
	      x2Counter[side-1][disk-1][ring-1]++;
	  
	      //if(side == 2 && disk == 4 && ring == 1) {    
	      m_trackerLayout2xZR->Fill(globalPosClu1.z(), globalPosClu1.perp());
	      m_trackerLayout2xYX->Fill(globalPosClu1.x(), globalPosClu1.y());
	      //}
	      
	    }
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
      //double dR = sqrt(pow(globalPosClu2.x() - globalPosClu1.x(), 2) + pow(globalPosClu2.y() - globalPosClu1.y(), 2));
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


DEFINE_FWK_MODULE(ITclusterAnalyzerCoincidences);
