// -*- C++ -*-
//
// Package:    BRIL_ITsim/ITclusterAnalyzer
// Class:      Ashish2xCoincidence
//
/**\class Ashish2xCoincidence Ashish2xCoincidence.cc BRIL_ITsim/ITclusterAnalyzer/plugins/Ashish2xCoincidence.cc
   Description: [one line class summary]
   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Ashish Sehrawat
//         Created:  Tue, 20 Sep 2020 13:41:06 GMT
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

// a struct to hold the residuals for each matched cluster
struct Residual {
  
  double dX;
  double dY;
  double dR;
  
  Residual(double X, double Y) : dX(X), dY(Y) {
    dR = sqrt(pow(dX, 2) + pow(dY, 2));
  }
  
  void print() {
    //std::cout << "(dX: " << dX << " dY: " << dY << " dR: " << dR << ") ";
  }
  
};

struct Residual1{

  double x_1;
  double x_2;
  double y_1;
  double y_2;
  
  double r_1;
  double r_2;
  double dr;

  Residual1(double x_1, double x_2, double y_1, double y_2) {
    
    r_1 = sqrt(pow(x_1, 2) + pow(y_1, 2));
    r_2 = sqrt(pow(x_2, 2) + pow(y_2, 2));
    dr = sqrt(pow(x_2, 2) + pow(y_2, 2)) - sqrt(pow(x_1, 2) + pow(y_1, 2));
 
    //std::cout <<" r_1 is: " << r_1 << " r_2 is: " << r_2 << " dr is: " << dr;
 
  }
};


struct deltaphiset{
  
  double phi_1;
  double phi_2;
  double deltaphi;
  
  
  deltaphiset(double x_1, double x_2, double y_1, double y_2, double_t phi_1, double_t phi_2){
    
    phi_1 = TMath::ATan2(y_2, x_2);
    phi_2 = TMath::ATan2(y_1, x_1);
    deltaphi = phi_2 - phi_1;
    
    //std::cout << " " << phi_1 << "  " << phi_2 << "  " << deltaphi; 
    
    
  }
  
};

class Ashish2xCoincidence : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Ashish2xCoincidence(const edm::ParameterSet&);
  ~Ashish2xCoincidence();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //bool findCoincidence(DetId, Global3DPoint, bool);
  const SiPixelCluster* findCoincidence2x(DetId, Global3DPoint, bool, unsigned int&);
  bool findCoincidenceInR2x(DetId, Global3DPoint, bool, std::vector<unsigned int>&, std::vector<edmNew::DetSet<SiPixelCluster>::const_iterator>&, std::vector<float>&);
  edm::DetSetVector<PixelDigiSimLink>::const_iterator findSimLinkDetSet(unsigned int thedetid);
  std::set<unsigned int> getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator, edmNew::DetSet<SiPixelCluster>::const_iterator, bool print);
  std::set<unsigned int> getSimTrackId1(edm::DetSetVector<PixelDigiSimLink>::const_iterator, edmNew::DetSet<SiPixelCluster>::const_iterator, bool print1);
  bool areSameSimTrackId(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>&);
  bool areSameSimTrackId1(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>&);
 
 uint32_t getModuleID(bool, unsigned int, unsigned int, unsigned int);
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
  TH2F* m_diskHistosCluster[8];
  
  //tracker maps for clusters
  TH2F* m_trackerLayoutClustersZR;
  TH2F* m_trackerLayoutClustersYX;
  
  
  //array of TH2F for 2xcoinc per disk per ring
  //first all coincidences
  TH2F* m_diskHistos2x[8];
  TH2F* m_diskHistos2xInR[8];
  //and the real ones
  TH2F* m_diskHistos2xreal[8];
  TH2F* m_diskHistos2xrealInR[8];
  //tracker maps for 2xcoinc
  TH2F* m_trackerLayout2xZR;
  TH2F* m_trackerLayout2xYX;
  TH2F* m_trackerLayout2xZR_InR;
  TH2F* m_trackerLayout2xYX_InR;
  
  //Residual histograms for all rings
  TH1F* m_residualX[5];
  TH1F* m_residualX_sametrack[5];
  TH1F* m_residualX_notsametrack[5];
  
  TH1F* m_residualY[5];
  TH1F* m_residualY_sametrack[5];
  TH1F* m_residualY_notsametrack[5];
  
  TH1F* m_residualR[5];
  TH1F* m_residualR_sametrack[5];
  TH1F* m_residualR_notsametrack[5];
  
  TH1F* m_residualR1[5];
  TH1F* m_residualR1_sametrack[5];
  TH1F* m_residualR1_notsametrack[5];
  
  
  TH1F* m_deltaphi[5];
  TH1F* m_deltaphi_sametrack[5];
  TH1F* m_deltaphi_notsametrack[5];
  
  
  
  
  //simple residual histograms for the cuts
  TH1F* m_residualX_Ring1;
  TH1F* m_residualY_Ring1;
  TH1F* m_residualR_Ring1;
  TH1F* m_residualR1_Ring1;
  TH1F* m_deltaphi_Ring1;
  
  
  TH1F* m_residualX_Ring2;
  TH1F* m_residualY_Ring2;
  TH1F* m_residualR_Ring2;
  TH1F* m_residualR1_Ring2;
  TH1F* m_deltaphi_Ring2;


  TH1F* m_residualX_Ring3;
  TH1F* m_residualY_Ring3;
  TH1F* m_residualR_Ring3;
  TH1F* m_residualR1_Ring3;
  TH1F* m_deltaphi_Ring3;

  TH1F* m_residualX_Ring4;
  TH1F* m_residualY_Ring4;
  TH1F* m_residualR_Ring4;
  TH1F* m_residualR1_Ring4;
  TH1F* m_deltaphi_Ring4;
  
  
  TH1F* m_residualX_Ring5;
  TH1F* m_residualY_Ring5;
  TH1F* m_residualR_Ring5;
  TH1F* m_residualR1_Ring5;
  TH1F* m_deltaphi_Ring5;
  
  
  
  TH1F* m_residualX_Ring1_sametrack;
  TH1F* m_residualY_Ring1_sametrack;
  TH1F* m_residualR_Ring1_sametrack;
  TH1F* m_residualR1_Ring1_sametrack;
  TH1F* m_deltaphi_Ring1_sametrack;
  

  TH1F* m_residualX_Ring2_sametrack;
  TH1F* m_residualY_Ring2_sametrack;
  TH1F* m_residualR_Ring2_sametrack;
  TH1F* m_residualR1_Ring2_sametrack;
  TH1F* m_deltaphi_Ring2_sametrack;
  
  TH1F* m_residualX_Ring3_sametrack;
  TH1F* m_residualY_Ring3_sametrack;
  TH1F* m_residualR_Ring3_sametrack;
  TH1F* m_residualR1_Ring3_sametrack;
  TH1F* m_deltaphi_Ring3_sametrack;
  
  TH1F* m_residualX_Ring4_sametrack;
  TH1F* m_residualY_Ring4_sametrack;
  TH1F* m_residualR_Ring4_sametrack;
  TH1F* m_residualR1_Ring4_sametrack;
  TH1F* m_deltaphi_Ring4_sametrack;
  
  TH1F* m_residualX_Ring5_sametrack;
  TH1F* m_residualY_Ring5_sametrack;
  TH1F* m_residualR_Ring5_sametrack;
  TH1F* m_residualR1_Ring5_sametrack;
  TH1F* m_deltaphi_Ring5_sametrack;
  
  TH1F* m_residualX_Ring1_notsametrack;
  TH1F* m_residualY_Ring1_notsametrack;
  TH1F* m_residualR_Ring1_notsametrack;
  TH1F* m_residualR1_Ring1_notsametrack;
  TH1F* m_deltaphi_Ring1_notsametrack;
  
  TH1F* m_residualX_Ring2_notsametrack;
  TH1F* m_residualY_Ring2_notsametrack;
  TH1F* m_residualR_Ring2_notsametrack;
  TH1F* m_residualR1_Ring2_notsametrack;
  TH1F* m_deltaphi_Ring2_notsametrack;

  TH1F* m_residualX_Ring3_notsametrack;
  TH1F* m_residualY_Ring3_notsametrack;
  TH1F* m_residualR_Ring3_notsametrack;
  TH1F* m_residualR1_Ring3_notsametrack;
  TH1F* m_deltaphi_Ring3_notsametrack;
  
  TH1F* m_residualX_Ring4_notsametrack;
  TH1F* m_residualY_Ring4_notsametrack;
  TH1F* m_residualR_Ring4_notsametrack;
  TH1F* m_residualR1_Ring4_notsametrack;
  TH1F* m_deltaphi_Ring4_notsametrack;
  
  TH1F* m_residualX_Ring5_notsametrack;
  TH1F* m_residualY_Ring5_notsametrack;
  TH1F* m_residualR_Ring5_notsametrack;
  TH1F* m_residualR1_Ring5_notsametrack;
  TH1F* m_deltaphi_Ring5_notsametrack;
  
  
  TH1F* m_residualX_InR;
  TH1F* m_residualY_InR;
  TH1F* m_residualR_InR;
  

  //the number of clusters per module
  TH1F* m_nClusters;
  
  //cuts for the coincidence
  
  double m_dr[5];
  double m_dphi[5];
  
  double m_dx;
  double m_dy;
  double m_dz;
  double m_dr_ring1;
  double m_dphi_ring1;
  double m_dr_ring2;
  double m_dphi_ring2;
  double m_dr_ring3;
  double m_dphi_ring3;
  double m_dr_ring4;
  double m_dphi_ring4;
  double m_dr_ring5;
  double m_dphi_ring5;
  
  
  //event counter
  uint32_t m_nevents;
  
  //coincidence counter
  uint32_t m_total2xcoincidences_ring1;
  uint32_t m_total2xcoincidences_ring2;
  uint32_t m_total2xcoincidences_ring3;
  uint32_t m_total2xcoincidences_ring4;
  uint32_t m_total2xcoincidences_ring5;
  
  uint32_t m_total2xcoincidencesInR;
  
  uint32_t m_true2xcoincidences_ring1;
  uint32_t m_true2xcoincidences_ring2;
  uint32_t m_true2xcoincidences_ring3;
  uint32_t m_true2xcoincidences_ring4;
  uint32_t m_true2xcoincidences_ring5;
  
  uint32_t m_fake2xcoincidences_ring1;
  uint32_t m_fake2xcoincidences_ring2;
  uint32_t m_fake2xcoincidences_ring3;
  uint32_t m_fake2xcoincidences_ring4;
  uint32_t m_fake2xcoincidences_ring5;
  
  uint32_t m_fake2xcoincidencesInR;
  
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
Ashish2xCoincidence::Ashish2xCoincidence(const edm::ParameterSet& iConfig)
  : //m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>> ("clusters"))
  m_tokenClusters(consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters")))
  , m_tokenSimLinks(consumes<edm::DetSetVector<PixelDigiSimLink>>(iConfig.getParameter<edm::InputTag>("simlinks")))
  , m_tokenDigis(consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("digis"))) //adding digis variable - COB 26.02.19
  , m_maxBin(iConfig.getUntrackedParameter<uint32_t>("maxBin"))
  , m_docoincidence(iConfig.getUntrackedParameter<bool>("docoincidence"))
  , m_dx(iConfig.getParameter<double>("dx_cut"))
  , m_dy(iConfig.getParameter<double>("dy_cut"))
  , m_dz(iConfig.getParameter<double>("dz_cut")) 
  , m_dr_ring1(iConfig.getParameter<double>("dr_cut_ring1"))
  , m_dphi_ring1(iConfig.getParameter<double>("dphi_cut_ring1"))
  , m_dr_ring2(iConfig.getParameter<double>("dr_cut_ring2"))
  , m_dphi_ring2(iConfig.getParameter<double>("dphi_cut_ring2"))
  , m_dr_ring3(iConfig.getParameter<double>("dr_cut_ring3"))
  , m_dphi_ring3(iConfig.getParameter<double>("dphi_cut_ring3"))
  , m_dr_ring4(iConfig.getParameter<double>("dr_cut_ring4"))
  , m_dphi_ring4(iConfig.getParameter<double>("dphi_cut_ring4"))
  , m_dr_ring5(iConfig.getParameter<double>("dr_cut_ring5"))
  , m_dphi_ring5(iConfig.getParameter<double>("dphi_cut_ring5")){
  
  
  
  //now do what ever initialization is needed
  m_nevents = 0;
  m_total2xcoincidences_ring1 = 0;
  m_total2xcoincidences_ring2 = 0;
  m_total2xcoincidences_ring3 = 0;
  m_total2xcoincidences_ring4 = 0;
  m_total2xcoincidences_ring5 = 0;
  
  m_total2xcoincidencesInR = 0;
  
  m_fake2xcoincidences_ring1 = 0;
  m_fake2xcoincidences_ring2 = 0;
  m_fake2xcoincidences_ring3 = 0;
  m_fake2xcoincidences_ring4 = 0;
  m_fake2xcoincidences_ring5 = 0;
  
  m_true2xcoincidences_ring1 = 0;
  m_true2xcoincidences_ring2 = 0;
  m_true2xcoincidences_ring3 = 0;
  m_true2xcoincidences_ring4 = 0;
  m_true2xcoincidences_ring5 = 0;
  
  m_fake2xcoincidencesInR = 0;
  
}
Ashish2xCoincidence::~Ashish2xCoincidence() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void Ashish2xCoincidence::beginJob() {
  
  edm::Service<TFileService> fs;

  fs->file().cd("/");
  TFileDirectory td = fs->mkdir("TEPX");
  
  td = fs->mkdir("TEPX/Residuals");
  
  
  m_residualX_Ring1 = td.make<TH1F>("ResidualsX_Ring1", "ResidualsX_Ring1;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring1 = td.make<TH1F>("ResidualsY_Ring1", "ResidualsY_Ring1;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring1 = td.make<TH1F>("ResidualsR_Ring1", "ResidualsR_Ring1;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring1 = td.make<TH1F>("ResidualsR1_Ring1", "ResidualsR1_Ring1;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring1 = td.make<TH1F>("Deltaphi_Ring1", "Deltaphi_Ring1;phi2-phi1;counts", 1000, -1, 1);
  
  
  
  m_residualX_Ring1_sametrack = td.make<TH1F>("ResidualsX_Ring1_sametrack", "ResidualsX_Ring1_sametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring1_sametrack = td.make<TH1F>("ResidualsY_Ring1_sametrack", "ResidualsY_Ring1_sametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring1_sametrack = td.make<TH1F>("ResidualsR_Ring1_sametrack", "ResidualsR_Ring1_sametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring1_sametrack = td.make<TH1F>("ResidualsR1_Ring1_sametrack", "ResidualsR1_Ring1_sametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring1_sametrack = td.make<TH1F>("Deltaphi_Ring1_sametrack", "Deltaphi_Ring1_sametrack;phi2-phi1;counts", 1000, -1, 1);
  
  
  m_residualX_Ring1_notsametrack = td.make<TH1F>("ResidualsX_Ring1_notsametrack", "ResidualsX_Ring1_notsametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring1_notsametrack = td.make<TH1F>("ResidualsY_Ring1_notsametrack", "ResidualsY_Ring1_notsametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring1_notsametrack = td.make<TH1F>("ResidualsR_Ring1_notsametrack", "ResidualsR_Ring1_notsametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring1_notsametrack = td.make<TH1F>("ResidualsR1_Ring1_notsametrack", "ResidualsR1_Ring1_notsametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring1_notsametrack = td.make<TH1F>("Deltaphi_Ring1_notsametrack", "Deltaphi_Ring1_notsametrack;phi2-phi1;counts", 1000, -1, 1);
  
  
   
  m_residualX_Ring2 = td.make<TH1F>("ResidualsX_Ring2", "ResidualsX_Ring2;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring2 = td.make<TH1F>("ResidualsY_Ring2", "ResidualsY_Ring2;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring2 = td.make<TH1F>("ResidualsR_Ring2", "ResidualsR_Ring2;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring2 = td.make<TH1F>("ResidualsR1_Ring2", "ResidualsR1_Ring2;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring2 = td.make<TH1F>("Deltaphi_Ring2", "Deltaphi_Ring2;phi2-phi1;counts", 1000, -1, 1);
  
  m_residualX_Ring2_sametrack = td.make<TH1F>("ResidualsX_Ring2_sametrack", "ResidualsX_Ring2_sametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring2_sametrack = td.make<TH1F>("ResidualsY_Ring2_sametrack", "ResidualsY_Ring2_sametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring2_sametrack = td.make<TH1F>("ResidualsR_Ring2_sametrack", "ResidualsR_Ring2_sametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring2_sametrack = td.make<TH1F>("ResidualsR1_Ring2_sametrack", "ResidualsR1_Ring2_sametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring2_sametrack = td.make<TH1F>("Deltaphi_Ring2_sametrack", "Deltaphi_Ring2_sametrack;phi2-phi1;counts", 1000, -1, 1);
  

  m_residualX_Ring2_notsametrack = td.make<TH1F>("ResidualsX_Ring2_notsametrack", "ResidualsX_Ring2_notsametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring2_notsametrack = td.make<TH1F>("ResidualsY_Ring2_notsametrack", "ResidualsY_Ring2_notsametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring2_notsametrack = td.make<TH1F>("ResidualsR_Ring2_notsametrack", "ResidualsR_Ring2_notsametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring2_notsametrack = td.make<TH1F>("ResidualsR1_Ring2_notsametrack", "ResidualsR1_Ring2_notsametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring2_notsametrack = td.make<TH1F>("Deltaphi_Ring2_notsametrack", "Deltaphi_Ring2_notsametrack;phi2-phi1;counts", 1000, -1, 1);
      
  
  m_residualX_Ring3 = td.make<TH1F>("ResidualsX_Ring3", "ResidualsX_Ring3;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring3 = td.make<TH1F>("ResidualsY_Ring3", "ResidualsY_Ring3;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring3 = td.make<TH1F>("ResidualsR_Ring3", "ResidualsR_Ring3;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring3 = td.make<TH1F>("ResidualsR1_Ring3", "ResidualsR1_Ring3;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring3 = td.make<TH1F>("Deltaphi_Ring3", "Deltaphi_Ring3;phi2-phi1;counts", 1000, -1, 1);
  
  m_residualX_Ring3_sametrack = td.make<TH1F>("ResidualsX_Ring3_sametrack", "ResidualsX_Ring3_sametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring3_sametrack = td.make<TH1F>("ResidualsY_Ring3_sametrack", "ResidualsY_Ring3_sametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring3_sametrack = td.make<TH1F>("ResidualsR_Ring3_sametrack", "ResidualsR_Ring3_sametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring3_sametrack = td.make<TH1F>("ResidualsR1_Ring3_sametrack", "ResidualsR1_Ring3_sametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring3_sametrack = td.make<TH1F>("Deltaphi_Ring3_sametrack", "Deltaphi_Ring3_sametrack;phi2-phi1;counts", 1000, -1, 1);
  
  m_residualX_Ring3_notsametrack = td.make<TH1F>("ResidualsX_Ring3_notsametrack", "ResidualsX_Ring3_notsametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring3_notsametrack = td.make<TH1F>("ResidualsY_Ring3_notsametrack", "ResidualsY_Ring3_notsametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring3_notsametrack = td.make<TH1F>("ResidualsR_Ring3_notsametrack", "ResidualsR_Ring3_notsametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring3_notsametrack = td.make<TH1F>("ResidualsR1_Ring3_notsametrack", "ResidualsR1_Ring3_notsametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring3_notsametrack = td.make<TH1F>("Deltaphi_Ring3_notsametrack", "Deltaphi_Ring3_notsametrack;phi2-phi1;counts", 1000, -1, 1);
  
  
  
  m_residualX_Ring4 = td.make<TH1F>("ResidualsX_Ring4", "ResidualsX_Ring4;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring4 = td.make<TH1F>("ResidualsY_Ring4", "ResidualsY_Ring4;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring4 = td.make<TH1F>("ResidualsR_Ring4", "ResidualsR_Ring4;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring4 = td.make<TH1F>("ResidualsR1_Ring4", "ResidualsR1_Ring4;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring4 = td.make<TH1F>("Deltaphi_Ring4", "Deltaphi_Ring4;phi2-phi1;counts", 1000, -1, 1);
  
  m_residualX_Ring4_sametrack = td.make<TH1F>("ResidualsX_Ring4_sametrack", "ResidualsX_Ring4_sametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring4_sametrack = td.make<TH1F>("ResidualsY_Ring4_sametrack", "ResidualsY_Ring4_sametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring4_sametrack = td.make<TH1F>("ResidualsR_Ring4_sametrack", "ResidualsR_Ring4_sametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring4_sametrack = td.make<TH1F>("ResidualsR1_Ring4_sametrack", "ResidualsR1_Ring4_sametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring4_sametrack = td.make<TH1F>("Deltaphi_Ring4_sametrack", "Deltaphi_Ring4_sametrack;phi2-phi1;counts", 1000, -1, 1);
  
  m_residualX_Ring4_notsametrack = td.make<TH1F>("ResidualsX_Ring4_notsametrack", "ResidualsX_Ring4_notsametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring4_notsametrack = td.make<TH1F>("ResidualsY_Ring4_notsametrack", "ResidualsY_Ring4_notsametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring4_notsametrack = td.make<TH1F>("ResidualsR_Ring4_notsametrack", "ResidualsR_Ring4_notsametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring4_notsametrack = td.make<TH1F>("ResidualsR1_Ring4_notsametrack", "ResidualsR1_Ring4_notsametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring4_notsametrack = td.make<TH1F>("Deltaphi_Ring4_notsametrack", "Deltaphi_Ring4_notsametrack;phi2-phi1;counts", 1000, -1, 1);
  
  
  m_residualX_Ring5 = td.make<TH1F>("ResidualsX_Ring5", "ResidualsX_Ring5;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring5 = td.make<TH1F>("ResidualsY_Ring5", "ResidualsY_Ring5;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring5 = td.make<TH1F>("ResidualsR_Ring5", "ResidualsR_Ring5;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring5 = td.make<TH1F>("ResidualsR1_Ring5", "ResidualsR1_Ring5;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring5 = td.make<TH1F>("Deltaphi_Ring5", "Deltaphi_Ring5;phi2-phi1;counts", 1000, -1, 1);
  
  m_residualX_Ring5_sametrack = td.make<TH1F>("ResidualsX_Ring5_sametrack", "ResidualsX_Ring5_sametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring5_sametrack = td.make<TH1F>("ResidualsY_Ring5_sametrack", "ResidualsY_Ring5_sametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring5_sametrack = td.make<TH1F>("ResidualsR_Ring5_sametrack", "ResidualsR_Ring5_sametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring5_sametrack = td.make<TH1F>("ResidualsR1_Ring5_sametrack", "ResidualsR1_Ring5_sametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring5_sametrack = td.make<TH1F>("Deltaphi_Ring5_sametrack", "Deltaphi_Ring5_sametrack;phi2-phi1;counts", 1000, -1, 1);
  
  m_residualX_Ring5_notsametrack = td.make<TH1F>("ResidualsX_Ring5_notsametrack", "ResidualsX_Ring5_notsametrack;x2-x1 (cm);counts", 1000, -1, 1);
  m_residualY_Ring5_notsametrack = td.make<TH1F>("ResidualsY_Ring5_notsametrack", "ResidualsY_Ring5_notsametrack;y2-y1 (cm);counts", 1000, -1, 1);
  m_residualR_Ring5_notsametrack = td.make<TH1F>("ResidualsR_Ring5_notsametrack", "ResidualsR_Ring5_notsametrack;deltaR (cm);counts", 1000, 0, 1);
  m_residualR1_Ring5_notsametrack = td.make<TH1F>("ResidualsR1_Ring5_notsametrack", "ResidualsR1_Ring5_notsametrack;r2-r1 (cm);counts", 1000, -1, 1);
  m_deltaphi_Ring5_notsametrack = td.make<TH1F>("Deltaphi_Ring5_notsametrack", "Deltaphi_Ring5_notsametrack;phi2-phi1;counts", 1000, -1, 1);    
  
  
  
  
  m_residualX_InR = td.make<TH1F>("ResidualsXInR", "ResidualsX;deltaX (cm);counts", 1000, -1, 1);
  m_residualY_InR = td.make<TH1F>("ResidualsYInR", "ResidualsY;deltaY (cm);counts", 1000, -1, 1);
  m_residualR_InR = td.make<TH1F>("ResidualsRInR", "ResidualsR;deltaR (cm);counts", 1000, 0, 1);
  
  fs->file().cd("/");
  td = fs->mkdir("TEPX/perModule");
  m_nClusters = td.make<TH1F>("Number of Clusters per module per event", "# of Clusters;# of Clusters; Occurence", 500, 0, 500);
  
  
  fs->file().cd("/");
  td = fs->mkdir("TEPX/Clusters");
  
  //now lets create the histograms
  for (unsigned int i = 0; i < 8; i++) {
    int disk = (i < 4) ? i - 4 : i - 3;
    std::stringstream histoname;
    histoname << "Number of clusters for Disk " << disk << ";Ring;# of Clusters per event";
    std::stringstream histotitle;
    histotitle << "Number of clusters for Disk " << disk;
    //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
    m_diskHistosCluster[i] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);
  }
  m_trackerLayoutClustersZR = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
  m_trackerLayoutClustersYX = td.make<TH2F>("XVsY", "x vs. y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    
  
  
  if (m_docoincidence) {
    fs->file().cd("/");
    td = fs->mkdir("TEPX/2xCoincidences");
    //now lets create the histograms
    for (unsigned int i = 0; i < 8; i++) {
      int disk = (i < 4) ? i - 4 : i - 3;
      std::stringstream histoname;
      histoname << "Number of 2x Coincidences for Disk " << disk << ";Ring;# of coincidences per event";
      std::stringstream histotitle;
      histotitle << "Number of 2x Coincidences for Disk " << disk;
      //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
      m_diskHistos2x[i] = td.make<TH2F>(histotitle.str().c_str(), histoname.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);
      
      std::stringstream histonamereal;
      histonamereal << "Number of real 2x Coincidences for Disk " << disk << ";Ring;# of real coincidences per event";
      std::stringstream histotitlereal;
      histotitlereal << "Number of real 2x Coincidences for Disk " << disk;
      //name, name, nbinX, Xlow, Xhigh, nbinY, Ylow, Yhigh
      m_diskHistos2xreal[i] = td.make<TH2F>(histotitlereal.str().c_str(), histonamereal.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);
      
      std::stringstream histonameInR;
      histonameInR << "Number of 2x Coincidences in R for Disk " << disk << ";Ring;# of coincidences per event";
      std::stringstream histotitleInR;
      histotitleInR << "Number of 2x Coincidences in R for Disk " << disk;
      m_diskHistos2xInR[i] = td.make<TH2F>(histotitleInR.str().c_str(), histonameInR.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);

      std::stringstream histonamerealInR;
      histonamerealInR << "Number of real 2x Coincidences in R for Disk " << disk << ";Ring;# of coincidences per event";
      std::stringstream histotitlerealInR;
      histotitlerealInR << "Number of real 2x Coincidences in R for Disk " << disk;
      m_diskHistos2xrealInR[i] = td.make<TH2F>(histotitlerealInR.str().c_str(), histonamerealInR.str().c_str(), 5, .5, 5.5, m_maxBin, 0, m_maxBin);
    }
    m_trackerLayout2xZR = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
    m_trackerLayout2xYX = td.make<TH2F>("XVsY", "X vs. Y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    m_trackerLayout2xZR_InR = td.make<TH2F>("RVsZ_InR", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
    m_trackerLayout2xYX_InR = td.make<TH2F>("XVsY_InR", "x vs. y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
  }
}






// ------------ method called for each event  ------------
void Ashish2xCoincidence::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
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
  unsigned int cluCounter[8][5];
  memset(cluCounter, 0, sizeof(cluCounter));
  
  //counter for 2x coincidences
  unsigned int x2Counter[8][5];
  memset(x2Counter, 0, sizeof(x2Counter));
  unsigned int x2CounterInR[8][5];
  memset(x2CounterInR, 0, sizeof(x2CounterInR));
  unsigned int x2Counterreal[8][5];
  memset(x2Counterreal, 0, sizeof(x2Counterreal));
  unsigned int x2CounterrealInR[8][5];
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
    //std::cout << "DetID " << std::hex << "0x" << rawid << std::dec << " " << detId.det() << " " 
    //       << detId.subdetId() << " " << ((rawid >> 23) & 0x3) << " " << ((rawid >> 18) & 0xF) << " " 
    //       << ((rawid >> 12) & 0x3F) << " " << ((rawid >> 2) & 0xFF) << std::endl;
    
    //find out which layer, side and ring
    unsigned int side = (tTopo->pxfSide(detId));  // values are 1 and 2 for -+Z
    unsigned int layer = (tTopo->pxfDisk(detId)); //values are 1 to 12 for disks TFPX1 to TFPX 8  and TEPX1 to TEPX 4
    unsigned int ring = (tTopo->pxfBlade(detId));
    
    if (layer > 8) { // TEPX modules
      
      //the index in my histogram map
      int hist_id = -1;
      unsigned int ring_id = ring - 1;
      
      if (side == 1) {
	//this is a TEPX- hit on side1
	hist_id = layer - 9;
      } else if (side == 2) {
	//this is a TEPX+ hit on side 2
	hist_id = 4 + layer - 9;
      }
      
      // Get the geomdet
      const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(detId));
      if (!geomDetUnit)
	continue;
      
      unsigned int nClu = 0;
	
      //fill the number of clusters for this module
      m_nClusters->Fill(DSVit->size());
      
      //now loop the clusters for each detector
      for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = DSVit->begin(); cluit != DSVit->end(); cluit++) {
	//increment the counters
	nClu++;
	cluCounter[hist_id][ring_id]++;
	
	// determine the position
	MeasurementPoint mpClu(cluit->x(), cluit->y());
	Local3DPoint localPosClu = geomDetUnit->topology().localPosition(mpClu);
	Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);
	
	
	m_trackerLayoutClustersZR->Fill(globalPosClu.z(), globalPosClu.perp());
	m_trackerLayoutClustersYX->Fill(globalPosClu.x(), globalPosClu.y());
	
	
	if (m_docoincidence) {
	  if(((rawid >> 2) & 0xFF) % 2 != 0){
	    
	    unsigned int coincidenceId;
	    const SiPixelCluster* found2xcoincidencecluster = this->findCoincidence2x(detId, globalPosClu, true, coincidenceId);
	    if (found2xcoincidencecluster) {
	      
	      x2Counter[hist_id][ring_id]++;
	      
	      
	      //const GeomDetUnit* geomDetUnit_2xcluster(tkGeom->idToDetUnit(coincidenceId));
	      
	      
	      //MeasurementPoint mpcoincidenceClu(found2xcoincidencecluster->x(), found2xcoincidencecluster->y());
	      //Local3DPoint thelocalPosClu = geomDetUnit_2xcluster->topology().localPosition(mpcoincidenceClu);
	      //Global3DPoint theglobalPosClu = geomDetUnit_2xcluster->surface().toGlobal(thelocalPosClu);
	      
	      //std::cout << " " << theglobalPosClu.x() << "  " << theglobalPosClu.y() << std::endl;
	      
	      
	      //now get the simlink detset
	      //edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = findSimLinkDetSet(rawid);
	      //std::set<unsigned int> simTrackId = this->getSimTrackId(simLinkDSViter, cluit, false);
	      //now get the simlink detset based on the coincidence hit detid
	      //simLinkDSViter = findSimLinkDetSet(coincidenceId);
	      //std::set<unsigned int> coincidencesimTrackId = this->getSimTrackId(simLinkDSViter, found2xcoincidencecluster, false);
	      //std::set<unsigned int> intersection;
	      //bool areSame = areSameSimTrackId(simTrackId, coincidencesimTrackId, intersection);
	      
	      
	      
	      //double delta_X = - globalPosClu.x() + theglobalPosClu.x();
	      //double delta_Y = - globalPosClu.y() + theglobalPosClu.y();
	      
	      
	      m_trackerLayout2xZR->Fill(globalPosClu.z(), globalPosClu.perp());
	      m_trackerLayout2xYX->Fill(globalPosClu.x(), globalPosClu.y());
	      
	      
	      
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
  for (unsigned int i = 0; i < 8; i++) {  // TEPX
    //loop the disks
    for (unsigned int j = 0; j < 5; j++) {
      //and the rings
      m_diskHistosCluster[i]->Fill(j + 1, cluCounter[i][j]);
    //  m_diskHistosHits[i]->Fill(j + 1, hitCounter[i][j]);
      if (m_docoincidence) {
	m_diskHistos2x[i]->Fill(j + 1, x2Counter[i][j]);
	m_diskHistos2xreal[i]->Fill(j + 1, x2Counterreal[i][j]);
	m_diskHistos2xInR[i]->Fill(j + 1, x2CounterInR[i][j]);
	m_diskHistos2xrealInR[i]->Fill(j + 1, x2CounterrealInR[i][j]);
	
      }
    }
 }
  
  m_nevents++;
  
}


// ------------ method called once each job just after ending the event loop  ------------
void Ashish2xCoincidence::endJob() {
  std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
  if (m_docoincidence) {
    //std::cout << "IT cluster Analyzer found " << m_true2xcoincidences_ring1 / (double)m_total2xcoincidences_ring1 * 100 
    //	      << "\% true double coincidences in TEPX Ring 1 modules." << std::endl;
    //std::cout << "IT cluster Analyzer found " << m_fake2xcoincidences_ring1 / (double)m_total2xcoincidences_ring1 * 100
    //        << "\% fake double coincidences in TEPX Ring 1 modules." << std::endl;
    
    //std::cout << "IT cluster Analyzer found " << m_true2xcoincidences_ring2 / (double)m_total2xcoincidences_ring2 * 100
    //        << "\% true double coincidences in TEPX Ring 2 modules." << std::endl;
    //std::cout << "IT cluster Analyzer found " << m_fake2xcoincidences_ring2 / (double)m_total2xcoincidences_ring2 * 100
    //        << "\% fake double coincidences in TEPX  Ring 2 modules." << std::endl;
    
    //std::cout << "IT cluster Analyzer found " << m_true2xcoincidences_ring3 / (double)m_total2xcoincidences_ring3 * 100
    //        << "\% true double coincidences in TEPX Ring 3 modules." << std::endl;
    //std::cout << "IT cluster Analyzer found " << m_fake2xcoincidences_ring3 / (double)m_total2xcoincidences_ring3 * 100
    //        << "\% fake double coincidences in TEPX Ring 3 modules." << std::endl;
    
    //std::cout << "IT cluster Analyzer found " << m_true2xcoincidences_ring4 / (double)m_total2xcoincidences_ring4 * 100
    //        << "\% true double coincidences in TEPX Ring 4 modules." << std::endl;
    //std::cout << "IT cluster Analyzer found " << m_fake2xcoincidences_ring4 / (double)m_total2xcoincidences_ring4 * 100
    //        << "\% fake double coincidences in TEPX Ring 4 modules." << std::endl;
    
    //std::cout << "IT cluster Analyzer found " << m_true2xcoincidences_ring5 / (double)m_total2xcoincidences_ring5 * 100
    //         << "\% true double coincidences in TEPX Ring 5 modules." << std::endl;
    //std::cout << "IT cluster Analyzer found " << m_fake2xcoincidences_ring5 / (double)m_total2xcoincidences_ring5 * 100
    //        << "\% fake double coincidences in TEPX Ring 5 modules." << std::endl;
    
    
  }
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Ashish2xCoincidence::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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


const SiPixelCluster* Ashish2xCoincidence::findCoincidence2x(DetId thedetid, Global3DPoint theglobalPosClu, bool isTEPX, unsigned int& foundDetId) {
  
  const SiPixelCluster* found2xcoincidencecluster = NULL; 
  //bool found = false;
  uint32_t rawid = thedetid.rawId();
  uint32_t newid = rawid;
  uint32_t newid1 = rawid;
  
  //now I have the raw ID and can mess with the bits
  //the side, layer and ring are the same and I just have to increment or decrement the module number
  unsigned int themodule = (tTopo->pxfModule(thedetid));
  unsigned int thering = (tTopo->pxfBlade(thedetid));
  //unsigned int thedisk = (tTopo->pxfDisk(thedetid));                                                                                    
  //unsigned int theside = (tTopo->pxfSide(thedetid));
  
  //const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(thedetid));                                                                         
  //std::pair<float,float> phiSpan = geomDetUnit->surface().phiSpan();                                                                     
  //std::pair<float,float> zSpan = geomDetUnit->surface().zSpan();                                                                         
  //std::pair<float,float> rSpan = geomDetUnit->surface().rSpan();  
  
  //in order to avoid duplicates, only look in the next module clockwise
  //depending on the ring, if I am already on the module with the highest module id in the ring, I need to go to the first one
  uint32_t newmodule = themodule + 1;
  uint32_t newmodule1 = themodule - 1;
  
  if (isTEPX) {
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
  
  if(isTEPX){
    if (thering == 1 && themodule == 1)
      newmodule1 = 20;
    if (thering == 2 && themodule == 1)
      newmodule1 = 28;
    if (thering == 3 && themodule == 1)
      newmodule1 = 36;
    if (thering == 4 && themodule == 1)
      newmodule1 = 44;
    if (thering == 5 && themodule == 1)
      newmodule1 = 48;
  }
  
  
  //now encode
  newid = (newid & 0xFFFFFC03) | ((newmodule & 0xFF) << 2);
  newid1 = (newid1 & 0xFFFFFC03) | ((newmodule1 & 0xFF) << 2);
  
  //now I have a raw id of the module I want to use
  DetId id(newid);
  DetId id1(newid1);
  
  //unsigned int ring = (tTopo->pxfBlade(id));
  //unsigned int module = (tTopo->pxfModule(id));
  
  edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(id);
  edmNew::DetSetVector<SiPixelCluster>::const_iterator theit1 = clusters->find(id1);
  
  if (theit == clusters->end()) {
    // return false;
    return found2xcoincidencecluster;   
  }
  
  if (theit1 == clusters->end()) {
    // return false;                                                                                                                       
    return found2xcoincidencecluster;
  }
  
  // Get the geomdet
  const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(id));
  const GeomDetUnit* geomDetUnit1(tkGeom->idToDetUnit(id1));
  //std::pair<float,float> phiSpan_new = geomDetUnit_new->surface().phiSpan();                                                            
  //std::pair<float,float> zSpan_new = geomDetUnit_new->surface().zSpan();                                                                 
  //std::pair<float,float> rSpan_new = geomDetUnit_new->surface().rSpan();
  
 
  unsigned int nClu = 0;
  //at the end of the day, need to find the closest coincidence hit, so store the minimum 2D distance in a temporary variable and a vector for all values
  double R_min = 1000.;
  //double r_min = 1000.;
  //double phi_min = 1000.;
  
  std::vector<Residual> r_vec;
  std::vector<Residual1>r_vec1;
  std::vector<deltaphiset>phiangle;
  
  //make the return value end();
  //found2xcoincidencecluster = theit->end();
  //found2xcoincidencecluster = theit1->end();
  
  
  for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit = theit->begin(); cluit != theit->end(); cluit++) {
    
    // determine the position
    MeasurementPoint mpClu(cluit->x(), cluit->y());
    Local3DPoint localPosClu = geomDetUnit->topology().localPosition(mpClu);
    Global3DPoint globalPosClu = geomDetUnit->surface().toGlobal(localPosClu);

    
    edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = findSimLinkDetSet(rawid);
    std::set<unsigned int> simTrackId = this->getSimTrackId(simLinkDSViter, cluit, false);
    //now get the simlink detset based on the coincidence hit detid                                                                             
    simLinkDSViter = findSimLinkDetSet(newid);
    std::set<unsigned int> coincidencesimTrackId = this->getSimTrackId(simLinkDSViter, found2xcoincidencecluster, false);
    std::set<unsigned int> intersection;
    bool areSame = areSameSimTrackId(simTrackId, coincidencesimTrackId, intersection);
    
    
    
    double dr = sqrt(pow(theglobalPosClu.x(), 2) + pow(theglobalPosClu.y(), 2)) - sqrt(pow(globalPosClu.x(), 2) + pow(globalPosClu.y(), 2));
    
    
    double phi1 = TMath::ATan2(globalPosClu.y(), globalPosClu.x());
    double phi2 = TMath::ATan2(theglobalPosClu.y(), theglobalPosClu.x());
    
     
    if((phi2-phi1) < m_dphi_ring1 && dr < m_dr_ring1 && thering == 1){
      
      nClu++;
      
      double delta_X = - globalPosClu.x() + theglobalPosClu.x();
      double delta_Y = - globalPosClu.y() + theglobalPosClu.y();
      
      double x_1 = theglobalPosClu.x();
      double x_2 = globalPosClu.x();
      double y_1 = theglobalPosClu.y();
      double y_2 = globalPosClu.y();
      
      double phi1 = TMath::ATan2(globalPosClu.y(), globalPosClu.x());                                                                   
      double phi2 = TMath::ATan2(theglobalPosClu.y(), theglobalPosClu.x());                                                             
      double deltaphi = phi2 - phi1;

      double dr = sqrt(pow(theglobalPosClu.x(), 2) + pow(theglobalPosClu.y(), 2)) - sqrt(pow(globalPosClu.x(), 2) + pow(globalPosClu.y(), 2));      

      double dR = sqrt(pow(theglobalPosClu.x() - globalPosClu.x(), 2) + pow(theglobalPosClu.y() - globalPosClu.y(), 2));

      deltaphiset clusterphiangle(x_1, x_2, y_1, y_2, phi1, phi2);
      phiangle.push_back(clusterphiangle);
      
      Residual r(delta_X, delta_Y);
      r_vec.push_back(r);
      
      Residual1 r1(x_1, x_2, y_1, y_2);
      r_vec1.push_back(r1);
      
      if (r.dR < R_min) {
	//if (r1.dr < r_min){
	//if(clusterphiangle.deltaphi < phi_min){
	R_min = r.dR;
	//   r_min = r1.dr;
	// phi_min = clusterphiangle.deltaphi;
	//found = true;
	// I assign this here to always have the closest cluster
	//foundCluster = cluit;
	
	foundDetId = newid;
	found2xcoincidencecluster = cluit; 
	//	}
	//	}
      }

      //std::cout << "Found matching cluster # " << nClu << std::endl;
      
      //std::cout << "Original x: " << theglobalPosClu.x() << " y: " << theglobalPosClu.y() << " z: " << theglobalPosClu.z() << " ring: " << thering << " module: " << themodule << std::endl;
      //std::cout << "New      x: " << globalPosClu.x() << " y: " << globalPosClu.y() << " z: " << globalPosClu.z() << " ring: " << ring << " module: " << module << std::endl;
      // }
      //  }
      // }
      //	}

    

        if (nClu > 1) {
      //if(areSame){
      //std::cout << "Warning, found " << nClu << "Clusters within the cuts - the minimum distance is " << r_min << "!" << std::endl;
      //std::cout << "All distances: ";
      //for (auto clusterphiangle : phiangle){
	//for (auto r : r_vec) {
	//for (auto r1 : r_vec1){
      //r.print();
      //if (r.dr == r_min) {
      //  if (isTEPX) {


      m_residualX_Ring1->Fill(delta_X);
      m_residualY_Ring1->Fill(delta_Y);
      m_residualR_Ring1->Fill(dR);
      m_residualR1_Ring1->Fill(dr);
      m_deltaphi_Ring1 -> Fill(deltaphi);
      
      // } 
      // }
      // }
      
      
      if(areSame){
	
	//for (auto clusterphiangle : phiangle){
	//for (auto r : r_vec) {
	// for (auto r1 : r_vec1){
	
	m_residualX_Ring1_sametrack->Fill(delta_X);
	m_residualY_Ring1_sametrack->Fill(delta_Y);
	m_residualR_Ring1_sametrack->Fill(dR);
	m_residualR1_Ring1_sametrack->Fill(dr);
	m_deltaphi_Ring1_sametrack->Fill(deltaphi);

      }
      // }
      //	}
      // }
      
      
      if(!areSame){
	
	// for (auto clusterphiangle : phiangle){
	// for (auto r : r_vec) {
	//for (auto r1 : r_vec1){
	
	 m_residualX_Ring1_notsametrack->Fill(delta_X);
	 m_residualY_Ring1_notsametrack->Fill(delta_Y);
         m_residualR_Ring1_notsametrack->Fill(dR);
         m_residualR1_Ring1_notsametrack->Fill(dr);
	 m_deltaphi_Ring1_notsametrack->Fill(deltaphi);
	      
      }
       }
      // }
      
      
      
      // }
      //std::cout << std::endl;
      // }
      // }
      // }
      //}
      
      
    }
    

 


  
    
    
    
  }
  
  
  
  
  
  
  // return found;
  return found2xcoincidencecluster;
}


edm::DetSetVector<PixelDigiSimLink>::const_iterator Ashish2xCoincidence::findSimLinkDetSet(unsigned int thedetid) {
  ////basic template
  edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDS = simlinks->find(thedetid);
  return simLinkDS;
}

std::set<unsigned int> Ashish2xCoincidence::getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter, edmNew::DetSet<SiPixelCluster>::const_iterator cluster, bool print) {
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
  //if(simTrackIds.size() != 1){
  //std::cout << "WARNING: have more than 1 simTrackId for this cluster! " << std::endl;
  //return 0;
  //}
  return simTrackIds;
}








bool Ashish2xCoincidence::areSameSimTrackId(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>& intersection) {
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










//----------
//Adding function to construct DetId from PXFDetId
//COB - 22.May.2019
uint32_t Ashish2xCoincidence::getModuleID(bool isTEPX, unsigned int side, unsigned int disk, unsigned int ring) {
  
  //std::cout << "isTEPX " << isTEPX << " side " << side << " disk " << disk << " ring " << ring << std::endl;

  uint32_t modid = -999;
  
  if (isTEPX) {
    if (side==1) {
      if (disk==9) {
	if (ring==2) {
	  modid = 346301444;
	} else if (ring==3) {
	  modid = 346305540;
	} else if (ring==4) {
	  modid = 346309636;
	} else if (ring==5) {
	  modid = 346313732;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << " side " << side << "!" << std::endl;
	  return modid;
	}
      } else if (disk==10) {
	if (ring==2) {
	  modid = 346563588;
	} else if (ring==3) {
	  modid = 346567684;
	} else if (ring==4) {
	  modid = 346571780;
	} else if (ring==5) {
	  modid = 346575876;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
	  return modid;
	}
      } else if (disk==11) {
	if (ring==2) {
	  modid = 346825732;
	} else if (ring==3) {
	  modid = 346829828;
	} else if (ring==4) {
	  modid = 346833924;
	} else if (ring==5) {
	  modid = 346838020;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
	  return modid;
	}
      } else if (disk==12) {
	if (ring==2) {
	  modid = 347087876;
	} else if (ring==3) {
	  modid = 347091972;
	} else if (ring==4) {
	  modid = 347096068;
	} else if (ring==5) {
	  modid = 347100164;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
	  return modid;
	}
      } else {
	std::cout << "Non-existent disk for TEPX!" << std::endl;
	return modid;
      }
    } else if (side==2) {
      if (disk==9) {
	if (ring==2) {
	  modid = 354690052;
	} else if (ring==3) {
	  modid = 354694148;
	} else if (ring==4) {
	  modid = 354698244;
	} else if (ring==5) {
	  modid = 354702340;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
	  return modid;
	}
      } else if (disk==10) {
	if (ring==2) {
	  modid = 354952196;
	} else if (ring==3) {
	  modid = 354956292;
	} else if (ring==4) {
	  modid = 354960388;
	} else if (ring==5) {
	  modid = 354964484;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
	  return modid;
	}
      } else if (disk==11) {
	if (ring==2) {
	  modid = 355214340;
	} else if (ring==3) {
	  modid = 355218436;
	} else if (ring==4) {
	  modid = 355222532;
	} else if (ring==5) {
	  modid = 355226628;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
	  return modid;
	}
      } else if (disk==12) {
	if (ring==2) {
	  modid = 355476484;
	} else if (ring==3) {
	  modid = 355480580;
	} else if (ring==4) {
	  modid = 355484676;
	} else if (ring==5) {
	  modid = 355488772;
	} else {
	  std::cout << "Non-existent ring number for disk " << disk << "!" << std::endl;
	  return modid;
	}
      } else {
	std::cout << "Non-existent disk for TEPX!" << std::endl;
	return modid;
      }
    } else {
      std::cout << "Non-existent side!" << std::endl;
      return modid;
    }
  }


  return modid;

}

DEFINE_FWK_MODULE(Ashish2xCoincidence);
