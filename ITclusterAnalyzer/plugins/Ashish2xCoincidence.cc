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
//         Created:  Tue, 5 Oct 2020 16:11:06 GMT
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



class Ashish2xCoincidence5 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Ashish2xCoincidence5(const edm::ParameterSet&);
  ~Ashish2xCoincidence5();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //bool findCoincidence(DetId, Global3DPoint, bool);
  const SiPixelCluster* findCoincidence2x(DetId, Global3DPoint, bool, unsigned int&, edmNew::DetSet<SiPixelCluster>::const_iterator, unsigned int);
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
  
  //simple residual histograms for the cuts
  TH1F* m_dX[8][5];
  TH1F* m_dY[8][5];
  TH1F* m_dR[8][5];
  TH1F* m_dr[8][5];
  TH1F* m_deltaphi[8][5];
  
  TH1F* m_dX_sametrack[8][5];
  TH1F* m_dY_sametrack[8][5];
  TH1F* m_dR_sametrack[8][5];
  TH1F* m_dr_sametrack[8][5];
  TH1F* m_deltaphi_sametrack[8][5];
  
  TH1F* m_dX_notsametrack[8][5];
  TH1F* m_dY_notsametrack[8][5];
  TH1F* m_dR_notsametrack[8][5];
  TH1F* m_dr_notsametrack[8][5];
  TH1F* m_deltaphi_notsametrack[8][5];
  
  //the number of clusters per module
  TH1F* m_nClusters;
  
  //cuts for the coincidence
  double m_dx;
  double m_dy;
  double m_dz;
  double m_dr_ring1;
  double m_dphi_ring1;
  
  double m_dr_cuts[8][5];
  double m_dphi_cuts[8][5]; 
  
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
Ashish2xCoincidence5::Ashish2xCoincidence5(const edm::ParameterSet& iConfig)
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
  , m_dphi_ring1(iConfig.getParameter<double>("dphi_cut_ring1")) {
  
  
  //now do what ever initialization is needed
  m_nevents = 0;
  
}
Ashish2xCoincidence5::~Ashish2xCoincidence5() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called once each job just before starting event loop  ------------
void Ashish2xCoincidence5::beginJob() {
  
  edm::Service<TFileService> fs;
  
  fs->file().cd("/");
  TFileDirectory td = fs->mkdir("TEPX");
  
  td = fs->mkdir("TEPX/Residuals");
  
  for (unsigned int i = 0; i < 8; i++) {
    for (unsigned j = 0; j < 5; j++) {
      int disk = (i < 4) ? i - 4 : i - 3;
      int ring = (j < 5) ? j+1 : j+1;
      
      std::stringstream histoname;
      histoname << "m_dX for Disk " << disk << " Ring " << ring << ";x2-x1 (cm);counts";
      std::stringstream histotitle;
      histotitle << "m_dX for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname1;
      histoname1 << "m_dY for Disk " << disk << " Ring " << ring << ";y2-y1 (cm);counts";
      std::stringstream histotitle1;
      histotitle1 << "m_dY for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname2;
      histoname2 << "m_dR for Disk " << disk << " Ring " << ring << ";deltaR (cm);counts";
      std::stringstream histotitle2;
      histotitle2 << "m_dR for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname3;
      histoname3 << "m_dr for Disk " << disk << " Ring " << ring << ";r2-r1 (cm);counts";
      std::stringstream histotitle3;
      histotitle3 << "m_dr for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname4;
      histoname4 << "m_deltaphi for Disk " << disk << " Ring " << ring << ";phi2-phi;counts";
      std::stringstream histotitle4;
      histotitle4 << "m_deltaphi for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname5;
      histoname5 << "m_dX_sametrack for Disk " << disk << " Ring " << ring << ";x2-x1 (cm);counts";
      std::stringstream histotitle5;
      histotitle5 << "m_dX_sametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname6;
      histoname6 << "m_dY_sametrack for Disk " << disk << " Ring " << ring << ";y2-y1 (cm);counts";
      std::stringstream histotitle6;
      histotitle6 << "m_dY_sametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname7;
      histoname7 << "m_dR_sametrack for Disk " << disk << " Ring " << ring << ";deltaR (cm);counts";
      std::stringstream histotitle7;
      histotitle7 << "m_dR_sametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname8;
      histoname8 << "m_dr_sametrack for Disk " << disk << " Ring " << ring << ";r2-r1 (cm);counts";
      std::stringstream histotitle8;
      histotitle8 << "m_dr_sametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname9;
      histoname9 << "m_deltaphi_sametrack for Disk " << disk << " Ring " << ring << ";phi2-phi1 ;counts";
      std::stringstream histotitle9;
      histotitle9 << "m_deltaphi_sametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname10;
      histoname10 << "m_dX_notsametrack for Disk " << disk << " Ring " << ring << ";x2-x1 (cm);counts";
      std::stringstream histotitle10;
      histotitle10 << "m_dX_notsametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname11;
      histoname11 << "m_dY_notsametrack for Disk " << disk << " Ring " << ring << ";y2-y1 (cm);counts";
      std::stringstream histotitle11;
      histotitle11 << "m_dY_notsametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname12;
      histoname12 << "m_dR_notsametrack for Disk " << disk << " Ring " << ring << ";deltaR (cm);counts";
      std::stringstream histotitle12;
      histotitle12 << "m_dR_notsametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname13;
      histoname13 << "m_dr_notsametrack for Disk " << disk << " Ring " << ring << ";r2-r1 (cm);counts";
      std::stringstream histotitle13;
      histotitle13 << "m_dr_notsametrack for Disk " << disk << " Ring " << ring;
      
      
      std::stringstream histoname14;
      histoname14 << "m_deltaphi_notsametrack for Disk " << disk << " Ring " << ring << ";phi2-phi1;counts";
      std::stringstream histotitle14;
      histotitle14 << "m_deltaphi_notsametrack for Disk " << disk << " Ring " << ring;
      
      
      
      m_dX[i][j] = td.make<TH1F>(histotitle.str().c_str(), histoname.str().c_str(), 1000, -1, 1);
      m_dY[i][j] = td.make<TH1F>(histotitle1.str().c_str(), histoname1.str().c_str(), 1000, -1, 1);
      m_dR[i][j] = td.make<TH1F>(histotitle2.str().c_str(), histoname2.str().c_str(), 1000, 0, 1);
      m_dr[i][j] = td.make<TH1F>(histotitle3.str().c_str(), histoname3.str().c_str(), 1000, -0.3, 0.3);
      m_deltaphi[i][j] = td.make<TH1F>(histotitle4.str().c_str(), histoname4.str().c_str(), 1000, -0.2, 0.2);
      
      m_dX_sametrack[i][j] = td.make<TH1F>(histotitle5.str().c_str(), histoname5.str().c_str(), 1000, -1, 1);
      m_dY_sametrack[i][j] = td.make<TH1F>(histotitle6.str().c_str(), histoname6.str().c_str(), 1000, -1, 1);
      m_dR_sametrack[i][j] = td.make<TH1F>(histotitle7.str().c_str(), histoname7.str().c_str(), 1000, 0, 1);
      m_dr_sametrack[i][j] = td.make<TH1F>(histotitle8.str().c_str(), histoname8.str().c_str(), 1000, -0.3, 0.3);
      m_deltaphi_sametrack[i][j] = td.make<TH1F>(histotitle9.str().c_str(), histoname9.str().c_str(), 1000, -0.2, 0.2);
      
      m_dX_notsametrack[i][j] = td.make<TH1F>(histotitle10.str().c_str(), histoname10.str().c_str(), 1000, -1, 1);
      m_dY_notsametrack[i][j] = td.make<TH1F>(histotitle11.str().c_str(), histoname11.str().c_str(), 1000, -1, 1);
      m_dR_notsametrack[i][j] = td.make<TH1F>(histotitle12.str().c_str(), histoname12.str().c_str(), 1000, 0, 1);
      m_dr_notsametrack[i][j] = td.make<TH1F>(histotitle13.str().c_str(), histoname13.str().c_str(), 1000, -0.3, 0.3);
      m_deltaphi_notsametrack[i][j] = td.make<TH1F>(histotitle14.str().c_str(), histoname14.str().c_str(), 1000, -0.2, 0.2);    
      
    }
  }
  
  
  
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
            
    }
    
    m_trackerLayout2xZR = td.make<TH2F>("RVsZ", "R vs. z position", 6000, -300.0, 300.0, 600, 0.0, 30.0);
    m_trackerLayout2xYX = td.make<TH2F>("XVsY", "X vs. Y position", 1000, -50.0, 50.0, 1000, -50.0, 50.0);
    
  }
}



// ------------ method called for each event  ------------
void Ashish2xCoincidence5::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
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
  
  unsigned int x2Counterreal[8][5];
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
      for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit1 = DSVit->begin(); cluit1 != DSVit->end(); cluit1++) {
	//increment the counters
	nClu++;
	cluCounter[hist_id][ring_id]++;
	
	//determine the position
	MeasurementPoint mpClu1(cluit1->x(), cluit1->y());
	Local3DPoint localPosClu1 = geomDetUnit->topology().localPosition(mpClu1);
	Global3DPoint globalPosClu1 = geomDetUnit->surface().toGlobal(localPosClu1);
	
	
	m_trackerLayoutClustersZR->Fill(globalPosClu1.z(), globalPosClu1.perp());
	m_trackerLayoutClustersYX->Fill(globalPosClu1.x(), globalPosClu1.y());
	
	
	if (m_docoincidence) {
	  if(((rawid >> 2) & 0xFF) % 2 != 0){
	    
	    unsigned int coincidenceId;
	    unsigned int overlappingmodule = 1;
	    
	    const SiPixelCluster* found2xcoincidencecluster = this->findCoincidence2x(detId, globalPosClu1, true, coincidenceId, cluit1, overlappingmodule);
	    
	    if (found2xcoincidencecluster) {
	      
	      x2Counter[hist_id][ring_id]++;
	      
	      m_trackerLayout2xZR->Fill(globalPosClu1.z(), globalPosClu1.perp());
	      m_trackerLayout2xYX->Fill(globalPosClu1.x(), globalPosClu1.y());
	      
	      
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
      if (m_docoincidence) {
	m_diskHistos2x[i]->Fill(j + 1, x2Counter[i][j]);
	m_diskHistos2xreal[i]->Fill(j + 1, x2Counterreal[i][j]);
       	
      }
    }
  }
  
  m_nevents++;
  
}


// ------------ method called once each job just after ending the event loop  ------------
void Ashish2xCoincidence5::endJob() {
  
  std::cout << "IT cluster Analyzer processed " << m_nevents << " events!" << std::endl;
  
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Ashish2xCoincidence5::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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


const SiPixelCluster* Ashish2xCoincidence5::findCoincidence2x(DetId thedetid, Global3DPoint globalPosClu1, bool isTEPX, unsigned int& foundDetId, edmNew::DetSet<SiPixelCluster>::const_iterator cluit1, unsigned int overlappingmodule) {
  
  const SiPixelCluster* found2xcoincidencecluster = NULL; 
  uint32_t rawid = thedetid.rawId();
  uint32_t newid = rawid;
  
  
  //now I have the raw ID and can mess with the bits
  //the side, layer and ring are the same and I just have to increment or decrement the module number
  
  unsigned int themodule = (tTopo->pxfModule(thedetid));
  unsigned int thering = (tTopo->pxfBlade(thedetid));
  unsigned int thelayer = (tTopo->pxfDisk(thedetid));
  unsigned int theside = (tTopo->pxfSide(thedetid));  
  
  //std::cout << thelayer << std::endl;
  //std::cout << theside << std::endl;
  
  
  //in order to avoid duplicates, only look in the next module clockwise
  //depending on the ring, if I am already on the module with the highest module id in the ring, I need to go to the first one
  uint32_t newmodule = themodule + 1;
  
  
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
  
  
  //now encode
  newid = (newid & 0xFFFFFC03) | ((newmodule & 0xFF) << 2);
  
  
  //now I have a raw id of the module I want to use
  DetId id(newid);
  
  edmNew::DetSetVector<SiPixelCluster>::const_iterator theit = clusters->find(id);
  
  if (theit == clusters->end()) {
    return found2xcoincidencecluster;   
  }
  
  
  // Get the geomdet
  const GeomDetUnit* geomDetUnit(tkGeom->idToDetUnit(id));
  
  
  //at the end of the day, need to find the closest coincidence hit, so store the minimum 2D distance in a temporary variable and a vector for all values
  double R_min = 1000.;
  
  
  for (edmNew::DetSet<SiPixelCluster>::const_iterator cluit2 = theit->begin(); cluit2 != theit->end(); cluit2++) {
    
    // determine the position
    MeasurementPoint mpClu2(cluit2->x(), cluit2->y());
    Local3DPoint localPosClu2 = geomDetUnit->topology().localPosition(mpClu2);
    Global3DPoint globalPosClu2 = geomDetUnit->surface().toGlobal(localPosClu2);    
    
    
    //double r1 = sqrt(pow(globalPosClu1.x(), 2) + pow(globalPosClu1.y(), 2));
    //double r2 = sqrt(pow(globalPosClu2.x(), 2) + pow(globalPosClu2.y(), 2));
    
    //double phi1 = TMath::ATan2(globalPosClu1.y(), globalPosClu1.x());
    //double phi2 = TMath::ATan2(globalPosClu2.y(), globalPosClu2.x());
    
    
    if (fabs(globalPosClu1.x() - globalPosClu2.x()) < m_dx
        && fabs(globalPosClu1.y() - globalPosClu2.y()) < m_dy
        && fabs(globalPosClu1.z() - globalPosClu2.z()) < m_dz) {
      
      
      //if((phi2-phi1) < m_dphi && (r2-r1) < m_dr && thering == 1) {
      
      
      double dX = - globalPosClu1.x() + globalPosClu2.x();
      double dY = - globalPosClu1.y() + globalPosClu2.y();
      
      double dr = sqrt(pow(globalPosClu2.x(), 2) + pow(globalPosClu2.y(), 2)) - sqrt(pow(globalPosClu1.x(), 2) + pow(globalPosClu1.y(), 2));
      double dR = sqrt(pow(globalPosClu2.x() - globalPosClu1.x(), 2) + pow(globalPosClu2.y() - globalPosClu1.y(), 2));
      
      double phi1 = TMath::ATan2(globalPosClu1.y(), globalPosClu1.x());        
      double phi2 = TMath::ATan2(globalPosClu2.y(), globalPosClu2.x()); 
      
      if (dR < R_min) {
	
	R_min = dR;
       	foundDetId = newid;
	found2xcoincidencecluster = cluit2; 
	
      }
      
      
      edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter = findSimLinkDetSet(rawid);                                       
      std::set<unsigned int> simTrackId = this->getSimTrackId(simLinkDSViter, cluit1, false);                                              
      
      //now get the simlink detset based on the coincidence hit detid	
      simLinkDSViter = findSimLinkDetSet(newid);                                                                                           
      std::set<unsigned int> coincidencesimTrackId = this->getSimTrackId(simLinkDSViter, cluit2, false);                
      std::set<unsigned int> intersection;                                                                                                 
      bool areSame = areSameSimTrackId(simTrackId, coincidencesimTrackId, intersection);   
      
      
      if (theside == 1) {
	
	
	m_dX[thelayer- 9][thering - 1] -> Fill(dX);
	m_dY[thelayer -9][thering - 1] -> Fill(dY);
	m_dR[thelayer - 9][thering - 1] -> Fill(dR);
	m_dr[thelayer - 9][thering - 1] -> Fill(dr);
	m_deltaphi[thelayer - 9][thering - 1] -> Fill(phi2-phi1);
	
	
	if(areSame) {
	  
	  
	  m_dX_sametrack[thelayer - 9][thering - 1] -> Fill(dX);
	  m_dY_sametrack[thelayer - 9][thering - 1] -> Fill(dY);
	  m_dR_sametrack[thelayer - 9][thering - 1] -> Fill(dR);
	  m_dr_sametrack[thelayer - 9][thering - 1] -> Fill(dr);
	  m_deltaphi_sametrack[thelayer - 9][thering - 1] -> Fill(phi2-phi1);
	  
	  
	}
	
	
	else if(!areSame) {
	  
	  
	  m_dX_notsametrack[thelayer - 9][thering - 1] -> Fill(dX);
	  m_dY_notsametrack[thelayer - 9][thering - 1] -> Fill(dY);
	  m_dR_notsametrack[thelayer - 9][thering - 1] -> Fill(dR);
	  m_dr_notsametrack[thelayer - 9][thering - 1] -> Fill(dr);
	  m_deltaphi_notsametrack[thelayer - 9][thering - 1] -> Fill(phi2-phi1);
	  
	  
	} 
      }
      
      
      
      if (theside == 2) {
	
	m_dX[thelayer- 5][thering - 1] -> Fill(dX);
	m_dY[thelayer -5][thering - 1] -> Fill(dY);
	m_dR[thelayer - 5][thering - 1] -> Fill(dR);
	m_dr[thelayer - 5][thering - 1] -> Fill(dr);
	m_deltaphi[thelayer - 5][thering - 1] -> Fill(phi2-phi1);
	
	
	if(areSame) {
	  
	  
	  m_dX_sametrack[thelayer - 5][thering - 1] -> Fill(dX);
	  m_dY_sametrack[thelayer - 5][thering - 1] -> Fill(dY);
	  m_dR_sametrack[thelayer - 5][thering - 1] -> Fill(dR);
	  m_dr_sametrack[thelayer - 5][thering - 1] -> Fill(dr);
	  m_deltaphi_sametrack[thelayer - 5][thering - 1] -> Fill(phi2-phi1);
	  
	  
	}
	
	else if(!areSame) {
	  
	  
	  m_dX_notsametrack[thelayer - 5][thering - 1] -> Fill(dX);
	  m_dY_notsametrack[thelayer - 5][thering - 1] -> Fill(dY);
	  m_dR_notsametrack[thelayer - 5][thering - 1] -> Fill(dR);
	  m_dr_notsametrack[thelayer - 5][thering - 1] -> Fill(dr);
	  m_deltaphi_notsametrack[thelayer - 5][thering - 1] -> Fill(phi2-phi1);
	  
	  
	}
      }
    }    
  }
  //  }  
  
  return found2xcoincidencecluster;
  
}


edm::DetSetVector<PixelDigiSimLink>::const_iterator Ashish2xCoincidence5::findSimLinkDetSet(unsigned int thedetid) {
  ////basic template
  edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDS = simlinks->find(thedetid);
  return simLinkDS;
}

std::set<unsigned int> Ashish2xCoincidence5::getSimTrackId(edm::DetSetVector<PixelDigiSimLink>::const_iterator simLinkDSViter, edmNew::DetSet<SiPixelCluster>::const_iterator cluster, bool print) {
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



bool Ashish2xCoincidence5::areSameSimTrackId(std::set<unsigned int> first, std::set<unsigned int> second, std::set<unsigned int>& intersection) {
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


DEFINE_FWK_MODULE(Ashish2xCoincidence5);
