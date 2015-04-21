// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Geometry & Topology
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"


// user includes
#include "HelperFunctions.h"

class ECALScaleFactorCalculator : public edm::EDAnalyzer {
  // This module is planned to be used for the energy response correction in fastsim.
  // Here, a tree dimensional histogram will be filled, which is used to calculate the
  // scale factors. Full/Fast sim are required to be run before this module.
  public:
    explicit ECALScaleFactorCalculator(const edm::ParameterSet&);
    ~ECALScaleFactorCalculator(){};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

    // "famosSimHits" for fastsim or "g4SimHits" for fullsim
    std::string module_;

    edm::Service<TFileService> file;
    TH2F* etaGaps;
    TH2F* simHitMap;
    TH3F* responseVsEVsEta;
    TH3F* responseVsEVsEta_woGaps;
    TH3F* responseVsEVsEta_woGaps3;
    TH3F* responseVsEVsEta_woGaps5;

    TH3F* sEtaEtaVsEVsEta_woGaps;
    TH3F* sPhiPhiVsEVsEta_woGaps;
    TH3F* sEtaPhiVsEVsEta_woGaps;
    TTree* tree;
    float tree_e, tree_eta, tree_response;
};

ECALScaleFactorCalculator::ECALScaleFactorCalculator(const edm::ParameterSet& iConfig):
  module_( iConfig.getUntrackedParameter<std::string>("module") )
{

  // Check for gaps in eta
  etaGaps = file->make<TH2F>("etaGaps", ";#eta_{gen};E/E_{gen}", 3000, 0, 3, 1000, 5, 1005 );

  simHitMap = file->make<TH2F>("hitMap", "hitmap;#eta;#phi", 1000, -3, 3, 1000, -3.5, 3.5 );

  // This 3D histogram can be used for calculating the scale factors
  responseVsEVsEta = file->make<TH3F>("responseVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 1.05 );

  // This 3D histogram is used for calculating the scale factors, the gaps are obmitted
  responseVsEVsEta_woGaps = file->make<TH3F>("responseVsEVsEta_woGaps", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 1.05 );
  responseVsEVsEta_woGaps3 = file->make<TH3F>("responseVsEVsEta_woGaps3", ";E_{gen};#eta_{gen};E^{3x3}/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 1.05 );
  responseVsEVsEta_woGaps5 = file->make<TH3F>("responseVsEVsEta_woGaps5", ";E_{gen};#eta_{gen};E^{5x5}/E_{gen}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 1.05 );


  sEtaEtaVsEVsEta_woGaps = file->make<TH3F>("sEtaEtaVsEVsEta_woGaps", ";E_{gen};#eta_{gen};#sigma_{#eta#eta}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 0.1 );
  sPhiPhiVsEVsEta_woGaps = file->make<TH3F>("sPhiPhiVsEVsEta_woGaps", ";E_{gen};#eta_{gen};#sigma_{#phi#phi}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 0.1 );
  sEtaPhiVsEVsEta_woGaps = file->make<TH3F>("sEtaPhiVsEVsEta_woGaps", ";E_{gen};#eta_{gen};#sigma_{#eta#phi}", 100, 5, 1005, 400, 0, 3.2, 2000, 0, 0.1 );

  tree = file->make<TTree>("responseTree", "same info as 3dhisto");
  tree->Branch( "e", &tree_e, "e/F");
  tree->Branch( "eta", &tree_eta, "eta/F");
  tree->Branch( "r", &tree_response, "r/F");
}

void
ECALScaleFactorCalculator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByLabel("genParticles", "", GenParticles);

  // As this module is intended for single particle guns, there should be exactly one
  // generated particle.
  float genE = -1, genEta = 100; // default values, which make no physical sensc
  if( GenParticles->size() ) {
    reco::GenParticle gen = GenParticles->at(0);
    genE = gen.energy();
    genEta = gen.eta(); // gen eta is positive for my photongun sample per definition
  }

  // set geometry for hits
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo = pG.product();

  // set topology for the calorimeter
  edm::ESHandle<CaloTopology> topo;
  iSetup.get<CaloTopologyRecord>().get(topo);

  edm::Handle<edm::PCaloHitContainer> SimHitsEB;
  iEvent.getByLabel( module_, "EcalHitsEB", SimHitsEB );
  edm::Handle<edm::PCaloHitContainer> SimHitsEE;
  iEvent.getByLabel( module_, "EcalHitsEE", SimHitsEE );
  edm::Handle<edm::PCaloHitContainer> SimHitsES;
  iEvent.getByLabel( module_, "EcalHitsES", SimHitsES );

  // merge them into one single vector
  auto SimHits = *SimHitsEB;
  SimHits.insert( SimHits.end(), SimHitsEE->begin(), SimHitsEE->end() );
  SimHits.insert( SimHits.end(), SimHitsES->begin(), SimHitsES->end() );

  float energyMaxCrystal = 0;
  float energy3x3 = 0;
  float energy5x5 = 0;
  float energyTotal = 0;

  PCaloHit maxHit;

  for( auto const& Hit : SimHits ) {
    if( Hit.energy() > energyMaxCrystal ) {
        energyMaxCrystal = Hit.energy();
        maxHit = Hit;
    }
    energyTotal += Hit.energy();
    auto pos = geo->getPosition( Hit.id() );
    simHitMap->Fill( pos.eta(), pos.phi(), Hit.energy() );
  }

  auto field3x3 = topo->getAllNeighbours( maxHit.id() );
  auto field5x5 = topo->getWindow( maxHit.id(), 5, 5 );
  for( auto const& Hit : SimHits ) {
    if( std::find( field3x3.begin(), field3x3.end(), Hit.id() ) != field3x3.end() ) {
      energy3x3 += Hit.energy();
    }
    if( std::find( field5x5.begin(), field5x5.end(), Hit.id() ) != field5x5.end() ) {
      energy5x5 += Hit.energy();
    }

  }

  double sigmaEtaEta, sigmaPhiPhi, sigmaEtaPhi;
  getCovariances( SimHits, *geo, *topo, sigmaEtaEta, sigmaPhiPhi, sigmaEtaPhi );

  responseVsEVsEta->Fill( genE, genEta, energyTotal/genE );
  etaGaps->Fill( genEta, genE, energyTotal/genE );

  // Gaps are obmitted. The values are estimated by looking at the 'etaGaps' histogram.
  // The scaling in the gaps could be studied in more detail, but as the fraction of gaps is low,
  // the influence should be small.
  if( (genEta>0.0035 && genEta<0.436)
      || (genEta>0.444 && genEta<0.783)
      || (genEta>0.794 && genEta<1.132)
      || (genEta>1.141 && genEta<1.475)
      || genEta>1.51 ) {
    responseVsEVsEta_woGaps->Fill( genE, genEta, energyTotal/genE, 1 );
    responseVsEVsEta_woGaps3->Fill( genE, genEta, energy3x3/genE, 1 );
    responseVsEVsEta_woGaps5->Fill( genE, genEta, energy5x5/genE, 1 );
    sEtaEtaVsEVsEta_woGaps->Fill( genE, genEta, sigmaEtaEta, 1 );
    sPhiPhiVsEVsEta_woGaps->Fill( genE, genEta, sigmaPhiPhi, 1 );
    sEtaPhiVsEVsEta_woGaps->Fill( genE, genEta, sigmaEtaPhi, 1 );
  }

  tree_e = genE;
  tree_eta = genEta;
  tree_response = energyTotal/genE;
  tree->Fill();

}

void
ECALScaleFactorCalculator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>( "module", "g4SimHits" );
  descriptions.add( "ecalScaleFactorCalculator", desc );
}

DEFINE_FWK_MODULE(ECALScaleFactorCalculator);
