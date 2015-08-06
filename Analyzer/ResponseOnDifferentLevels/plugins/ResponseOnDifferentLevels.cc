// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"


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
#include "TProfile.h"

class ResponseOnDifferentLevels : public edm::EDAnalyzer {
  // This module is planned to be used for the energy response correction in fastsim.
  // Here, a tree dimensional histogram will be filled, which is used to calculate the
  // scale factors. Full/Fast sim are required to be run before this module.
  public:
    explicit ResponseOnDifferentLevels(const edm::ParameterSet&);
    ~ResponseOnDifferentLevels(){};

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob();

    // "famosSimHits" for fastsim or "g4SimHits" for fullsim
    std::string module_;

    edm::Service<TFileService> file;

    std::vector<std::string> sclStrings;
    std::vector<std::string> recHitStrings;
    std::vector<std::string> pfRecHitStrings;
    std::vector<std::string> caloHitStrings;
};

ResponseOnDifferentLevels::ResponseOnDifferentLevels(const edm::ParameterSet& iConfig):
  module_( iConfig.getUntrackedParameter<std::string>("module") )
{

  file->make<TH1F>( "gedPhotons:_eb", "", 100, 0.8, 1.1 );
  file->make<TH1F>( "gedPhotons:_ee", "", 100, 0.8, 1.1 );
  file->make<TProfile>( "gedPhotons:_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );

  // Superclusters
  sclStrings.push_back( "hybridSuperClusters:uncleanOnlyHybridSuperClusters" );
  sclStrings.push_back( "multi5x5SuperClusters:uncleanOnlyMulti5x5EndcapSuperClusters" );

  sclStrings.push_back( "hybridSuperClusters:" );
  sclStrings.push_back( "multi5x5SuperClustersWithPreshower:" );

  sclStrings.push_back( "correctedHybridSuperClusters:" );
  sclStrings.push_back( "correctedMulti5x5SuperClustersWithPreshower:" );

  sclStrings.push_back( "particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel" );
  sclStrings.push_back( "particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower" );

  sclStrings.push_back( "multi5x5SuperClusters:multi5x5EndcapSuperClusters" );
  sclStrings.push_back( "particleFlowEGamma:" );
  //sclStrings.push_back( "hfEMClusters:" ); not filled?

  // RecHits
  recHitStrings.push_back( "ecalRecHit:EcalRecHitsEB" );
  recHitStrings.push_back( "ecalRecHit:EcalRecHitsEE" );
  recHitStrings.push_back( "ecalPreshowerRecHit:EcalRecHitsES" );
  recHitStrings.push_back( "reducedEcalRecHitsEB:" );
  recHitStrings.push_back( "reducedEcalRecHitsEE:" );
  recHitStrings.push_back( "reducedEcalRecHitsES:" );

  // PFRecHits
  pfRecHitStrings.push_back( "particleFlowRecHitECAL:" );
  pfRecHitStrings.push_back( "particleFlowRecHitECAL:Cleaned" );

  // Calo Hits
  caloHitStrings.push_back( "EcalHitsEB" );
  caloHitStrings.push_back( "EcalHitsEE" );
  caloHitStrings.push_back( "EcalHitsES" );


  // make all histograms
  std::vector<std::string> allStrings;
  allStrings.reserve( sclStrings.size() + recHitStrings.size() + pfRecHitStrings.size() + caloHitStrings.size() );
  allStrings.insert( allStrings.end(), sclStrings.begin(), sclStrings.end() );
  allStrings.insert( allStrings.end(), recHitStrings.begin(), recHitStrings.end() );
  allStrings.insert( allStrings.end(), pfRecHitStrings.begin(), pfRecHitStrings.end() );
  allStrings.insert( allStrings.end(), caloHitStrings.begin(), caloHitStrings.end() );

  for( auto s : allStrings ) {
      file->make<TH1F>( (s+"_eb").c_str(), "", 100, 0.8, 1.1 );
      file->make<TH1F>( (s+"_ee").c_str(), "", 100, 0.8, 1.1 );
      file->make<TProfile>( (s+"_VsEta").c_str(), "", 100, 0., 3.2, 0.8, 1.1 );
      file->make<TProfile>( (s+"_VsEta_cone01").c_str(), "", 100, 0., 3.2, 0.8, 1.1 );
  }

  // extra for raw energy of superclusters
  for( auto s : sclStrings ) {
      file->make<TH1F>( (s+"_eb_raw").c_str(), "", 100, 0.8, 1.1 );
      file->make<TH1F>( (s+"_ee_raw").c_str(), "", 100, 0.8, 1.1 );
      file->make<TProfile>( (s+"_VsEta_raw").c_str(), "", 100, 0., 3.2, 0.8, 1.1 );
  }

}

template< typename T >
float energyInCone( const T& hits, const CaloGeometry& geo, float eta_, float phi_, float dr=0.3 ) {
    float sumE = 0;
    for( auto const& hit : hits ) {
        auto pos = geo.getPosition( hit.id() );
        if( reco::deltaR( pos.eta(), pos.phi(), eta_, phi_ ) < dr ) {
            sumE += hit.energy();
        }
    }
    return sumE;
}

template< typename T >
float energyInConeDet( const T& hits, const CaloGeometry& geo, float eta_, float phi_, float dr=0.3 ) {
    float sumE = 0;
    for( auto const& hit : hits ) {
        auto pos = geo.getPosition( hit.detId() );
        if( reco::deltaR( pos.eta(), pos.phi(), eta_, phi_ ) < dr ) {
            sumE += hit.energy();
        }
    }
    return sumE;
}

void
ResponseOnDifferentLevels::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByLabel("genParticles", "", GenParticles);

  reco::GenParticle g1, g2;
  for( auto const& gen : *GenParticles ) {
    if( gen.pdgId() != 22 || gen.status() != 3 || !gen.numberOfMothers() || gen.mother()->pdgId() != 25) continue;
      if( !g1.pdgId() ){
        g1 = gen;
      } else{
        g2 = gen;
      }
  }
  if( reco::deltaR( g1.eta(), g1.phi(), g2.eta(), g2.phi() ) < 0.6 ) {
      return;
      // photons are too close to each other, so no matching on sim-level can be done
  }

  bool isG1eb = std::abs( g1.eta() ) < 1.4442;
  bool isG2eb = std::abs( g2.eta() ) < 1.4442;

  bool isG1ee = std::abs( g1.eta() ) > 1.566 && std::abs( g1.eta() ) < 2.5;
  bool isG2ee = std::abs( g2.eta() ) > 1.566 && std::abs( g2.eta() ) < 2.5;


  // reco level
  edm::Handle<std::vector<reco::Photon> > gedPhotons;
  if( iEvent.getByLabel("gedPhotons", "", gedPhotons) ) {

      for( auto const& g : *gedPhotons ) {
          if( reco::deltaR( g.eta(), g.phi(), g1.eta(), g1.phi() ) < 0.5 ) {
              if( isG1eb ) file->getObject<TH1F>("gedPhotons:_eb")->Fill( g.energy()/g1.energy() );
              if( isG1ee ) file->getObject<TH1F>("gedPhotons:_ee")->Fill( g.energy()/g1.energy() );
              file->getObject<TProfile>("gedPhotons:_VsEta")->Fill( std::abs(g1.eta()), g.energy()/g1.energy() );
          }
          if( reco::deltaR( g.eta(), g.phi(), g2.eta(), g2.phi() ) < 0.5 ) {
              if( isG2eb ) file->getObject<TH1F>("gedPhotons:_eb")->Fill( g.energy()/g2.energy() );
              if( isG2ee ) file->getObject<TH1F>("gedPhotons:_ee")->Fill( g.energy()/g2.energy() );
              file->getObject<TProfile>("gedPhotons:_VsEta")->Fill( std::abs(g2.eta()), g.energy()/g2.energy() );
          }

      }
  }

  // sim level

  // set geometry for hits
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo = pG.product();

  // Calo hits
  for( auto s : caloHitStrings ) {
    edm::Handle< edm::PCaloHitContainer > caloHits;
    if( iEvent.getByLabel( "famosSimHits", s, caloHits ) || iEvent.getByLabel( "g4SimHits", s, caloHits ) ) {

        float r1 = energyInCone( *caloHits, *geo, g1.eta(), g1.phi() );
        float r2 = energyInCone( *caloHits, *geo, g2.eta(), g2.phi() );

        if( isG1eb ) file->getObject<TH1F>(s+"_eb")->Fill( r1/g1.energy() );
        if( isG1ee ) file->getObject<TH1F>(s+"_ee")->Fill( r1/g1.energy() );
        if( isG2eb ) file->getObject<TH1F>(s+"_eb")->Fill( r2/g2.energy() );
        if( isG2ee ) file->getObject<TH1F>(s+"_ee")->Fill( r2/g2.energy() );
        file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g1.eta()), r1/g1.energy() );
        file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g2.eta()), r2/g2.energy() );
        file->getObject<TProfile>(s+"_VsEta_cone01")->Fill( std::abs(g1.eta()), energyInCone( *caloHits, *geo, g1.eta(), g1.phi(), 0.1 )/g1.energy() );
        file->getObject<TProfile>(s+"_VsEta_cone01")->Fill( std::abs(g2.eta()), energyInCone( *caloHits, *geo, g2.eta(), g2.phi(), 0.1 )/g2.energy() );
    }
  }

  // pfRecHits
  for( auto s : pfRecHitStrings ) {
    edm::Handle< std::vector<reco::PFRecHit> > pfRecHits;
    if( iEvent.getByLabel( module_, s, pfRecHits ) ) {

        float r1 = energyInConeDet( *pfRecHits, *geo, g1.eta(), g1.phi() );
        float r2 = energyInConeDet( *pfRecHits, *geo, g2.eta(), g2.phi() );

        if( isG1eb ) file->getObject<TH1F>(s+"_eb")->Fill( r1/g1.energy() );
        if( isG1ee ) file->getObject<TH1F>(s+"_ee")->Fill( r1/g1.energy() );
        if( isG2eb ) file->getObject<TH1F>(s+"_eb")->Fill( r2/g2.energy() );
        if( isG2ee ) file->getObject<TH1F>(s+"_ee")->Fill( r2/g2.energy() );
        file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g1.eta()), r1/g1.energy() );
        file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g2.eta()), r2/g2.energy() );
        file->getObject<TProfile>(s+"_VsEta_cone01")->Fill( std::abs(g1.eta()), energyInConeDet( *pfRecHits, *geo, g1.eta(), g1.phi(), .1 )/g1.energy() );
        file->getObject<TProfile>(s+"_VsEta_cone01")->Fill( std::abs(g2.eta()), energyInConeDet( *pfRecHits, *geo, g2.eta(), g2.phi(), .1 )/g2.energy() );
    }
  }

  ////////////////////////////////////////////////////////////////
  // ECALRecHits
  ////////////////////////////////////////////////////////////////

  for( auto s : recHitStrings ) {
    edm::Handle< edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHits;
    if( iEvent.getByLabel( s.substr(0, s.find(":")), s.substr( s.find(":")+1 ), recHits ) ) {

        float r1 = energyInCone( *recHits, *geo, g1.eta(), g1.phi() );
        float r2 = energyInCone( *recHits, *geo, g2.eta(), g2.phi() );

        if( isG1eb ) file->getObject<TH1F>(s+"_eb")->Fill( r1/g1.energy() );
        if( isG1ee ) file->getObject<TH1F>(s+"_ee")->Fill( r1/g1.energy() );
        if( isG2eb ) file->getObject<TH1F>(s+"_eb")->Fill( r2/g2.energy() );
        if( isG2ee ) file->getObject<TH1F>(s+"_ee")->Fill( r2/g2.energy() );
        file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g1.eta()), r1/g1.energy() );
        file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g2.eta()), r2/g2.energy() );
        file->getObject<TProfile>(s+"_VsEta_cone01")->Fill( std::abs(g1.eta()), energyInCone( *recHits, *geo, g1.eta(), g1.phi(), .1 )/g1.energy() );
        file->getObject<TProfile>(s+"_VsEta_cone01")->Fill( std::abs(g2.eta()), energyInCone( *recHits, *geo, g2.eta(), g2.phi(), .1 )/g2.energy() );
    }
  }


  ////////////////////////////////////////////////////////////////
  // SuperCluster
  ////////////////////////////////////////////////////////////////
  for( auto s : sclStrings ) {
    edm::Handle< std::vector<reco::SuperCluster> > cls;
    if( iEvent.getByLabel( s.substr(0, s.find(":")), s.substr( s.find(":")+1 ), cls ) ) {
        for( auto const& cl : *cls ) {

          if( reco::deltaR( cl.eta(), cl.phi(), g1.eta(), g1.phi() ) < 0.5 ) {
            if( isG1eb ) file->getObject<TH1F>(s+"_eb")->Fill( cl.energy()/g1.energy() );
            if( isG1ee ) file->getObject<TH1F>(s+"_ee")->Fill( cl.energy()/g1.energy() );
            file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g1.eta()), cl.energy()/g1.energy() );

            if( isG1eb ) file->getObject<TH1F>(s+"_eb_raw")->Fill( cl.rawEnergy()/g1.energy() );
            if( isG1ee ) file->getObject<TH1F>(s+"_ee_raw")->Fill( cl.rawEnergy()/g1.energy() );
            file->getObject<TProfile>(s+"_VsEta_raw")->Fill( std::abs(g1.eta()), cl.rawEnergy()/g1.energy() );
          }

          if( reco::deltaR( cl.eta(), cl.phi(), g2.eta(), g2.phi() ) < 0.5 ) {
            if( isG2eb ) file->getObject<TH1F>(s+"_eb")->Fill( cl.energy()/g2.energy() );
            if( isG2ee ) file->getObject<TH1F>(s+"_ee")->Fill( cl.energy()/g2.energy() );
            file->getObject<TProfile>(s+"_VsEta")->Fill( std::abs(g2.eta()), cl.energy()/g2.energy() );

            if( isG2eb ) file->getObject<TH1F>(s+"_eb_raw")->Fill( cl.rawEnergy()/g2.energy() );
            if( isG2ee ) file->getObject<TH1F>(s+"_ee_raw")->Fill( cl.rawEnergy()/g2.energy() );
            file->getObject<TProfile>(s+"_VsEta_raw")->Fill( std::abs(g2.eta()), cl.rawEnergy()/g2.energy() );
          }

        }
    }
  }



}

void
ResponseOnDifferentLevels::endJob() {

  TProfile* tmp;

  // merge TPRofilesVsEta
  tmp = file->make<TProfile>( "EcalHits_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "EcalHitsEB_VsEta") );
  tmp->Add( file->getObject<TProfile>( "EcalHitsEE_VsEta") );

  tmp = file->make<TProfile>( "ecalRecHit:EcalRecHits_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "ecalRecHit:EcalRecHitsEB_VsEta") );
  tmp->Add( file->getObject<TProfile>( "ecalRecHit:EcalRecHitsEE_VsEta") );

  tmp = file->make<TProfile>( "reducedEcalRecHits:_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "reducedEcalRecHitsEB:_VsEta") );
  tmp->Add( file->getObject<TProfile>( "reducedEcalRecHitsEE:_VsEta") );

  tmp = file->make<TProfile>( "uncleanSuperClusters_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "hybridSuperClusters:uncleanOnlyHybridSuperClusters_VsEta") );
  tmp->Add( file->getObject<TProfile>( "multi5x5SuperClusters:uncleanOnlyMulti5x5EndcapSuperClusters_VsEta") );

  tmp = file->make<TProfile>( "SuperClusters_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "hybridSuperClusters:_VsEta") );
  tmp->Add( file->getObject<TProfile>( "multi5x5SuperClustersWithPreshower:_VsEta") );

  tmp = file->make<TProfile>( "correctedSuperClusters_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "correctedHybridSuperClusters:_VsEta") );
  tmp->Add( file->getObject<TProfile>( "correctedMulti5x5SuperClustersWithPreshower:_VsEta") );

  tmp = file->make<TProfile>( "particleFlowSuperClusterECAL_VsEta", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel_VsEta") );
  tmp->Add( file->getObject<TProfile>( "particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower_VsEta") );

  tmp = file->make<TProfile>( "uncleanSuperClusters_VsEta_raw", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "hybridSuperClusters:uncleanOnlyHybridSuperClusters_VsEta_raw") );
  tmp->Add( file->getObject<TProfile>( "multi5x5SuperClusters:uncleanOnlyMulti5x5EndcapSuperClusters_VsEta_raw") );

  tmp = file->make<TProfile>( "SuperClusters_VsEta_raw", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "hybridSuperClusters:_VsEta_raw") );
  tmp->Add( file->getObject<TProfile>( "multi5x5SuperClustersWithPreshower:_VsEta_raw") );

  tmp = file->make<TProfile>( "correctedSuperClusters_VsEta_raw", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "correctedHybridSuperClusters:_VsEta_raw") );
  tmp->Add( file->getObject<TProfile>( "correctedMulti5x5SuperClustersWithPreshower:_VsEta_raw") );

  tmp = file->make<TProfile>( "particleFlowSuperClusterECAL_VsEta_raw", "", 100, 0., 3.2, 0.8, 1.1 );
  tmp->Add( file->getObject<TProfile>( "particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel_VsEta_raw") );
  tmp->Add( file->getObject<TProfile>( "particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower_VsEta_raw") );

}

void
ResponseOnDifferentLevels::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>( "module", "g4SimHits" );
  descriptions.add( "ecalScaleFactorCalculator", desc );
}

DEFINE_FWK_MODULE(ResponseOnDifferentLevels);
/
