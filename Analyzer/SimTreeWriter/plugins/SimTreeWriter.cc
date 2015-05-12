// -*- C++ -*-
//
// Package:    FastSim/SimTreeWriter
// Class:      SimTreeWriter
//
/**\class SimTreeWriter SimTreeWriter.cc FastSim/SimTreeWriter/plugins/SimTreeWriter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Maximilian Knut Kiesel
//         Created:  Fri, 17 Oct 2014 08:04:04 GMT
//
//

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

//Geometry
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "TProfile.h"
#include "TTree.h"
#include "TVector3.h"

//
// class declaration
//


// helper wrapper
std::ostream& operator << ( std::ostream& os, const TVector3& vec ) {
    os << "pt, eta, phi = " << std::setprecision(3) << vec.Pt() << "\t" << vec.Eta() << "\t" << vec.Phi();
    return os;
}

void computeShowerWidth( const std::vector<TVector3>& hits, float& sigmaEtaEta, float& sigmaPhiPhi );
void energyInConeAroundSeed( const std::vector<TVector3>& hits, TVector3& coneVector, float deltaR=0.1 );



class SimTreeWriter : public edm::EDAnalyzer {
   public:
      explicit SimTreeWriter(const edm::ParameterSet&);
      ~SimTreeWriter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;




      // ----------member data ---------------------------
      TTree* outTree_;
      std::string Module_;


      // variables added to tree
      TVector3 _genVec;
      TVector3 _hitVec;
      TVector3 _coneVector;
      std::vector< TVector3 > _hitVecs;
      float _sigmaEtaEta;
      float _sigmaPhiPhi;
      float _showerStdPhi;
      float _showerStdEta;
      float _time;
      float _timeStd;
      TH2F* simHitMap;
      edm::Service<TFileService> file;
};

SimTreeWriter::SimTreeWriter(const edm::ParameterSet& iConfig) :
    Module_( iConfig.getParameter<std::string>("ProducerModule") )
{
}


SimTreeWriter::~SimTreeWriter()
{
}

void
SimTreeWriter::beginJob()
{
    // create dictionary
    gInterpreter->GenerateDictionary("std::vector<TVector3>","TVector3.h;vector");

    outTree_ = file->make<TTree>("SimTree", "Used for FastSim-FullSim comparisons");
    outTree_->Branch( "genVec", &_genVec );
    outTree_->Branch( "hitVec", &_hitVec );
    outTree_->Branch( "coneVector", &_coneVector );
    outTree_->Branch( "sigmaEtaEta", &_sigmaEtaEta, "sigmaEtaEta/F" );
    outTree_->Branch( "sigmaPhiPhi", &_sigmaPhiPhi, "sigmaPhiPhi/F" );
    outTree_->Branch( "showerStdEta", &_showerStdEta, "showerStdEta/F" );
    outTree_->Branch( "showerStdPhi", &_showerStdPhi, "showerStdPhi/F" );
    outTree_->Branch( "time", &_time, "time/F" );
    outTree_->Branch( "timeStd", &_timeStd, "timeStd/F" );

    simHitMap = file->make<TH2F>("simHitMap", "title", 1000, -3, 3, 1000, -3.5, 3.5 );
    // set split level to 0 to spress warning
    //outTree_->Branch( "hitVecs", &_hitVecs, 32000, 0 );

    file->make<TH1F>( "e_eb", ";E (GeV);", 50, 0, 130 );
    file->make<TH1F>( "e_ee", ";E (GeV);", 50, 0, 130 );
    file->make<TH1F>( "eta", ";|#eta|;", 50, 0, 3 );
    file->make<TH1F>( "response", ";E_{sim}/E_{gen}; Normalized Entries", 100, 0.8, 1.1 );
    file->make<TH1F>( "response_eb", ";Barrel E_{sim}/E_{gen}; Normalized Entries", 100, 0.8, 1.1 );
    file->make<TH1F>( "response_ee", ";Endcap E_{sim}/E_{gen}; Normalized Entries", 100, 0.8, 1.1 );
    file->make<TH1F>( "time", "t", 50, 4, 14 );
    file->make<TProfile>( "responseVsEta", ";|#eta|;E_{sim}/E_{gen}", 50, 0, 3, 0.8, 1.1 );
    file->make<TProfile>( "responseVsE", ";E_{sim};E_{sim}/E_{gen}", 50, 0, 130, 0.8, 1.1 );
    file->make<TProfile>( "responseVsEgen", ";E_{gen};E_{sim}/E_{gen}", 10, 5, 105, 0.8, 1.1 );
}



// ------------ method called for each event  ------------
void
SimTreeWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::PCaloHitContainer> ECALEBSimHits;
    edm::Handle<edm::PCaloHitContainer> ECALEESimHits;
    edm::Handle<edm::PCaloHitContainer> ECALESSimHits;
    edm::Handle<reco::GenParticleCollection> GenParticles;

    iEvent.getByLabel(Module_,"EcalHitsEB",ECALEBSimHits);
    iEvent.getByLabel(Module_,"EcalHitsEE",ECALEESimHits);
    iEvent.getByLabel(Module_,"EcalHitsES",ECALESSimHits);
    iEvent.getByLabel("genParticles","",GenParticles);

    _genVec.SetXYZ(0,0,0);
    _sigmaEtaEta = 0;
    _sigmaPhiPhi = 0;
    _showerStdPhi = 0;
    _showerStdEta = 0;
    _time = 0;
    _timeStd = 0;

    if( GenParticles->size() ) {
        reco::GenParticle gen = GenParticles->at(0);
        _genVec.SetPtEtaPhi( gen.pt(), gen.eta(), gen.phi() );
    }

    // set geometry for hits
    edm::ESHandle<CaloGeometry> pG;
    iSetup.get<CaloGeometryRecord>().get(pG);
    const CaloGeometry* geo = pG.product();

    PCaloHit Hit;
    GlobalPoint pos;

    std::map< unsigned int, TVector3 > hitMap;
    TVector3 tmpVector;

    // Barrel
    for(unsigned int i = 0; i<ECALEBSimHits->size(); ++i){
        Hit = ECALEBSimHits->at(i);

        int hitId = Hit.id();
        pos = geo->getPosition( hitId );
        simHitMap->Fill( pos.eta(), pos.phi(), Hit.energy() );

        tmpVector.SetPtEtaPhi( 1, pos.eta(), pos.phi() );
        tmpVector.SetMag( Hit.energy() );

        // hits with the same id have same eta and phi, and their energy/pt is added up
        if( hitMap.find( hitId ) == hitMap.end() ) {
            hitMap[ hitId ] = tmpVector;
        } else {
            hitMap[ hitId ] += tmpVector;
        }
    }

    // Endcap
    for(unsigned int i = 0; i<ECALEESimHits->size(); ++i){
        Hit = ECALEESimHits->at(i);

        int hitId = Hit.id();
        pos = geo->getPosition( hitId );
        simHitMap->Fill( pos.eta(), pos.phi(), Hit.energy() );

        tmpVector.SetPtEtaPhi( 1, pos.eta(), pos.phi() );
        tmpVector.SetMag( Hit.energy() );

        // hits with the same id have same eta and phi, and their energy/pt is added up
        if( hitMap.find( hitId ) == hitMap.end() ) {
            hitMap[ hitId ] = tmpVector;
        } else {
            hitMap[ hitId ] += tmpVector;
        }
    }

    // Preshower detector
    for(unsigned int i = 0; i<ECALESSimHits->size(); ++i){
        Hit = ECALESSimHits->at(i);

        int hitId = Hit.id();
        pos = geo->getPosition( hitId );
        simHitMap->Fill( pos.eta(), pos.phi(), Hit.energy() );

        tmpVector.SetPtEtaPhi( 1, pos.eta(), pos.phi() );
        tmpVector.SetMag( Hit.energy() );

        // hits with the same id have same eta and phi, and their energy/pt is added up
        if( hitMap.find( hitId ) == hitMap.end() ) {
            hitMap[ hitId ] = tmpVector;
        } else {
            hitMap[ hitId ] += tmpVector;
        }
    }

    // now fill the hitVecs vector
    _hitVecs.clear();
    _hitVec.SetPtEtaPhi(0,0,0);
    float sum_E = 0;
    for( std::map< unsigned int, TVector3 >::iterator it = hitMap.begin(); it != hitMap.end(); ++it ) {
        _hitVec += it->second;
        _hitVecs.push_back( it->second );
        float e = _hitVec.Pt() / sin( 2*atan( exp( -_hitVec.Eta() ) ) );
        _showerStdPhi += e* pow( _hitVec.DeltaPhi( _genVec ), 2 );
        _showerStdEta += e* pow( _hitVec.Eta() - _genVec.Eta(), 2 );
        sum_E += e;
    }
    computeShowerWidth( _hitVecs, _sigmaEtaEta, _sigmaPhiPhi );
    _showerStdPhi = sqrt( _showerStdPhi ) / sum_E;
    _showerStdEta = sqrt( _showerStdEta ) / sum_E;

    energyInConeAroundSeed( _hitVecs, _coneVector );

    bool isEB = abs(_genVec.Eta()) < 1.4442;

    isEB ? file->getObject<TH1F>("e_eb")->Fill( _hitVec.Mag() ) :
        file->getObject<TH1F>("e_ee")->Fill( _hitVec.Mag() );

    file->getObject<TH1F>("eta")->Fill( _hitVec.Eta() );

    isEB ? file->getObject<TH1F>("response_eb")->Fill( _hitVec.Mag() / _genVec.Mag() ) :
        file->getObject<TH1F>("response_ee")->Fill( _hitVec.Mag() / _genVec.Mag() );

    file->getObject<TH1F>("response")->Fill( _hitVec.Mag() / _genVec.Mag() );
    file->getObject<TProfile>("responseVsEta")->Fill( std::abs(_genVec.Eta()), _hitVec.Mag() / _genVec.Mag() );
    file->getObject<TProfile>("responseVsE")->Fill( _hitVec.Mag(), _hitVec.Mag() / _genVec.Mag() );
    file->getObject<TProfile>("responseVsEgen")->Fill( _genVec.Mag(), _hitVec.Mag() / _genVec.Mag() );

    file->getObject<TH1F>("time")->Fill( _time );

    outTree_->Fill();
}

void energyInConeAroundSeed( const std::vector<TVector3>& hits, TVector3& coneVector, float deltaR )
{
    TVector3 maxEnergyVec = hits[0];
    for( std::vector<TVector3>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit ) {
        if( hit->Mag() > maxEnergyVec.Mag() ) {
            maxEnergyVec = *hit;
        }
    }
    coneVector.SetXYZ(0,0,0);
    for( std::vector<TVector3>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit ) {
        if( maxEnergyVec.DeltaR( *hit ) < deltaR ) {
            coneVector += *hit;
        }
    }
}

void computeShowerWidth( const std::vector<TVector3>& hits, float& sigmaEtaEta, float& sigmaPhiPhi )
{
    sigmaEtaEta = 0;
    sigmaPhiPhi = 0;
    if( !hits.size() ) return;

    // find the crystal with most energy
    float etaMaximum = 0;
    float phiMaximum = 0;
    float maxEnergy = 0;

    for( std::vector<TVector3>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit ) {
        if( hit->Perp() > maxEnergy ) {
            maxEnergy = hit->Perp();
            etaMaximum = hit->Eta();
            phiMaximum = hit->Phi();
        }
    }

    float ecalCrystalSize = 0;
    if( fabs(etaMaximum) < 1.4442 ) { //barrel
        ecalCrystalSize = 0.0174;
    } else if ( fabs(etaMaximum) > 1.52 ) { // endcap
        ecalCrystalSize = 0.0226358;
    } else {
        return;
    }

    // compute mean
    float meanEta = 0;
    float meanPhi = 0;
    float energy = 0;
    for( std::vector<TVector3>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit ) {
        if( abs( hit->Eta() - etaMaximum ) > 2*ecalCrystalSize || abs( hit->Phi() - phiMaximum ) > 2*ecalCrystalSize ) continue; // use 5x5 matrix
        meanEta += hit->Eta();
        meanPhi += hit->Phi();
        energy += hit->Perp();
    }
    meanEta /= hits.size();
    meanPhi /= hits.size();

    float dEta = 0;
    float dPhi = 0;
    float sum_weight = 0;
    for( std::vector<TVector3>::const_iterator hit = hits.begin(); hit != hits.end(); ++hit ) {
        if( abs( hit->Eta() - etaMaximum ) > 2*ecalCrystalSize || abs( hit->Phi() - phiMaximum ) > 2*ecalCrystalSize ) continue; // use 5x5 matrix

        float w = 4.7+log( hit->Perp() / energy );
        if( w < 0 ) w = 0;
        dEta += w* pow(meanEta - hit->Eta(), 2 );
        dPhi += w* pow(meanPhi - hit->Phi(), 2 );
        sum_weight += w;
    }

    if( sum_weight > 0 ) {
        sigmaEtaEta = sqrt( dEta / sum_weight );
        sigmaPhiPhi = sqrt( dPhi / sum_weight );
    }
}


void
SimTreeWriter::endJob()
{
}

void
SimTreeWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimTreeWriter);
