
void getCovariances( const edm::PCaloHitContainer& SimHits, const CaloGeometry& geometry, const CaloTopology& topo, double& sigmaEtaEta, double& sigmaPhiPhi, double& sigmaEtaPhi )
{

  // Find crystal with most energy
  float energyMaxCrystal = 0;
  PCaloHit maxHit;

  for( auto const& Hit : SimHits ) {
    if( Hit.energy() > energyMaxCrystal ) {
        energyMaxCrystal = Hit.energy();
        maxHit = Hit;
    }
  }

  // Copy all interesting crrystal enregyies to temporary map
  auto field5x5 = topo.getWindow( maxHit.id(), 5, 5 );
  std::vector< std::pair<DetId, double> > energyMap;
  for( auto const& Hit : SimHits ) {
    if( std::find( field5x5.begin(), field5x5.end(), Hit.id() ) != field5x5.end() ) {
      energyMap.push_back( std::pair<DetId, double>( Hit.id(), Hit.energy() ) );
    }
  }


  // first find energy-weighted mean position - doing it when filling the energy map might save time
  double e5x5_ = 0;
  math::XYZVector meanPosition(0.0, 0.0, 0.0);
  for ( auto energyMapEntry : energyMap)
  {
      DetId id = energyMapEntry.first;
      if (id != DetId(0))
      {
        GlobalPoint positionGP = geometry.getPosition(id);
        math::XYZVector position(positionGP.x(),positionGP.y(),positionGP.z());
        meanPosition = meanPosition + energyMapEntry.second * position;
      }
      e5x5_ += energyMapEntry.second;
  }

  meanPosition /= e5x5_;

  // now we can calculate the covariances
  double numeratorEtaEta = 0;
  double numeratorEtaPhi = 0;
  double numeratorPhiPhi = 0;
  double denominator     = 0;

  for ( auto energyMapEntry : energyMap )
  {
      DetId id = energyMapEntry.first;
      if (id != DetId(0))
      {
        GlobalPoint position = geometry.getPosition(id);

        double dPhi = position.phi() - meanPosition.phi();
        if (dPhi > + Geom::pi()) { dPhi = Geom::twoPi() - dPhi; }
        if (dPhi < - Geom::pi()) { dPhi = Geom::twoPi() + dPhi; }

        double dEta = position.eta() - meanPosition.eta();
        double w = 0.;
        if ( energyMapEntry.second > 0.)
          w = std::max(0.0, 4.7 + log( energyMapEntry.second / e5x5_));

        denominator += w;
        numeratorEtaEta += w * dEta * dEta;
        numeratorEtaPhi += w * dEta * dPhi;
        numeratorPhiPhi += w * dPhi * dPhi;
      }
  }

  sigmaEtaEta = sqrt(numeratorEtaEta / denominator);
  sigmaEtaPhi = sqrt(numeratorEtaPhi / denominator);
  sigmaPhiPhi = sqrt(numeratorPhiPhi / denominator);
}

