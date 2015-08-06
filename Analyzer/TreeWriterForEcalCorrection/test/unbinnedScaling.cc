#include<iomanip> // provides setprecision
#include<iostream>
#include<sstream>
#include<string>

// ROOT
#include<TChain.h>
#include<TFile.h>
#include<TH3F.h>
#include<TMath.h>
#include<TProfile.h>
#include<TROOT.h>

using namespace std;

std::map<int,std::vector<float>> getVectorsFromTree( TChain& tree, const TH3F& h, float minR ) {
  /* map.first:  global bin in TH3F h, belonging to e,eta,-1
   * map.second: vector of r in this bin (with this e,eta)
   */

  std::map<int,std::vector<float>> vectors;

  float r,e,eta;
  tree.SetBranchAddress("r",&r);
  tree.SetBranchAddress("e",&e);
  tree.SetBranchAddress("eta",&eta);
  for( int i=0; i<tree.GetEntries(); i++ ) {
    tree.GetEntry(i);
    if( r < minR ) continue;
    int bin = h.FindFixBin( e, eta, -1 );

    // check if this bin already exists in map
    if( !vectors.count( bin ) ) {
      vectors[bin] = std::vector<float>();
    }

    vectors.at(bin).push_back( r );
  }

  for( auto& m : vectors ) {
    std::sort( m.second.begin(), m.second.end() );
  }

  return vectors;
}



TH3F calculateResponseUnbinned( TChain& fasttree, TChain& fulltree, TH3F h, float minR ) {

  auto fastvectors = getVectorsFromTree( fasttree, h, minR );
  auto fullvectors = getVectorsFromTree( fulltree, h, minR );

  for( auto m : fastvectors ) {

    if( ! fullvectors.count( m.first ) ) {
      std::cout << "no compatible vector found for fullsim" << std::endl;
      continue;
    }

    auto fastvector = m.second;
    auto fullvector = fullvectors.at( m.first );

    TProfile profile("profile", "title", h.GetZaxis()->GetNbins(), h.GetZaxis()->GetXmin(), h.GetZaxis()->GetXmax(), "s" );
    for( unsigned i=0; i<fastvector.size(); i++ ) {
      auto fa = fastvector[i];
      auto jRel = 1.*i/fastvector.size()*fullvector.size();
      int j = (int) jRel;
      auto fu = fullvector[j]*(jRel-j)+fullvector[j+1]*(j-jRel+1);
      if( fa ) profile.Fill( fa, fu/fa );
    }

    int xbin, ybin, zbin;
    h.GetBinXYZ( m.first, xbin, ybin, zbin );

    for( int i=0; i<profile.GetNbinsX()+2; i++ ) {
      h.SetBinContent( xbin, ybin, i, profile.GetBinContent(i) );
    }

  }

  return h;

}

void write( const TObject& ob, const std::string& filename ) {
  TFile f( filename.c_str(), "recreate" );
  ob.Write();
  f.Close();
  cout << "Info: Wrote " << ob.GetName() << " to file " << filename << "." << std::endl;
}

int main( int argc, char** argv ) {
  string fastname = "../fast_tree.root";
  string fullname = "../full_tree.root";
  string treeName = "treeWriterForEcalCorrection/responseTree";

  TChain fasttree( treeName.c_str() );
  fasttree.AddFile( fastname.c_str() );
  TChain fulltree( treeName.c_str() );
  fulltree.AddFile( fullname.c_str() );

  TH3F h3default("scaleVsEVsEta", ";E_{gen};#eta_{gen};E/E_{gen}", 100, 5, 1005, 40, 0, 3.2, 1000, 0.3, 1.05 );

  auto scales3d = calculateResponseUnbinned( fasttree, fulltree, h3default, /*minR*/0.3 );
  write( scales3d, "scaleECALFastsim.root" );
}




