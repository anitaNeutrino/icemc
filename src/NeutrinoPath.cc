#include "NeutrinoPath.h"
#include "WorldModel.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanFreqsGenerator.h"
#include "Crust2.h"

// ClassImp(icemc::NeutrinoPath);

icemc::NeutrinoPath::NeutrinoPath(){

  reset();
}


icemc::NeutrinoPath::NeutrinoPath(const Geoid::Position& interaction, const TVector3& rfDir, const WorldModel* m)
  : fW(m), fInteractionPos(interaction)

{
  
  reset();

  LocalCoordinateSystem lc(interaction, -interaction, rfDir);

  // const double n = AskaryanFreqsGenerator::NICE;
  ///@todo put this somewhere better?
  
  const double theta_cherenkov = AskaryanFreqsGenerator::CHANGLE_ICE;
  ///@todo add a theta integration!!!
  ///@todo make integrator / random number generator
  const double phi = gRandom->Uniform(0, TMath::TwoPi());
  fNeutrinoDir.SetMagThetaPhi(1.0, theta_cherenkov, phi);

  fNeutrinoDir = lc.localTranslationToGlobal(fNeutrinoDir);
  
  ((icemc::Crust2*)fW)->integratePath(fInteractionPos, fNeutrinoDir);
}









void icemc::NeutrinoPath::reset(){
  // zero all members
  theta_in = 0;
  lat_in = 0;
  nearthlayers = 0;
  weight_prob = 0.;
  weight1 = 0;
  weight = 0.;
  logweight = 0.;
  len_int = 0;
  pieceofkm2sr = 0;
}




