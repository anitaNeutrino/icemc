#include "DiffuseFlux.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanFreqsGenerator.h"


TVector3 icemc::Source::DiffuseFlux::pickNeutrinoDirection(const Geoid::Position& interaction, const TVector3& rfDir) {

  LocalCoordinateSystem lc(interaction, -interaction, rfDir);

  // const double n = AskaryanFreqsGenerator::NICE;
  ///@todo put this somewhere better?

  //@todo return WEIGHT!
  const double theta_cherenkov = AskaryanFreqsGenerator::CHANGLE_ICE;
  ///@todo add a theta integration!!!  
  const double phi = pickUniform(0, TMath::TwoPi());

  TVector3 nuDir;
  nuDir.SetMagThetaPhi(1.0, theta_cherenkov, phi);

  nuDir = lc.localTranslationToGlobal(nuDir);

  return nuDir;
  // fColumnDensity = ((icemc::Crust2*)fW)->integratePath(fInteractionPos, fNeutrinoDir);
}
