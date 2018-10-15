#include "DiffuseFlux.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanRadiationModel.h"


TVector3 icemc::Source::DiffuseFlux::pickNeutrinoDirection(const Geoid::Position& interaction, const TVector3& rfDir) {

  LocalCoordinateSystem lc(interaction, -interaction, rfDir);

  // const double n = AskaryanRadiationModel::NICE;
  ///@todo put this somewhere better?

  //@todo return WEIGHT!
  const double theta_cherenkov = AskaryanRadiationModel::CHANGLE_ICE;
  ///@todo add a theta integration!!!  
  const double phi = pickUniform(0, TMath::TwoPi());

  TVector3 nuDir;
  nuDir.SetMagThetaPhi(1.0, theta_cherenkov, phi);

  nuDir = lc.localTranslationToGlobal(nuDir);

  return nuDir;

}
