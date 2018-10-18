#include "DiffuseFlux.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanRadiationModel.h"
#include "Geoid.h"

TVector3 icemc::Source::DiffuseFlux::pickNeutrinoDirection(const OpticalPath& opticalPath){


  /**
   * We have the 
   * 
   */
  const Geoid::Position& start = opticalPath.steps.at(0).start;
  const Geoid::Position pointLocalXTowards(Geoid::Pole::South);
  const TVector3& rfDir = opticalPath.steps.at(0).direction.Unit(); ///@todo better interface

  // here we construct a local coordinate system to pick the neutrino direction
  // the z-axis is along the RF direction.
  LocalCoordinateSystem lc(start, pointLocalXTowards, rfDir);
  
  //@todo return WEIGHT this with the 1./cos_theta_range factor
  const double thetaCherenkov = AskaryanRadiationModel::CHANGLE_ICE;
  // const double theta = pickUniform(0, TMath::Pi());
  const double phi = pickUniform(0, TMath::TwoPi());

  TVector3 nuDir;
  nuDir.SetMagThetaPhi(1.0, thetaCherenkov, phi);

  nuDir = lc.localTranslationToGlobal(nuDir);
  
  return nuDir.Unit();

}
