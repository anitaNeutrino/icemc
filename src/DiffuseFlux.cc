#include "DiffuseFlux.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanRadiationModel.h"
#include "Geoid.h"

TVector3 icemc::Source::DiffuseFlux::pickNeutrinoDirection(const OpticalPath& opticalPath, double dtheta){

  const Geoid::Position& start = opticalPath.steps.at(0).start;
  const Geoid::Position pointLocalXTowards(Geoid::Pole::South);
  const TVector3& rfDir = opticalPath.steps.at(0).direction().Unit(); ///@todo better interface

  // here we construct a local coordinate system to pick the neutrino direction
  // the z-axis is along the RF direction.
  LocalCoordinateSystem lc(start, pointLocalXTowards, rfDir);
  
  //@todo return WEIGHT this with the 1./cos_theta_range factor
  //const double thetaCherenkov = AskaryanRadiationModel::CHANGLE_ICE;
  
  const double theta = pickUniform(AskaryanRadiationModel::CHANGLE_ICE-dtheta, AskaryanRadiationModel::CHANGLE_ICE+dtheta);
  const double phi = pickUniform(0, TMath::TwoPi());
  //const double phi = 0; ///@attention test condition, remove me!

  TVector3 nuDir;
  //nuDir.SetMagThetaPhi(1.0, thetaCherenkov, phi);
  nuDir.SetMagThetaPhi(1.0, theta, phi);

  nuDir = lc.localTranslationToGlobal(nuDir);

  // std::cout << nuDir.Angle(rfDir) << std::endl;
  
  return nuDir.Unit();

}
