#include "DiffuseFlux.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanRadiationModel.h"
#include "Geoid.h"

TVector3 icemc::Source::DiffuseFlux::pickNeutrinoDirection(const OpticalPath& opticalPath, const double dtheta){

  const Geoid::Position& start = opticalPath.steps.at(0).start;
  const Geoid::Position pointLocalXTowards(Geoid::Pole::North);
  const TVector3& rfDir = opticalPath.steps.at(0).direction().Unit(); ///@todo better interface

  // here we construct a local coordinate system to pick the neutrino direction
  // the z-axis is along the RF direction.
  LocalCoordinateSystem lc(start, pointLocalXTowards, rfDir);

  if (dtheta <= 0){
    directionWeight = DBL_MAX; // Since we divide by weight this weight out this neutrino
  }
  else{
    double theta2 = AskaryanRadiationModel::CHANGLE_ICE - dtheta > 0 ? AskaryanRadiationModel::CHANGLE_ICE - dtheta : 0;
    double theta1 = AskaryanRadiationModel::CHANGLE_ICE + dtheta < TMath::Pi() ? AskaryanRadiationModel::CHANGLE_ICE + dtheta : TMath::Pi();

    directionWeight = 2/(cos(theta2)-cos(theta1));
  }  

  const double theta = pickUniform(AskaryanRadiationModel::CHANGLE_ICE-dtheta, AskaryanRadiationModel::CHANGLE_ICE+dtheta);
  const double phi = pickUniform(0, TMath::TwoPi());

  TVector3 nuDir;
  //nuDir.SetMagThetaPhi(1.0, AskaryanRadiationModel::CHANGLE_ICE, 0);
  //nuDir.SetMagThetaPhi(1.0, AskaryanRadiationModel::CHANGLE_ICE, phi);
  nuDir.SetMagThetaPhi(1.0, theta, phi);

  nuDir = lc.localTranslationToGlobal(nuDir);

  // std::cout << nuDir.Angle(rfDir) << std::endl;
  
  return nuDir.Unit();
}

double icemc::Source::DiffuseFlux::getDirectionWeight(){
  return directionWeight;
}
