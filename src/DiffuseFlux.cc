#include "DiffuseFlux.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanRadiationModel.h"
#include "Geoid.h"

TVector3 icemc::Source::DiffuseFlux::pickNeutrinoDirection(const OpticalPath& opticalPath, const LoopInfo& loop){


  const double dtheta = loop.dTheta;
  double Ch_angle = loop.inFirn? AskaryanRadiationModel::CHANGLE_FIRN : AskaryanRadiationModel::CHANGLE_ICE;

  const Geoid::Position& start = opticalPath.steps.at(0).start;
  const Geoid::Position pointLocalXTowards(Geoid::Pole::North);
  const TVector3& rfDir = opticalPath.steps.at(0).direction().Unit(); ///@todo better interface

  // here we construct a local coordinate system to pick the neutrino direction
  // the z-axis is along the RF direction.
  LocalCoordinateSystem lc(start, pointLocalXTowards, rfDir);

  double theta2 = Ch_angle - dtheta > 0 ? Ch_angle - dtheta : 0;
  double theta1 = Ch_angle + dtheta < TMath::Pi() ? Ch_angle + dtheta : TMath::Pi();

  // Since we divide by weight this weight out neutrino that can't be detected
  directionWeight = dtheta > 0 ? 2/(cos(theta2)-cos(theta1)) : DBL_MAX; 

  const double theta = pickUniform(theta2, theta1);
  const double phi = pickUniform(0, TMath::TwoPi());
    
  TVector3 nuDir;
  //nuDir.SetMagThetaPhi(1.0, Ch_angle-dtheta, 0);
  //nuDir.SetMagThetaPhi(1.0, Ch_angle, phi); // Force signal to be on-cone
  nuDir.SetMagThetaPhi(1.0, theta, phi);
  
  nuDir = lc.localTranslationToGlobal(nuDir);
  //std::cout << "nudir=" << nuDir << std::endl;

  return nuDir.Unit();
}

double icemc::Source::DiffuseFlux::getDirectionWeight(){
  return directionWeight;
}
