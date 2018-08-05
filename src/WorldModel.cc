#include "Detector.h"
#include "IcemcLog.h"
#include "WorldModel.h"
#include "LocalCoordinateSystem.h"

inline TVector3 normalToPlane(const Geoid::Position& p1,  const Geoid::Position& p2,  const Geoid::Position& p3){
  TVector3 d3 = p3 - p1;
  TVector3 d2 = p2 - p1;

  return d3.Cross(d2).Unit();
}






TVector3 icemc::WorldModel::GetSurfaceNormal(const Geoid::Position& p) const {
  const double d = 1; // delta meters
  // which coordinates do we want to vary?
  TVector3 delta1;
  TVector3 delta2;

  if(p.Z() >= p.X() && p.Z() >= p.Y()){
    delta1.SetXYZ(d, 0, 0);
    delta2.SetXYZ(0, d, 0);
  }
  else if(p.Y() >= p.Z() && p.Y() >= p.X()){
    delta1.SetXYZ(0, 0, d);
    delta2.SetXYZ(d, 0, 0);
  }
  else if(p.X() >= p.Y() && p.X() >= p.Z()){
    delta1.SetXYZ(0, d, 0);
    delta2.SetXYZ(0, 0, d);
  }
  
  Geoid::Position pPlusDelta1 = p + delta1;
  double alt1 = SurfaceAboveGeoid(pPlusDelta1);
  pPlusDelta1.SetAltitude(alt1);

  Geoid::Position pPlusDelta2 = p + delta2;
  double alt2 = SurfaceAboveGeoid(pPlusDelta2);
  pPlusDelta2.SetAltitude(alt2);
      
  TVector3 dv1 = (pPlusDelta1 - p).Unit();
  TVector3 dv2 = (pPlusDelta2 - p).Unit();
  TVector3 surfaceNorm = dv1.Cross(dv2);

  // make the vector radial
  int sign = p.Dot(surfaceNorm) > 0 ? 1 : -1;
  surfaceNorm *= sign;
  
  return surfaceNorm.Unit();
}





double icemc::WorldModel::CreateHorizons(Detector* detector,  double horizonDistanceMeters, double timeStepSeconds){

  if(fHorizons.GetN()==0){
  
    if(timeStepSeconds <= 0){
      icemcLog() << icemc::error << "got timeStepSeconds = " << timeStepSeconds << ", setting to default " << icemc::WorldModel::defaultTimeStep << std::endl;
      timeStepSeconds = defaultTimeStep;
    }

    double startTime = detector->getStartTime();
    double endTime = detector->getEndTime();

    
    const int nSteps = (endTime - startTime)/timeStepSeconds;

    icemcLog() << icemc::info << "Creating horizons for " << nSteps << " time samples" << std::endl;
    
    for(int step=0; step <= nSteps; step++){
      double sampleTime = startTime + step*timeStepSeconds;

      Geoid::Position pos = detector->getPosition(sampleTime);

      double volume = IceVolumeWithinHorizon(pos, horizonDistanceMeters);

      fHorizons.SetPoint(fHorizons.GetN(), sampleTime, volume);

      std::cout << "\r" << step << " / " << nSteps << std::flush;
    }
    std::cout << std::endl;
  }

  std::cout << fHorizons.GetN() << std::endl;
  exit(1);
  
  return 0;
}


double icemc::WorldModel::IceVolumeWithinHorizon(const Geoid::Position& p,  double horizonDistanceMeters, double stepSizeMeters) const {

  const LocalCoordinateSystem localCoords(p);
  double localVolume = 0;
  const double drSq = horizonDistanceMeters*horizonDistanceMeters;
  const double kilometersToMeters = 1e3;
  for(double dy=-horizonDistanceMeters; dy <= horizonDistanceMeters; dy += stepSizeMeters){
    for(double dx=-horizonDistanceMeters; dx <= horizonDistanceMeters; dx += stepSizeMeters){
      if(dx*dx + dy*dy <= drSq){
	TVector3 deltaLocal(dx, dy, 0);
	TVector3 deltaGlobal = localCoords.localTranslationToGlobal(deltaLocal);
	double iceThickness = kilometersToMeters*this->IceThickness(p + deltaGlobal);
	localVolume += stepSizeMeters*stepSizeMeters*iceThickness;
      }
    }
  }
  return localVolume;
  
}
