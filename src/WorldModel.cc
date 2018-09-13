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











































icemc::Mesh::Mesh(){

}

icemc::Mesh::~Mesh(){
  if(fSurfUp){
    delete fSurfUp;
  }
  if(fSurfDown){
    delete fSurfDown;
  }  
}

size_t icemc::Mesh::addPoint(const Geoid::Position& p, double val){
  
  if(fDoneInit){
    icemcLog() << icemc::warning << "Can't add a point after calling distanceZ!" << std::endl;
    return N();
  }
  Geoid::Position p2 = p;
  p2.SetMag(p.Surface());
  
  if(p.Z() >= 0){
    fXUps.push_back(p.X());
    fYUps.push_back(p.Y());
    fZUps.push_back(val);
  }
  else{
    fXDowns.push_back(p.X());
    fYDowns.push_back(p.Y());
    fZDowns.push_back(val);
  }
  return N();
}

double icemc::Mesh::eval(const Geoid::Position& p) const {
  
  
  if(N() == 0){
    icemcLog() << icemc::warning << "Can't interpolate with 0 points" << std::endl;
    return TMath::QuietNaN();
  }

  if(!fDoneInit){
    init();
  }
  
  Geoid::Position p2 = p;
  p2.SetMag(p.Surface());
  if(p.Z() >= 0){
    return fSurfUp->Interpolate(p2.X(), p2.Y());
  }
  else{
    fSurfDown->Interpolate(p2.X(), p2.Y());
  }
  return 0;
}

void icemc::Mesh::init() const {

  
  fSurfUp = new ROOT::Math::Delaunay2D(fXUps.size(),   &fXUps[0],   &fYUps[0],   &fZUps[0]);
  fSurfDown = new ROOT::Math::Delaunay2D(fXDowns.size(), &fXDowns[0], &fYDowns[0], &fZDowns[0]);  
  

  fDoneInit = true;
}



// double icemc::Surface::distanceZ(const TVector3& p) const {

//   if(N() == 0){
//     icemcLog() << icemc::warning << "Can't interpolate with 0 points" << std::endl;
//     return TMath::QuietNaN();
//   }

//   if(!fDoneInit){
//     init();
//   }
  
//   if(p.Z() >= 0){
//     return p.Z() - fSurfUp->Interpolate(p.X(), p.Y());
//   }
//   else{
//     return p.Z() - fSurfDown->Interpolate(p.X(), p.Y());
//   }
// }


// double icemc::Surface::distanceR(const TVector3& p) const {

//   /**
//    * Here we try to get the radial distance above the surface
//    * We interpolate in z(x,y) but we want 3D deltaR.
//    * 
//    * Since the surface covers the globe, there will be a value of dr which cuts the surface.
//    * At the point we cut the surface, the position, p_c, will be (x_c, y_c, distanceZ(x_c, y_c))
//    * Therefore we need to scale the p until this is true (or close enough).
//    * 
//    * The key to make this efficient is to make intelligent guesses as to the scale factor in the loop iteration.
//    */

//   const double epsilon = 0.001; // meters, i.e. 1mm

//   TVector3 p2 = p;

//   double dz = distanceZ(p2);
//   const double cosTheta = p2.CosTheta();
  
//   int iter = 0;
//   while(TMath::Abs(dz) > epsilon){

//     // if we're close to the surface near the south pole, dz will be close to dr.
//     double scaleFactor = 1 - (distanceZ(p2)/p2.Z()*cosTheta);
//     p2 *= scaleFactor;
//     dz = distanceZ(p2);

//     iter++;
//     if(iter >= 1000){
//       printf("After %d iterations, ", iter);
//       printf("sf = %8.7lf, ", scaleFactor);
//       printf("dz = %8.7lf, ", dz);
//       printf("p2.Z() = %8.7lf, ", p2.Z());
//       printf("\n");
//     }
//     if(iter >= 2000){
//       std::cout <<  "Something went wrong in " << __PRETTY_FUNCTION__ << " check this!" << std::endl;
//       exit(1);
//     }
//   }
//   // std::cout << "Converged after " << iter << " iterations" << std::endl;
//   return p.Mag() - p2.Mag();
// }
