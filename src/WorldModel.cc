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

  std::vector<ROOT::Math::Delaunay2D*> meshes {fXUp, fXDown, fYUp, fYDown, fZUp, fZDown};
  for(auto m : meshes){
    if(m){
      delete m;
    }
  }
}

size_t icemc::Mesh::addPoint(const Geoid::Position& p, double val){
  
  if(fDoneInit){
    icemcLog() << icemc::warning << "Can't add a point after calling eval!" << std::endl;
    return N();
  }
  Geoid::Position p2 = p;
  p2.SetMag(p.Surface());
  fPositions.emplace_back(p2);
  fValues.emplace_back(val);  
  return N();
}


inline void nudgeInsideBounds(Geoid::Position& p, ROOT::Math::Delaunay2D* d){

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

  double interpolatedMeshVal = 0;

  double X = TMath::Abs(p2.X());
  double Y = TMath::Abs(p2.Y());
  double Z = TMath::Abs(p2.Z());  
  
  if(Z >= X && Z >= Y){
    if(p2.Z() >= 0){
      interpolatedMeshVal = fZUp->Interpolate(p2.X(),  p2.Y());
    }
    else{
      interpolatedMeshVal = fZDown->Interpolate(p2.X(),  p2.Y());
    }
  }
  else if(X >= Y && X >= Z){
    if(p2.X() >= 0){
      interpolatedMeshVal = fXUp->Interpolate(p2.Y(),  p2.Z());
    }
    else{
      interpolatedMeshVal = fXDown->Interpolate(p2.Y(),  p2.Z());
    }
  }
  else if(Y >= X && Y >= Z){
    if(p2.Y() >= 0){
      interpolatedMeshVal = fYUp->Interpolate(p2.X(),  p2.Z());
    }
    else{
      interpolatedMeshVal = fYDown->Interpolate(p2.X(),  p2.Z());
    }    
  }
  else{
    icemcLog() << icemc::error << "You shouldn't get here!" << std::endl;
  }

  return interpolatedMeshVal;
}

void icemc::Mesh::init() const {

  if(!fDoneInit){

    for(int i=0; i < N(); i++){
      const Geoid::Position& p = fPositions.at(i);
      const double val = fValues.at(i);
      
      if(p.X() >= 0){
	xUpX.push_back(p.Y());
	xUpY.push_back(p.Z());
	xUpZ.push_back(val);
      }
      else{
	xDownX.push_back(p.Y());
	xDownY.push_back(p.Z());
	xDownZ.push_back(val);
      }

      if(p.Y() >= 0){
	yUpX.push_back(p.X());
	yUpY.push_back(p.Z());	
        yUpZ.push_back(val);
      }
      else{
	yDownX.push_back(p.X());
	yDownY.push_back(p.Z());
	yDownZ.push_back(val);
      }
      
      if(p.Z() >= 0){
	zUpX.push_back(p.X());
	zUpY.push_back(p.Y());	
        zUpZ.push_back(val);
      }
      else{
	zDownX.push_back(p.X());
	zDownY.push_back(p.Y());
	zDownZ.push_back(val);
      }            
    }

    fXUp = new ROOT::Math::Delaunay2D(xUpX.size(), &xUpX[0], &xUpY[0], &xUpZ[0]);
    fXDown = new ROOT::Math::Delaunay2D(xDownX.size(), &xDownX[0], &xDownY[0], &xDownZ[0]);    

    fYUp = new ROOT::Math::Delaunay2D(yUpX.size(), &yUpX[0], &yUpY[0], &yUpZ[0]);
    fYDown = new ROOT::Math::Delaunay2D(yDownX.size(), &yDownX[0], &yDownY[0], &yDownZ[0]);    

    fZUp = new ROOT::Math::Delaunay2D(zUpX.size(), &zUpX[0], &zUpY[0], &zUpZ[0]);
    fZDown = new ROOT::Math::Delaunay2D(zDownX.size(), &zDownX[0], &zDownY[0], &zDownZ[0]);    

    fDoneInit = true;
  }
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
