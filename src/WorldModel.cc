#include "Detector.h"
#include "Report.h"
#include "WorldModel.h"
#include "LocalCoordinateSystem.h"
#include <memory>


/** 
 * Plane is defined by three points
 * 
 * @param p1 is the first point
 * @param p2 is the second points
 * @param p3 is the third point
 * 
 * @return vector normal to the plane surface
 */
inline TVector3 normalToPlane(const Geoid::Position& p1,  const Geoid::Position& p2,  const Geoid::Position& p3){
  TVector3 d3 = p3 - p1;
  TVector3 d2 = p2 - p1;

  return d3.Cross(d2).Unit();
}

/** 
 * Get the area of the triangle defined by the three points.
 * 
 * @param p1 
 * @param p2 
 * @param p3 
 * 
 * @return 
 */
inline double area(const Geoid::Position& p1, const Geoid::Position& p2, const Geoid::Position& p3){
  TVector3 normal = normalToPlane(p1, p2, p3);
  // length of normal gives area of parallelogram of vectors, therefore divide by 2
  return 0.5*normal.Mag();
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





// double icemc::WorldModel::CreateHorizons(Detector* detector,  double horizonDistanceMeters, double timeStepSeconds){

//   if(fHorizons.GetN()==0){
  
//     if(timeStepSeconds <= 0){
//       icemc::report() << severity::error << "got timeStepSeconds = " << timeStepSeconds << ", setting to default " << icemc::WorldModel::defaultTimeStep << std::endl;
//       timeStepSeconds = defaultTimeStep;
//     }

//     double startTime = detector->getStartTime();
//     double endTime = detector->getEndTime();
    
//     const int nSteps = (endTime - startTime)/timeStepSeconds;

//     icemc::report() << severity::info << "Creating horizons for " << nSteps << " time samples" << std::endl;
    
//     for(int step=0; step <= nSteps; step++){
//       double sampleTime = startTime + step*timeStepSeconds;

//       Geoid::Position pos = detector->getPosition(sampleTime);

//       double volume = IceVolumeWithinHorizon(pos, horizonDistanceMeters);

//       fHorizons.SetPoint(fHorizons.GetN(), sampleTime, volume);

//       std::cout << "\r" << step << " / " << nSteps << std::flush;
//     }
//     std::cout << std::endl;
//   }

//   std::cout << fHorizons.GetN() << std::endl;
//   exit(1);
  
//   return 0;
// }






































size_t icemc::Mesh::addPoint(const Geoid::Position& p, double val){
  
  if(fDoneInit){
    icemc::report() << severity::warning << "Can't add a point after calling eval!" << std::endl;
    return N();
  }
  Geoid::Position p2 = p;
  p2.SetMag(p.EllipsoidSurface());
  fPoints.emplace_back(Point(p, val));
  return N();
}


double icemc::Mesh::eval(const Geoid::Position& p) const {

  if(N() == 0){
    icemc::report() << severity::warning << "Can't interpolate with 0 points" << std::endl;
    return TMath::QuietNaN();
  }

  if(!fDoneInit){
    init();
  }
  
  Geoid::Position p2 = p;
  p2.SetMag(p.EllipsoidSurface());

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
    icemc::report() << severity::error << "You shouldn't get here!" << std::endl;
  }

  return interpolatedMeshVal;
}

void icemc::Mesh::init() const {

  if(!fDoneInit){

    for(int i=0; i < N(); i++){
      const Geoid::Position& p = fPoints.at(i).position;
      const double val = fPoints.at(i).value;
      
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

    const int minPoints = 2;
    if(xUpX.size() > minPoints){
      fXUp = std::unique_ptr<Delaunay2D>(new Delaunay2D(xUpX.size(), &xUpX[0], &xUpY[0], &xUpZ[0]));
    }

    if(xDownX.size() > minPoints){
    fXDown = std::unique_ptr<Delaunay2D>(new Delaunay2D(xDownX.size(), &xDownX[0], &xDownY[0], &xDownZ[0]));
    }

    if(yUpX.size() > minPoints){
      fYUp = std::unique_ptr<Delaunay2D>(new Delaunay2D(yUpX.size(), &yUpX[0], &yUpY[0], &yUpZ[0]));
    }

    if(yDownX.size() > minPoints){
      fYDown = std::unique_ptr<Delaunay2D>(new Delaunay2D(yDownX.size(), &yDownX[0], &yDownY[0], &yDownZ[0]));
    }

    if(zUpX.size() > minPoints){
    fZUp = std::unique_ptr<Delaunay2D>(new Delaunay2D(zUpX.size(), &zUpX[0], &zUpY[0], &zUpZ[0]));
    }

    if(zDownX.size() > minPoints){
      fZDown = std::unique_ptr<Delaunay2D>(new Delaunay2D(zDownX.size(), &zDownX[0], &zDownY[0], &zDownZ[0]));
    }

    
    fDoneInit = true;
  }
}
