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

double icemc::WorldModel::areaLonLat(const Geoid::Position& p, double d){

  // Simple conversion of lat into theta for use in math
  double minTheta = TMath::DegToRad()*(90 - (p.Latitude()+d/2));
  double maxTheta = minTheta + d*TMath::DegToRad();
  
  double geoidRadius = Geoid::getGeoidRadiusAtLatitude(p.Latitude());
  return  geoidRadius*geoidRadius*d*TMath::DegToRad()*(TMath::Cos(minTheta)-TMath::Cos(maxTheta));
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


size_t icemc::Mesh::addPoint(const Geoid::Position& p2, double val){
  
  if(fDoneBuild){
    icemc::report() << severity::warning << "Can't add a point after calling build!" << std::endl;
    return N();
  }
  Geoid::Position p = p2;
  p.SetMag(p2.EllipsoidSurface());


  /**
   * To follow what's going on here you need to image a sphere cut in half.
   * There's only one sphere of data but we're going to represent it with 6 hemispheres (Delaunays).
   * We take the same sphere and cut it in half 3 ways:
   * along the x-y plane, the x-z plane, and the y-z plane.
   * Annoyingly, this does triple the amount of memory we need, but it ensures no edge effects.
   * 
   * i.e. we get a point (x,y,z) and a value v.
   * If z>0 we put it in the fUpZ category, otherwise fDownZ
   * If y>0 we put it in the fUpY category, otherwise fDownY
   * If x>0 we put it in the fUpX category, otherwise fDownx
   */


  if(p.X() >= 0){
    fUpX.fPoints.add(p, val);
  }
  else{
    fDownX.fPoints.add(p, val);
  }

  if(p.Y() >= 0){
    fUpY.fPoints.add(p, val);
  }
  else{
    fDownY.fPoints.add(p, val);    
  }
  

  if(p.Z() >= 0){
    fUpZ.fPoints.add(p, val);
  }
  else{
    fDownZ.fPoints.add(p, val);
  }

  fNumPoints++;
  
  return N();
}


double icemc::Mesh::eval(const Geoid::Position& p2) const {

  if(fDoneBuild==false){
    return buildWarning(__PRETTY_FUNCTION__);
  }
  
  Geoid::Position p = p2;
  p.SetMag(p2.EllipsoidSurface());

  double interpolatedMeshVal = 0;

  double X = TMath::Abs(p.X());
  double Y = TMath::Abs(p.Y());
  double Z = TMath::Abs(p.Z());
  
  if(Z >= X && Z >= Y){
    // for z we interpolate v in x-y
    if(p.Z() >= 0){
      interpolatedMeshVal = fUpZ.fDelaunay->Interpolate(p.X(),  p.Y());
    }
    else{
      interpolatedMeshVal = fDownZ.fDelaunay->Interpolate(p.X(),  p.Y());
    }
  }
  else if(X >= Y && X >= Z){
    // for x we interpolate v in y-z
    if(p.X() >= 0){
      interpolatedMeshVal = fUpX.fDelaunay->Interpolate(p.Y(),  p.Z());
    }
    else{
      interpolatedMeshVal = fDownX.fDelaunay->Interpolate(p.Y(),  p.Z());
    }
  }
  else if(Y >= X && Y >= Z){
    // for y we interpolate v in z-x
    if(p.Y() >= 0){
      interpolatedMeshVal = fUpY.fDelaunay->Interpolate(p.Z(),  p.X());
    }
    else{
      interpolatedMeshVal = fDownY.fDelaunay->Interpolate(p.Z(),  p.X());
    }
  }
  else{
    icemc::report() << severity::error << "lon=" << p.Longitude() << " lat=" << p.Latitude() << "  outside interpolation bounds!" << std::endl;
    exit(1);
  }

  return interpolatedMeshVal;
}


double icemc::Mesh::buildWarning(const char* funcName) const {
  icemc::report() << severity::warning << "Can't call " << funcName << "without first calling build! returning 0." << std::endl;
  return 0;
}

void icemc::Mesh::build() {

  if(!fDoneBuild){

    fUpX.build();
    fDownX.build();
    fUpY.build();
    fDownY.build();
    fUpZ.build();
    fDownZ.build();    
    
    fDoneBuild = true;
  }
}


bool icemc::Mesh::Hemisphere::build(){

  //@todo remove these redeclarations from WorldModel.h, needed these to get it to compile
  const int minPoints = 2;
  static const int nDim=3;
  static const UInt_t kdTreeBinSize = 1e6;
  if(fPoints.v.size() > minPoints){

    fKDTree = std::make_shared<TKDTreeID>(fPoints.v.size(), nDim, kdTreeBinSize);
    fKDTree->SetData(0, &fPoints.x[0]);
    fKDTree->SetData(1, &fPoints.y[0]);
    fKDTree->SetData(2, &fPoints.z[0]);
    fKDTree->Build();    

    fDelaunay = std::make_shared<Delaunay2D>(fPoints.v.size(),
					     fAxis==Axis::x ? &fPoints.y[0] : fAxis==Axis::y ? &fPoints.z[0] : &fPoints.x[0],
					     fAxis==Axis::x ? &fPoints.z[0] : fAxis==Axis::y ? &fPoints.x[0] : &fPoints.y[0],
					     &fPoints.v[0]);

    return true;    
  }
  else{
    icemc::report() << severity::error << "Can't initialize with only " << fPoints.v.size() << " points!" << std::endl;
    return false;
  }
}


double icemc::Mesh::findMaxWithinDistance(const Geoid::Position &p, double distanceMeters) const {

  if(fDoneBuild==false){
    return buildWarning(__PRETTY_FUNCTION__);
  }

  double maxVal = -DBL_MAX;
  bool foundAny = false;
  for(auto h : {&fUpX, &fDownX, &fUpY, &fDownY, &fUpZ, &fDownZ}){
    if(h->fKDTree!=nullptr){
      double xyz[3];
      p.GetXYZ(xyz);
      std::vector<int> indices;
      h->fKDTree->FindInRange(xyz, distanceMeters, indices);

      for(int i : indices){
	foundAny = true;	
	double val = h->fPoints.v[i];
	if(val > maxVal){
	  maxVal = val;
	}
      }
    }
  }

  if(!foundAny){
    icemc::report() << severity::warning << "Found no points within " << distanceMeters << "m of " << p << ", returning 0" << std::endl;
    maxVal = 0;
  }
  return maxVal;
}
