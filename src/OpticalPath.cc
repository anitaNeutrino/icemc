#include "OpticalPath.h"
#include "LocalCoordinateSystem.h"

double icemc::OpticalPath::distance() const {
  double totalDistance = 0;
  for(const auto& step : steps){
    totalDistance += step.distance();
  }
  return totalDistance;
}

double icemc::OpticalPath::attenuation() const {
  double totalAttenuation = 1;
  for(const auto& step : steps){
    totalAttenuation *= step.attenuation();
  }
  return totalAttenuation;
}


///@todo make a unified fresnel namespace set of functions
///@todo test me properly
TVector3 updatePolarization(const TVector3& pol,
			    const TVector3 &surfaceNormal,
			    const TVector3 &incoming,
			    double n_incoming,
			    const TVector3 &outgoing,
			    double n_outgoing) {

  // Here we construct a new coordinate system based on the orientation of the refraction
  // the bending plane is the y-z plane.
  // z-axis is normal to the surface, with z increasing in the direction of the outgoing ray
  const TVector3 localZ = surfaceNormal.Dot(outgoing) > 0 ? surfaceNormal.Unit() : -surfaceNormal.Unit();
  const TVector3 localX = outgoing.Cross(localZ).Unit(); // perpendicular to the bending plane, parallel to surface
  // const TVector3 localY = localZ.Cross(localX).Unit(); // would complete a RH coordinate system

  const double cosTheta_i = localZ.Dot(incoming.Unit());
  const double cosTheta_t = localZ.Dot(outgoing.Unit());

  const double n1_cosTheta_i = n_incoming*cosTheta_i;
  const double n2_cosTheta_t = n_outgoing*cosTheta_t;

  const double n2_cosTheta_i = n_outgoing*cosTheta_i;
  const double n1_cosTheta_t = n_incoming*cosTheta_t;

  // s is for E-field perpendicular to plane of incidence (a.k.a bending plane)
  // p is for E-field parallel to plane of incidence (a.k.a bending plane)

  const double r_s = (n1_cosTheta_i - n2_cosTheta_t)/(n1_cosTheta_i + n2_cosTheta_t);
  // const double t_s = r_s + 1;
  
  const double r_p = (n2_cosTheta_i - n1_cosTheta_t)/(n2_cosTheta_i + n1_cosTheta_t);
  // const double t_p = (r_p + 1)*n_incoming/n_outgoing;

  // std::cout << "************************************************" << "\n"; 
  // std::cout << "*\t" << surfaceNormal << "\n";
  // std::cout << "*\t" << localZ << "\n";
  // std::cout << "*\t" << incoming << "\n";
  // std::cout << "*\t" << outgoing << "\n";  
  // std::cout << "*\t" << cosTheta_i << "\t" << cosTheta_t << "\n";
  // std::cout << "*\t" << n_incoming <<  "\t" << n_outgoing << "\n";
  // std::cout << "*\t" << n1_cosTheta_i << "\t" << n2_cosTheta_t << "\n";
  // std::cout << "*\t" << n2_cosTheta_i << "\t" << n1_cosTheta_t << "\n";    
  // std::cout << "*\t" << r_s << "\t" << r_p  <<  "\n";
  // std::cout << "*\t" << t_s << "\t" << t_p  <<  "\n";
  // std::cout << "*\t" << TMath::Sqrt(1 - r_s*r_s) << "\t" << TMath::Sqrt(1 - r_p*r_p)  <<  "\n";  
  // std::cout << "************************************************" << std::endl;  

  const double polPerp = pol.Dot(localX); //component of incoming pol perpendicular to plane of incidence
  // Mag2() is the square of the length of the TVector3
  const double polPara = TMath::Sqrt(pol.Mag2() - polPerp*polPerp);
  
  const double polPerpOutgoing = TMath::Sqrt(1 - r_s*r_s) * polPerp;
  const double polParaOutgoing = TMath::Sqrt(1 - r_p*r_p) * polPara;
  // const double polPerpOutgoing = t_s * polPerp;
  // const double polParaOutgoing = t_p * polPara;

  // in the limit that the incoming is totally parallel to the surfaceNormal, then newPara would be the same as localY
  // therefore do z.cross(x) to get our new y.
  TVector3 newParaVector = outgoing.Unit().Cross(localX).Unit();

  TVector3 transmittedPol = polPerpOutgoing*localX + polParaOutgoing*newParaVector;

  return transmittedPol;
}
//end GetFresnel()



TVector3 icemc::OpticalPath::polarization(const TVector3& initialPolarization) const {

  TVector3 finalPolarization = initialPolarization;
  
  for(int i=1; i < steps.size(); i++){
    if(i < steps.size()){
      const Step& s1 = steps.at(i-1);
      const Step& s2 = steps.at(i);
      finalPolarization = updatePolarization(finalPolarization, s1.boundaryNormal,
					     s1.direction(), s1.n,
					     s2.direction(), s2.n);
    }

  }
  return finalPolarization;  
};


