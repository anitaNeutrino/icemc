#include "SourceModel.h"
#include "EarthModel.h"

TVector3 icemc::DiffuseFlux::pickDirection(const Geoid::Position& interaction, double refractiveIndex, const TVector3* rfToDetector) {

  TVector3 nuDir(0,  0, 1);

  if(rfToDetector && rfToDetector->Mag() > 0){

    nuDir = *rfToDetector;

    // right now the nuDir points directly where the RF needs to go...
    // we must pick the neutrino direction such that 

    // so we need to nudge it a certain angle away
    ///@todo randomly choose off cone angle in theta, 
    // double thetaOffConeRad = 0;

    // this gives us an orthogonal unit vector
    const TVector3 orthogonalToRfDir = (*rfToDetector).Orthogonal();
    const double cherenkovAngle = TMath::ACos(1/refractiveIndex);

    const double deltaTheta = cherenkovAngle;

    nuDir.Rotate(deltaTheta, orthogonalToRfDir);

    // just a check I can picture 3D vector things in my head.
    const double epsilon = 1e-10;
    if(TMath::Abs(nuDir.Angle(*rfToDetector) - deltaTheta) > epsilon){
      std::cerr << "You can't do 3D vector things in your head!" << std::endl;
    }

    double phi = fRandom.Uniform(0, TMath::TwoPi());

    nuDir.Rotate(phi, *rfToDetector);    
    
    
  }
  else{
    double cosTheta = fRandom.Uniform(-1, 1);
    double phi = fRandom.Uniform(0, TMath::TwoPi());
    
    nuDir.SetMagThetaPhi(1, TMath::ACos(cosTheta), phi);
    
    // pick unbiasaed
  }
  return nuDir;  
}
