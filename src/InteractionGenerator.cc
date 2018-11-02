#include "InteractionGenerator.h"
#include "Constants.h"
#include "Settings.h"
#include "Antarctica.h"
#include "Neutrino.h"
#include "LocalCoordinateSystem.h"
// #include "TGraphAntarctica.h" ///@todo remove after debugging
// #include "TFile.h" ///@todo remove after debugging

icemc::InteractionGenerator::InteractionGenerator(const Settings *settings, std::shared_ptr<WorldModel> worldModel) :
  fSettings(settings), fWorldModel(worldModel)
{

}






Geoid::Position icemc::InteractionGenerator::pickInteractionPosition(const Geoid::Position &detector) {

  Geoid::Position nuPos;

  const double MAX_HORIZON_DIST = 800e3; ///@todo get me from settings?
  double localMaxIceThickness = fWorldModel->maxIceThicknessWithinDistance(detector, MAX_HORIZON_DIST);

  // This is a pretty important statement about icemc, we sample the ICE uniformly.
  // To do that we first sample the x/y plane uniformly within the horizon radius,
  // then we weight randomly chosen positions by ice thickness so that x/y positions
  // where the ice is twice as thick are twice as likely to be chosen.

  LocalCoordinateSystem lc(detector);  

  // int numTries = 0;
  bool iceThickEnough = false;
  while(!iceThickEnough){
    double dx=1, dy=1;
    while(dx*dx + dy*dy > 1){
      dx = pickUniform(-1, 1);
      dy = pickUniform(-1, 1);
    }

    dx*=MAX_HORIZON_DIST;
    dy*=MAX_HORIZON_DIST;
    
    TVector3 deltaPos(dx, dy, 0);
    nuPos = lc.localPositionToGlobal(deltaPos);
    
    // then roll a dice to see if it's deep enough
    double iceThicknessHere = fWorldModel->IceThickness(nuPos);
    double randomThickness = pickUniform()*localMaxIceThickness;

    // numTries++;

    if(iceThicknessHere >= randomThickness){
      iceThickEnough = true;
      double surfaceElevation = fWorldModel->SurfaceAboveGeoid(nuPos);
      nuPos.SetAltitude(surfaceElevation - randomThickness);

      // std::cout << "Picked an interaction position after " <<  numTries << " tries..." << std::endl;
      std::cout << "nuPos = " << nuPos.Longitude() << ", " << nuPos.Latitude() << ", " << nuPos.Altitude() << std::endl;
      // std::cout << surfaceElevation << ", " << randomThickness << std::endl;
    }
  }

  // static TGraphAntarctica* grD = new TGraphAntarctica();
  // static TGraphAntarctica* grN = new TGraphAntarctica();
  
  // if(grN){grN->SetPoint(grN->GetN(), nuPos.Longitude(), nuPos.Latitude());}
  // if(grD){grD->SetPoint(grD->GetN(), detector.Longitude(), detector.Latitude());}

  // static int nCount = 0;
  // nCount++;

  // if(nCount==5000){
  //   TFile* grTest = new TFile("pickPosition.root", "recreate");
  //   grD->SetName("grDetector");
  //   grN->SetName("grNeutrino");    
  //   grD->Write();
  //   grN->Write();
  //   delete grD; grD = nullptr;
  //   delete grN; grN = nullptr;
  //   grTest->Write();
  //   grTest->Close();
  //   delete grTest;
  // }
  
  
  return nuPos;
}






// int icemc::InteractionGenerator::PickDownwardInteractionPoint(const Geoid::Position&r_bn, const Settings *settings1, const Antarctica *antarctica1) {

//   if(settings1->UNBIASED_SELECTION==1){
//     if(antarctica1->PickUnbiased(this)){ // pick neutrino direction and interaction point
//       dtryingdirection=1.;
//       iceinteraction=1;
//     }
//     else{
//       iceinteraction=0;
//     }
//   }
//   else{
//     iceinteraction=1;

//     // If we require neutrinos from a particular position
//     // we generate that cartesian position here

//     static Geoid::Position specific_position;

//     if (settings1->SPECIFIC_NU_POSITION) 
//       {
// 	specific_position.SetLonLatAlt(settings1->SPECIFIC_NU_POSITION_LONGITUDE,
// 				       settings1->SPECIFIC_NU_POSITION_LATITUDE,
// 				       settings1->SPECIFIC_NU_POSITION_ALTITUDE);
//       }

//     do{
//       ///@todo ibnposition
//       posnu = antarctica1->pickInteractionPosition(r_bn);
//     } while(settings1->SPECIFIC_NU_POSITION &&  (posnu - specific_position).Mag() > settings1->SPECIFIC_NU_POSITION_DISTANCE);
//   }
  
//   return 0;
// }//PickDownwardInteractionPoint





/**
 * Choose CC or NC: get from ratios in Ghandi etal paper,
 * updated for the CTEQ6-DIS parton distribution functions (M.H. Reno, personal communication).
 * Need to add capability of using ratios from Connolly et al.
 */

icemc::Neutrino::Interaction::Current icemc::InteractionGenerator::pickCurrent() {
  double rnd = pickUniform();
  if (rnd<=0.6865254){ // 10^18 eV - 10^21 eV (use this one for ANITA)
//if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    return Neutrino::Interaction::Current::Charged;//"cc";
  }
  else{
    return Neutrino::Interaction::Current::Neutral;//"nc";  
  }
} //GetCurrent
