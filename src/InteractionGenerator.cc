#include "InteractionGenerator.h"
#include "Constants.h"
#include "Settings.h"
#include "Antarctica.h"
#include "Neutrino.h"


icemc::InteractionGenerator::InteractionGenerator(const Settings *settings, std::shared_ptr<WorldModel> worldModel) :
  fSettings(settings), fWorldModel(worldModel)
{

}



Geoid::Position icemc::InteractionGenerator::pickInteractionPosition(const Geoid::Position &detector) {
  //@todo restore me!
  Geoid::Position nuPos;

  // ///@todo unbreak this!
  // const double MAX_HORIZON_DIST = 800e3; ///@todo Get this from settings or something?
  // bool iceThickEnough = false;

  // const int numDimensions = 3;
  // double pCart[numDimensions];
  // detector.GetXYZ(pCart);
  // std::vector<int> indicesOfPointsWithinRange;
  // fKDTree->FindInRange(pCart, MAX_HORIZON_DIST, indicesOfPointsWithinRange);


  // ///@todo put this inside mesh?
  // std::sort(indicesOfPointsWithinRange.begin(), indicesOfPointsWithinRange.end());
  // double localMax = -999999;
  // {
  //   const auto& thickness = fThicknesses.at(Layer::Ice);
  //   auto point = thickness.begin();
  //   int i=0;
  //   for(int index : indicesOfPointsWithinRange){
  //     while (i != index){
  // 	i++;
  // 	++point;
  //     }

  //     if(point->value > localMax){
  // 	localMax = point->value;
  //     }
  //   }
  // }

  // This is a pretty important statement about icemc, we sample the ICE uniformly.
  // To do that we first sample the x/y plane uniformly within the horizon radius,
  // then we weight randomly chosen positions by ice thickness so that x/y positions
  // where the ice is twice as thick are twice as likely to be chosen.
  // int numTries = 0;
  // while(!iceThickEnough){
  //   double dx=1, dy=1;
  //   while(dx*dx + dy*dy > 1){
  //     dx = pickUniform(-1, 1);
  //     dy = pickUniform(-1, 1);
  //   }

  //   dx*=MAX_HORIZON_DIST;
  //   dy*=MAX_HORIZON_DIST;
    
  //   TVector3 deltaPos(dx, dy, 0);
  //   nuPos = detector + deltaPos;
    
  //   // then roll a dice to see if it's deep enough
  //   double iceThicknessHere = fWorldModel->IceThickness(nuPos);
  //   double randomHeight = pickUniform()*localMax;

  //   numTries++;
  //   // std::cout << numTries << "\t" << iceThicknessHere << std::endl;

  //   if(randomHeight < iceThicknessHere){
  //     iceThickEnough = true;
  //     double surfaceElevation = fWorldModel->SurfaceAboveGeoid(nuPos);
  //     nuPos.SetAltitude(surfaceElevation - randomHeight);
  //   }
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
