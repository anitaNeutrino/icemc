#include "InteractionGenerator.h"
#include "Constants.h"
#include "Settings.h"
#include "Antarctica.h"
#include "Neutrino.h"
#include "LocalCoordinateSystem.h"
#include "CrossSectionModel.h"
#include "Inelasticity.h"
#include "Report.h"

icemc::InteractionGenerator::InteractionGenerator(const Settings *settings,
						  std::shared_ptr<WorldModel> worldModel,
						  std::shared_ptr<CrossSectionModel> crossSectionModel,
						  std::shared_ptr<YGenerator> yGenerator						  
						  ) :
  fSettings(settings),
  fWorldModel(worldModel),
  fCrossSectionModel(crossSectionModel),
  fYGenerator(yGenerator)
{

  if(fSettings->SPECIFIC_NU_POSITION > 0){
    fSpecificInteractionCenter.SetLonLatAlt(fSettings->SPECIFIC_NU_POSITION_LONGITUDE,
					    fSettings->SPECIFIC_NU_POSITION_LATITUDE,
					    fSettings->SPECIFIC_NU_POSITION_ALTITUDE);
    std::cout << "SPECIFIC ALTITUDE " << fSpecificInteractionCenter << std::endl;    
  }
}



Geoid::Position icemc::InteractionGenerator::pickInteractionPosition(const Geoid::Position &detector) {

  if(fSettings->SPECIFIC_NU_POSITION > 0){
    // if we need a specific place, then don't pick with reference to the detector, instead pick with reference to the specific place
    return pickInteractionPosition(fSpecificInteractionCenter, fSettings->SPECIFIC_NU_POSITION_DISTANCE);
  }
  else{
    return pickInteractionPosition(detector, fSettings->MAX_HORIZON_DISTANCE);
  }
}


Geoid::Position icemc::InteractionGenerator::pickInteractionPosition(const Geoid::Position &center, double rangeMeters) {

  Geoid::Position interactionPosition;

  double localMaxIceThickness = fWorldModel->maxIceThicknessWithinDistance(center, fSettings->MAX_HORIZON_DISTANCE);

  // This is a pretty important statement about icemc, we sample the ICE uniformly.
  // To do that we first sample the x/y plane uniformly within the horizon radius,
  // then we weight randomly chosen positions by ice thickness so that x/y positions
  // where the ice is twice as thick are twice as likely to be chosen.
  LocalCoordinateSystem lc(center);  

  int numTries = 0;
  const int maxTries = 10000;
  while(numTries < maxTries){
    double dx=1, dy=1;
    while(dx*dx + dy*dy > 1){
      dx = pickUniform(-1, 1);
      dy = pickUniform(-1, 1);
    }

    dx*=rangeMeters;
    dy*=rangeMeters;
    
    TVector3 deltaPos(dx, dy, 0);
    interactionPosition = lc.localPositionToGlobal(deltaPos);
    
    // then roll a dice to see if it's deep enough
    double iceThicknessHere = fWorldModel->IceThickness(interactionPosition);
    double randomThickness = pickUniform()*localMaxIceThickness;

    numTries++;

    if(iceThicknessHere >= randomThickness){
      double surfaceElevation = fWorldModel->SurfaceAboveGeoid(interactionPosition);
      interactionPosition.SetAltitude(surfaceElevation - randomThickness);
      std::cout << "surf/thick\t"  <<  surfaceElevation << "\t" << randomThickness << "\t" << localMaxIceThickness << std::endl;

      // std::cout << "Picked an interaction position after " <<  numTries << " tries..." << std::endl;
      // std::cout << "interaction = " << interactionPosition.Longitude() << ", " << interactionPosition.Latitude() << ", " << interactionPosition.Altitude() << std::endl;
      break;
    }
  }

  if(numTries >= maxTries){
    icemc::report() << severity::error << "Failed to pick a neutrino position after " << numTries << "!" << std::endl;
  }   
  return interactionPosition;
}




icemc::Interaction icemc::InteractionGenerator::generate(const Neutrino& n, const Geoid::Position& detector){

  Interaction interaction;
  interaction.current = pickCurrent();
  interaction.crossSection = fCrossSectionModel->getSigma(n.energy, n.leptonNumber, interaction.current);
  interaction.length = CrossSectionModel::getInteractionLength(interaction.crossSection);
  interaction.y = fYGenerator->pickY(n.energy, n.leptonNumber, interaction.current);
  interaction.position = pickInteractionPosition(detector);
  // std::cout << interaction.position << std::endl;
  return interaction;
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

icemc::Interaction::Current icemc::InteractionGenerator::pickCurrent() {
  double rnd = pickUniform();
  if (rnd<=0.6865254){ // 10^18 eV - 10^21 eV (use this one for ANITA)
//if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    return Interaction::Current::Charged;//"cc";
  }
  else{
    return Interaction::Current::Neutral;//"nc";  
  }
} //GetCurrent
