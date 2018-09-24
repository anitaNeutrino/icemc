#include "NeutrinoPath.h"
#include "WorldModel.h"
#include "LocalCoordinateSystem.h"
#include "AskaryanFactory.h"
#include "Crust2.h"

// ClassImp(icemc::NeutrinoPath);

icemc::NeutrinoPath::NeutrinoPath(){

  reset();
}


icemc::NeutrinoPath::NeutrinoPath(const Geoid::Position& interaction, const TVector3& rfDir, const WorldModel* m)
  : fW(m), fInteractionPos(interaction)

{
  
  // fColumnDensity = ((icemc::Crust2*)fW)->integratePath(fInteractionPos, fNeutrinoDir);
}









void icemc::NeutrinoPath::reset(){
  // zero all members
  theta_in = 0;
  lat_in = 0;
  nearthlayers = 0;
  weight_prob = 0.;
  weight1 = 0;
  weight = 0.;
  logweight = 0.;
  len_int = 0;
  pieceofkm2sr = 0;
}




