#include "NeutrinoPath.h"

ClassImp(icemc::NeutrinoPath);


icemc::NeutrinoPath::NeutrinoPath(){

  reset();
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
