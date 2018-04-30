#include "GeneratedNeutrino.h"

ClassImp(icemc::GeneratedNeutrino)
ClassImp(icemc::PassingNeutrino)


icemc::GeneratedNeutrino::GeneratedNeutrino(int ithLoop)
: inu(ithLoop)
{
  weight = -1;
  passCutNoWay = -1;
  passCut2 = -1;
  passCut3 = -1;
  passCutWithinHorizon = -1;
}

icemc::GeneratedNeutrino::~GeneratedNeutrino(){

}




icemc::PassingNeutrino::PassingNeutrino(const icemc::GeneratedNeutrino& gNu, const icemc::RadioSignal& radioSignal)
  : GeneratedNeutrino(gNu), signal(radioSignal)
{
  
}




icemc::PassingNeutrino::~PassingNeutrino(){

}
