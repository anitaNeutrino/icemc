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




icemc::PassingNeutrino::PassingNeutrino(const icemc::GeneratedNeutrino& gNu, const icemc::AskaryanFreqs& askFreqs)
  : GeneratedNeutrino(gNu), askaryanFreqs(askFreqs)
{
  
}




icemc::PassingNeutrino::~PassingNeutrino(){

}
