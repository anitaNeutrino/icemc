#include "Integrator.h"
#include "TString.h"

icemc::Integrator::Integrator(int run){

  TString hashStr = GetName();
  hashStr += TString::Format("_%d");
  hash = hashStr.Hash();

  std::cout << hashStr <<  "\t"  << hash << std::endl;
  
  fRandom.SetSeed(run);  
}

icemc::Integrator::reset(){
  hash++;
  fRandom.SetSeed(hash);
}
