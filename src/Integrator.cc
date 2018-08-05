#include "Integrator.h"
#include "TString.h"
#include <iostream>

icemc::Integrator::Integrator(int run){

  TString hashStr = GetName();
  hashStr += TString::Format("_%d", run);
  hash = hashStr.Hash();

  std::cout << hashStr <<  "\t"  << hash << std::endl;
  
  fRandom.SetSeed(run);  
}


void icemc::Integrator::reset(){
  hash++;
  fRandom.SetSeed(hash);
}
