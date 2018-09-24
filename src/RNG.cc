#include "RNG.h"
#include "Report.h"

#include "TRandom3.h"
#include <typeinfo>

/**
 * static member initialization
 */
int icemc::RNG::gRun(-1);
UInt_t icemc::RNG::gEventNumber(0);
std::map<std::string, unsigned int> icemc::RNG::fNameCount;

void icemc::RNG::newSeeds(int run, UInt_t eventNumber){
  gRun = run;
  gEventNumber = eventNumber;
}



icemc::RNG::RNG() :
  fRandom(std::unique_ptr<TRandom3>(new TRandom3())){  
}

icemc::RNG::~RNG(){

}


void icemc::RNG::updateSeed(){

  if(fID==-1){
    std::string name = typeid(*this).name();
    fID = fNameCount[name]++;
    icemc::report() << severity::info << "Initilized RNG: " << name << " with ID " << fID << std::endl;
  }

  if(fRun != gRun || fEventNumber != gEventNumber){
    fRun = gRun;
    fEventNumber = gEventNumber;
    std::string hashMe = typeid(*this).name() + std::to_string(fID)  + std::to_string(fRun) + std::to_string(fEventNumber);
    static std::hash<std::string> hashFcn;
    size_t hash = hashFcn(hashMe);
    fRandom->SetSeed((unsigned int)hash);
  }
}

double icemc::RNG::pickUniform(double x1){
  updateSeed();
  return fRandom->Uniform(x1);
}

double icemc::RNG::pickUniform(double x1, double x2){
  updateSeed();  
  return fRandom->Uniform(x1, x2);
}

double icemc::RNG::pickGaus(double mean, double sigma){
  updateSeed();  
  return fRandom->Gaus(mean, sigma);
}


int icemc::RNG::pickPoisson(double mean){
  updateSeed();  
  return fRandom->Poisson(mean);
}
