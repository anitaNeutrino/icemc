#include "LoopInfo.h"
#include "Report.h"

ClassImp(icemc::LoopInfo);
ClassImp(icemc::LoopInfo::Step);

bool operator==(bool b, icemc::LoopInfo::Step s){
  return s==b;
}
bool operator!=(bool b, icemc::LoopInfo::Step s){
  return s!=b;
}

std::ostream& operator<<(std::ostream& os, icemc::LoopInfo::Step s){
  if(s.wasEvaluated()==true){
    if(s.getResult()==true){
      os << "true";
    }
    else{
      os << "false";
    }
  }
  else{
    os << "unknown";
  }
  return os;
}


bool icemc::LoopInfo::Step::operator==(bool b) const {
  return evaluated ? result == b : false;
}

bool icemc::LoopInfo::Step::operator!=(bool b) const {
  return !(*this==b);
}

icemc::LoopInfo::Step& icemc::LoopInfo::Step::operator=(bool b){
  result = b;
  evaluated = true;
  return *this;
}


bool icemc::LoopInfo::Step::wasEvaluated() const {return evaluated;}
bool icemc::LoopInfo::Step::getResult() const {
  if(wasEvaluated()==false){
    report() << severity::warning << "Getting result for Step that was not evaluated!" << std::endl;
  }
  return result;
}

void icemc::LoopInfo::setWeights(double volumeFraction, double dtheta){

  positionWeight = volumeFraction > 0 ? 1/volumeFraction : DBL_MAX; //@reasonable if no ice visible?
  
  if (dtheta > TMath::Pi()/2)
    dtheta = TMath::Pi()/2;

  if (dtheta == -1 || dtheta == 0) //@todo is this reasonable?
    directionWeight = 0;
  else
    directionWeight = 2/(1-cos(dtheta));
}

double icemc::LoopInfo::dTheta() const {
  if(directionWeight == 0)
    return 0;
  else
    return acos(1-2/directionWeight);
}
