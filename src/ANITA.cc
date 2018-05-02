#include "ANITA.h"

icemc::ANITA::ANITA(){
  TVector v;
  testVecNotRealYet.push_back(v);
}

icemc::ANITA::~ANITA(){
  // does nothing
}


/** 
 * @todo Make the match the current icemc implementation
 * 
 * @param n ANITA's desired n (for trigger and digitizer paths both!)
 * @param dt ANITA's desired dt (for trigger and digitizer paths both!)
 */
void icemc::ANITA::getNDtForTimeDomainAskaryanSignals(int& n, double& dt) const {
  n = 10000;
  dt = 0.01;
}



icemc::GeographicCoordinate icemc::ANITA::getCenterOfDetector(){
  GeographicCoordinate gc;
  return gc;
}

const std::vector<TVector>&  icemc::ANITA::getFieldCalcLocationsRelativeToAveDetPos(){
  return testVecNotRealYet;
}

bool icemc::ANITA::applyTrigger(const std::vector<TGraph>& pureSignalVoltageTimeGraphs, const TVector& poyntingVector, const TVector& polarizationVector){
  return true;
}

