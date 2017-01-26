/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             This is a test program for the new input.yaml file format
***********************************************************************************************************/

#include <iostream>


#include "Settings.h"

int main(){

  Settings s("inputs.yaml");
  s.printAllKeyValuePairStrings();

  std::cout << std::endl << std::endl;

  const char* whichPayloadKey = "Which payload";
  int whichPayload = 0;
  s.getSetting(whichPayloadKey, whichPayload);
  std::cout << whichPayloadKey << " -> " << whichPayload << std::endl;



  const char* bandThresholdsKey = "Thresholds for each band";
  std::vector<double> bandThresholds;
  s.getSetting(bandThresholdsKey, bandThresholds);

  std::cout << bandThresholdsKey << " -> ";
  for(unsigned int i=0; i < bandThresholds.size(); i++){
    std::cout << bandThresholds.at(i) << " ";
  }
  std::cout << std::endl;






  return 0;



}
