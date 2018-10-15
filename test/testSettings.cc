/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             This is a test program for the new input.conf file format
***********************************************************************************************************/

#include <iostream>

#include "Constants.h"
#include "Settings.h"

#include "Crust2.h"
#include "Tools.h"
#include "RayTracer.h"
#include "anita.hh"
#include "balloon.hh"
#include "Antarctica.h"
#include "Spectra.h"
#include "AskaryanRadiationModel.h"
#include "ShowerModel.h"
#include "RayTracer.h"
#include "ConnollyEtAl2011.h"
#include "Taumodel.hh"


int main(){

  icemc::Settings s;

  std::ofstream outputsFile("/tmp/outputs.txt");
  s.ReadInputs("inputs.anita3.conf",  outputsFile);
  

  icemc::Balloon *bn1 = new icemc::Balloon(&s); // instance of the balloon
  icemc::Anita *anita1 = new icemc::Anita(&s, "/tmp/", bn1);// right now this constructor gets banding info
  icemc::ShowerModel *sec1 = new icemc::ShowerModel(&s);
  icemc::AskaryanRadiationModel *sig1 = new icemc::AskaryanRadiationModel(&s, 1024, 1e-9/2.6);
  // input parameters
  s.ApplyInputs(anita1);



  //  s.printAllKeyValuePairStrings();

  std::cout << std::endl << std::endl;

  const char* whichPayloadKey = "Which payload";
  int whichPayload = 0;
  s.getSetting(whichPayloadKey, whichPayload);
  std::cout << whichPayloadKey << " -> " << whichPayload << std::endl;


  const char* bandThresholdsKey = "Band thresholds";
  std::vector<double> bandThresholds;
  s.getSetting(bandThresholdsKey, bandThresholds);

  std::cout << bandThresholdsKey << " -> ";
  for(unsigned int i=0; i < bandThresholds.size(); i++){
    std::cout << bandThresholds.at(i) << " ";
  }
  std::cout << std::endl;


  // std::cout << "Now about to try and get a setting that doesn't exist" << std::endl;
  const char* nonExistentKey = "adsfasdfasdfasdf";
  std::cout << "Now about to try and get a setting, called " << nonExistentKey << ", that doesn't exist..." << std::endl;
  std::vector<double> shouldBeEmpty;
  s.getSetting(nonExistentKey, shouldBeEmpty);
  std::cout << std::endl;

  delete sig1;
  delete sec1;
  delete anita1;
  delete bn1;


  return 0;



}
