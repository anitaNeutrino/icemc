/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: strutt@physics.ucla.edu

 Description:
             This is a test program for the new input.yaml file format
***********************************************************************************************************/

#include <iostream>


#include "Constants.h"
#include "Settings.h"
#include "position.hh"

#include "earthmodel.hh"
#include "Tools.h"
#include "vector.hh"
#include "roughness.hh"
#include "anita.hh"
#include "balloon.hh"
#include "icemodel.hh"
#include "Spectra.h"
#include "signal.hh"
#include "secondaries.hh"
#include "ray.hh"
#include "counting.hh"
#include "Primaries.h"
#include "Taumodel.hh"


int main(){

  Settings s;


  std::ofstream outputsFile("/tmp/outputs.txt");

  Balloon *bn1=new Balloon(); // instance of the balloon
  Anita *anita1=new Anita();// right now this constructor gets banding info
  Secondaries *sec1=new Secondaries();
  Signal *sig1=new Signal();
  Ray *ray1=new Ray(); // create new instance of the ray class
  // input parameters
  int NNU;
  double RANDOMISEPOL;

  s.ReadInputs("inputs.anita3.conf",  outputsFile, NNU, RANDOMISEPOL);
  s.ApplyInputs(anita1,  sec1,  sig1,  bn1,  ray1);



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

  delete ray1;
  delete sig1;
  delete sec1;
  delete anita1;
  delete bn1;


  return 0;



}
