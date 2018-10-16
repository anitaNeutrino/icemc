#include "EventGenerator.h"

// system includes
#include "signal.h"

// ROOT includes
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TVector3.h"

// libAntarcticaRoot includes
#include "Geoid.h"

// icemc includes
#include "RootOutput.h"
#include "Constants.h"
#include "Settings.h"
#include "Crust2.h"
#include "Antarctica.h"
#include "Spectra.h"
#include "AskaryanRadiationModel.h"
#include "ShowerModel.h"
#include "RayTracer.h"
#include "ConnollyEtAl2011.h"
#include "Report.h"
#include "Neutrino.h"
#include "NeutrinoFactory.h"
#include "DiffuseFlux.h"
#include "RNG.h"
#include "TRandom3.h" ///@todo remove when gRandom is excized

#include <string>
#include <sstream>

#if __cplusplus > 199711L
#define isnan std::isnan 
#include <type_traits>
#endif

#include <typeinfo>
#include <fenv.h> 

bool ABORT_EARLY = false;    // This flag is set to true when interrupt_signal_handler() is called

void icemc::EventGenerator::interrupt_signal_handler(int sig){
  signal(sig,  SIG_IGN);
  ABORT_EARLY = true;
  static int numSignals = 0;
  numSignals++;
  const int maxSignals = 3;
  if(numSignals >= maxSignals){
    std::cerr << "Got " <<  numSignals << " interrupts, quitting more aggressively!" << std::endl;
    raise(SIGINT);
  }
  return;
}




/** 
 * Constructor 
 */
icemc::EventGenerator::EventGenerator(const Settings* settings) :  
  fSettings(settings)
{
}


/** 
 * Destructor
 * 
 */
icemc::EventGenerator::~EventGenerator()
{

}


void icemc::EventGenerator::generate(Detector& detector){

  icemc::report() << severity::info << "Seed is " << fSettings->SEED << std::endl;

  ShowerModel showerModel(fSettings);
  ConnollyEtAl2011 crossSectionModel(fSettings);
  auto antarctica = std::make_shared<Antarctica>();

  icemc::report() << "Area of the earth's surface covered by antarctic ice is " << antarctica->ice_area << std::endl;
  
  icemc::RootOutput ro(this, fSettings, fSettings->getOutputDir(), fSettings->getRun());

  int n;
  double dt;
  detector.getDesiredNDt(n, dt);
  AskaryanRadiationModel askaryanModel(fSettings, n, dt);

  int NNU = fSettings->NNU;
  // int RANDOMISEPOL = fSettings->RANDOMISEPOL;
  
  gRandom->SetSeed(fSettings->SEED); // settings seed is now updated with run_no to uniquify it elsewhere

  Source::Spectra nuSpectra(fSettings);
  NeutrinoFactory nuFactory(fSettings);
  Source::DiffuseFlux directionModel;

  /**
   * Main loop over generated neutrinos
   */
  signal(SIGINT, icemc::EventGenerator::interrupt_signal_handler);  
  int run = fSettings->getRun();
  for (int inu = fSettings->getStartNu(); inu < NNU && !ABORT_EARLY; inu++) {
    if (NNU >= 100 && ((inu % (NNU / 100)) == 0)){
      std::cout << inu << " neutrinos. " << (double(inu)/double(NNU)) * 100 << "% complete.\n";
    }
    else if(NNU < 100){
      std::cout << inu << " / " << NNU << "neutrinos" << std::endl;
    }
    
    UInt_t eventNumber=(UInt_t)(run)*NNU+inu;

    RNG::newSeeds(run, eventNumber); // updates all RNG seeds
    gRandom->SetSeed(eventNumber+6e7);
    
    double eventTime = pickUniform(detector.getStartTime(), detector.getEndTime()); 
    Geoid::Position detectorPos = detector.getPosition(eventTime);
    Geoid::Position interactionPos = antarctica->pickInteractionPosition(detectorPos);

    RayTracer rayTracer(antarctica.get(), detectorPos);
    OpticalPath opticalPath = rayTracer.findPathToDetector(interactionPos);
    
    if(opticalPath.residual > 1){
      ///@todo better interface for no good solution?
      // std::cout << "bad" << std::endl;
      // no solution was found withing the fitter tolerance
      ///@todo Fill some tree indicating this!
      continue;
    }
    else{
      // std::cout << "good" << std::endl;
    }
    
    // if we get here, then there's a ray tracing solution
    // to get from our interaction point to the payload
    // now we have that, we can calculate the neutrino path

    const TVector3& rfDirFromInteraction = opticalPath.steps.at(0).direction;
    TVector3 v = directionModel.pickNeutrinoDirection(interactionPos, rfDirFromInteraction);
    
    // make a neutrino, we've picked energy, flavor, interaction current
    Neutrino nu = nuFactory.makeNeutrino();

    nu.path.direction = v.Unit();
    
    Shower shower = showerModel.generate(nu);
    PropagatingSignal signal = askaryanModel.generate(nu, shower, rfDirFromInteraction);
    signal.propagate(opticalPath);
    
    ///@todo put this in some useful function tucked away somewhere...
    double tofDetector = opticalPath.steps.back().distance()/constants::CLIGHT;
    Geoid::Position rfExit = detectorPos - opticalPath.steps.back().direction;


    int closestRX = -1;
    double minTof= DBL_MAX;
    for(int rx=0; rx < detector.getNumRX(); rx++){
      double tofRX = (detector.getPositionRX(rx) - rfExit).Mag()/constants::CLIGHT;
      if(tofRX < minTof){
	minTof = tofRX ;
	closestRX = rx;
      }
    }

    std::cout << "closestRX = " << closestRX << std::endl;
    
    for(int rx = 0; rx < detector.getNumRX(); rx++){      
      double tofRX = (detector.getPositionRX(rx) - rfExit).Mag()/constants::CLIGHT;
      PropagatingSignal signalRX = signal;
      signal.waveform.applyConstantGroupDelay(tofRX - tofDetector, false);
      std::cout << rx << "\t" << 1e9*(tofRX - tofDetector) << std::endl;
      detector.addSignalToRX(signalRX, rx);
    }

    bool triggered =  detector.applyTrigger();

    if(triggered){
      std::cout << inu << "\tPASSED!" << std::endl;
      break;
    }
    // else{
    //   std::cout << "FAILED!" << std::endl;
    // }

  }
  signal(SIGINT, SIG_DFL); /// unset signal handler
}
