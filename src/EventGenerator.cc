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
#include "MonoEnergetic.h"
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
  initialize();
}


void icemc::EventGenerator::initialize(){

  ///@todo can the details of this be hidden?
  if(fSettings->EXPONENT > 10 && fSettings->EXPONENT < 10){
    fSourceEnergyModel = std::make_shared<Source::MonoEnergetic>(1e20*Energy::Unit::eV);
  }
  else{
    fSourceEnergyModel = std::make_shared<Source::Spectra>(fSettings);
  }

  fSourceDirectionModel = std::make_shared<Source::DiffuseFlux>();
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

  NeutrinoFactory nuFactory(fSettings,  fSourceEnergyModel, fSourceDirectionModel);
  RayTracer rayTracer(antarctica);
  
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
    
    fEventTime = pickUniform(detector.getStartTime(), detector.getEndTime()); 
    fDetectorPos = detector.getPosition(fEventTime);
    Geoid::Position interactionPos = antarctica->pickInteractionPosition(fDetectorPos);

    OpticalPath opticalPath = rayTracer.findPath(fDetectorPos, interactionPos);
    
    if(opticalPath.residual > 1){
      ///@todo better interface for no good solution?
      // std::cout << "bad" << std::endl;
      // no solution was found withing the fitter tolerance
      ///@todo Fill some tree indicating this!      
      ro.allTree.Fill();
      continue;
    }
    
    fNeutrino = nuFactory.makeNeutrino(opticalPath);
    
    fShower = showerModel.generate(fNeutrino);
    PropagatingSignal signal = askaryanModel.generate(fNeutrino, fShower, opticalPath);
    signal.propagate(opticalPath);
    
    ///@todo put this in some useful function tucked away somewhere...
    double tofDetector = opticalPath.steps.back().distance()/constants::CLIGHT;
    Geoid::Position rfExit = fDetectorPos - opticalPath.steps.back().direction;

    std::vector<double> delays(detector.getNumRX());
    for(int rx = 0; rx < detector.getNumRX(); rx++){
      delays.at(rx) = (detector.getPositionRX(rx) - rfExit).Mag()/constants::CLIGHT - tofDetector;
    }
    double minDelay = *std::min_element(delays.begin(), delays.end());
    for(int rx=0; rx < detector.getNumRX(); rx++){
      PropagatingSignal signalRX = signal;
      signalRX.waveform.applyConstantGroupDelay(delays.at(rx) - minDelay, false);
      detector.addSignalToRX(signalRX, rx);
    }

    bool triggered =  detector.applyTrigger();

    if(triggered){
      ro.passTree.Fill();
      std::cout << inu << "\tPASSED!" << std::endl;
    }
  }
  signal(SIGINT, SIG_DFL); /// unset signal handler
}
