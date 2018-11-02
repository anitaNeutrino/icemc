#include "EventGenerator.h"

// system includes
#include "signal.h"

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
#include "InteractionGenerator.h"



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


  report() << severity::info << "EXPONENT = " << fSettings->EXPONENT << std::endl;

  if(fSettings->EXPONENT > 10 && fSettings->EXPONENT < 10){ ///@todo fix this
    fSourceEnergyModel = std::make_shared<Source::MonoEnergetic>(1e20*Energy::Unit::eV);
  }
  else{
    fSourceEnergyModel = std::make_shared<Source::Spectra>(fSettings);
  }

  ///@todo add support for directional sources
  fSourceDirectionModel = std::make_shared<Source::DiffuseFlux>();

  if(fSettings->YPARAM==0){ ///@todo enumify YPARAM
    fYGenerator = std::make_shared<GhandiEtAl::YGenerator>();
  }
  else if(fSettings->YPARAM==1){
    fYGenerator = std::make_shared<ConnollyEtAl2011::YGenerator>(fSettings);            
  }
  else{
    icemc::report() << severity::error << "Unknown YPARAM!" << std::endl;
  }

  
  if(fSettings->SIGMAPARAM==0){ ///@todo enumify SIGMAPARAM
    fCrossSectionModel = std::make_shared<MHReno::CrossSectionModel>(fSettings);
  }
  else if(fSettings->SIGMAPARAM==1){
    fCrossSectionModel = std::make_shared<ConnollyEtAl2011::CrossSectionModel>(fSettings);
  }
  else{
    icemc::report() << severity::error << "Unknown SIGMAPARAM!" << std::endl;
  }


  


  
}


/** 
 * Destructor
 * 
 */
icemc::EventGenerator::~EventGenerator()
{

  
  
}



void icemc::EventGenerator::delayAndAddSignalToEachRX(const PropagatingSignal& signal, const OpticalPath& opticalPath, Detector& detector) const {
  double tofDetector = opticalPath.steps.back().distance()/constants::CLIGHT;
  ///@todo this interface needs to be improved
  Geoid::Position rfExit = fEvent.detector - opticalPath.steps.back().direction();

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
}


void icemc::EventGenerator::generate(Detector& detector){

  icemc::report() << severity::info << "Seed is " << fSettings->SEED << std::endl;


  
  ShowerModel showerModel(fSettings);
  
  auto antarctica = std::make_shared<Antarctica>();
  auto interactionGenerator = std::make_shared<InteractionGenerator>(fSettings, antarctica);

  icemc::report() << "Area of the earth's surface covered by antarctic ice is " << antarctica->ice_area << std::endl;
  
  icemc::RootOutput output(this, fSettings, fSettings->getOutputDir(), fSettings->getRun());

  int n;
  double dt;
  detector.getDesiredNDt(n, dt);
  AskaryanRadiationModel askaryanModel(fSettings, n, dt);

  gRandom->SetSeed(fSettings->SEED); // settings seed is now updated with run_no to uniquify it elsewhere

  NeutrinoFactory nuFactory(fSettings, fSourceEnergyModel, fSourceDirectionModel, fCrossSectionModel, fYGenerator, antarctica, interactionGenerator);
  RayTracer rayTracer(antarctica);

  int numEvents = fSettings->NNU - fSettings->getStartNu();
  std::vector<double> eventTimes;
  eventTimes.reserve(numEvents);
  for(int e=0; e < numEvents; e++){
    eventTimes.push_back(pickUniform(detector.getStartTime(),  detector.getEndTime()));
  }
  if(fOrderedEventTimes==true){
    std::sort(eventTimes.begin(), eventTimes.end());
  }
  
  
  
  /**
   * Main loop over generated neutrinos
   */
  signal(SIGINT, icemc::EventGenerator::interrupt_signal_handler);  

  int run = fSettings->getRun();
  for(int entry=0; entry < eventTimes.size() && ABORT_EARLY==false; entry++){
    
    UInt_t eventNumber = (UInt_t)(fEvent.loop.run)*fSettings->NNU+entry;
    fEvent = Event(run, eventNumber,  eventTimes.at(entry));
    
    RNG::newSeeds(run, eventNumber); // updates all RNG seeds
    gRandom->SetSeed(eventNumber+6e7); ///@todo kill me
    
    fEvent.detector = detector.getPosition(fEvent.loop.eventTime);
    Geoid::Position interactionPos = interactionGenerator->pickInteractionPosition(fEvent.detector);

    OpticalPath opticalPath = rayTracer.findPath(interactionPos, fEvent.detector);
    fEvent.loop.rayTracingSolution = opticalPath.residual < 1; // meters
    
    if(fEvent.loop.rayTracingSolution==false){
      output.allTree.Fill();
      continue;
    }
    
    fEvent.neutrino = nuFactory.makeNeutrino(opticalPath);
    
    fEvent.shower = showerModel.generate(fEvent.neutrino);
    
    PropagatingSignal signal = askaryanModel.generate(fEvent.neutrino, fEvent.shower, opticalPath);

    // std::cout << signal.maxEField() << "\n";
    
    signal.propagate(opticalPath);

    // std::cout << signal.maxEField() << "\n"  << std::endl;

    fEvent.loop.chanceInHell = false; //detector.chanceInHell(signal);

    if(fEvent.loop.chanceInHell==false){
      output.allTree.Fill();
      continue;
    }
    
    delayAndAddSignalToEachRX(signal, opticalPath, detector);

    fEvent.loop.passesTrigger = detector.applyTrigger();

    if(fEvent.loop.passesTrigger==true){
      output.passTree.Fill();
    }
    output.allTree.Fill();
  }
  signal(SIGINT, SIG_DFL); /// unset signal handler
}
