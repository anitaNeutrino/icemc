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
#include "NeutrinoGenerator.h"
#include "DiffuseFlux.h"
#include "RNG.h"
#include "TRandom3.h" ///@todo remove when gRandom is excized
#include "MonoEnergetic.h"
#include "InteractionGenerator.h"



bool ABORT_EARLY = false;    // This flag is set to true when interrupt_signal_handler() is called

void icemc::EventGenerator::interrupt_signal_handler(int sig){
  signal(sig,  SIG_IGN);
  if(!ABORT_EARLY){
    icemc::report() << severity::info << "Aborting early" << std::endl;
  }
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

  ///@todo I don't like selecting things based on numbers, make this more human readable?
  if(fSettings->EXPONENT > 4 && fSettings->EXPONENT < 30){
    Energy energy(pow(10, fSettings->EXPONENT), Energy::Unit::eV);
    fSourceEnergyModel = std::make_shared<Source::MonoEnergetic>(energy);
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
    icemc::report() << severity::error << "Unknown YPARAM " << fSettings->YPARAM << std::endl;
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



void icemc::EventGenerator::printProgress(int entry, size_t n){

  const int maxPrints = 100;
  const int printEvery = n/maxPrints > 0 ? n/maxPrints : 1;
  if((entry%printEvery)==0){
    std::cout << entry << " of " << n << std::endl;
  }
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
  auto interactionGenerator = std::make_shared<InteractionGenerator>(fSettings,
								     antarctica,
								     fCrossSectionModel,
								     fYGenerator);

  icemc::report() << "Area of the earth's surface covered by antarctic ice is " << antarctica->ice_area << std::endl;

  icemc::RootOutput output(this, fSettings, fSettings->getOutputDir(), fSettings->getRun());

  int n;
  double dt;
  detector.getDesiredNDt(n, dt);
  AskaryanRadiationModel askaryanModel(fSettings, n, dt);

  NeutrinoGenerator nuGenerator(fSettings, fSourceEnergyModel, fSourceDirectionModel, antarctica);
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

    printProgress(entry, eventTimes.size());

    UInt_t eventNumber = (UInt_t)(fEvent.loop.run)*fSettings->NNU+entry;

    RNG::newSeeds(run, eventNumber); // updates all RNG seeds
    gRandom->SetSeed(eventNumber+6e7); ///@todo kill me

    fEvent = Event(run, eventNumber,  eventTimes.at(entry));
    fEvent.neutrino = nuGenerator.generate();
    fEvent.detector = detector.getPosition(fEvent.loop.eventTime);

    fEvent.interaction = interactionGenerator->generate(fEvent.neutrino, fEvent.detector);
    // std::cout << fEvent.interaction.position << std::endl;

    OpticalPath opticalPath = rayTracer.findPath(fEvent.interaction.position, fEvent.detector);
    fEvent.loop.rayTracingSolution = opticalPath.residual < 1; // meters

    if(fEvent.loop.rayTracingSolution==false){
      std::cout << "No ray tracing solution\t" << fEvent.interaction.position << std::endl;
      std::cout << (fEvent.interaction.position - fEvent.detector).Mag() << std::endl;
      output.allTree.Fill();
      continue;
    }

    fEvent.neutrino.path.direction = fSourceDirectionModel->pickNeutrinoDirection(opticalPath);
    fEvent.shower = showerModel.generate(fEvent.neutrino, fEvent.interaction);

    PropagatingSignal signal = askaryanModel.generate(fEvent.neutrino, fEvent.shower, opticalPath);

    fEvent.signalSummary = signal.propagate(opticalPath);

    static double bestMaxE = 0;
    if(signal.maxEField() > bestMaxE){
      bestMaxE = signal.maxEField();
      std::cout << bestMaxE << std::endl;
    }
    // std::cout  << signal.maxEField() << std::endl;

    fEvent.loop.chanceInHell = detector.chanceInHell(signal);

    if(fEvent.loop.chanceInHell==false){
      // std::cout << "No chance in hell\t" << fEvent.interaction.position << std::endl;
      output.allTree.Fill();
      continue;
    }

    delayAndAddSignalToEachRX(signal, opticalPath, detector);

    fEvent.loop.passesTrigger = detector.applyTrigger();

    if(fEvent.loop.passesTrigger==true){
      std::cout << "PASSED!" << std::endl;

      fEvent.neutrino.path.integrate(fEvent.interaction.position, antarctica);
      // std::cout << "PASSED\t" << fEvent.interaction.position << std::endl;

      output.passTree.Fill();
    }
    // std::cout << std::endl;
    output.allTree.Fill();
  }
  signal(SIGINT, SIG_DFL); /// unset signal handler
}
