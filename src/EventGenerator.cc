#include "EventGenerator.h"

// system includes
#include "signal.h"

// libAntarcticaRoot includes
#include "Geoid.h"

// icemc includes
#include "RootOutput.h"
#include "Constants.h"
#include "Settings.h"
#include "Crust.h"
#include "Antarctica.h"
#include "Spectra.h"
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

#include "Summary.h"

#include "ImpulsiveRadioGenerator.h"
#include "AskaryanRadiationModel.h"
#include "ShowerModel.h"

#include "InteractionGenerator.h"
#include "NeutrinoInteractionGenerator.h"
#include "WaisPulser.h"

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

  if(fSettings->CONSTANTY > 0){
    fYGenerator = std::make_shared<ConstantY>(0.2);
  }
  else if(fSettings->YPARAM==0){ ///@todo enumify YPARAM
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
    std::cout << "~~~~~~~~~~~~~~~" << (int)(100*entry/n) << "% done:  " << entry << " of " << n << " neutrinos ~~~~~~~~~~~~~~~" << std::endl;
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

  int n;
  double dt;
  detector.getDesiredNDt(n, dt);
  fRadioModel = std::make_shared<AskaryanRadiationModel>(fSettings, n, dt);

  ShowerModel showerModel(fSettings);

  auto antarctica = std::make_shared<Antarctica>(fSettings->ICE_MODEL);
  std::shared_ptr<InteractionGenerator> interactionGenerator = std::make_shared<NeutrinoInteractionGenerator>(fSettings,
													      antarctica,
													      fCrossSectionModel,
													      fYGenerator);  
  // auto waisModel = std::make_shared<WaisPulser>(fSettings, antarctica);
  // interactionGenerator = waisModel;
  // fRadioModel = waisModel;

  //icemc::report() << "Volume of antarctic ice is " << antarctica->GetTotalIceVolume() << " m^3" << std::endl;
  //icemc::report() << "Area of the earth's surface covered by antarctic ice is " << antarctica->GetTotalIceArea() << " m^2" << std::endl;

  icemc::RootOutput output(this, fSettings, fSettings->getOutputDir(), fSettings->getRun());

  Summary summary(fSettings->EXPONENT, antarctica->GetTotalIceVolume(), CrossSectionModel::getInteractionLength(fCrossSectionModel->getSigma(fSourceEnergyModel->pickNeutrinoEnergy(), Neutrino::L::Matter, Interaction::Current::Charged)));

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

    UInt_t eventNumber = (UInt_t)(fEvent.loop.run)*fSettings->NNU+entry; ///@todo fix me

    RNG::newSeeds(run, eventNumber); // updates all RNG seeds
    gRandom->SetSeed(eventNumber+6e7); ///@todo kill me

    fEventSummary = EventSummary(run, eventNumber,  eventTimes.at(entry));
    fEventSummary.neutrino = nuGenerator.generate();
    fEventSummary.detector = detector.getPosition(fEventSummary.loop.eventTime);
    fEventSummary.interaction = interactionGenerator->generateInteraction(fEventSummary.neutrino, fEventSummary.detector);
    OpticalPath opticalPath = rayTracer.findPath(fEventSummary.interaction.position, fEventSummary.detector);
    fEventSummary.loop.rayTracingSolution = opticalPath.residual < 1; // meters

    if(false){
      std::cout << "Interaction position: " << fEventSummary.interaction.position << "\n";
      std::cout << "Interaction altitude: " << fEventSummary.interaction.position.Altitude() << "\n";
      std::cout << "Detector: " << fEventSummary.detector << "\n";
      std::cout << "Direction from interaction to balloon:" << (fEventSummary.detector - fEventSummary.interaction.position).Unit() << std::endl;
      std::cout << "Separation = " << fEventSummary.interaction.position.Distance(fEventSummary.detector) << std::endl;
      std::cout << "Optical path distance = " << opticalPath.distance() << ", residual = " << opticalPath.residual << std::endl;
    }

    if(fEventSummary.loop.rayTracingSolution==false){
      //std::cout << "No ray tracing solution between source " << fEventSummary.interaction.position << " and detector " << fEventSummary.detector << ": " << (fEventSummary.interaction.position - fEventSummary.detector).Mag() << std::endl;
      output.allTree().Fill();
      summary.addEvent(fEventSummary);
      continue;
    }

    fEventSummary.loop.dTheta = fSettings->USEDIRECTIONWEIGHTS ? fRadioModel->getThetaRange(detector.signalThreshold(), fEventSummary.neutrino, showerModel.generate(fEventSummary.neutrino, fEventSummary.interaction), opticalPath) : TMath::Pi();
    fEventSummary.neutrino.path.direction = fSourceDirectionModel->pickNeutrinoDirection(opticalPath, fEventSummary.loop.dTheta);

    fEventSummary.shower = showerModel.generate(fEventSummary.neutrino, fEventSummary.interaction);
    PropagatingSignal signal = fRadioModel->generateImpulse(opticalPath, fEventSummary.neutrino, fEventSummary.shower);

    fEventSummary.loop.viewAngle =  (fEventSummary.shower.axis.Angle(opticalPath.steps.at(0).direction().Unit()));
    
    fEvent.signalAt1m = signal.waveform.getTimeDomain();
    fEventSummary.signalSummaryAt1m = signal.summarize();
    
    fEvent.loop.fresnel = signal.propagate(opticalPath);

    fEvent.loop.magnification = opticalPath.magnification();
    fEvent.loop.attenuation = opticalPath.attenuation();

    
    fEvent.signalAtDetector = signal.waveform.getTimeDomain();
    fEventSummary.signalSummaryAtDetector = signal.summarize();
  
    //fEventSummary.loop.chanceInHell = true;
    fEventSummary.loop.chanceInHell = detector.chanceInHell(signal);    

    if(fEventSummary.loop.chanceInHell==false){
      //std::cout << "No chance in hell\t" << fEventSummary.interaction.position << std::endl;
      output.allTree().Fill();      
      summary.addEvent(fEventSummary);
      continue;
    }
    
    delayAndAddSignalToEachRX(signal, opticalPath, detector);
    fEventSummary.loop.passesTrigger = detector.applyTrigger();

    if(fEventSummary.loop.passesTrigger==true){
      std::cout << "Passed trigger!\t" << fEventSummary.interaction.position << std::endl;
      fEventSummary.neutrino.path.integrate(fEventSummary.interaction, antarctica);
      fEventSummary.loop.setPositionWeight(antarctica->IceVolumeWithinHorizon(fEventSummary.detector)/antarctica->GetTotalIceVolume());
      fEventSummary.loop.directionWeight = fSourceDirectionModel->getDirectionWeight();

      fEvent.copy(fEventSummary);
      detector.write(fEvent);
      output.passTree().Fill();
    }
    // Add event to summary, with whether or not it passes trigger
    summary.addEvent(fEventSummary, fEventSummary.loop.passesTrigger.getResult());
    output.allTree().Fill();
  }
  signal(SIGINT, SIG_DFL); /// unset signal handler

  //@todo only works for monoenergetic energy
  summary.summarize();
}
