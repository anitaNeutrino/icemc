#ifndef ICEMC_EVENT_GENERATOR_H
#define ICEMC_EVENT_GENERATOR_H

#include "TVector3.h"

#include "Settings.h"

#include "Constants.h"
#include "CommandLineOptions.h"
#include "ShowerModel.h"
#include "RNG.h"

namespace icemc {

  namespace Source {
    class EnergyModel;
    class DirectionModel;
  }

  class Detector;
  class Neutrino;

  /**
   * @class EventGenerator
   * @brief The glue that brings the simulation together
   * 
   * Contains the main neutrino generation loop in generate()
   */
  class EventGenerator : public RNG {
  private:
    friend class RootOutput;
    const Settings* fSettings;
  public:

    EventGenerator(const Settings* settings);
    virtual ~EventGenerator();
    static void interrupt_signal_handler(int);  // This catches ctrl-C interrupt (SIGINT)

    /** 
     * @brief Run the neutrino generation
     * 
     * This function does the stuff that used to be the main in the icemc executable
     * @todo needs some love to constify... probably not doable...
     */    
    void generate(Detector& detector);

  private:
    void initialize();
    std::shared_ptr<Source::EnergyModel> fSourceEnergyModel = nullptr;
    std::shared_ptr<Source::DirectionModel> fSourceDirectionModel = nullptr;


    UInt_t fEventNumber;
    Int_t fRun;
    Neutrino fNeutrino;
    Shower fShower;
    Double_t fEventTime;
    Geoid::Position fDetectorPos;
  };
}



#endif //ICEMC_EVENT_GENERATOR_H
