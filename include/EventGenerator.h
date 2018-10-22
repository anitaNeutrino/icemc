#ifndef ICEMC_EVENT_GENERATOR_H
#define ICEMC_EVENT_GENERATOR_H

#include "TVector3.h"

#include "Settings.h"

#include "Constants.h"
#include "CommandLineOptions.h"
#include "RNG.h"
#include "Event.h"
#include "PropagatingSignal.h"

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

    void generate(Detector& detector);
    // bool nextEvent();

  private:
    void initialize();
    std::shared_ptr<Source::EnergyModel> fSourceEnergyModel = nullptr;
    std::shared_ptr<Source::DirectionModel> fSourceDirectionModel = nullptr;

    void delayAndAddSignalToEachRX(const PropagatingSignal& signal, const OpticalPath& opticalPath, Detector& detector) const;


    
    bool fOrderedEventTimes = true;
    Event fEvent;
  };
}


#endif //ICEMC_EVENT_GENERATOR_H
