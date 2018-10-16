#ifndef ICEMC_EVENT_GENERATOR_H
#define ICEMC_EVENT_GENERATOR_H

#include "TVector3.h"

#include "Settings.h"

#include "Constants.h"
#include "CommandLineOptions.h"
#include "ShowerModel.h"
#include "RNG.h"

namespace icemc {

  class Taumodel;
  class AskaryanRadiationModel;
  class Interaction;
  class RayTracer;
  class Roughness;
  class Screen;
  class ConnollyEtAl2011;
  class ANITA;

  namespace Source {
    class Spectra;
  }

  class Detector;

  /**
   * @class EventGenerator
   * @brief The glue that brings the simulation together
   * 
   * Contains the main neutrino generation loop in generateNeutrinos()
   */
  class EventGenerator : public RNG {
  private:
    const Settings* fSettings;
  public:

    // enum RayDirection {
    //   direct = 0,
    //   downgoing = 1,
    // };

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


    //do a threshold scan
    double threshold_start=-1.;
    double threshold_end=-6.;
    static const int NTHRESHOLDS=20;
    double threshold_step=(threshold_end-threshold_start)/(double)NTHRESHOLDS;

    double npass_v_thresh[NTHRESHOLDS]={0.};
    double denom_v_thresh[NTHRESHOLDS]={0.};
    double npass_h_thresh[NTHRESHOLDS]={0.};
    double denom_h_thresh[NTHRESHOLDS]={0.};
    double thresholds[NTHRESHOLDS];

  private:

    

  };
}



#endif //ICEMC_EVENT_GENERATOR_H
