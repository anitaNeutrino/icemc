#ifndef ICEMC_GENERATED_NEUTRINOS_H
#define ICEMC_GENERATED_NEUTRINOS_H

#include "TObject.h" ///< ClassDef/ClassImp
#include "balloon.hh" /// < For BalloonInfo
#include "AskaryanFreqs.h" ///< For AskaryanFreqs
#include "secondaries.hh" ///< For AskaryanFreqs

namespace icemc {

  /**
   * @class GeneratedNeutrino
   * @brief Each EventGenerator loop makes one of these
   * 
   * This class stores information for each walk though the main loop in EventGenerator.
   * And should store enough information as to why and event did/did not trigger.
   * And not much more than that as there will be MANY in the output trees.
   * (ROOT is good at compression so keep things mostly -1, 0, 1 and you'll be fine).
   * 
   * Conventions:
   * You want some record of whether or not a generated neutrino passed cut.
   * Cuts are applied in sequence and you can leave the EventGenerator loop early 
   * so not all variables will be filled.
   * Therefore "passesCutX" type variables should be filled with -1 on construction.
   * If a cut is failed, set "passesCutX = 0 " and exit the loop.
   * If a cut succeeds, set "passesCut X = 1 " and continue.
   * 
   * 
   * Lives in Rootoutput::allTree;
   */
  class GeneratedNeutrino {
  public:
    /** 
     * Constructor, set all cut tracking variables to -1
     * 
     * @param ithLoop 
     */
    GeneratedNeutrino(int ithLoop = -1); ///< Set everything to -1
    virtual ~GeneratedNeutrino();
    BalloonInfo balloon;
    double weight;
    int inu;
    int passCutNoWay;
    int passCut2;
    int passCut3;
    int passCutWithinHorizon;    

    ClassDef(GeneratedNeutrino, 1)
  };
  

  /**
   * @class PassingNeutrino
   * @brief Generated neutrino, which passes the trigger
   * 
   * If you pass the detector trigger you will want to know a lot more stuff.
   * What does the signal look like? Etc. etc...
   * Lives in Rootoutput::passTree.
   */

  class PassingNeutrino : public GeneratedNeutrino {
  public:    
    PassingNeutrino(const GeneratedNeutrino& genNu, const AskaryanFreqs& askFreqs, const ShowerProperties& sp);
    virtual ~PassingNeutrino();

    AskaryanFreqs askaryanFreqs;
    ShowerProperties showerProps;

    ClassDef(PassingNeutrino, 1)    
  };
}

#endif // ICEMC_GENERATED_NEUTRINOS_H
