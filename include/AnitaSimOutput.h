#ifndef ANITA_SIM_OUTPUT
#define ANITA_SIM_OUTPUT

#include "TTree.h"
#include "VoltsRX.h"

class TFile;
class UsefulAnitaEvent;
class Adu5Pat;
class RawAnitaHeader;
class TruthAnitaEvent;

namespace icemc {
  class ANITA;
  class Settings;
  class RayTracer;
  class Screen;

  /**
   * @class AnitaSimOutput
   * @brief Manages the creation of fake ANITA data from the simulation of the ANITA instrument
   */

  class AnitaSimOutput {
  public:
    AnitaSimOutput(const ANITA* detector, const Settings* settings, const char* outputDir, int run);
    virtual ~AnitaSimOutput();

    void fillRootifiedAnitaDataTrees();

  private:
    const ANITA* fDetector; ///< The ANITA detector, parent and owner of this object
    const Settings* fSettings; ///< The simulation settings
    TString fOutputDir; ///< The output directory
    int fRun; ///< The simulated run number (used to uniquely name output files)
    
    UsefulAnitaEvent* fEvent;
    RawAnitaHeader* fHeader;
    Adu5Pat* fGps;
    TruthAnitaEvent* fTruth;

    TFile* fHeadFile;
    TFile* fGpsFile;
    TFile* fEventFile;
    TFile* fTruthFile;
    
    TTree eventTree;
    TTree headTree;
    TTree adu5PatTree;
    
    TTree triggerSettingsTree;
    TTree configTree;
    TTree truthTree;


    void initRootifiedAnitaDataFiles();
    int getIceMCAntfromUsefulEventAnt(int UsefulEventAnt);
    
  };


}


#endif
