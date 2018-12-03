#ifndef ICEMC_ROOT_OUTPUT
#define ICEMC_ROOT_OUTPUT

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Event.h"


// need some ifdefs?
class TruthAnitaEvent;

class UsefulAnitaEvent;
class RawAnitaHeader;
class Adu5Pat;

namespace icemc {

  class EventGenerator;
  class Settings;
  class Anita;
  class RayTracer;
  class Screen;

  /**
   * @class RootOutput
   * @brief Looks after all ROOT output from running icemc.
   *
   * If stuff needs to be added, make sure it's a member of the 
   * event generator class.
   */
  class RootOutput {


  public:
    RootOutput(const EventGenerator* uhen = NULL, const Settings* settings = NULL, const char* outputDir = ".", int run = 0);
    RootOutput(const char* fileName,  int run);    
    virtual ~RootOutput();

    /**
     * The goal is to have all info in these trees.
     * 
     */
    TTree& allTree(){return  *fAllTree;}  ///< For every neutrino generated,  not too detailed.
    TTree& passTree(){return *fPassTree;} ///< For every neutrino that passes the detector trigger, very detailed.

    static void initHist(TH1* h, const char* name, const char* title, int nx, double xMin, double xMax, TFile* f);
    static void initHist(TH2* h, const char* name, const char* title, int nx, double xMin, double xMax, int ny, double yMin, double yMax, TFile* f);
    static void initTree(TTree* t, const char* name, const char* title, TFile* f);

    const TString& getOutputDir() const {return fOutputDir;}
    int getRun() const {return fRun;}
    
    EventSummary& eventSummary() const {return *fEventSummary;}
    Event& event() const {return *fEvent;}

  private:
    void initIceFinal(const EventGenerator* uhen, const Settings* settings);
    void readIceFinal();


    TString fOutputDir; ///< The output directory
    int fRun; ///< The simulated run number (used to uniquely name output files)
    TFile* fIceFinal = nullptr;
    TTree* fAllTree = nullptr;
    TTree* fPassTree = nullptr;
    EventSummary* fEventSummary = nullptr;
    Event* fEvent = nullptr;

    

  };

}


#endif // ICEMC_ROOT_OUTPUT_MANAGER
