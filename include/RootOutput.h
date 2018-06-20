#ifndef ICEMC_ROOT_OUTPUT
#define ICEMC_ROOT_OUTPUT

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"


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
    virtual ~RootOutput();

    /**
     * The goal is to have all info in these trees.
     * 
     */
    TTree allTree;  ///< Every neutrino generated
    TTree passTree; ///< Everything that passes the instrument trigger, very detailed


    /**
     * Keep for now while we move to new tree format
     */
    TTree ytree;		///<To record y distributions
    TTree summarytree;		///< finaltree filled for all events that pass
    TTree mytaus_tree;
    TTree finaltree;		///< finaltree filled for all events that pass // move to pass tree?
    TTree nupathtree;
    // TTree viewangletree;	///< signal as it is produced at the interaction
    // TTree balloontree;		///< filled for all events

    TH1D h1mybeta;
    TH1D h1mytheta;             ///< 90-incidentangle when neutrinos enter the Earth.
    TH1F hundogaintoheight_e;
    TH1F hundogaintoheight_h;
    TH1F rec_diff;
    TH1F recsum_diff;
    TH1F rec_diff0;
    TH1F rec_diff1;
    TH1F rec_diff2;
    TH1F rec_diff3;
    TH1F sampleweights;         ///< we sample weights for early events and save them in this histogram, to determine where the cut should be.
    TH2F ref_int_coord;
    TH2F dir_int_coord;
    TH1F prob_eachphi_bn;
    TH1F prob_eachilon_bn;
    TH2F h6;
    TH1F h10;
    TH1F hy;
    TH1F fraction_sec_muons ;
    TH1F fraction_sec_taus;
    TH1F n_sec_muons;
    TH1F n_sec_taus;

    static void initHist(TH1* h, const char* name, const char* title, int nx, double xMin, double xMax, TFile* f);
    static void initHist(TH2* h, const char* name, const char* title, int nx, double xMin, double xMax, int ny, double yMin, double yMax, TFile* f);
    static void initTree(TTree* t, const char* name, const char* title, TFile* f);

    const TString& getOutputDir() const {return fOutputDir;}
    int getRun() const {return fRun;}

  private:
    TString fOutputDir; ///< The output directory
    int fRun; ///< The simulated run number (used to uniquely name output files)
    TFile* fIceFinal;

    void initIceFinal(const EventGenerator* uhen, const Settings* settings);

  };

}


#endif // ICEMC_ROOT_OUTPUT_MANAGER
