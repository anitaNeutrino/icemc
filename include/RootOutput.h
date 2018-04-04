#ifndef ICEMC_ROOT_OUTPUT
#define ICEMC_ROOT_OUTPUT

#include "IcemcLog.h"

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
  class Ray;
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
    RootOutput(Log& iLog, const EventGenerator* uhen = NULL, const Settings* settings = NULL, const char* outputDir = ".", int run = 0);
    virtual ~RootOutput();

    TTree tree2;		///< Filled for each event that is beyond the horizon.
    TTree tree18;
    TTree tree16;
    TTree tree11;		///< tree11
    TTree groundtree;
    TTree icetree;
    TTree ytree;		///<To record y distributions
    TTree banana_tree;		///<To record banana plot info - Stephen
    TTree summarytree;		///< finaltree filled for all events that pass
    TTree mytaus_tree;
    TTree finaltree;		///< finaltree filled for all events that pass
    TTree nupathtree;
    TTree neutrino_positiontree;
    TTree viewangletree;	///< signal as it is produced at the interaction
    TTree jaimetree;		///< signal as it is produced at the interaction
    TTree tree7;		///< tree6 filled just after flavor is set
    TTree tree6b;		///< tree6b filled for the closest antenna to the interaction
    TTree tree6;		///< tree6 filled for neutrinos that enter S of 60 deg S latitude.
    TTree tree5;		///< tree5 filled for each nutau.
    TTree tree3;		///< tree3 if signal is detectable.
    TTree balloontree;		///< filled for all events
    TTree vmmhz_tree;           ///< To record frequency spread at point where it is first filled
    TTree tree1;                ///< tree1 filled for each neutrino
    
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

    void fillRootifiedAnitaDataTrees(const EventGenerator* uhen, const Settings& settings1, const Ray* ray1, const Screen* panel1);

  private:
    TString fOutputDir; ///< The output directory
    int fRun; ///< The simulated run number (used to uniquely name output files)
    TFile* fIceFinal;
    UsefulAnitaEvent* realEvPtr;
    RawAnitaHeader* rawHeaderPtr;
    Adu5Pat* Adu5PatPtr;
    TruthAnitaEvent* truthEvPtr;
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

    void initIceFinal(const EventGenerator* uhen, const Settings* settings);
    void initHist(TH1* h, const char* name, const char* title, int nx, double xMin, double xMax);
    void initHist(TH2* h, const char* name, const char* title, int nx, double xMin, double xMax, int ny, double yMin, double yMax);
    void initTree(TTree* t, const char* name, const char* title, TFile* f);
    void initRootifiedAnitaDataFiles(Log& iLog, const EventGenerator* uhen, const Settings* settings1);
    int getIceMCAntfromUsefulEventAnt(const Settings *settings1,  int UsefulEventAnt);

  };

}


#endif // ICEMC_ROOT_OUTPUT_MANAGER
