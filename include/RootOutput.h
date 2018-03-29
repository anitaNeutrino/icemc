#ifndef ICEMC_ROOT_OUTPUT
#define ICEMC_ROOT_OUTPUT

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

namespace icemc {

  /**
   * @class RootOutput
   * @brief Looks after all ROOT output from running icemc
   */


  class RootOutput {

    
  public:
    RootOutput();
    RootOutput(const char* outputDir, int run);
    virtual ~RootOutput();
    
    /** 
     * Generate the icefinal file and contained trees
     */
    void make_icefinal();

    class Tree2Output  {
      int inu;
      double horizcoord;
      double vertcoord;
      double scalefactor_distance;
      double scalefactor_attenuation;
      ClassDefNV(Tree2Output, 1);
    } fTree2Output;
    // struct Tree2Output {
    //   int inu;
    //   double horizcoord;
    //   double vertcoord;
    //   double scalefactor_distance;
    //   double scalefactor_attenuation;
    //   ClassDefNV(Tree2Output, 1);
    // };

    TFile* fIceFinalFile; ///< icefinal output file
    TTree* fTree2;
    //Tree2Output* fTree2Output;
    
    
  private:
    TString fOutputDir; ///< The output directory
    int fRun; ///< The simulated run number (used to uniquely name output files)
    


    void zeroFilePointers();
    
  };

  



}


#endif // ICEMC_ROOT_OUTPUT_MANAGER
