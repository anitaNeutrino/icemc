#ifndef ICEMC_ENVIRONMENT_VARIABLES
#define ICEMC_ENVIRONMENT_VARIABLES

#include "TString.h"


namespace icemc{

  /** 
   * @namespace EnvironmentVariable
   * @brief Access required environmental variables in one place.
   * 
   * Little namespace to provide sanitized access to any environmental variables 
   * that are necessary to run icemc
   */
  namespace EnvironmentVariable {

    const char* ICEMC_SRC_DIR();
  
    const char* ICEMC_VERSION(TString outputdir);

  }
}
  
#endif
