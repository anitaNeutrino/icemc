#include "EnvironmentVariable.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "TString.h"

/** 
 * Get an environment variable with a nice wrapper function, that tells you what went wrong
 * 
 * @param envName is the environment varible name (e.g. ICEMC_SRC_DIR)
 * @param description is what should get printed if the program can't find it (does nothing if NULL)
 * @param failHard is an optional boolian (default true) that terminates the program if the environment variable isn't found.
 * 
 * @return a c-string containing the value of the environment variable (e.g. ~/ANITA/icemc/)
 */
const char* getEnv(const char* envName, const char* description, bool failHard=true){

  const char* envVar = std::getenv(envName);
  if(!envVar){
    const char* RED = "\x1b[31m";
    const char* RESET = "\x1b[0m";

    // print fatal error or warning
    if(failHard){
      std::cerr << RED << "Fatal error! " << RESET;
    }
    else{
      std::cerr << RED << "Warning! " << RESET;
    }

    // what we couldn't find
    std::cerr << "Could not find environment variable "
              << RED << envName << RESET << std::endl;

    // helpful message, if you gave me one.
    if(description){
      std::cerr << description << std::endl;
    }

    // maybe actually give up
    if(failHard){
      std::cerr << "Giving up." << std::endl;
      exit(1);
    }
  }
  return envVar;
}


/** 
 * Get the source directory of icemc.
 * Will cause the program to terminte if ICEMC_SRC_DIR is not defined.
 * 
 * @return the value of the environment variable (e.g. ~/ANITA/anitaBuildTool/components/icemc)
 */
const char* EnvironmentVariable::ICEMC_SRC_DIR(){

#ifdef ANITA_BUILD_TOOL  
  const char* icemc_src_dir = getEnv("ICEMC_SRC_DIR", "Without this environment variable I can't find input data or config files!", true);
  return icemc_src_dir;
#else
  const char* icemc_src_dir = getEnv("ICEMC_SRC_DIR", "Will guess icemc source directory is present working directory", false);
  if(!icemc_src_dir){
    icemc_src_dir = ".";
  }
  return icemc_src_dir;
  
#endif
}

const char* EnvironmentVariable::ICEMC_VERSION(TString outputdir){

  system(Form("git rev-parse HEAD > %s/gitversion.txt", outputdir.Data()));
  static std::string gitversion;
  std::ifstream gitversionfile (Form("%s/gitversion.txt", outputdir.Data()));
  if (gitversionfile.is_open())
    {
      getline (gitversionfile,gitversion);
      gitversionfile.close();
    }

  return gitversion.c_str();
  
}
