#ifndef ICEMC_COMMAND_LINE_OPTIONS_H
#define ICEMC_COMMAND_LINE_OPTIONS_H

#include <string>
#include "Settings.h"
#include "IcemcLog.h"

namespace icemc {

  /** 
   * @class CommandLineOpts 
   * @brief A simple command line option parser, updates icemc::Settings
   * 
   * 
   * The job of this class is to parse the command line arguments to an icemc 
   * executable and update the icemc::Settings appropriately.
   * Also populates its public members with the things it found.
   */
  class CommandLineOpts {

  public:
    CommandLineOpts(int argc, char* argv[], icemc::Settings& settings, icemc::Log& log);
    std::string executable; ///< derived from argv[0]
    int startNu;
    int nnu_tmp;
    double exp_tmp;
    std::string outputdir;
    std::string input;
    std::string run_num;
    int run_no;
    double trig_thresh;


    bool are_good; ///< True if the options are coherent, otherwise false


  private:
    void setExecName(char* argv0);
  };



}

#endif


