#ifndef ICEMC_COMMAND_LINE_OPTIONS_H
#define ICEMC_COMMAND_LINE_OPTIONS_H

#include <string>

namespace icemc {

  /** 
   * @class CommandLineOpts 
   * @brief A simple command line option parser.
   * 
   * Parses a default set of command line arguments for icemc executables
   * Useage: In your executable, do
   * 
   * int main(int argc, char* argv[]){
   *   icemc::CommandLineOpts args(argc, argv);
   *   args.run; // this is the run that was passed 
   * 
   */
  class CommandLineOpts {

  public:
    CommandLineOpts(int argc, char* argv[]);
    std::string executable; ///< derived from argv[0]
    int startNu;
    int nnu_tmp=0;
    double exp_tmp=0;    
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


