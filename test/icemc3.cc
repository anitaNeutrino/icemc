#include "EventGenerator.h"
#include "CommandLineOpts.h"
#include "IcemcLog.h"

int main(int argc,  char **argv) {
  //--------------------------------------------------------------
  //  MC Anita
  //
  // 12/01/03
  //
  //--------------------------------------------------------------

  //floating point exceptions 
  //
  //
 #ifdef ICEMC_FEEXCEPT
  feenableexcept(FE_INVALID | FE_DIVBYZERO); 
#endif

  icemc::Log log("v3/run1/", 1);
  log << "This is normal text output.\n";
  log << icemc::Log::info << "This should stand out slightly" << std::endl;
  log << icemc::Log::warning << "This should stand out slightly more" << std::endl;
  log << icemc::Log::error << "This should stand out loads!" << std::endl;


  sleep(2);
  return 0;

  
  icemc::Settings settings;
  icemc::CommandLineOpts clOpts(argc, argv, settings);

  if(clOpts.are_good){
    icemc::EventGenerator uhen;
    uhen.generateNeutrinos(settings, clOpts);
  }
  
}
