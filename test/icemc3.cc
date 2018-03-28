#include "EventGenerator.h"
#include "CommandLineOpts.h"

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

  icemc::Settings settings;
  icemc::CommandLineOpts clOpts(argc, argv, settings);

  if(clOpts.are_good){
    icemc::EventGenerator uhen;
    uhen.generateNeutrinos(settings, clOpts);
  }
  
}
