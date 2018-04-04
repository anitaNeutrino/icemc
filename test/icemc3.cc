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

  icemc::Log iLog;
  icemc::Settings settings;
  icemc::CommandLineOpts clOpts(argc, argv, settings, iLog);

  if(clOpts.are_good){
    icemc::EventGenerator uhen;
    uhen.generateNeutrinos(settings, clOpts, iLog);
  }
  
  return 0;
}
