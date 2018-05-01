#include "Settings.h"
#include "CommandLineOpts.h"
#include "EventGenerator.h"

int main(int argc,  char **argv) {

  //--------------------------------------------------------------
  //  MC Anita
  //
  // 12/01/03
  //
  //--------------------------------------------------------------

  icemc::Settings settings;
  icemc::CommandLineOpts clOpts(argc, argv, settings);

  if(clOpts.are_good){
    icemc::EventGenerator uhen;
    uhen.generateNeutrinos(settings, clOpts);
  }
  
  return 0;
}
