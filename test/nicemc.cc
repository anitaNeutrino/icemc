#include "Settings.h"
#include "CommandLineOptions.h"
#include "EventGenerator.h"
#include "ANITA.h"

int main(int argc,  char **argv) {

  //--------------------------------------------------------------
  //  MC Anita
  //
  // 12/01/03
  //
  //--------------------------------------------------------------

  icemc::Settings settings;
  icemc::CommandLineOptions clOpts(argc, argv, settings);

  if(clOpts.are_good){
    icemc::ANITA anita(&settings);
    icemc::EventGenerator uhen(&anita);
    uhen.generateNeutrinos(settings);
  }
  
  return 0;
}
