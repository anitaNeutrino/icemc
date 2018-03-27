#include "EventGenerator.h"


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

  icemc::EventGenerator uhen;
  uhen.generateNeutrinos(argc, argv);
  
}
