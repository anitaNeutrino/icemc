#include "AskaryanFreqs.h"
#include "IcemcLog.h"
#include "anita.hh"

icemc::AskaryanFreqs::AskaryanFreqs()
  :  vmmhz(Anita::NFREQ, 0)
{
  
  
}


icemc::AskaryanFreqs::AskaryanFreqs(int nf, const double* vmmhz_input)
  :  vmmhz(vmmhz_input, vmmhz_input + nf)
{
  
  
}




double icemc::AskaryanFreqs::operator[](int i) const {
  if(i >= 0 && i < vmmhz.size()){
    return vmmhz[i];
  }
  else {
    Log() << icemc::error << "Attempt to access icemc::AskaryanFreqs at index " << i << ", which is out of bounds. (min=0, max=" << vmmhz.size() <<"). Returning 0.\n";
    return 0;
  }
}
