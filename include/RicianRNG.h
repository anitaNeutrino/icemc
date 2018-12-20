#ifndef ICEMC_RICIAN_RNG_H
#define ICEMC_RICIAN_RNG_H

#include "RNG.h"


namespace icemc {

  class RicianRNG : public RNG {
  public:
    virtual ~RicianRNG();
    void pickRician(double signal, double signal_phase, double sigma, double& amplitude, double &phase);

  private:
    void get_circular_bivariate_normal_random_variable(double& rand_gauss_a, double& rand_gauss_b);
    
  };
  

}


#endif //ICEMC_RICIAN_RNG_H
