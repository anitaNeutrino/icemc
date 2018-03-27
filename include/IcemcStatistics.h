#ifndef ICEMC_STATISTICS_H
#define ICEMC_STATISTICS_H


namespace icemc {

  /**
   * @namespace statistics
   * @brief Various statistics related functions should go in here
   */

  namespace statistics {

    /** 
     * Get the Poisson error associated from n events, currently valid up to 20
     * 
     * @param n the number of observed events
     * @param poissonErrorPlus the error in the positive direction
     * @param poissonErrorMinus the error in the negative direction
     */
    void getPoissonError(int n, double& poissonErrorPlus, double& poissonErrorMinus);
  }
}

#endif
