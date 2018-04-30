#ifndef ICEMC_RADIO_SIGNAL_H
#define ICEMC_RADIO_SIGNAL_H

#include <vector>
#include <algorithm>


namespace icemc {

  class Signal; ///< Let the linker worry about finding Signal

  /**
   * @class RadioSignal
   * @brief The full signal with all factors accounted for (1/r,  atten. etc.)
   * 
   * Wraps the stupidly name vmmhz array with somethings a bit more humanly readable.
   * And stops people passing raw c-style arrays around.
   */

  class RadioSignal {

    friend Signal; ///< Allow this generating class to manipulate the private members of the RadioSignal class

  public:

    /** 
     * Default constructor
     */    
    RadioSignal();

    /** 
     * Constructor from c-style array
     * 
     * @param nf is the number of frequencies
     * @param vmmhz_input points to the first element of an array of length nf
     */
    RadioSignal(int nf, const double* vmmhz_input);


    /** 
     * Access the i-th element of the frequency represenation of the signal. Does a bounds check.
     * 
     * @param i is the element to access
     * 
     * @return the i_th element of the frequency representation
     */
    double operator[](int i) const;

    /** 
     * @brief Get the largest value in the frequency array 
     * @return the maximum value
     */
    double max() const {
      return *std::max_element(vmmhz.begin(), vmmhz.end());
    }


    /** 
     * @brief Get the smallest value in the frequency array 
     * @return the minimum value
     */
    double min() const {
      return *std::min_element(vmmhz.begin(), vmmhz.end());
    }
    
  private:

    std::vector<double> vmmhz;

    
    // double vmmhz[Anita::NFREQ]; // never again
  };




}


#endif // ICEMC_RADIO_SIGNAL_H
