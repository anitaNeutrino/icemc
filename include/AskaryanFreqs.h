#ifndef ICEMC_ASKARYAN_FREQS_H
#define ICEMC_ASKARYAN_FREQS_H

#include <vector>
#include <algorithm>


namespace icemc {

  class AskaryanFreqsGenerator; ///< The generator class can access the private elements

  /**
   * @class AskaryanFreqs
   * @brief The frequencies of the signal with all factors accounted for (1/r,  atten. etc.)
   * 
   * Currently this class does not contain any phase information although, 
   * spoiler alert, the frequencies will be peak aligned to make an impulse.
   * Wraps the stupidly name vmmhz array with somethings a bit more humanly readable.
   * And stops people passing raw c-style arrays around.
   */

  class AskaryanFreqs {

    friend AskaryanFreqsGenerator; ///< Allow this generating class to manipulate the private members of the AskaryanFreqs class

  public:

    /** 
     * Default constructor
     */    
    AskaryanFreqs();

    /** 
     * Constructor from c-style array(s)
     * 
     * @param nf is the number of frequencies
     * @param vmmhz_input points to the first element of an array of length nf
     * @param vmmz_em_input optional pointer to an array of length nf, just the electromagnetic shower contribution to the askaryan frequencies
     */
    AskaryanFreqs(int nf, const double* vmmhz_input, const double* vmmz_em_input = NULL);

    /** 
     * Access the i-th element of the frequency magnitudes of the signal. Does a bounds check.
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
    double max_element() const {
      return *std::max_element(vmmhz.begin(), vmmhz.end());
    }    

    /** 
     * @brief Get the smallest value in the frequency array 
     * @return the minimum value
     */
    double min_element() const {
      return *std::min_element(vmmhz.begin(), vmmhz.end());
    }
    
  private:

    std::vector<double> vmmhz; ///< Binned frequencies in V/m/MHz  (Volts per meter per MHz)
    std::vector<double> vmmhz_em; ///< Just from the EM component of the shower,  also binned frequencies in V/m/MHz  (Volts per meter per MHz)    

  };




}


#endif // ICEMC_ASKARYAN_FREQS_H
