#ifndef ICEMC_ASKARYAN_FREQS_H
#define ICEMC_ASKARYAN_FREQS_H

#include <vector>
#include <algorithm>
#include "TObject.h"
#include "TGraph.h"

namespace icemc {

  class AskaryanFreqsGenerator;
  class ShowerProperties;

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
     * @param minFreqHz lower bound frequency (Hz)
     * @param deltaFreqHz change between bins
     * @param vmmhz_input points to the first element of an array of length nf
     */
    AskaryanFreqs(int nf, double minFreqHz, double deltaFreqHz, double cherenkovAngle, const ShowerProperties* sp, const double* vmmhz_input);



    
    /** 
     * @brief When you view the Askaryan frequencies off cone, how much of each frequency you see will depend on the view angle.
     * This function takes in a viewAngle relative to the shower axis (@todo check it is relative to that!) and reduces
     * the amplitude in each frequency bin accordingly.
     * 
     * @param viewAngleRadians view angle relative to shower axis
     * @param deltaThetaEmTest 
     * @param deltaThetaHadTest 
     * @param testFreqHz 
     */
    void taperAmplitudesForOffConeViewing(double viewAngleRadians, double deltaThetaEmTest=0, double deltaThetaHadTest=0, double testFreqHz=0);

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
    inline double maxElement() const {
      return *std::max_element(vmmhz.begin(), vmmhz.end());
    }    

    /** 
     * @brief Get the smallest value in the frequency array 
     * @return the minimum value
     */
    inline double minElement() const {
      return *std::min_element(vmmhz.begin(), vmmhz.end());
    }

    /** 
     * @brief Creates a TGraph of the stored Askaryan E-field (Volts/M/MHz) vs. frequency
     * @return the created TGraph
     */
    TGraph makeGraph() const;


    /** 
     * Get the total power in the frequency array
     * 
     * @return 
     */
    double power() const {
      if(fTotalPowerDirty){
	fTotalPower = 0;
	for(auto amp : vmmhz){
	  //@todo confirm a power normalization scheme
	  fTotalPower += amp*amp*fDeltaFreqHz;
	}
	fTotalPowerDirty = false;
      }
      return fTotalPower;
    }


  private:

    std::vector<double> vmmhz; ///< Binned frequencies in V/m/MHz  (Volts per meter per MHz)
    double fMinFreqHz = 0; ///< Frequency of vmmhz[0] (Hz)
    double fDeltaFreqHz = 0; ///< Space between frequency bins (Hz)

    /// for tapering... (a.k.a viewing the frequencies of the shower off-cone)
    
    double fCherenkovAngleRad = 0; ///< Cherenkov angle of frequencies
    
    double fEmFrac = 0; // fraction of the shower from EM component 
    double fHadFrac = 0; // fraction of the shower from hadronic component

    // the tapering goes as 1/frequency, therefore only need to find a single frequency
    // to know the tapering at all frequencies...
    double fSpreadTestFreqHz = 0; // frequency at which to do full tapering calculation
    double fDeltaThetaEmTest = 0; // angular spread of test frequency from EM component
    double fDeltaThetaHadTest = 0;// angular spread of test frequency from Had component


    
    mutable double fTotalPower = 0;
    mutable bool fTotalPowerDirty = true;

    ClassDef(AskaryanFreqs, 1)
  };




  

}


#endif // ICEMC_ASKARYAN_FREQS_H
