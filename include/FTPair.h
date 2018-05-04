#ifndef ICEMC_FT_PAIR_H
#define ICEMC_FT_PAIR_H

#include <complex>
#include "TGraph.h"
#include "FFTtools.h"

namespace icemc {


  /**
   * @class FTPair
   * @brief A coupled time domain and frequency domain representation of a waveform.
   *  
   * This class is inevitably a bit complicated because it keeps track of time and frequency 
   * domain representations of waveforms and swaps between them on demand as either representation
   * changes.
   * 
   * It's probably not a good class for writing to disk as it has a lot of internal state.
   */
  class FTPair {
  public:
    
    /** 
     * Time domain c-style array constructor
     * 
     * @param n number of time domain points
     * @param dt time between samples 
     * @param t0 time of first sample
     * @param timeDomainAmplitude array of length n of amplitudes
     */
    FTPair(int n, const double* timeDomainAmplitude,  double dt, double t0 = 0);

    /** 
     * Time domain c-style array constructor
     * 
     * @param n number of time domain points
     * @param dt time between samples 
     * @param t0 time of first sample
     * @param timeDomainAmplitude array of length n of amplitudes
     */
    FTPair(const std::vector<double>& timeDomainAmplitude, double dt, double t0 = 0);

    /** 
     * TGraph constructor
     * 
     * @param grTimeDomain 
     */
    FTPair(const TGraph& grTimeDomain);


    /** 
     * Use this to view, but not change the signal in the time domain.
     * (TGraph::Draw is actually non-const so not sure how useful this is)
     * 
     * @return non-const reference to the time domain waveform.
     */    
    const TGraph& getTimeDomain() const;

    
    /** 
     * Use this to alter the signal in the time domain.
     * e.g. to window a waveform (zero) the first M points by
     * for(int i=0; i < M; i++){
     *   yourFtPair.changeTimeDomain().GetY()[i] = 0;
     * }
     * Then the next time you ask for the frequency domain,
     * 
     *   std::vector<std::complex<double> > yourFtPair.getFreqDomain();
     * 
     * all the FTs are done internally and the frequency bins
     * are updated automagically.
     * 
     * @return non-const reference to the binned frequencies.
     */    
    TGraph& changeTimeDomain();

    
    /** 
     * Use this to view, but not change the signal in the frequency domain.
     * @return const reference to the binned frequencies.
     */    
    const std::vector<std::complex<double> >& getFreqDomain() const;

    /** 
     * Use this to alter the signal in the frequency domain.
     * e.g. to change (zero) the jth frequency bin do:
     *   yourFtPair.changeFreqDomain()[j].real(0);
     *   yourFtPair.changeFreqDomain()[j].imag(0);
     * Then the next time you ask for the time domain it should get updated automagically.
     * 
     * @return non-const reference to the binned frequencies.
     */
    std::vector<std::complex<double> >& changeFreqDomain();


    /** 
     * Force a regeneration of the time domain from the frequency domain.
     * This shouldn't be necessary if you use the accessors correctly.
     */
    void forceUpdateTimeDomain() const;

    /** 
     * Force a regeneration of the frequency domain from the time domain.
     * This shouldn't be necessary if you use the accessors correctly.
     */
    void forceUpdateFreqDomain() const;
    
    
    
    /** 
     * @brief Make a power spectrum graph from the complex frequency domain
     * 
     * @todo make sure units make sense here!
     * 
     * @return new TGraph showing power per unit frequency
     */
    TGraph makePowerSpectrumGraph() const;
    



    
    /** 
     * @brief Use the numerical recipes FT directly.
     * 
     * NOTE: nsize MUST be a power of 2 otherwise this function 
     * will give you nonsense and cause undefined behaviour!
     * This is not checked inside this function.
     * 
     * Also the outputted data is formatted a little differently to 
     * the FFTW convention.
     * 
     * Perhaps you should use the FTPair class to do all the hard work for you?
     * 
     * @param data array to transform 
     * @param isign +1 for real to complex, -1 for inverse
     * @param nsize length of data array, MUST BE A POWER OF 2.
     */
    static void realft(double *data, const int isign, int nsize);








    /** 
     * @brief return true if input is a power of 2
     * 
     * @param n number to test
     * 
     * @return true if a power of 2, false otherwise...
     */    
    inline static bool isPowerOf2(int n){
      /**
       * for n > 0, this is a cute bit trick.
       * In binary, if n is a power of 2...
       *    n     = 100000... (some number of 0s)
       *    n - 1 = 011111.... (number of 1s)
       * therefore (n & n-1) should give 0 for a power of 2
       * (where & is the bitwise AND operator)
       */
      return  (n > 0) && ((n & (n - 1)) == 0);
    }

    /** 
     * @brief Find the next power of 2 greater than or equal to n.
     * 
     * @param n find the next power of 2 greater than or equal to this number
     * 
     * @return the next power of 2 (or n if n is a power of 2)
     */
    inline static int nextPowerOf2(int n){
      int nOut = 1; 
      while(nOut < n){
	nOut = nOut << 1;
      }
      return nOut;
    }

    static inline int getFreqInfo(int n, double dt, double& df) {
      // once and only once...
      df = 1./(n*dt);
      int nf = (n/2)+1;
      return nf;
    }    


    /** 
     * Turn on/off the debug output.
     * 
     * @param db desired state of debug statement switch.
     */
    void setDebug(bool db = true){fDebug = true;}

    
  private:
    
    mutable TGraph fTimeDomainGraph; ///< The time domain representation of the waveform
    mutable std::vector<std::complex<double> > fFreqDomain; ///< The frequency domain representation
    mutable bool fNeedToUpdateTimeDomain; ///< If you call changeFreqDomain, this is set to true, set to false after maybeUpdateFreqDomain
    mutable bool fNeedToUpdateFreqDomain; ///< If you call changeTimeDomain, this is set to true, set to false after maybeUpdateTimeDomain
    bool fDebug; ///< Toggle debug output, probably not useful unless you're developing and there's problem.
    void maybeUpdateFreqDomain() const; ///< If fNeedToUpdateFreqDomain is true, does the appropriate forward FTs
    void maybeUpdateTimeDomain() const; ///< If fNeedToUpdateTimeDomain is true, does the appropriate inverse FTs
    int zeroPadTimeDomainToNearestPowerOf2() const; ///< Because numerical recipes is less good than FFTW


    
  };
}




#endif //ICEMC_FT_PAIR_H
