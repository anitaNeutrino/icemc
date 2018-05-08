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
   * It's probably not a good class to write to disk, as it has a lot of internal state.
   * Maybe just store the time domain graph?
   */
  class FTPair {
  public:
    
    /** 
     * Time domain c-style array constructor
     * 
     * @param n number of time domain points
     * @param timeDomainAmplitude is an array of length n, containing the time domain amplitudes
     * @param dt is the time between samples
     * @param t0 is the time of first sample
     */
    FTPair(int n, const double* timeDomainAmplitude,  double dt, double t0 = 0);

    /** 
     * Time domain vector constructor
     * 
     * @param timeDomainAmplitude vector containing the time domain amplitudes
     * @param dt is the time between samples
     * @param t0 is the time of the first sample
     */
    FTPair(const std::vector<double>& timeDomainAmplitude, double dt, double t0 = 0);

    /** 
     * Time domain TGraph constructor
     * 
     * @param grTimeDomain 
     */
    FTPair(const TGraph& grTimeDomain);




    /** 
     * Freq domain array constructor
     * 
     * @param nf is the length of the freqDomainPhasors array
     * @param freqDomainPhasors is the array of complex numbers characterizing the signal in the frequency domain (length nf)
     * @param deltaF is the time between frequency bins
     * @param t0 time of first sample (in the time domain)
     */
    FTPair(int nf,  const std::complex<double>* freqDomainPhasors, double deltaF, double t0 = 0);


    /** 
     * Freq domain std::vector constructor
     * 
     * 
     * @param freqDomainPhasors is the vector of complex numbers characterizing the signal in the frequency domain
     * @param deltaF is the time between frequency bins
     * @param t0 time of first sample (in the time domain)
     */    
    FTPair(const std::vector<std::complex<double> >& freqDomainPhasors, double deltaF, double t0 = 0);
    
    





    
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
     *   const std::vector<std::complex<double> >&  constFreqs = yourFtPair.getFreqDomain();
     * 
     * or
     * 
     *   std::vector<std::complex<double> >&  nonConstFreqs = yourFtPair.changeFreqDomain();
     * 
     * The frequency  array is updated automagically and all the FTs are done internally.
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
     * @todo allow Ryan style normalization options
     * 
     * @return new TGraph showing power per unit frequency
     */
    TGraph makePowerSpectrumGraph() const;
    

    /** 
     * Get the number of times from the current state of the FTPair
     * @return the length of the time domain (after doing internal updates)
     */
    inline int getNumTimes() const {return getTimeDomain().GetN();}

    /** 
     * Get the deltaT from the current state of the FTPair
     * @return the time domain sample separation (after doing internal updates)
     */
    inline double getDeltaT() const {return getTimeDomain().GetX()[1] - getTimeDomain().GetX()[0];}


    /** 
     * Get the number of frequency bins from the current state of the FTPair
     * @return the length of the frequency domain (after doing internal updates)
     */
    inline int getNumFreqs() const {return getFreqDomain().size();}

    /** 
     * Get the deltaF from the current state of the FTPair
     * @return the frequency domain bin separation (after doing internal updates)
     */
    inline double getDeltaF() const {return getDeltaF(getNumTimes(), getDeltaT());}


    
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
     * @todo hide this function inside the cc file when the refactor is complete.
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
       * (where & is the bitwise AND operator).
       * 
       * Put another way, in general subtracting 1 will 
       * set the lowest 1 to 0 and all bits below that will 
       * be set to 1. In a power of two, that lowest 1 bit is
       * the only 1 bit.
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


    /** 
     * How many frequency bins will there be given this many time bins?
     * 
     * @param nt is the number of time samples
     * 
     * @return the number of frequency bins
     */
    inline static int getNumFreqs(int nt){return 1 + nt/2;}

    /** 
     * How many time points will there be given this many frequency bins?
     * 
     * @param nf is the number of frequency bins
     * 
     * @return the number of time points
     */
    inline static int getNumTimes(int nf){return 2*(nf - 1);}


    /** 
     * What is the frequency bin spacing for a generic DFT?
     * 
     * @param nt is the number of time points
     * @param dt is the separation between time points
     * 
     * @return the frequency bin spacing
     */
    inline static double getDeltaF(int nt, double dt){return 1./(nt*dt);}



    /** 
     * What is the time point spacing for a generic inverse DFT?
     * 
     * @param nf is the number of frequency bins
     * @param df is the separation between frequency bins
     * 
     * @return the time bins spacing
     */
    inline static double getDeltaT(int nf, double df){
      int nt = getNumTimes(nf);      
      return 1./(nt*df);
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
    bool fDebug; ///< To toggle debug output, probably not useful unless you're developing and there's problem.
    void maybeUpdateFreqDomain() const; ///< If fNeedToUpdateFreqDomain is true, does the appropriate forward FTs
    void maybeUpdateTimeDomain() const; ///< If fNeedToUpdateTimeDomain is true, does the appropriate inverse FTs
    int zeroPadTimeDomainLengthToPowerOf2() const; ///< Because numerical recipes is less good than FFTW
    int zeroPadFreqDomainSoTimeDomainLengthIsPowerOf2(double df) const; ///< Because numerical recipes is less good than FFTW
  };
}




#endif //ICEMC_FT_PAIR_H
