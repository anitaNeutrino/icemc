#ifndef CHANTRIGGER_HH_
#define CHANTRIGGER_HH_

////////////////////////////////////////////////////////////////////////////////////////////////
//class ChanTrigger:
////////////////////////////////////////////////////////////////////////////////////////////////

#include "TRandom3.h"
#include "vector.hh"
#include "anita.hh"

// using std::vector;

namespace icemc {

  class GlobalTrigger;
  class Balloon;
  class Screen;
  class IceModel;
  class Settings;
  class RX;
  

  //! Class that handles the channel trigger
  class ChanTrigger {
    
  private:

    /**
     * ADC counts to relative power threshold
     * 
     * Function to convert adc thresholds into
     * relative power thresholds that can be handled by icemc
     * So far this only works for ANITA-3 and full band trigger
     *
     * @param  anita1 :: Anita - anita object
     * @param  ipol :: int - which polarisation 
     * @param  iant :: int - which antenna number
     * @return relative power threshold (double)
     */
    double ADCCountstoPowerThreshold(Anita *anita1, int ipol, int iant);

    const static int NSURF=9;                                     ///< Number of surfs
    const static int NSURFPLUSONE=10;                             ///< Number of surfs plus one
    const static int NSURFMINUSONE=8;                             ///< Number of surfs minus on1
    const static int NCHANNELS=32;                                ///< Number of channells on each surf
    const static int NPOINTS=4073;                                ///< Max number of points from surf measurements
    static const unsigned NFOUR = 1024;                           ///< Number of points in Fourier space
    static const unsigned HALFNFOUR = 512;                        ///< Half of the number of points in the Fourier space

    TRandom3 Rand3;                                               ///< Random number generator instance
    double thisrate;                                              ///< Rate in MHz
    double thispowerthresh;                                       ///< Relative power threshold
    double e_component;                                           ///< E comp along polarization
    double h_component;                                           ///< H comp along polarization
    double n_component;                                           ///< normal comp along polarization
    double e_component_kvector;                                   ///< component of e-field along the rx e-plane
    double h_component_kvector;                                   ///< component of the e-field along the rx h-plane
    double n_component_kvector;                                   ///< component of the e-field along the normal  
    double hitangle_e;                                            ///< angle the ray hits the antenna wrt e-plane
    double hitangle_h;                                            ///< angle the ray hits the antenna wrt h-plane
      
  public:

    //!  Channel trigger constructur
    ChanTrigger(); 

  
    //! Initialize trigger bands
    /**
     * Initialize trigger bands
     * 
     * @param  anita1 :: Anita - anita object
     */
    void InitializeEachBand(Anita *anita1);

    //! Apply the antenna gain
    /**
     *
     * Loop through all points in the screen, and apply the antenna gain to each point
     * All the waveforms are then summed with screen point weights and delay
     * volts_rx_forfft = time domain waveforms after antenna gain
     * vhz_rx          = amplitude in Fourier domain after antenna gain
     * 
     * @param  settings1 :: Settings - the simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  bn1 :: Balloon - balloon object
     * @param  panel1 :: Screen - screen object
     * @param  ant :: int- antenna number
     * @param  n_eplane :: Vector 
     * @param  n_hplane :: Vector 
     * @param  n_normal :: Vector
     */  
    void ApplyAntennaGain(Settings *settings1, Anita *anita1, Balloon *bn1, Screen *panel1, int ant, Vector &n_eplane, Vector &n_hplane, Vector &n_normal);

    //! Apply trigger path
    /**
     *
     * Apply trigger path to each channel
     * v_banding_rfcm_forfft = time domain waveforms after trigger path
     * v_banding_rfcm         = amplitude in Fourier domain after trigger path
     *
     * @param  settings1 :: Settings - the simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  ant :: int - antenna number
     */ 
    void TriggerPath(Settings *settings1, Anita *anita1, int ant, Balloon *bn1);

    //! Apply digitizer path
    /**
     *
     * Apply digitizer path to each channel
     * volts_rx_rfcm_lab  = time domain waveforms after digitizer path
     * vhz_rx_rfcm_lab_e  = amplitude in Fourier domain after digitizer path
     *
     * @param  settings1 :: Settings - the simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  ant :: int - antenna number
     */ 
    void DigitizerPath(Settings *settings1, Anita *anita1, int ant, Balloon *bn1);

    //! Time shift and fluctuate signal
    /**
     *
     * Apply time delays to each channel and add thermal noise
     *
     * @param  settings1 :: Settings - the simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  ilayer :: int - payload layer
     * @param  ifold :: int - payload phi sector
     * @param  volts_rx_rfcm_lab_e_all :: double [48][512] - time domain waveform for each channel (VPOL)
     * @param  volts_rx_rfcm_lab_h_all :: double [48][512] - time domain waveform for each channel (HPOL)
     */ 
    void TimeShiftAndSignalFluct(Settings *settings1, Anita *anita1, int ilayer, int ifold, double volts_rx_rfcm_lab_e_all[48][512], double volts_rx_rfcm_lab_h_all[48][512]);
  
    //!  Convert E and H to left and right e field
    /**
     * Convert E and H to left and right e field
     * 
     * @param  e_component :: double - 
     * @param  h_component :: double - 
     * @param  lcp_component :: double - referenced lcp component
     * @param  rcp_component :: double - referenced rcp component
     */ 
    static void ConvertEHtoLREfield(double,double,double&,double&);

    //!Convert E and H to left and right energy
    /**
     * Convert E and H to left and right energy
     * 
     * @param  e_component :: double - 
     * @param  h_component :: double - 
     * @param  lcp_component :: double - referenced lcp component
     * @param  rcp_component :: double - referenced rcp component
     */ 
    static void ConvertEHtoLREnergy(double,double,double&,double&);

    //! Convert L1 trigger of the Anita trigger scheme
    /**
     * Convert H and V components to left and right circular polarization
     * in time domain
     * 
     * @param nfour :: const int - number of fourier points
     * @param vvolts :: double * - array of voltage values (VPOL)
     * @param hvolts :: double * - array of voltage values (HPOL)
     * @param left :: double * - array of voltage values (LCP)
     * @param right :: double * - array of voltage values (RCP)
     */ 
    void ConvertHVtoLRTimedomain(const int nfour,double *vvolts,
				 double *hvolts,
				 double *left,double *right);


    //!	The L1 trigger of the Anita trigger scheme
    /**
     * 
     *  @param anita1 :: Anita - anita payload object
     *  @param timedomain_output_1 :: double [4][NFOUR] - time domain output for each band e
     *  @param timedomain_output_2 :: double [4][NFOUR] - time domain output for each band h 
     *  @param powerthreshold :: double[2][5] - relative power thresholds for each pol and each band
     *  @param channels_passing_e_forglob :: double* - array of channels passing the L1 trigger for e (used in the GlobalTrigger class)
     *  @param channels_passing_h_forglob :: double* - array of channels passing the L1 trigger for h (used in the GlobalTrigger class)
     *  @param npass :: &int - number of bands passing the L1 trigger
     */
    void L1Trigger(Anita *anita1,double timedomain_output_1[5][Anita::NFOUR],double timedomain_output_2[5][Anita::NFOUR],double powerthreshold[2][5],int *channels_passing_e_forglob,int *channels_passing_h_forglob,int &npass);    


    //!	Returns the thisrate variable value (in MHz)
    /**
     *	\todo	Refactor the class so that this function is either deprecated or "thisrate"
     *			is implemented differently.
     */
    double getRate();

    //!	Calculates the trigger threshold for an antenna via a fit from the "singles" rate and band identifier
    /**
     *	This code is stolen from Stephen's AnitaHardwareTrigger in his Aesop software
     *	@param rate is the desired singles rate, in MHz.
     *	@param band ranges from 0 to 3
     *    @return  trigger threshold, in units of Diode Output / Average Diode Output
     *	The output is good for all antennas, and the curves are good fits in the
     *	region ~1 MHz to 25 MHz.  The curves look like they should continue beyond
     *	that range as well.
     *
     *	\todo	The hard-coded values should be explained and the band dependence should be more explicit,
     *			because the actual bands can be modified elsewhere in the code and this function won't be
     *			changed unless the user specifically knows to.
     */
    double rateToThreshold(double rate, int band);

    //!	Get Noise
    /**
     *  
     * @param  settings1 :: Settings - the simulation settings
     * @param  altitude_bn :: double - altitude of the balloon
     * @param  geoid :: double -
     * @param  bw :: double - 
     * @param  theta :: double -
     * @param  temp :: double -
     *
     */
    static double GetNoise(Settings *settings1,double altitude_bn,double geoid,double theta,double bw,double temp);

    //! Which bands passes the trigger
    /**
     *
     * @param  settings1 :: Settings - simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  globaltrig1 :: GlobalTrigger - global trigger object
     * @param  bn1 :: Balloon - balloon object
     * @param  ilayer :: int - layer number
     * @param  ifold :: int - phi sector
     * @param  dangle :: int - 
     * @param  emfrac :: double - electromagnetic fraction of the shower
     * @param  hadfrac :: double - hadronic fraction of the shower
     * @param  thresholds :: double [2][5] - relative power thresholds for each pol and band
     */
    void WhichBandsPass(Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double dangle, double emfrac, double hadfrac, double thresholds[2][5]);

    //! Which bands passes the trigger (for trigger scheme 0 and 1)
    /**
     *
     * @param  settings1 :: Settings - simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  globaltrig1 :: GlobalTrigger - global trigger object
     * @param  bn1 :: Balloon - balloon object
     * @param  ilayer :: int - layer number
     * @param  ifold :: int - phi sector
     * @param  thresholds :: double [2][5] - relative power thresholds for each pol and band
     */
    void WhichBandsPassTrigger1(Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double thresholds[2][5]);

    //! Which bands passes the trigger (for trigger scheme larger than 2)
    /**
     *
     * @param  settings1 :: Settings - simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  globaltrig1 :: GlobalTrigger - global trigger object
     * @param  bn1 :: Balloon - balloon object
     * @param  ilayer :: int - layer number
     * @param  ifold :: int - phi sector
     * @param  dangle :: int - 
     * @param  emfrac :: double - electromagnetic fraction of the shower
     * @param  hadfrac :: double - hadronic fraction of the shower
     * @param  thresholds :: double [2][5] - relative power thresholds for each pol and band
     */
    void WhichBandsPassTrigger2(Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double dangle, double emfrac, double hadfrac, double thresholds[2][5]);

    //! Find peak voltage of a waveform
    /**
     * @param waveform :: double* - waveform
     * @param n :: int - number of points in waveform
     */
    static double FindPeak(double *waveform,int n);

  
    //!	Sets the threshold values based on which payload and where the antenna is located physically
    /**
     *	The nadir antennas had a separate threshold from the other antennas due to the way that they
     *	were "OR"'d between their two neighbors.
     *
     *	\todo	Deprecate in favor of the more robust boost::multi_array or the more specialized
     *			PayloadArray class. Both have multi-index access to the same items.
     */
    void GetThresholds(Settings *settings1,Anita *anita1,int ilayer,double thresholds[2][5]); // get thresholds for this layer

    //! Apply the diode convolution
    /**
     * @param  settings1 :: Settings - simulation settings
     * @param  anita1 :: Anita - anita object
     * @param  globaltrig1 :: GlobalTrigger - global trigger object
     * @param  ilayer :: int - layer number
     * @param  ifold :: int - phi sector
     * @param  mindiodeconvl :: double[5]
     * @param  onediodeconvl :: double[5]
     * @param  psignal :: double[5][NFOUR]
     * @param  timedomain_output :: double[5][NFOUR]
     * @param  ibinshift :: int
     * @param  ipol :: int - which polarization
     * @param  thresholds :: double[2][5] - relative power thresholds for each pol and band
     */
    void DiodeConvolution(Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, int ilayer, int ifold, double mindiodeconvl[5], double onediodeconvl[5], double psignal[5][Anita::NFOUR],  double timedomain_output[5][Anita::NFOUR], int ibinshift, int ipol, double thresholds[2][5]);


    //! Increment the volts in each band 
    /**
     * This increments the bwslice_volts_...'s at each frequency
     * @param  settings1 :: Settings - simulation settings
     * @param  anita1 :: Anita - anita payload object
     * @param  ibw :: int - which band
     * @param  k :: int - frequency bin
     */
    void addToChannelSums(Settings *settings1,Anita *anita1,int ibw,int k);

    //!	Returns whether the indicated antenna and band are "masked"
    /**
     *  Only works for Anita-1 and 2	
     * 
     * @param  surfTrigBandMask :: unsigned short [9][2] - surf masks 
     * @param  ibw :: int - which band
     * @param  ilayer :: int - layer number
     * @param  ifold :: int - phi sector
     * @param  ipol :: int - pol number
     * @return whether the channel and band are masked or not
     */
    static int IsItUnmasked(unsigned short surfTrigBandMask[9][2],int ibw,int ilayer, int ifold, int ipol);

    //! Apply impulse response to digitizer path
    /**
     * This can only work when FFTtools is also linked
     *
     * @param  settings1 :: Settings - simulation settings
     * @param  anita1 :: Anita - anita payload object
     * @param  nPoints :: int - number of points in time domain
     * @param  ant :: int - antenna number
     * @param  x :: double* - time values
     * @param  y :: double[512] - output voltages
     * @param  pol :: bool - which polarization
     */
    void applyImpulseResponseDigitizer(Settings *settings1, Anita *anita1, int nPoints, int ant, double *x, double y[512], bool pol);

    //! Apply impulse response to trigger path
    /**
     * This can only work when FFTtools is also linked
     *
     * @param  settings1 :: Settings - simulation settings
     * @param  anita1 :: Anita - anita payload object
     * @param  ant :: int - antenna number
     * @param  y :: double[512] - output voltages
     * @param  vhz :: double* - amplitude in Fourier domain
     * @param  pol :: bool - which polarization
     */
    void applyImpulseResponseTrigger(Settings *settings1, Anita *anita1, int ant, double y[512], double *vhz, bool pol);

    //! Add noise from ANITA-3 flight to the time domain waveforms
    /**
     * Rayleigh distribution parameters of ANITA-3 thermal noise are read in anita.cc
     * This function generates random noise using amplitude of Rayleigh distributions.
     * This can only work when FFTtools is also linked
     *
     * @param  anita1 :: Anita - anita payload object
     * @param  pol :: int - which polarization
     * @param  ant :: int - which antennta
     * @param  also_digi :: bool - also fill in digitizer 
     */
    void getNoiseFromFlight(Anita* anita1, int ant, bool also_digi = true);

    //! Inject pulse after the antenna (used for trigger efficiency scans)
    /**
     * Pulser waveform is read in anita.cc
     * THIS FUNCTION IS ONLY USED WITH THE OLD TRIGGER IMPULSE RESPONSE
     * 
     * @param  anita1 :: Anita - anita payload object
     * @param  ant :: int - which antennta
     */
    void injectImpulseAfterAntenna(Anita *anita1, int ant);

    //! Get time domain graph of pulse at AMPA (used for trigger efficiency scans)
    /**
     * Pulser waveform is read in anita.cc
     * 
     * @param  anita1 :: Anita - anita payload object
     * @param  ant :: int - which antennta
     */
    TGraph *getPulserAtAMPA(Anita *anita1, int ant);

    //! Save signal and noise waveforms at trigger
    /**
     * @param anita1 :: Anita - anita payload object
     * @param sig    :: double[2][NFOUR/2] - output 
     * @param noise  :: double[2][NFOUR/2] - output
     */
    void saveTriggerWaveforms(Anita *anita1, double sig0[48], double sig1[48], double noise0[48], double noise1[48]);

    //! Save signal and noise waveforms at digitizer
    /**
     * @param anita1 :: Anita - anita payload object
     * @param sig    :: double[2][NFOUR/2] - output 
     * @param noise  :: double[2][NFOUR/2] - output
     */
    void saveDigitizerWaveforms(Anita *anita1, double sig0[48], double sig1[48], double noise0[48], double noise1[48]);

    //! Inject pulse at the surf (used for trigger efficiency scans)
    /**
     * Pulser waveforms are read in anita.cc
     * volts_triggerPath_e/h is substituded with pulser waveform + noise
     * This can only work when FFTtools is also linked
     * 
     * @param  anita1 :: Anita - anita payload object
     * @param  ant :: int - which antennta
     * @param  volts_triggerPath_e :: double[NFOUR] - time domain waveform at surf
     * @param  volts_triggerPath_h :: double[NFOUR] - time domain waveform at surf
     */
    void injectImpulseAtSurf(Anita *anita1, double volts_triggerPath_e[Anita::HALFNFOUR], double volts_triggerPath_h[Anita::HALFNFOUR], int ant);

    //! Add CW
    /**
     * Add CW in time domain
     *
     * @param  anita1 :: Anita - anita payload object
     * @param frequency at which to simulate cw
     * @param phase
     * @param amplitude
     */  
    void calculateCW(Anita *anita1, double frequency, double phase, double amplitude);
    //! Apply Butterworth Filter
    /**
     * This is an approximation of the notch filters flown during ANITA4
     * 
     * @param ff :: frequenzy in Hz
     * @param ampl :: amplitude
     * @param filters :: array of three integers indicating notch status (on/off)
     * @return frequency domain amplitude scaled by the filter
     */
    double applyButterworthFilter(double ff, double ampl, int notchStatus[3]);

  
    double vhz_rx[2][5][Anita::NFREQ];                          ///< Array of amplitudes in the Fourier domain (V/Hz) after the antenna gain. Indeces stand for [ipol][iband][ifreq] 
    double volts_rx_forfft[2][5][Anita::HALFNFOUR];             ///< Array of time domain after the antenna gain. Indeces stand for [ipol][iband][itime] 
    std::vector<int> flag_e[5];                                 ///< Which bands pass trigger e
    std::vector<int> flag_h[5];                                 ///< Which bands pass trigger h
    double bwslice_volts_pol0[5];                               ///< Sum voltage for each slice in bandwidth for the lcp polarization
    double bwslice_volts_pol1[5];                               ///< Sum voltage for each slice in bandwidth for the rcp polarization
    double bwslice_energy_pol0[5];                              ///< Square the sum of voltage for each slice in bandwidth for the 0th polarization
    double bwslice_energy_pol1[5];                              ///< Square the sum of voltage for each slice in bandwidth for the 1st polarization
    double bwslice_volts_pol0_em[5];                            ///< Component of the voltage that comes from the em shower for 0th polarization
    double bwslice_volts_pol1_em[5];                            ///< Component of the voltage that comes from the em shower for 1st polarization
    double bwslice_volts_pole[5];                               ///< Sum voltage for each slice in bandwidth for the e polarization
    double bwslice_energy_pole[5];                              ///< Square the sum of voltage for each slice in bandwidth for e polarization.  The 5th element is the full band
    double bwslice_volts_polh[5];                               ///< Sum voltage for each slice in bandwidth for the h polarization
    double bwslice_energy_polh[5];                              ///< Square the sum of voltage for each slice in bandwidth for h polarization.  The 5th element is the full band
    double volts_rx_rfcm_lab[2][Anita::HALFNFOUR];              ///< For digitizer path, time domain voltage vs. time after rx, rfcm's and lab
    double volts_rx_rfcm_lab_all[2][48][Anita::HALFNFOUR];      ///< For digitizer path, time domain voltage vs. time after rx, rfcm's and lab
    double volts_rx_rfcm[2][Anita::HALFNFOUR];                  ///< For digitizer path, time domain voltage vs. time after rx, rfcm's
    double justNoise_digPath[2][Anita::HALFNFOUR];              ///< For digitizer path, time domain noise from flight
    double justNoise_trigPath[2][Anita::HALFNFOUR];             ///< For trigger path, time domain noise from flight
    double cw_digPath[2][Anita::HALFNFOUR];                     ///< For digitizer path, time domain cw
    double justSig_trigPath[2][Anita::HALFNFOUR];               ///< Just signal in trigger path
    double justSig_digPath[2][Anita::HALFNFOUR];                ///< Just signal in trigger path
  
    // these are filled for triggerscheme==0 and triggerscheme==1
    // frequency domain voltage and energy based
    double signal_eachband[2][Anita::NBANDS_MAX];               ///< Signal in each band			     
    double threshold_eachband[2][Anita::NBANDS_MAX];            ///< Threshold in each band		     
    double noise_eachband[2][Anita::NBANDS_MAX];                ///< Noise in each band			     
    int passes_eachband[2][Anita::NBANDS_MAX];                  ///< Whether the signal passes or not each band
    std::vector<double> vsignal_eachband[2];                    ///< Signal in each band			     
    std::vector<double> vthreshold_eachband[2];                 ///< Threshold in each band		     
    std::vector<double> vnoise_eachband[2];                     ///< Noise in each band			     
    std::vector<int>    vpasses_eachband[2];                    ///< Whether the signal passes or not each band
    
    double v_banding_rfcm[2][5][Anita::NFREQ];                  ///< This is Volts/m as a function of frequency after rfcm's and banding
    double v_banding_rfcm_forfft[2][5][HALFNFOUR];              ///< Starts out as V/s vs. freq after banding, rfcm, after fft it is V vs. t
    double vm_banding_rfcm_forfft[2][5][HALFNFOUR];             ///< Starts out as V/s vs. freq after banding, rfcm, after fft it is V vs. t
    double v_banding_rfcm_forfft_temp[2][5][HALFNFOUR];         ///< Use for the averaging over 10 neighboring bins
    double integral_vmmhz;                                      ///< Electric field integral    
    int unwarned;                                               ///< Whether we have warned the user about resetting thresholds when they are beyond the measured bounds
  }; //class ChanTrigger

}
  
#endif

