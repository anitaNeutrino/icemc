#ifndef TRIGGER_HH_
#define TRIGGER_HH_

////////////////////////////////////////////////////////////////////////////////////////////////
//class Trigger:
////////////////////////////////////////////////////////////////////////////////////////////////

class RX;
class Anita;
using std::vector;
class Balloon;
class IceModel;
class Settings;
class TRandom3;
class Screen;
class GlobalTrigger;

//! Handles L0 and L1 Triggers for an antenna
class AntTrigger {
    
private:
    
    TRandom3 Rand3;
    const static int NSURF=9; // number of surfs
    const static int NSURFPLUSONE=10; // for some reason the threshold array is 1 longer
    const static int NSURFMINUSONE=8; // for some reason the threshold array is 1 longer
    const static int NCHANNELS=32; // number of channels on each surf
    const static int NPOINTS=4073; // max number of points from surf measurements
    
    static const unsigned NFOUR = 1024;
    static const unsigned HALFNFOUR = 512;
    
    // // these are from Ryan's threshold scans
    // int threshold[NSURFPLUSONE][NCHANNELS][NPOINTS];
    // int rate[NSURFMINUSONE][NCHANNELS][NPOINTS];    
    // int minadcthresh[NSURFMINUSONE][NCHANNELS];
    // int maxadcthresh[NSURFMINUSONE][NCHANNELS];
    
    
    //int Getiangle(double viewangle);
    
  //    double ADCCountstoPowerThreshold(int threshadc,int isurf,int ichan);
    double ADCCountstoPowerThreshold(Anita *anita1, int ipol, int iant);
    double thisrate;// set when getFunction is called for each channel. this is in MHz
    double thispowerthresh;// set when getFunction is called for each channel

    double e_component; // E comp along polarization
    double h_component; // H comp along polarization
    double n_component; // normal comp along polarization
    double e_component_kvector; // component of e-field along the rx e-plane
    double h_component_kvector; // component of the e-field along the rx h-plane
    double n_component_kvector; // component of the e-field along the normal 
    double hitangle_e, hitangle_h;       // angle the ray hits the antenna wrt e-plane, h-plane
    
    double vhz_rx_e[Anita::NFREQ]; // V/Hz after antenna gains
    double vhz_rx_h[Anita::NFREQ];
    // same but with binning for fft
    double volts_rx_e_forfft[Anita::HALFNFOUR];
    double volts_rx_h_forfft[Anita::HALFNFOUR];

public:
    
    AntTrigger(); // constructor
    void InitializeEachBand(Anita *anita1);
    void ConvertInputWFtoAntennaWF(Settings *settings1, Anita *anita1, Balloon *bn1, Screen *panel1, Vector n_eplane, Vector n_hplane, Vector n_normal, int ilayer, int ifold);
    void DigitizerPath(Settings *settings1, Anita *anita1, int ilayer, int ifold);
    void TimeShiftAndSignalFluct(Settings *settings1, Anita *anita1, int ilayer, int ifold, double volts_rx_rfcm_lab_e_all[48][512], double volts_rx_rfcm_lab_h_all[48][512]);
  void PrepareTriggerPath(Settings *settings1, Anita *anita1, Screen *panel1, int ilayer, int ifold, double hitangle_e, double hitangle_h, double e_component, double h_component);
  
    // just for historical reference, this function:
    //void AntTrigger::RunTrigger(Settings *settings1,int ilayer,int ifold,double *vmmhz, Screen *panel1, Anita *anita1,double hitangle_e,double hitangle_h,double e_component,double h_component,double *arrival_times,double volts_rx_rfcm_lab_e_all[48][512],double volts_rx_rfcm_lab_h_all[48][512])
    // was split into InitializeEachBand, ConvertInputWFtoAntennaWF, ImpulseResponse, TimeShiftAndSignalFluct, Banding

    static void ConvertEHtoLREfield(double,double,double&,double&);
    static void ConvertEHtoLREnergy(double,double,double&,double&);
    void ConvertHVtoLRTimedomain(const int nfour,double *vvolts,
				 double *hvolts,
				 double *left,double *right);
    //int Passes(double strength,double angle,int trigger_band); // whether a particular channel passes or not
  void L1Trigger(Anita *anita1,double timedomain_output_1[5][Anita::NFOUR],double timedomain_output_2[5][Anita::NFOUR],double powerthreshold[2][5],int *channels_passing_e_forglob,int *channels_passing_h_forglob,int &npass);    

    vector<int> flag_e[5];
    vector<int> flag_h[5];
    
    // inputs are:  signal strength, angle off cone, and the trigger band
    //void getFunctions(int thresh,int isurf,int ichan,double viewangle);
    double getRate();
    double rateToThreshold(double rate, int band); // converts a single channel rate to threshold in p/<p>
    static double GetNoise(Settings *settings1,double altitude_bn,double geoid,double theta,double bw,double temp);
  void WhichBandsPass(int inu,Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double dangle, double emfrac, double hadfrac);
  void WhichBandsPassTrigger1(Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double thresholds[2][5]);
  void WhichBandsPassTrigger2(int inu,Settings *settings1, Anita *anita1, GlobalTrigger *globaltrig1, Balloon *bn1, int ilayer, int ifold, double dangle, double emfrac, double hadfrac, double thresholds[2][5]);
  static double FindPeak(double *waveform,int n); // find peak voltage of a waveform
    void GetThresholds(Settings *settings1,Anita *anita1,int ilayer,double thresholds[2][5]); // get thresholds for this layer
    
    double bwslice_volts_pol0[5];  // sum voltage for each slice in bandwidth for the lcp polarization
    double bwslice_volts_pol1[5]; // same, for rcp polarization
    
    double bwslice_energy_pol0[5];  // square the sum of voltage for each slice in bandwidth for the 0th polarization
    double bwslice_energy_pol1[5]; // same, for 1st polarization
    
    double bwslice_volts_pol0_em[5];  // component of the voltage that comes from the em shower
    double bwslice_volts_pol1_em[5]; // same, for 1st polarization
    
    double bwslice_volts_pole[5];
    double bwslice_energy_pole[5]; // square the sum of voltage for each slice in bandwidth.  The 5th element is the full band
    double bwslice_volts_polh[5];
    double bwslice_energy_polh[5];
    
        
    double v_banding_rfcm_e[5][Anita::NFREQ];// this is Volts/m as a function of frequency after rfcm's and banding
    double v_banding_rfcm_h[5][Anita::NFREQ];
    
    //static const double bwslice_center[4]; // center frequencies
    //static const double bwslice_width[4]; // 3 dB bandwidths, without overlap
    
    // these are filled for triggerscheme==0 and triggerscheme==1
    // frequency domain voltage and energy based
    double signal_eachband[2][Anita::NBANDS_MAX];
    double threshold_eachband[2][Anita::NBANDS_MAX];
    double noise_eachband[2][Anita::NBANDS_MAX];
    int passes_eachband[2][Anita::NBANDS_MAX];
    
    vector<double> vsignal_eachband[2];
    vector<double> vthreshold_eachband[2];
    vector<double> vnoise_eachband[2];
    vector<int> vpasses_eachband[2];
    
    // Used for AntTrigger::PrepareBandWaveforms(...) and AntTrigger::WhichBandsPass(...)
    double v_banding_rfcm_e_forfft[5][HALFNFOUR]; // starts out as V/s vs. freq after banding, rfcm, after fft it is V vs. t
    double v_banding_rfcm_h_forfft[5][HALFNFOUR];
    double vm_banding_rfcm_1_forfft[5][HALFNFOUR];
    double vm_banding_rfcm_2_forfft[5][HALFNFOUR];
    double v_banding_rfcm_e_forfft_temp[5][HALFNFOUR];
    double v_banding_rfcm_h_forfft_temp[5][HALFNFOUR];
    // End of band waveform triggering arrays
    
    double integral_vmmhz;
    double integral_vmmhz_r[Anita::NFREQ];
    
    // this increments the bwslice_volts_...'s at each frequency
    void addToChannelSums(Settings *settings1,Anita *anita1,int ibw,int k);
    
    static int IsItUnmasked(unsigned short surfTrigBandMask[9][2],int ibw,int ilayer, int ifold, int ipol);
#ifdef ANITA_UTIL_EXISTS
  void applyImpulseResponseDigitizer(Settings *settings1, Anita *anita1, int nPoints, int ant, double *x, double y[512], bool pol);
  void applyImpulseResponseTrigger(Settings *settings1, Anita *anita1, int nPoints, int ant, double *x, double y[512], double *vhz, bool pol);
  void applyImpulseResponseTrigger(Settings *settings1, Anita *anita1, int ant, double y[512], bool pol);
  double *getNoiseFromFlight(Anita* anita1, int pol, int ant);
  void injectImpulseAfterAntenna(Anita *anita1, double volts_triggerPath_e[Anita::HALFNFOUR], double volts_triggerPath_h[Anita::HALFNFOUR], int ant);
  void injectImpulseAmplitudeAfterAntenna(Anita *anita1, double vhz_triggerPath_e[Anita::NFREQ], double vhz_triggerPath_h[Anita::NFREQ], int ant);
  void injectImpulseAtSurf(Anita *anita1, double volts_triggerPath_e[Anita::HALFNFOUR], double volts_triggerPath_h[Anita::HALFNFOUR], int ant);
#endif
    
    int unwarned;  // whether we have warned the user about resetting thresholds when they are beyond the measured bounds
}; //class AntTrigger

//namespace Bands {
    //   double bwslice_max[4]=; //maximum of each bandwidth slice  
    //   for (int i=0;i<4;i++) {
    //     bwslice_min[i]=bwslice_center[i]-bwslice_width[i]/2; // get low edge of bandwidth slices
    //     bwslice_max[i]=bwslice_center[i]+bwslice_width[i]/2;  // get upper edge of bandwidth slices
    //   } //for
//}


// This is just a place holder
//double power_thresholds[5]={3.,3.,3.,3.,3.};

#endif

