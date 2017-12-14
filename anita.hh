////////////////////////////////////////////////////////////////////////////////////////////////
//class Anita:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef ICEMC_ANITA_HH
#define ICEMC_ANITA_HH


#include "rx.hpp"
#include <array>

#ifdef ANITA_UTIL_EXISTS
#include "FFTtools.h"
#endif

#include "TF1.h"

class RX;
class TGraph;
class TFile;
class TTree;
class TH1F;
class Vector;
class Position;
class Balloon;
class IceModel;
class Settings;
class TRandom3;

using std::ofstream;
using std::vector;
using std::array;
//! Contains everything about positions within payload and signals it sees for each event, in both the trigger and signal paths.
class Anita {

private:

  std::string stemp;

  double GaintoHeight(double gain,double freq,double nmedium_receiver);
  TGraph *gshort[4];
  void setTrigRequirement(int WHICH);


public:

  
  int number_all_antennas;                                                                                       ///< this keeps count of the number of antennas for use with timing calculations, etc.

  static const int NBANDS_MAX=100;                                                                               ///< max number of bands
  static const int NPOL=2;                                                                                       ///< number of polarizations
  //static const int NFREQ=128;                                                                                  ///< number of frequency bins
  static const int NFREQ=128;
  //const int NFREQ=4096;


  static const int NTRIG=5;
  static const int NANTENNAS_MAX=2000;
  static const int NLAYERS_MAX=5;                                                                                ///< max number of layers (in smex design, it's 4)
  static const int NTRIGGERLAYERS_MAX=3;
  static const int NPHI_MAX=400;                                                                                 ///< max number of antennas around in phi (in smex, 16)
  Vector ANTENNA_POSITION_START[2][NLAYERS_MAX][NPHI_MAX];                                                          ///< antenna positions from Kurt's measurements
  double ANTENNA_DOWN[NLAYERS_MAX][NPHI_MAX];                                                                    ///< down angles of antennas from Kurt's measurements
  double SIMON_DELTA_R[NLAYERS_MAX][NPHI_MAX];                                                                   ///< measurements by Simon used in analysis ANITA-2
  double SIMON_DELTA_PHI[NLAYERS_MAX][NPHI_MAX];                                                                 ///< measurements by Simon used in analysis ANITA-2

  Vector antenna_positions[2][NLAYERS_MAX * NPHI_MAX];                                                              ///< these are the antenna positions in space in a coordinate system where x=north and y=west and the origin is at the center of the payload

  int NRX_PHI[NLAYERS_MAX];                                                                                      ///< number of antennas around in each layer. (radians)
  double PHI_EACHLAYER[NLAYERS_MAX][NPHI_MAX];                                                                   ///< phi of the center of each antenna on each layer
  
   //before correcting for offset for the layer.
  //only used if it is cylindrically symmetric (radians)
  double PHI_OFFSET[NLAYERS_MAX];                                                                               ///< antenna offset in phi for each layer (radians)
  double THETA_ZENITH[NLAYERS_MAX];                                                                             ///< how the antenna is tilted in theta (in radians with 0=up)
  // 0=horizontal,+90=down

  int inu;            ///< Neutrino number
  // what the payload looks like

  double LAYER_VPOSITION[Anita::NLAYERS_MAX];                                                                   ///< position of layers in z relative to vertical center of the payload
  // anita proposal "says that the separation between upper and lower
  // 2 layers of antennas is just under 4m.
  // for anita hill, consider the positions of the "layers" of the "payload" (the stations) to be sitting on the horizontal grid defined by polar coordinates
  double LAYER_HPOSITION[Anita::NLAYERS_MAX];                                                                   ///< distance in horizontal plane between center axis of the "payload" and each "layer".
  double LAYER_PHIPOSITION[Anita::NLAYERS_MAX];                                                                 ///< phi corresponding to the position of each "layer" on the "payload"
  double RRX[Anita::NLAYERS_MAX];                                                                               ///< radius that the antenna sits from the axis of the payload (feedpoint)
  Double_t deltaTPhaseCentre[2][NLAYERS_MAX][NPHI_MAX];                                                         ///< Relative to photogrammetry + ring offset

  double THERMALNOISE_FACTOR;                                                                                   ///< factor to multiply thermal noise for error analysis


  Anita(); // constructor
  ~Anita();
  void Initialize(Settings *settings1,ofstream &foutput,int inu, TString outputdir);                                                ///< initialize a bunch of stuff
  void initializeFixedPowerThresholds(ofstream &foutput);
  void readVariableThresholds(Settings *settings1);
  void readAmplification();
  void getDiodeDataAndAttenuation(Settings *settings1, TString outputdir);
  void getPulserData();
  
  // takes arrays that span NFREQ and turn them into arrays that span HALFNFOUR
  void MakeArrayforFFT(double *vsignalarray_e,double *vsignal_e_forfft, double phasedelay, bool useconstantdelay);
  
  void GetArrayFromFFT(double *tmp_fftvhz, double *vhz_rx);
  
  int Match(int ilayer,int ifold,int rx_minarrivaltime);
  int getLabAttn(int NPOINTS_LAB,double *freqlab,double *labattn);

  void labAttn(double *vhz);
  void SetNoise(Settings *settings1,Balloon *bn1,IceModel *antarctica);
  void calculate_antenna_positions(Settings *settings1,double pitch, double roll, double phi_spin,Vector n_north,Vector n_east);// this calculates the above

  TFile *fnoise;
  TTree *tdiode;

  static const int NFOUR=1024;
  static const int HALFNFOUR=512;

  // these are used for the satellite thing
  int NBANDS;                                                                                        ///< number of frequency sub-bands (not counting full band)

  int PERCENTBW;                                                                                     ///< percent bandwidth

  // these variables are for filling the tsignals tree
  double signal_vpol_inanita[5][HALFNFOUR];                                                          ///< this is the signal waveform in the vertical polarization, before converting to LCP, RCP where applicable
  //double noise_vpol_inanita[5][HALFNFOUR];                                                         ///< this is the noise waveform in the vertical polarization, before converting to LCP, RCP where applicable
  double total_vpol_inanita[5][HALFNFOUR];                                                           ///< this is the sum of the signal and noise in the vertical polarization, before converting to LCP, RCP where applicable

  double timedomainsignal_rfcm[HALFNFOUR];
  double timedomainsignal_lab[HALFNFOUR];

  TTree *turfratechain;
  TTree *surfchain;
  TFile *fturf;
  TFile *fsurf;

  UShort_t phiTrigMask;
  UShort_t phiTrigMaskH;
  UShort_t l1TrigMask;
  UShort_t l1TrigMaskH;
  Double_t deadTime;                                                                                 ///< fractional deadTime
  unsigned int realTime_turfrate;                                                                    ///< realtime from the turf rate file
  unsigned int realTime_tr_min;                                                                      ///< min realtime from the turf rate file
  unsigned int realTime_tr_max;                                                                      ///< max realtime from the turf rate file
  unsigned int realTime_surf;                                                                        ///< realtime from the surf file
  unsigned int realTime_surf_min;                                                                    ///< min realtime from the surf file
  unsigned int realTime_surf_max;                                                                    ///< max realtime from the surf file
  UShort_t thresholds[2][48];                                                                        ///< thresholds as read from the surf file: first index is pol, second is antenna number (only working for Anita3)
  UShort_t scalers[2][48];                                                                           ///< scalers as read from the surf file: first index is pol, second is antenna number (only working for Anita3)
  Double_t fakeThresholds[2][48];                                                                    ///< Fake thresholds (coming from converting fake scalers to thresholds)
  Double_t fakeThresholds2[2][48];                                                                   ///< Fake thresholds 2 (coming from converting flight scalers to thresholds)
  Double_t fakeScalers[2][48];                                                                       ///< Fake scalers (coming from converting threhsolds during flight to scalers using threshold scan)

  int iturf;// for indexing
  int isurf;
  int iturfevent;

  static const int npointThresh = 1640;
  Float_t threshScanThresh[2][48][npointThresh];                                                     ///< adc thresholds from threshold scan
  Float_t threshScanScaler[2][48][npointThresh];                                                     ///< scalers from threshold scan
  Float_t minadcthresh[2][48];
  Float_t maxadcthresh[2][48];

  void setphiTrigMaskAnita3(UInt_t realTime_flightdata);
  void setphiTrigMask(UInt_t realTime_flightdata);
  void setTimeDependentThresholds(UInt_t realTime_flightdata);

  double total_diodeinput_1_inanita[5][HALFNFOUR];                                                   ///< this is the waveform that is input to the tunnel diode in the first (LCP or vertical) polarization
  double total_diodeinput_2_inanita[5][HALFNFOUR];                                                   ///< this is the waveform that is input to the tunnel diode in the second (RCP or horizontal) polarization

  double total_diodeinput_1_allantennas[48][HALFNFOUR];                                              ///< this is across all antennas, just the full band
  double total_diodeinput_2_allantennas[48][HALFNFOUR];                                              ///< needs comment

  // these are just for the antenna that receives the signal first
  double timedomain_output_inanita[2][5][HALFNFOUR];                                                 ///< this is just for writing out to the following tree

  double time_trig[HALFNFOUR];
  double weight_inanita; // weight of the event
  int arrayofhits_inanita[3][16][2][HALFNFOUR];
    
  //std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>  arrayofhits_inanita;



  // same as arrayofhits_inanita but it's time reversed
  int arrayofhits_forgaryanderic[3][16][2][HALFNFOUR];
  //std::array< std::array< std::array< std::array<std::vector<int>,5>, 2>, 16>, 3>  arrayofhits_inanita;

  int l1trig_anita3and4_inanita[2][16][HALFNFOUR];



  int l1trig_anita4lr_inanita[3][16][HALFNFOUR];

  int l1trig_anita4lr_forgaryanderic[3][16][HALFNFOUR];


  int l2trig_anita4lr_inanita[16][3][HALFNFOUR];

  int l2trig_anita4lr_forgaryanderic[16][HALFNFOUR];                                                           ///< when it passes 2/3

  int l3type0trig_anita4lr_inanita[16][HALFNFOUR];
  int l3trig_anita4lr_inanita[16][HALFNFOUR];

  int l3type0trig_anita4lr_forgaryanderic[16][HALFNFOUR];
  int l3type1trig_anita4lr_forgaryanderic[16][HALFNFOUR];




  double timedomain_output_corrected_forplotting[2][6][HALFNFOUR];                                            ///< this is just for writing out to the following tree


  double timedomain_output_allantennas[2][48][HALFNFOUR];                                                     ///< this is across all antennas, just the full band



  int flag_e_inanita[5][HALFNFOUR];
  int flag_h_inanita[5][HALFNFOUR];
  double dangle_inanita,emfrac_inanita,hadfrac_inanita;
  double ston[5];                                                                                             ///< signal to noise;

  int iminbin[5];                                                                                             ///< this is the minimum bin to start
  int imaxbin[5];
  int maxbin_fortotal[5];                                                                                     ///< when it sums the noise and signal together it shortens the waveform

  double peak_v_banding_rfcm[2][5];                                                                           ///< peak V in e/h polarization after rfcm's and banding
  double peak_rx_signalonly[2];                                                                               ///< peak voltage in e/h polarization received by the antenna
  double peak_rx_rfcm[2];                                                                                     ///< peak voltage in e/h polarization received by the antenna

  double peak_rx_rfcm_signalonly[2];                                                                          ///< peak voltage in e/h polarization received by the antenna

  double peak_rx_rfcm_lab[2];                                                                                 ///< peaks of the previous arrays
  //I think this is the numerator of the vertical axis on Matt's plot




  int channels_passing[2][5];                                                                                 ///< channels passing.  This is reset for every antenna for every event
  int l1_passing; // l1 passing
  int l1_passing_allantennas[48]; // l1 passing
    
  int irx;
  void BoxAverageComplex(double *array,const int n,int navg);
  void BoxAverage(double *array,const int n,int navg);
  int GetRx(int ilayer, int ifold);                                                                           ///< get antenna number based on which layer and position it is
  int GetRxTriggerNumbering(int ilayer, int ifold);                                                           ///< get antenna number based on which layer and position it is


  double avgfreq_rfcm[NFREQ];
  double avgfreq_rfcm_lab[NFREQ];



  double vmmhz_banding[NFREQ];                                                                                ///< V/m/MHz after banding
  double vmmhz_banding_rfcm[NFREQ];                                                                           ///< V/m/MHz after banding and rfcms

  double rms_rfcm_e_single_event;                                                                             ///< This is in Volts, not mV!

  // Note: The following 4 RMS noise variables are for all antennas of all events.
  // In fact, they don't represent RMS until after all events are finished!
  double rms_rfcm[2];                                                                                         ///< rms noise just after rfcm's
  double rms_lab[2];                                                                                          ///< rms noise at lab chip


  TFile *fsignals;
  TTree *tsignals;

  TFile *fdata;
  TTree *tdata;                                                                                               ///< writing data out for the analysers
  TTree *tgaryanderic;                                                                                        ///< writing data out for the analysers

  TTree *tglob;

  TH1F *hsignals[5];                                                                                          ///< s/n (max diode output/mean diode output) for vertical polarization in each band

  double f_pulser[NFOUR/4];
  double f_phases[NFOUR/4];
  double f_noise[NFOUR/4];
  double v_pulser[NFOUR/4];
  double v_phases[NFOUR/4];
  double v_noise[NFOUR/4];



  double cumulat_prob[9];
  double cumulat_prob_plus1[9];


  // for filling tsignals tree
  double timedomainnoise_rfcm_banding[2][5][HALFNFOUR];
  double timedomainnoise_rfcm_banding_long[2][5][HALFNFOUR];
  double timedomainnoise_rfcm[2][HALFNFOUR];
  double timedomainnoise_lab[2][HALFNFOUR];
  double timedomainnoise_rfcm_long[2][HALFNFOUR];
  double timedomainnoise_lab_long[2][HALFNFOUR];

  double phases[5][HALFNFOUR];

  // for filling tglob
  // for each polarization
  int passglobtrig[2];
  double integral_vmmhz_foranita;


  int nnoiseevents;                                                                                          ///< total number of noise events we're choosing from
  int noiseeventcounter;                                                                                     ///< counts which event we're on so we go in order

  double FREQ_LOW;                                                                                           ///< lowest frequency
  double FREQ_HIGH;                                                                                          ///< highest frequency

  double NOTCH_MIN;                                                                                          ///< low edge of notch filter.  This is set in the input file
  double NOTCH_MAX;                                                                                          // upper edge of notch filter

  int BANDING;// set in the input file
  // whether or not you set your own banding (1)
  // or use anita-1 banding

  double freq[NFREQ];  // frequency for each bin
  double freq_forfft[NFOUR]; // frequencies for taking fft of signal
  double freq_forplotting[NFOUR/4]; // just one entry for frequency, unlike the above.
  double freq_forfft_long[2*NFOUR]; // frequencies for taking fft of signal
  double freq_forplotting_long[NFOUR/2]; // just one entry for frequency, unlike the above.
  double time[NFOUR/2];
  double time_long[NFOUR];

  double time_centered[NFOUR/2];
  double freqdomain_rfcm_banding[5][HALFNFOUR/2]; // average noise in frequency domain
  double freqdomain_rfcm_banding_long[5][HALFNFOUR]; // average noise in frequency domain

  double freqdomain_rfcm[HALFNFOUR/2]; // average noise in frequency domain
  double freqdomain_rfcm_long[HALFNFOUR]; // average noise in frequency domain

  double freqdomain_rfcm_theory[HALFNFOUR/2]; // average noise in frequency domain
  double avgfreqdomain_lab[HALFNFOUR/2]; // average noise in frequency domain
  double avgfreqdomain_lab_long[HALFNFOUR]; // average noise in frequency domain

  double phases_rfcm_banding[2][5][HALFNFOUR/2];
  double phases_rfcm_banding_long[2][5][HALFNFOUR];
  double phases_rfcm[2][HALFNFOUR/2];
  double phases_rfcm_long[2][HALFNFOUR];
  double phases_lab[2][HALFNFOUR];
  double phases_lab_long[2][HALFNFOUR];


  // this goes from 0 to Fmax, and represents both real and imaginary components




  void getDiodeModel();
  void setDiodeRMS(Settings *settings1, TString outputdir);
  
  TF1 fdiode;
  double maxt_diode;
  int idelaybeforepeak[5];
  int iwindow[5];
  double diode_real[5][NFOUR]; // This is the time domain of the diode response. (actually NFOUR/2 array is used.)
  double fdiode_real[5][NFOUR]; // This is the fft of the diode response. (use all NFOUR array. This is for doubling array size for zero padding)


  void myconvlv(double *timedomain_forconvl,const int NFOUR,double *fdiode,double &maxdiodeconvl,double &onediodeconvl,double *power_noise,double *diodeconv);

  void GetArrivalTimes(const Vector& rf_direction,Balloon *bn1,Settings *settings1);
  void GetArrivalTimesBoresights(const Vector rf_direction[NLAYERS_MAX][NPHI_MAX]);

  void GetArrivalTimesBoresights(const Vector rf_direction[NLAYERS_MAX][NPHI_MAX],Balloon *bn1, Settings *settings1);

  int rx_minarrivaltime;
  double arrival_times[2][NLAYERS_MAX*NPHI_MAX];

  static int SurfChanneltoBand(int isurf);
  int AntennaWaveformtoSurf(int ilayer,int ifold); // find surf that generates this antenna's waveform
  static int AntennaNumbertoSurfNumber(int ilayer,int ifold); // find surf where this antenna is triggered
  static int GetAntennaNumber(int ilayer,int ifold); // given icemc indices ilayer, ifold, find antenna number as defined officially on anita
  static int GetLayer(int rx);
  static int GetIfold(int rx);
  static int GetSurfChannel(int antenna, int ibw,int ipol); // which channel on the surf this channel on this antenna corresponds to.
  static int WhichBand(int ibw,int ipol); // which band, 1-8, in order as they are on the surf
  void Banding(int j,double *freq_noise,double *powerperfreq,int NPOINTS_NOISE);
  void Banding(int iband,double *vmmhz);
  void RFCMs(int ilayer,int ifold,double *vmmhz);
  void normalize_for_nsamples(double *spectrum, double nsamples, double nsamp);
  void convert_power_spectrum_to_voltage_spectrum_for_fft(int length,double *spectrum, double domain[], double phase[]);
  void GetNoiseWaveforms(); // make time domain noise waveform based on avgnoise being the v^2
  //void GetNoiseWaveform(int iband); // make time domain noise waveform based on avgnoise being the v^2
  void GetPhases();
  int count_getnoisewaveforms; //count how many times we're in GetNoiseWaveforms for calculating rms voltages

  // each of the above graphs has 601 bins in it
  static const int NPOINTS_BANDS=601;

  double freq_bands[5][NPOINTS_BANDS]; // a frequency array for each of the four bands
  double attn_bands[5][NPOINTS_BANDS]; // attn array for each of the four bands in dB
  double bandsattn[5][NPOINTS_BANDS]; // as a fraction
  //double correl[4][NPOINTS_BANDS]; // correlation between each band and the fullband
  double correl_banding[5][NPOINTS_BANDS]; // correlation between each band and the fullband
  double correl_lab[NPOINTS_BANDS]; // correlation between each band and the fullband
  //double correl[5][NPOINTS_BANDS]; // correlation between each band and the fullband


  static const int NPOINTS_AMPL=58;// bins in amplification
  double freq_ampl[NANTENNAS_MAX][NPOINTS_AMPL]; // frequencies for each amp bin
  double ampl[NANTENNAS_MAX][NPOINTS_AMPL]; // amplification
  double ampl_notdb[NANTENNAS_MAX][NPOINTS_AMPL];// amplification again, but as a fraction
  double noisetemp[NANTENNAS_MAX][NPOINTS_AMPL]; // noise temp each amp bin


  static const int NPOINTS_NOISE=2000;


  //double bwslice_thresholds[4]={2.319,2.308,2.300,2.290}; // this allows you to set different thresholds for each band
  double bwslice_vnoise[NLAYERS_MAX][5]; // expected noise voltage for antenna layer and
  //each slice in bandwidth

  double probability[5];
  double bwslice_enoise[5]; // average integrated power in each band
  double bwslice_fwhmnoise[5]; // 1/2 of fwhm of hnoise
  double bwslice_rmsdiode[5]; // average rms diode output across noise waveforms in each band
  double bwslice_meandiode[5]; // mean diode output across all samples in a sample of noise waveforms generated for each band
  double bwslice_vrms[5]; // rms noise voltage for this bandwidth slice
  double bwslice_dioderms_fullband_allchan[2][48]; // diode rms for noise read from flight
  double bwslice_diodemean_fullband_allchan[2][48]; // diode rms for noise read from flight
  double freq_noise[5][NPOINTS_NOISE]; // frequency array that goes with vnoise array


  double impedence;
  double phase;
  double powerthreshold[5];
  double powerthreshold_nadir[5];
  int NCH_PASS; // for ANITA 3 trigger - requires some number of channels pass

  double l1window; // time window where we require coincidences at L1

  double minsignalstrength; // minimum signal strength (measured as output of the diode) that a signal has to be for it to be worth adding to noise and performing the diode integration (each time we do this is uses up a noise waveform)

  double INTEGRATIONTIME; // integration time of the tunnel diode
  static const int nsamp=100; // number of samples that were used to measure the noise data
  double TIMESTEP; // time step between samples for digitization


  double DEADTIME;

  double TRIG_TIMESTEP; // this is the l1 trigger window for the anita 3 trigger.
  unsigned N_STEPS_PHI;
  unsigned N_STEPS_THETA;

  static const unsigned N_SUMMED_PHI_SECTORS = 4;
  static const unsigned N_SUMMED_LAYERS = 3;



  double energythreshold; // relative to expected energy from noise
  double MIN_PHI_HYPOTHESIS;
  double MAX_PHI_HYPOTHESIS;
  double MIN_THETA_HYPOTHESIS;
  double MAX_THETA_HYPOTHESIS;

  int USEPHASES;

  int NTRIGGERLAYERS; // number of layers considered by the trigger.  may be different from nlayers, the number of physical layers on the payload.
  // In Anita 1 and Anita 2, the number of physical layers were 3 while the number of trigger layers were 2.
  int PHITRIG[NLAYERS_MAX]; // number of positions in phi for each trigger layer
  int REQUIRE_CENTRE; // require centre antenna in clump to be one of those hit
  static const int NTRIGPHISECTORS=16; // number of phi sectors in the trigger

  int GAINS;// whether to use constant gains as entered in GetBeamWidths (0) or to use Ped's measurements as entered in ReadGains (1)
  static const int NPOINTS_GAIN =131; // number of pointqs within bandwidth that gain is measured. from 200 MHz to 1.5 GHz with step size of 10 MHz
  double gainv_measured[NPOINTS_GAIN]; // may way of making the program use 3994760 fewer bytes than if these five arrays had still had 100000 elements
  double gainh_measured[NPOINTS_GAIN];
  double gainhv_measured[NPOINTS_GAIN];
  double gainvh_measured[NPOINTS_GAIN];
  double frequency_forgain_measured[NPOINTS_GAIN];

  double gain_angle[4][NPOINTS_GAIN][7]; /* first term: 0 = v polarization channel, a angle
					    1 = h polarization channel, a angle
					    2 = h polarization channel, e angle
					    3 = v polarization channel, e angle
					    second term: frequency
					    third term: angle */

  // frequency binning
  // anita proposal says frequency range is 0.2-1.2 GHz.
  // specs for the quad ridge antenna Model 0312-810 say 0.3-1.5 GHz
  double flare[4][NFREQ];  // for coarse antenna models:  beam width: e-plane: vp/hp, h-plane: vp/hp
  double gain[2][NFREQ];   // for coarse antenna models:  gain vert pol,h pol

  int GetBeamWidths(Settings *settings1); // for getting beam widths using coarse models (horn specs or simple model for EeVA)
  void Set_gain_angle(Settings *settings1,double nmedium_receiver);
  double Get_gain_angle(int gain_type, int k, double hitangle);
  void ReadGains();
  void AntennaGain(Settings *settings1,double hitangle_e,double hitangle_h,double e_component,double h_component,int k,double &vsignalarray_e,double &vsignalarray_h);


  double reference_angle[7]; // reference angles for finding gains of antenna

  double inv_angle_bin_size[6];
  int whichbin[NFREQ]; // these are for finding gains as a function of frequency
  double scalef2[NFREQ], scalef1[NFREQ]; // they are set in Set_gain_angle
  double vvGaintoHeight[NFREQ], hhGaintoHeight[NFREQ], hvGaintoHeight[NFREQ], vhGaintoHeight[NFREQ]; // holds results of the function double GaintoHeight

  double diffraction[2][89][NFREQ];
  void SetDiffraction();
  double GetDiffraction(int ilayer, double zenith_angle, int ifreq);



  static const int NPOINTS_LAB=272; // from note 137

  double freqlab[NPOINTS_LAB]; // frequency for each lab attn. bin

  double labattn[NPOINTS_LAB]; // lab attenuation


  double VNOISE[NLAYERS_MAX]; // noise calculated for each antenna layer depending on cant angle- this is only used right now for the chance in hell cuts


  int trigRequirements[NLAYERS_MAX];//  0th element - L1 - how many channels per antenna should pass
  // 1st element- L2 - how many antennas on a layer
  // 2nd element - L3 - how many L2 triggers should be coincident

  int antennatosurf[32];

  double maxthreshold;
  double bwslice_thresholds[5]; // thresholds for each band -- this is just an initialization- this is set in the input file
  int bwslice_allowed[5]; // these bands are allowed to contribute to the trigger sum -- this is set in the input file
  int bwslice_required[5]; // these bands are required to be among the channels that pass -- this is set in the input file
  int pol_allowed[2];// which polarisations are allowed to have channels that fire (V,H)
  int pol_required[2];// which polarisations are required to have channels that fire (V,H)



  double bwslice_center[5]; // center frequencies
  double bwslice_width[5]; // 3 dB bandwidths, without overlap


  double bwslice_min[5]; //minimum of each bandwidth slice

  double bwslice_max[5]; //minimum of each bandwidth slice
  double bwmin; // minimum width of any allowed bandwidth slice

  TRandom3* summed_power_trigger_digitizer_zero_random;
  TFile* coherent_datafile;
  TTree* coherent_waveform_sum_tree;
  static const unsigned int NUM_COHERENT_ANTENNAS = 9;
  unsigned hypothesis_offsets[16][200][200][4][3]; // Time bin offsets for each hypothesis - [center_phi_sector_index][phi_angle_index][theta_angle_index][phi_sector_index][layer_index]
  vector< vector< vector <double> > > hypothesis_angles; // Time bin offsets for each hypothesis - [center_phi_sector_index][phi_angle_index][theta_angle_index][phi_sector_index][layer_index]
  //unsigned antenna_indices[16][9];	// These are the indices of the antennas used for a given hypothesis' center phi center index - [center_phi_sector-index][which of the nine]
  vector< vector <int> > vdifferent_offsets;
  vector< vector <double> > vdifferent_angles;

  void calculate_all_offsets(void);	// This function creates offsets for coherent sum trigger
  void getDifferentOffsets();
  void printDifferentOffsets();
  void calculate_single_offset(const unsigned center_phi_sector_index, const double angle_phi, const double angle_theta, double hypothesis_offset[][3]);
  void calculate_single_offset(const unsigned center_phi_sector_index, const unsigned index_phi, const unsigned index_theta, double hypothesis_offset[][3]);
  unsigned cwst_event_number;
  unsigned cwst_center_phi_sector;
  double cwst_rms_noise;
  double cwst_actual_rms;
  double cwst_threshold;
  unsigned cwst_window_start;
  unsigned cwst_window_end;
  double cwst_deg_theta;
  double cwst_deg_phi;
  double cwst_actual_deg_theta;
  double cwst_actual_deg_phi;
  Vector cwst_rf_direction;
  Vector cwst_0th_sector_position;
  double cwst_timesteps[HALFNFOUR];
  RX cwst_RXs[48];
  RX cwst_aligned_wfms[9];
  //vector <double>* cwst_whole_wfms[NUM_COHERENT_ANTENNAS];
  //vector <double>* cwst_wfms[NUM_COHERENT_ANTENNAS];
  //vector <double>* cwst_aligned_wfms[NUM_COHERENT_ANTENNAS];
  vector <double> cwst_summed_wfm;
  vector <double> cwst_power_of_summed_wfm;
  double cwst_power;
  void fill_coherent_waveform_sum_tree(unsigned inu, unsigned center_phi_sector, Settings* settings1, double rms_noise, double actual_rms, unsigned window_start, unsigned window_end, double deg_theta, double deg_phi, double actual_deg_theta, double actual_deg_phi, vector <double>& summed_wfm, vector <double>& power_of_summed_wfm, double power);
  void GetPayload(Settings*, Balloon*);
  double VNOISE_ANITALITE[NPHI_MAX]; // noise for each antenna, for the anita-lite trigger configuration.
  double INCLINE_TOPTHREE; // cant angle of top three layers of antennas
  double INCLINE_NADIR; // cant angle of nadir (bottom) layer of antennas
  double LIVETIME;

  double SIGMA_THETA; // resolution on the polar angle of the signal

  double extraCableDelays[2][48];
  TRandom3 *fRand;
#ifdef ANITA_UTIL_EXISTS
  void readImpulseResponseDigitizer(Settings *settings1);
  void readImpulseResponseTrigger(Settings *settings1);
  void readTriggerEfficiencyScanPulser(Settings *settings1);
  void readNoiseFromFlight(Settings *settings1);
  void getQuickTrigNoiseFromFlight(double justNoise[HALFNFOUR], int ipol, int iant);
  TGraph *RayleighFits[2][48];
  Int_t numFreqs;
  Double_t *freqs;
  TGraph *gPulseAtAmpa;
  RFSignal *fSignalChainResponseDigitizer[2][3][16]; // 0:VPOL, 1:HPOL ---- 0:TOP, 1:MIDDLE, 2:BOTTOM
  RFSignal *fSignalChainResponseTrigger[2][3][16]; // 0:VPOL, 1:HPOL ---- 0:TOP, 1:MIDDLE, 2:BOTTOM
#endif
  void calculateDelaysForEfficiencyScan();

  void GetPhasesFromFFT(double *tmp_fftvhz, double *phases);
  void FromTimeDomainToIcemcArray(double *vsignalarray, double vhz[NFREQ]);

  
  Double_t fTimes[HALFNFOUR];
  Double_t fSignalChainResponseDigitizerFreqDomain[2][3][16][400];
  Double_t fSignalChainResponseTriggerFreqDomain[2][3][16][400];
  Double_t fRatioTriggerDigitizerFreqDomain[2][3][16][400];
  Double_t deltaT;

  // Trigger efficiency scan parameters
  Int_t trigEffScanPhi;                      // central phi sector of trigger efficiency scan
  Double_t trigEffScanAtt[5];                // attenuations to apply to central and adjecent antennas
  Double_t trigEffScanPhiDelay[5];           // delays between phi sectors
  Double_t trigEffScanRingDelay[3];          // delays between rings
  Int_t    trigEffScanApplyRingDelay[5];     // to which phi sectors apply ring delays 
  Int_t    trigEffScanRingsUsed[3];          // to which rings apply scan
  Double_t trigEffScanPulseAtAmpa[HALFNFOUR];
  Double_t trigEffScanPulseAtAmpaUpsampled[NFOUR];
  Double_t trigEffScanAmplitudeAtAmpa[NFREQ];
  Double_t trigEffScanPulseAtSurf[250][HALFNFOUR];
  int TUFFstatus[3];

}; //class Anita




#endif //ICEMC_ANITA_HH
