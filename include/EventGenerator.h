#ifndef ICEMC_EVENT_GENERATOR_H
#define ICEMC_EVENT_GENERATOR_H

#include "anita.hh"
#include "TVector3.h"

#include "Settings.h"

#include "NeutrinoPath.h"
#include "Constants.h"
#include "CommandLineOpts.h"
#include "secondaries.hh"


namespace icemc {

  class Taumodel;
  class AskaryanFreqsGenerator;
  class Interaction;
  class RayTracer;
  class Roughness;
  class Screen;
  class Counting;
  class Spectra;
  class Primaries;
  class GeneratedNeutrino;
  class PassingNeutrino;
  class ANITA;


  /**
   * @class EventGenerator
   * @brief The glue that brings the simulation together
   * 
   * Contains the main neutrino generation loop in generateNeutrinos()
   */

  class EventGenerator {
  public:

    enum RayDirection {
      direct = 0,
      downgoing = 1,
    };

    EventGenerator();
    virtual ~EventGenerator();
    static void interrupt_signal_handler(int);  // This catches ctrl-C interrupt (SIGINT)

    static const int NVIEWANGLE=100; // number of viewing angles to look at the signal, on Seckel's request
    // int inu; // counts neutrinos as they are generated
    double eventsfound_beforetrigger=0.;
    double eventsfound_crust=0; //number of events that only traverse the crust
    double eventsfound_weightgt01=0; // summing weights > 0.1
    double eventsfound_belowhorizon=0; // how many are below horizon
    double eventsfound=0;   // how many events found
    double eventsfound_prob=0;   // how many events found,  probabilities
    double sum[3]; // sum of weight for events found for 3 flavors
    static const int NBINS=10; // keep track of the number of events found,  binned by weights
    double MIN_LOGWEIGHT=-3;
    double MAX_LOGWEIGHT=-1;
    int index_weights=0; // which bin the weight falls into
    double eventsfound_binned[NBINS];
    double eventsfound_binned_e[NBINS];
    double eventsfound_binned_mu[NBINS];
    double eventsfound_binned_tau[NBINS];
    double km3sr = 0;                      // total km3sr
    double km3sr_e=0;                       // to calculate km3sr for electrons
    double km3sr_mu=0;                       // to calculate km3sr for muons
    double km3sr_tau=0;                       // to calculate km3sr for taus
    double error_plus=0;          // keeping track of asymmetric error bars
    double error_e_plus=0;
    double error_mu_plus=0;
    double error_tau_plus=0;
    double error_minus=0;
    double error_e_minus=0;
    double error_mu_minus=0;
    double error_tau_minus=0;
    int ierr=0; // integer for returning error from functions.
    double gain_dipole=2.15;  // antenna gain (nominal value for dipole is 2.15)
    double changle_deg=0; // same,  in degrees.

    // inputs
    int NNU;        // number of neutrinos
    // int whichray=0; // indexes the rays that we look at (actually just used for ice but we use it in GetNuFlavor so keep it here)
    double RANDOMISEPOL=0.;

    
    double volume_thishorizon; // for plotting volume within the horizon of the balloon
    // int realtime_this;  // for plotting real unix time
    // double longitude_this; // for plotting longitude
    // double latitude_this; // for plotting latitude
    // double altitude_this; // for plotting altitude
    // double heading_this=0.;// for plotting heading
    // double gps_offset=0;
    double pnu;   ///< energy of neutrinos

    double MEANX=0;
    double MEANY=0.;
    double SIGNALRADIUS=2.; // in degrees

    // frequency binning
    const double FREQ_LOW_DISCONES=120.E6; // min frequency for discones
    const double FREQ_HIGH_DISCONES=1000.E6; // max frequency for discones

    int passes_thisevent=0; // this event passes
    int unmasked_thisevent=0; // this event is unmasked


    int NDISCONES=8;
    double heff_discone=0; // effective height of a discone antenna
    double polarfactor_discone=0.;// factor to account for antenna beam pattern.
    double thislambda=0;// for finding wavelength at each frequency
    // double volts_discone=0.;// for finding voltage at each discone
    double vnoise_discone=0.; // noise on each discone


    double BW_DISCONES=300.E6-120.E6; // bandwidth of the discones

    // ray tracing
    double fresnel1=0;  ///< net fresnel factor on field at ice-firn interface
    double fresnel1_eachboresight[icemc::Anita::NLAYERS_MAX][icemc::Anita::NPHI_MAX];  // for slac simulation
    double fresnel2=0;  ///< net fresnel factor on field at firn-air
    double mag1=0;  // magnification factor on field at ice-firn interface
    double mag1_eachboresight[icemc::Anita::NLAYERS_MAX][icemc::Anita::NPHI_MAX];// for slac simulation
    double mag2=0;  // magnification factor on field in firn-air
    double t_coeff_pokey, t_coeff_slappy;
    double rflength=0;  // distance from interaction to ice-air exit point

    double diffexit=0;  // checking exit point between_MAX interations
    double diff_3tries=0;
    double diffnorm=0;  // checking angle of surf normal between iterations
    double diffrefr=0;  // checking angle of refr between iterations
    // double costheta_inc=0;  // cos angle of incidence wrt surface normal
    // double costheta_exit=0; // theta of exit point wrt earth (costheta=1 at south pole)
    double theta_rf_atbn; // polar angle of the signal as seen by perfect eyes at the balloon.
    double theta_rf_atbn_measured; //polar angle of the signal as measured at the balloon (just the above variable smeared by 0.5 degrees)

    double costhetanu=-1000; // costheta of neutrino direction wrt earth (costheta=1 at south pole)

    // neutrino path
    // double theta_in=0; // theta where neutrino enters earth (radians, south pole=0)
    // double lat_in=0; // latitude where neutrino enters earth (degrees, south pole=-90)
    // double nearthlayers=0; // how many layers (core, mantle, crust) does nnu traverse
    // double weight_prob=0.;//event weight,  including probability it interacts somewhere in ice along its path
    // double weight1=0; // event weight,  just based on absorption in earth,  see note
    // double weight=0.; // total event weight (product of the previous 2)
    // double logweight=0.;// log of the previous number
    // double len_int=0;// interaction length in m
    // double pieceofkm2sr=0; // Use this for making plots comparing different cross sections.  The integral of a plot from a run will be the total Area*sr of the detector.  That way it is proportional to the cross section and the integral is something meaningful to people.
    NeutrinoPath* fNeutrinoPath;

    double CUTONWEIGHTS=1.E-10; // cut out events with small enough weight that they don't matter,  to save time

    // counting variables
    int count_inthisloop1=0; // for debugging
    int count_inthisloop2=0;
    int count_inthisloop3=0;
    double averaging_thetas1=0; // for debugging
    double averaging_thetas2=0;
    double averaging_thetas3=0;
    int count_total=0; // number of neutrinos looped over.  should always equal inu

    // i.e.,  is theoretically visible by balloon

    int count_asktrigger=0;
    int count_asktrigger_nfb=0;

    double count_passestrigger_w=0; // same as above,  but sum weights
    int passestrigger=0; // 1=this event passes trigger, 0=does not
    // so that we know to increment count_passestrigger_w at the end of the event
    int allcuts[2]={0, 0}; // index is which ray (upward or downward)
    // 1=this ray for this event passes all cuts,  0=does not
    double allcuts_weighted[2]={0, 0}; // same as above but weighted

    //signal has a chance to pass after accounting for 1/r
    int count_chanceofsurviving=0; // based on neutrino direction,  has a chance of making it through earth.

    int count_chanceinhell0=0; // based on 1/r,  and best case attenuation
    // (after 2nd guess at rf exit point)
    // signal has a chance of passing

    double count_chanceinhell2_w=0; //same as above,  but sum weights
    int chanceinhell2=0; // 1=this event has chance to pass, 0=does not-
    // so that we know to increment count_chanceinhell2_w by weight1 at the end.

    int count_chordgoodlength=0; // Incremented if neutrino path through earth is more than 1m
    int count_d2goodlength=0; // neutrino path through ice is more than 1m
    // int count_rx=0; // counting antennas we loop through them

    double sum_frac[3] = {0}; // fraction of passing events that are e, mu, tau adding weights
    double sum_frac_db[3] = {0}; // same for double bangs

    static const int NBINS_DISTANCE=28; // keep track of number that pass as a function of distance to make Peter's plot
    double eventsfound_binned_distance[NBINS_DISTANCE] = {0.};  // binning cumulative events found vs. distance
    int index_distance=0; // index for filling array above
    double km3sr_distance[NBINS_DISTANCE] = {0.}; // result of conversion of events to sensitivity
    double error_distance_plus[NBINS_DISTANCE] = {0.}; // errors on above
    double error_distance_minus[NBINS_DISTANCE] = {0.};
    int eventsfound_binned_distance_forerror[NBINS_DISTANCE][NBINS] = {{0}}; // for propagation of errors

    //taus
    double km3sr_db=0;
    double km3sr_nfb=0;
    double ptau=0;
    int count_passestrigger_nfb=0;
    double percent_increase_db=0;
    double percent_increase_nfb=0;
    double percent_increase_total=0;
    double error_nfb=0;
    double error_km3sr_nfb=0;
    double error_percent_increase_nfb=0;

    TVector3 n_exit2bn_db[5];
    TVector3 nrf_iceside_db[5];  // direction of rf [tries][3d]
    double n_exit_phi;  //phi angle of the ray from the surface to the balloon
    int count_dbexitsice=0;

    // int count_interacts_in_earth=0;
    double eventsfound_nfb_binned[NBINS]; // counting events without first bang

    // rf parameters
    double heff_max=0.62639; // maximum value of the effective height based on antenna specs

    // event geometry
    // double scalefactor_distance=0; // 1/r scalefactor
    double scalefactor_attenuation=0; //scalefactor due to attenuation in ice
    double MAX_ATTENLENGTH=1671;
    // double maxtaper=0; // this is just for plotting - maximum you are ever off cerenkov cone while
    //an event is detectable
    double dviewangle_deg=0; ///< deviation from the cherenkov angle

    double forseckel[NVIEWANGLE][Anita::NFREQ];// Per Seckel's request,  get strength of signal across frequencies for different viewing angles.
    double viewangles[NVIEWANGLE];

    //Input files
    int err=0; // errors from GetDirection function

    double djunk; // junk variable





    // double vmmhz1m_visible = 0; //Actual V/m/Mhz at 1m
    int freq_bins = Anita::NFREQ; //Because the compiler objected to using the const directly
    
    double total_kgm2 = 0; // output of Getchord
    double nnu_array[3];
    double r_in_array[3];
    double nsurf_rfexit_array[3];
    double nsurf_rfexit_db_array[3];
    double r_bn_array[3];
    double n_bn_array[3];
    double posnu_array[3];
    double nrf_iceside_array[5][3];
    double nrf_iceside_db_array[5][3];
    double ant_max_normal0_array[3];
    double ant_max_normal1_array[3];
    double ant_max_normal2_array[3];
    double n_pol_array[3];
    double n_exit2bn_array[5][3];
    double r_enterice_array[3];
    double n_exit2bn_db_array[5][3];
    double rfexit_array[5][3];
    double rfexit_db_array[5][3];

    int times_crust_entered_det=0;  //Counter for total times each Earth layer is entered for detected neutrinos only
    int times_mantle_entered_det=0;
    int times_core_entered_det=0;

    int crust_entered=0;
    int mantle_entered=0;
    int core_entered=0;

    int neutrinos_passing_all_cuts=0;
    double sum_weights=0;
    //End verification plot block

    int xsecParam_nutype = 0; // neutrino = 0, antineutrino = 1;
    CurrentType xsecParam_nuint  = CurrentType::Neutral;



    // ray tracing
    double viewangle=0;
    double viewangle_temp=0; // angle of ray in ice relative to neutrino direction
    double viewangle_deg=0; // viewing angle in degrees.
    double cosviewangle=0; // cosine of viewing angle
    double offaxis=0; // viewangle-changle,  for plotting
    double nsigma_offaxis=0;// offaxis,  relative to deltheta_had,  for plotting
    double theta_threshold=0; // maximum angle you can be away from the cerenkov angle and still have a chance of seeing the event.
    double theta_threshold_deg=0; // maximum angle you can be away from the cerenkov angle and still have a chance of seeing the event.

    double nsigma_em_threshold=0; //number of sigma away theta_threshold is.
    double nsigma_had_threshold=0;// just for plotting
    double slopeyangle=0; // angle between nominal surface normal and the surface normal after slopeyness is applied
    //-------------------
    int beyondhorizon = 0;  //Switch tells if neutrino interacts beyond the balloon's horizon. (0: inside horizon,  1: outside horizon)

    // frequency binning
    double vmmhz1m_max=0; // maximum V/m/MHz at 1m from Jaime (highest frequency)
    // double vmmhz_lowfreq=0.; // V/m/MHz after 1/r,  attenuation at the lowest freq.
    double bestcase_atten=0;// attenuation factor,  best case
    double vmmhz1m_fresneledonce=0; // above,  after fresnel factor applied for ice-air interface
    double vmmhz1m_fresneledtwice=0; // above,  after fresnel factor applied for firn

    // given the angle you are off the Cerenkov cone,  the fraction of the observed e field that comes from the em shower
    // double vmmhz_em[Anita::NFREQ];
    double vmmhz_min=0;   // minimum of the above array
    double vmmhz_max=0;                        // maximum of the above array
    double deltheta_em[Anita::NFREQ], deltheta_had[Anita::NFREQ];     // for ch angular distribution
    double deltheta_em_max, deltheta_had_max;     // maximum value of above array angular distribution
    double deltheta_em_mid2, deltheta_had_mid2;     // widths of cones for the mid2 band

    // shower properties
    // double emfrac, hadfrac, sumfrac;               // em and had fractions
    // int n_interactions=1;           // count number of interactions for this event,  including secondaries.
    // double emfrac_db, hadfrac_db; //db stands for "double bang"
    // int nuflavorint2;
    double costheta_nutraject2;
    double phi_nutraject2;
    double altitude_int2;
    int currentint2;
    double d12;
    double d22;
    double dtryingdirection2;
    double logchord2;
    double r_fromballoon2;
    double chord_kgm2_bestcase2;
    double chord_kgm2_ice2;
    double weight_bestcase2;
    double r_exit2bn2;
    double r_exit2bn_measured2;
    std::string taudecay;                   // tau decay type: e, m, h

    double elast_y=0;                   // inelasticity
    // double volts_rx_0=0;              // voltage on an individual antenna,  lc polarization
    // double volts_rx_1=0;              // voltage on an individual antenna,  rc polarization
    // double volts_rx_max=0; // max voltage seen on an antenna - just for debugging purposes
    // double volts_rx_ave=0; // ave voltage seen on an antenna,  among hit antennas
    // double volts_rx_sum=0; // ave voltage seen on an antenna,  among hit antennas

    // double volts_rx_max_highband; // max voltage seen on an antenna - just for debugging purposes
    // double volts_rx_max_lowband; // max voltage seen on an antenna - just for debugging purposes
    // double volts_rx_rfcm_lab_e_all[48][512];
    // double volts_rx_rfcm_lab_h_all[48][512];



    double chengji = 0;


    double sigma = 0;                       // for cross section
    double len_int_kgm2=0;              // interaction length in kg/m^2
    double eventsfound_db=0; // same,  for double bang
    double eventsfound_nfb=0; // for taus

    // positions in earth frame
    double horizcoord=0; // x component of interaction position
    double vertcoord=0; // y component of interaction position

    double dist_int_bn_2d=0; // 2-d (surface) distance between interaction and balloon
    double dist_int_bn_2d_chord=0; //chord distance between interaction and balloon (between surface points)

    double viewangle_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // viewing angle for each antenna

    double cosalpha; // angle between nu momentum and surface normal at earth entrance point
    double mytheta; ///< alpha minus 90 degrees
    double cosbeta0; // angle between nu momentum and surface normal at interaction point.
    double mybeta; ///< beta minus 90 degrees
    double nuexitlength=0; // distance from interaction to where neutrino would leave
    double nuexitice=0;
    double nuentrancelength=0; // for taus
    double taulength=0;  // distance tau travels in ice before decaying
    double icethickness=0; // for taus
    double theta_pol_measured; // theta of the polarization as measured at the payload (for now we don't correct for the 10 degree cant)

    double ptaui=0;
    double ptauf=0;
    double tauweight=0;
    double nutauweightsum=0;
    double tauweightsum=0;
    double nutauweight=0;
    int tautrigger=0;
    int tauweighttrigger=0;


    double sourceLon;
    double sourceAlt;
    double sourceLat;
    double sourceMag;

    // TVector3 n_nutraject_ontheground; //direction of the neutrino from the person standing on the ground just below the balloon.
    // TVector3 n_pol; // direction of polarization
    // TVector3 n_pol_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // direction of polarization of signal seen at each antenna
    // TVector3 n_pol_db; // same,  double bangs

    // int l3trig[Anita::NPOL];  // 16 bit number which says which phi sectors pass L3 V-POL
    // // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs
    // int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];
    // //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs
    // int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];



    double icethck;
    double elev;
    double lon_ground;
    double lat_ground;
    double lon_ice;
    double lat_ice;
    double h20_depth;
    double lon_water;
    double lat_water;

    int pdgcode;

    double rms_rfcm_e;
    double rms_rfcm_h;
    double rms_lab_e;
    double rms_lab_h;

    double avgfreq_rfcm[Anita::NFREQ];
    double avgfreq_rfcm_lab[Anita::NFREQ];
    double freq[Anita::NFREQ];

    UInt_t eventNumber;






    
    // @todo constify... if I can const this, then we're probably near the end of the refactor...
    void applyRoughness(const Settings& settings1, const int& inu, Interaction* interaction1,  RayTracer* ray1, Screen* panel1, Antarctica* antarctica1, Balloon* bn1, const AskaryanFreqsGenerator* askFreqGen, Anita* anita1, const ShowerProperties& showerProps);


    
    void   GetSmearedIncidentAngle(TVector3 &specular, TVector3 &nrf_iceside, TVector3 &n_exit2bn, double SMEARINCIDENTANGLE) const;
    void   GetAir(double *col1) const;
    double GetThisAirColumn(const Settings*,  Geoid::Position r_in,  TVector3 nnu, Geoid::Position posnu,  double *col1,  double& cosalpha, double& mytheta,  double& cosbeta0, double& mybeta) const;
    double IsItDoubleBang(double exitlength,  double plepton) const;
    int WhereIsSecondBang(const Geoid::Position& posnu,  const TVector3& nnu,  double nuexitlength,  double pnu,  Antarctica *antarctica1,
			  const Geoid::Position& r_bn, Geoid::Position &posnu2,  Geoid::Position &rfexit_db,  TVector3 &n_exit2bn_db) const;

    TVector3 GetPolarization(const TVector3 &nnu,  const TVector3 &nrf2_iceside, int inu) const;
    void Attenuate(const Antarctica *antartica1, const Settings *settings1,  double& vmmhz_max,  double rflength,  const Geoid::Position &posnu) const ;
    void Attenuate_down(Antarctica *antarctica1,  const Settings *settings1,  double& vmmhz_max,  const Geoid::Position &rfexit2,  const Geoid::Position &posnu,  const Geoid::Position &posnu_down) const ;
    void IsAbsorbed(double chord_kgm2,  double len_int_kgm2,  double& weight) const;
    int GetRayIceSide(const TVector3 &n_exit2rx,  const TVector3 &nsurf_rfexit,  double nexit,  double nenter,  TVector3 &nrf2_iceside) const;

    // @todo constify... needs some love to constify
    int GetDirection(const Settings *settings1,  Interaction *interaction1,  const TVector3 &refr,  double deltheta_em,  double deltheta_had,  const ShowerProperties& sp,  double vmmhz1m_max,  double r_fromballoon,  RayTracer *ray1,  const AskaryanFreqsGenerator* askFreqGen,  Geoid::Position posnu,  Anita *anita1,  Balloon *bn1,  TVector3 &nnu,  double& costhetanu,  double& theta_threshold) ;
    void GetFresnel(Roughness *rough1,  int ROUGHNESS_SETTING,  const TVector3 &nsurf_rfexit,  const TVector3 &n_exit2rx,  TVector3 &n_pol,  const TVector3 &nrf2_iceside,  double efield, double deltheta_em, double deltheta_had,  double &t_coeff_pokey,  double &t_coeff_slappy,  double &fresnel,  double &mag) const;
    // void GetFresnel(Roughness *rough1,  int ROUGHNESS_SETTING,  const TVector3 &nsurf_rfexit,  const TVector3 &n_exit2rx,  TVector3 &n_pol,  const TVector3 &nrf2_iceside,  double efield,  const ShowerProperties& ,  double deltheta_em, double deltheta_had,  double &t_coeff_pokey,  double &t_coeff_slappy,  double &fresnel,  double &mag) const;    
    double GetViewAngle(const TVector3 &nrf2_iceside,  const TVector3 &nnu) const;
    int TIR(const TVector3 &n_surf,  const TVector3 &nrf2_iceside,  double N_IN,  double N_OUT) const;
    // void IntegrateBands(Anita *anita1,  int k,  Screen *panel1,  double *freq,  double scalefactor,  double *sumsignal) const;

    // @todo constify... needs some love to constify
    void Summarize(const Settings *settings1,  Anita* anita1,  Counting *count1,  Spectra *spectra1, const AskaryanFreqsGenerator* askFreqGen,  Primaries *primary1,  double,  double eventsfound,  double,  double,  double,  double*,  double,  double,  double&,  double&,  double&,  double&, TString);
    void WriteNeutrinoInfo(const int& inu, const Geoid::Position&,  const TVector3&,  const Geoid::Position&,  double,  NuFlavor,  CurrentType,  double,  std::ofstream &nu_out) const;

    /** 
     * @brief Run the neutrino generation
     * 
     * This function does the stuff that used to be the main in the icemc executable
     * @todo needs some love to constify... probably not doable...
     */    
    void generateNeutrinos(const Settings& settings1, const CommandLineOpts& clOpts);


    //do a threshold scan
    double threshold_start=-1.;
    double threshold_end=-6.;
    static const int NTHRESHOLDS=20;
    double threshold_step=(threshold_end-threshold_start)/(double)NTHRESHOLDS;

    double npass_v_thresh[NTHRESHOLDS]={0.};
    double denom_v_thresh[NTHRESHOLDS]={0.};
    double npass_h_thresh[NTHRESHOLDS]={0.};
    double denom_h_thresh[NTHRESHOLDS]={0.};
    double thresholds[NTHRESHOLDS];

    Interaction* interaction1 = nullptr;
    // Balloon* bn1;
    // Anita* anita1;
    ANITA* fDetector = nullptr;
    Taumodel* fTauPtr = nullptr;

    GeneratedNeutrino* fGenNu = nullptr;
    PassingNeutrino* fPassNu = nullptr;

  private:

    

  };
}



#endif //ICEMC_EVENT_GENERATOR_H
