#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <ctype.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <time.h>
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"
#include "signal.h"

//#include "rx.hpp"
#include "Constants.h"
#include "Settings.h"
#include "position.hh"

#include "earthmodel.hh"
#include "Tools.h"
#include "vector.hh"
#include "roughness.hh"
#include "anita.hh"
#include "balloon.hh"
#include "icemodel.hh"
#include "trigger.hh"
#include "Spectra.h"
#include "signal.hh"
#include "secondaries.hh"
#include "ray.hh"
#include "counting.hh"
#include "Primaries.h"
#include "Taumodel.hh"
#include "screen.hh"

#include <string>
#include <sstream>

#if __cplusplus > 199711L
	#include <type_traits>
#endif

#include <typeinfo>

//#define ANITA3_EVENTREADER

#ifdef ANITA_UTIL_EXISTS
#include "UsefulAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "RawAnitaHeader.h"
#include "Adu5Pat.h"
#include "FFTtools.h"
UsefulAnitaEvent*     realEvPtr    = NULL;
RawAnitaHeader*       rawHeaderPtr = NULL;
Adu5Pat* 	      Adu5PatPtr   = NULL;
#ifdef ANITA3_EVENTREADER
#include "TruthAnitaEvent.h"
TruthAnitaEvent*      truthEvPtr   = NULL;
#endif
#endif

Taumodel* TauPtr = NULL;

ClassImp(RX);

using namespace std;

class EarthModel;
class Position;

/************MOVED FROM shared.hh and shared.cc*****************/
// These need to be moved elsewhere.
// They were all global variables inside of shared.hh.
// Most of these should not be global variables.
const int NVIEWANGLE=100; // number of viewing angles to look at the signal,  on Seckel's request
// int irays; // counts rays (for roughness) // LC: commented out because not used

int inu; // counts neutrinos as they are generated
double eventsfound_beforetrigger=0.;
double eventsfound_crust=0; //number of events that only traverse the crust
double eventsfound_weightgt01=0; // summing weights > 0.1
double eventsfound_belowhorizon=0; // how many are below horizon
double eventsfound=0;   // how many events found
double eventsfound_prob=0;   // how many events found,  probabilities
double sum[3]; // sum of weight for events found for 3 flavors
// These numbers are from Feldman and Cousins wonderful paper:  physics/9711021
double poissonerror_minus[21] = {0.-0.00, 1.-0.37, 2.-0.74, 3.-1.10, 4.-2.34, 5.-2.75, 6.-3.82, 7.-4.25, 8.-5.30, 9.-6.33, 10.-6.78, 11.-7.81, 12.-8.83, 13.-9.28, 14.-10.30, 15.-11.32, 16.-12.33, 17.-12.79, 18.-13.81, 19.-14.82, 20.-15.83};
double poissonerror_plus[21] = {1.29-0., 2.75-1., 4.25-2., 5.30-3., 6.78-4., 7.81-5., 9.28-6., 10.30-7., 11.32-8., 12.79-9., 13.81-10., 14.82-11., 16.29-12., 17.30-13., 18.32-14., 19.32-15., 20.80-16., 21.81-17., 22.82-18., 23.82-19., 25.30-20.};
const int NBINS=10; // keep track of the number of events found,  binned
// by weights
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
int whichray=0; // indexes the rays that we look at (actually just used for ice but we use it in GetNuFlavor so keep it here)
double RANDOMISEPOL=0.;
/************MOVED FROM shared.hh and shared.cc*****************/


double volume_thishorizon; // for plotting volume within the horizon of the balloon
int realtime_this;  // for plotting real unix time
double longitude_this; // for plotting longitude
double latitude_this; // for plotting latitude
double altitude_this; // for plotting altitude
double heading_this=0.;// for plotting heading
double gps_offset=0;

// inputs

double pnu=pow(10., 20);   //!< energy of neutrinos

double MEANX=0;
double MEANY=0.;

double SIGNALRADIUS=2.; // in degrees

// frequency binning
const double FREQ_LOW_DISCONES=120.E6; // min frequency for discones
const double FREQ_HIGH_DISCONES=1000.E6; // max frequency for discones

double bwslice_vnoise_thislayer[4];// for filling tree6b,  noise for each bandwidth on each layer
int passes_thisevent=0; // this event passes
int unmasked_thisevent=0; // this event is unmasked

int discones_passing;  // number of discones that pass
int NDISCONES=8;
double heff_discone=0; // effective height of a discone antenna
double polarfactor_discone=0.;// factor to account for antenna beam pattern.
double thislambda=0;// for finding wavelength at each frequency
double volts_discone=0.;// for finding voltage at each discone
double vnoise_discone=0.; // noise on each discone


double BW_DISCONES=300.E6-120.E6; // bandwidth of the discones

// ray tracing
double fresnel1=0;  //!< net fresnel factor on field at ice-firn interface
double fresnel1_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX];  // for slac simulation
double fresnel2=0;  //!< net fresnel factor on field at firn-air
double mag1=0;  // magnification factor on field at ice-firn interface
double mag1_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX];// for slac simulation
double mag2=0;  // magnification factor on field in firn-air
double t_coeff_pokey, t_coeff_slappy;
double rflength=0;  // distance from interaction to ice-air exit point

double e_comp_max1=0;
double h_comp_max1=0;
double e_comp_max2=0;
double h_comp_max2=0;
double e_comp_max3=0;
double h_comp_max3=0;

double diffexit=0;  // checking exit point between_MAX interations
double diff_3tries=0;
double diffnorm=0;  // checking angle of surf normal between iterations
double diffrefr=0;  // checking angle of refr between iterations
double costheta_inc=0;  // cos angle of incidence wrt surface normal
double costheta_exit=0; // theta of exit point wrt earth (costheta=1 at south pole)
double theta_rf_atbn; // polar angle of the signal as seen by perfect eyes at the balloon.
double theta_rf_atbn_measured; //polar angle of the signal as measured at the balloon (just the above variable smeared by 0.5 degrees)

double costhetanu=-1000; // costheta of neutrino direction wrt earth (costheta=1 at south pole)

// neutrino path
double theta_in=0; // theta where neutrino enters earth (radians, south pole=0)
double lat_in=0; // latitude where neutrino enters earth (degrees, south pole=-90)
double nearthlayers=0; // how many layers (core, mantle, crust) does nnu traverse
double weight_prob=0.;//event weight,  including probability it interacts somewhere in ice along its path
double weight1=0; // event weight,  just based on absorption in earth,  see note
double weight=0.; // total event weight (product of the previous 2)
double logweight=0.;// log of the previous number
double len_int=0;// interaction length in m
double pieceofkm2sr=0; // Use this for making plots comparing different cross sections.  The integral of a plot from a run will be the total Area*sr of the detector.  That way it is proportional to the cross section and the integral is something meaningful to people.

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
int count_pass=0;  // how many total trigger channels pass (4 bandwidth slices*2 pol * nrx)

double count_passestrigger_w=0; // same as above,  but sum weights
int passestrigger=0; // 1=this event passes trigger, 0=does not
// so that we know to increment count_passestrigger_w at the end of the event
int allcuts[2]={0, 0}; // index is which ray (upward or downward)
// 1=this ray for this event passes all cuts,  0=does not
double allcuts_weighted[2]={0, 0}; // same as above but weighted
double allcuts_weighted_polarization[3]={0, 0, 0}; // same as above but divided into [0] vpol, [1] hpol, [2] both vpol and hpol


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
int count_rx=0; // counting antennas we loop through them

double sum_frac[3]; // fraction of passing events that are e, mu, tau adding weights
double sum_frac_db[3]; // same for double bangs

const int NBINS_DISTANCE=28; // keep track of number that pass as a function of distance to make Peter's plot
double eventsfound_binned_distance[NBINS_DISTANCE] = {0.};  // binning cumulative events found vs. distance
int index_distance=0; // index for filling array above
double km3sr_distance[NBINS_DISTANCE] = {0.}; // result of conversion of events to sensitivity
double error_distance_plus[NBINS_DISTANCE] = {0.}; // errors on above
double error_distance_minus[NBINS_DISTANCE] = {0.};
int eventsfound_binned_distance_forerror[NBINS_DISTANCE][NBINS] = {{0}}; // for propagation of errors



//taus
double km3sr_db = 0;
double km3sr_nfb=0;
double ptau=0;
int count_passestrigger_nfb=0;
double percent_increase_db=0;
double percent_increase_nfb=0;
double percent_increase_total=0;
double error_nfb=0;
double error_km3sr_nfb=0;
double error_percent_increase_nfb=0;

Vector n_exit2bn_db[5];
Vector nrf_iceside_db[5];  // direction of rf [tries][3d]
double n_exit_phi;  //phi angle of the ray from the surface to the balloon
int count_dbexitsice=0;

// int count_interacts_in_earth=0;
double eventsfound_nfb_binned[NBINS]; // counting events without first bang

// rf parameters
double heff_max=0.62639; // maximum value of the effective height based on antenna specs

// event geometry
double scalefactor_distance=0; // 1/r scalefactor
double scalefactor_attenuation=0; //scalefactor due to attenuation in ice
double MAX_ATTENLENGTH=1671;
double maxtaper=0; // this is just for plotting - maximum you are ever off cerenkov cone while
//an event is detectable
double dviewangle_deg=0; //!< deviation from the cherenkov angle

double forseckel[NVIEWANGLE][Anita::NFREQ];// Per Seckel's request,  get strength of signal across frequencies for different viewing angles.
double viewangles[NVIEWANGLE];
double GetAirDistance(double altitude_bn,  double beta); // given beta=angle wrt horizontal that the ray hits the balloon,  calculate distance that the ray traveled in air,  including curvature of earth

//Input files
int err=0; // errors from GetDirection function

double djunk; // junk variable

//For verification plots - added by Stephen
int max_antenna0 = -1;  //antenna with the peak voltage,  top layer
int max_antenna1 = -1;  //antenna with the peak voltage,  middle layer
int max_antenna2 = -1;  //antenna with the peak voltage,  bottom layer
double max_antenna_volts0 = 0; //Voltage on the antenna with maximum signal,  top layer
double max_antenna_volts0_em = 0; //Component of voltage from em shower on the antenna with maximum signal,  top layer

double max_antenna_volts1 = 0; //Voltage on the antenna with maximum signal,  middle layer
double max_antenna_volts2 = 0; //Voltage on the antenna with maximum signal,  bottom layer

double rx0_signal_eachband[2][5];
double rx0_threshold_eachband[2][5];
double rx0_noise_eachband[2][5];
int rx0_passes_eachband[2][5];

double voltagearray[Anita::NLAYERS_MAX*Anita::NPHI_MAX]; //Records max voltages on each antenna for one neutrino

Vector ant_max_normal0; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  top layer
Vector ant_max_normal1; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  middle layer
Vector ant_max_normal2; //Vector normal to the face of the antenna with the maximum signal for a single neutrino,  bottom layer
double vmmhz1m_visible = 0; //Actual V/m/Mhz at 1m
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


// functions

// set up array of viewing angles for making plots for seckel
void SetupViewangles(Signal *sig1);

void GetAir(double *col1); // get air column as a function of theta- only important for black hole studies
double GetThisAirColumn(Settings*,  Position r_in,  Vector nnu, Position posnu,  double *col1,  double& cosalpha, double& mytheta,  double& cosbeta0, double& mybeta);

double ScaleVmMHz(double vmmhz1m_max,  const Position &posnu,  const Position &r_bn);


double IsItDoubleBang(double exitlength,  double plepton);
int WhereIsSecondBang(const Position& posnu,  const Vector& nnu,  double nuexitlength,  double pnu,  IceModel *antarctica1,  const Position& r_bn, Position &posnu2,  Position &rfexit_db,  Vector &n_exit2bn_db);
void GetCurrent(string& current);


double GetAverageVoltageFromAntennasHit(Settings *settings1,  int *nchannels_perrx_triggered,  double *voltagearray,  double& volts_rx_sum);


Vector GetPolarization(const Vector &nnu,  const Vector &nrf2_iceside);

void Attenuate(IceModel *antartica1, Settings *settings1,  double& vmmhz_max,  double rflength,  const Position &posnu);
void Attenuate_down(IceModel *antarctica1,  Settings *settings1,  double& vmmhz_max,  const Position &rfexit2,  const Position &posnu,  const Position &posnu_down);
void IsAbsorbed(double chord_kgm2,  double len_int_kgm2,  double& weight);

void GetBalloonLocation(Interaction *interaction1,Ray *ray1,Balloon *bn1,IceModel *antarctica);

int GetRayIceSide(const Vector &n_exit2rx,  const Vector &nsurf_rfexit,  double nexit,  double nenter,  Vector &nrf2_iceside);

int GetDirection(Settings *settings1,  Interaction *interaction1,  const Vector &refr,  double deltheta_em,  double deltheta_had,  double emfrac,  double hadfrac,  double vmmhz1m_max,  double r_fromballoon,  Ray *ray1,  Signal *sig1,  Position posnu,  Anita *anita1,  Balloon *bn1,  Vector &nnu,  double& costhetanu,  double& theta_threshold);

// int GetDirectionRough(Settings *settings1, const Vector &refr,  double deltheta_em,  double deltheta_had, 
// 					  double emfrac,  double hadfrac,  double pnu,  double vmmhz1m_max, 
// 					  double r_fromballoon,  Ray *ray1, 
// 					  Signal *sig1,  Position posnu,  Anita *anita1,  Balloon *bn1, 
// 					  Vector &nnu,  double& costhetanu,  double& theta_threshold);

void GetFresnel(Roughness *rough1,  int ROUGHNESS_SETTING,  const Vector &nsurf_rfexit,  const Vector &n_exit2rx,  Vector &n_pol,  const Vector &nrf2_iceside,  double efield,  double emfrac,  double hadfrac,  double deltheta_em, double deltheta_had,  double &t_coeff_pokey,  double &t_coeff_slappy,  double &fresnel,  double &mag);

double GetViewAngle(const Vector &nrf2_iceside,  const Vector &nnu);
int TIR(const Vector &n_surf,  const Vector &nrf2_iceside,  double N_IN,  double N_OUT);

void IntegrateBands(Anita *anita1,  int k,  Screen *panel1,  double *freq,  double scalefactor,  double *sumsignal);
void Integrate(Anita *anita1,  int j,  int k,  double *vmmhz,  double *freq,  double scalefactor,  double sumsignal);
void interrupt_signal_handler(int);	// This catches the Control-C interrupt,  SIGINT
bool ABORT_EARLY = false;		// This flag is set to true when interrupt_signal_handler() is called

void Summarize(Settings *settings1,  Anita* anita1,  Counting *count1,  Spectra *spectra1, Signal *sig1,  Primaries *primary1,  double,  double eventsfound,  double,  double,  double,  double*,  double,  double,  double&,  double&,  double&,  double&,  ofstream&,  ofstream&);
    
void WriteNeutrinoInfo(Position&,  Vector&,  Position&,  double,  string,  string,  double,  ofstream &nu_out);
void CloseTFile(TFile *hfile);

int Getmine(double*,  double*,  double*,  double*);
void Getearth(double*,  double*,  double*,  double*);


#ifdef ANITA_UTIL_EXISTS
int GetIceMCAntfromUsefulEventAnt(Settings *settings1,  int UsefulEventAnt);
#ifdef R_EARTH
#undef R_EARTH
#endif
#endif



//do a threshold scan
double threshold_start=-1.;
double threshold_end=-6.;
const int NTHRESHOLDS=20;
double threshold_step=(threshold_end-threshold_start)/(double)NTHRESHOLDS;

double npass_v_thresh[NTHRESHOLDS]={0.};
double denom_v_thresh[NTHRESHOLDS]={0.};
double npass_h_thresh[NTHRESHOLDS]={0.};
double denom_h_thresh[NTHRESHOLDS]={0.};
double thresholds[NTHRESHOLDS];

int main(int argc,  char **argv) {
  //--------------------------------------------------------------
  //  MC Anita
  //
  // 12/01/03
  //
  //--------------------------------------------------------------


  // for comparing with peter
  double sumsignal[5]={0.};
  double sumsignal_aftertaper[5]={0.};
    
  string stemp;
  
  Settings* settings1 = new Settings();
  
  string input="inputs.txt";
  string run_num;//current run number as string
  int run_no = 0;//current run number as integer

  ifstream inputsfile(input.c_str());
  if (!inputsfile.good()) {
    cout << "Can't find inputs file.\n";
    exit(1);
  }
  
  if( (argc%3!=1)&&(argc%2!=1)) {
    cout << "Syntax for options: -i inputfile -o outputdir -r run_number\n";
    return 0;
  }
  int nnu_tmp=0;
  double trig_thresh=0.;
  char clswitch; // command line switch
  if (argc>1) {
    while ((clswitch = getopt(argc, argv, "t:i:o:r:n:")) != EOF) {
      switch(clswitch) {
      case 'n':
	nnu_tmp=atoi(optarg);
        break;
      case 't':
	trig_thresh=atof(optarg);
        break;
      case 'i':
        input=optarg;
        cout << "Changed input file to: " << input << endl;
        inputsfile.close();
        inputsfile.open(input.c_str());
        break;
      case 'o':
        settings1->outputdir=optarg;
        cout << "Changed output directory to: " << settings1->outputdir << endl;
        stemp="mkdir " + settings1->outputdir;
        system(stemp.c_str());
        break;
      case 'r':
        run_num=optarg;
        stringstream convert(run_num);
        convert>>run_no;
        break;
      } // end switch
    } // end while
  } // end if arg>1
  

  settings1->SEED=settings1->SEED +run_no;
  cout <<"seed is " << settings1->SEED << endl;

  TRandom *rsave = gRandom;
  TRandom3 *Rand3 = new TRandom3(settings1->SEED);//for generating random numbers
  gRandom=Rand3;

  stemp=settings1->outputdir+"/nu_position"+run_num+".txt";
  ofstream nu_out(stemp.c_str(),  ios::app); //Positions,  direction of momentum,  and neutrino type for Ryan.
  
  stemp=settings1->outputdir+"/veff"+run_num+".txt";
  ofstream veff_out(stemp.c_str(),  ios::app);//to output only energy and effective volume to veff.txt
  
  stemp=settings1->outputdir+"/distance"+run_num+".txt";
  ofstream distanceout(stemp.c_str());
  
  stemp=settings1->outputdir+"/debug"+run_num+".txt";
  fstream outfile(stemp.c_str(), ios::out);
  
  stemp=settings1->outputdir+"/forbrian"+run_num+".txt";
  fstream forbrian(stemp.c_str(), ios::out);
  
  stemp=settings1->outputdir+"/al_voltages_direct"+run_num+".dat";
  fstream al_voltages_direct(stemp.c_str(), ios::out); //added djg ------provide anita-lite voltages and noise from MC for anita-lite analysis
  
  stemp=settings1->outputdir+"/events"+run_num+".txt";
  ofstream eventsthatpassfile(stemp.c_str());
  
  stemp=settings1->outputdir+"/numbers"+run_num+".txt";
  ofstream fnumbers(stemp.c_str()); // debugging
  
  stemp=settings1->outputdir+"/output"+run_num+".txt";
  ofstream foutput(stemp.c_str(),  ios::app);
  
  stemp=settings1->outputdir+"/slacviewangles"+run_num+".dat";
  ofstream fslac_viewangles(stemp.c_str()); // this outputs numbers that we need for analyzing slac data
  
  stemp=settings1->outputdir+"/slac_hitangles"+run_num+".dat";
  ofstream fslac_hitangles(stemp.c_str()); // this outputs numbers that we need for analyzing slac data
    
  Balloon *bn1=new Balloon(); // instance of the balloon
  Anita *anita1=new Anita();// right now this constructor gets banding info
  Secondaries *sec1=new Secondaries();
  Primaries *primary1=new Primaries();
  Signal *sig1=new Signal();
  Ray *ray1=new Ray(); // create new instance of the ray class
  Counting *count1=new Counting();
  GlobalTrigger *globaltrig1;
  Taumodel *taus1 = new Taumodel();
  // input parameters
  settings1->ReadInputs(inputsfile,  foutput,  anita1,  sec1,  sig1,  bn1,  ray1);
  


  settings1->SEED=settings1->SEED + run_no;
  gRandom->SetSeed(settings1->SEED);

  bn1->InitializeBalloon();
  anita1->Initialize(settings1, foutput, inu);

  if (trig_thresh!=0)
    anita1->powerthreshold[4]=trig_thresh;
  if (nnu_tmp!=0)
    NNU=nnu_tmp;

  Spectra *spectra1 = new Spectra((int)settings1->EXPONENT);
  Interaction *interaction1=new Interaction("nu", primary1, settings1, 0, count1);
  Interaction *int_banana=new Interaction("banana", primary1, settings1, 0, count1);
  
  Roughness *rough1=new Roughness(settings1->ROUGHSIZE); // create new instance of the roughness class
  int fSCREEN_NUMPOINTS_EDGE = settings1->ROUGHNESS;
  Screen *panel1 = new Screen(fSCREEN_NUMPOINTS_EDGE);              // create new instance of the screen class

  // roughness testing
  //for (int ii=0; ii<90; ii++){
  //  std::cerr<<ii<<"  "<<rough1->ConvertTheta0AirGlass_to_GlassAir(ii)<<"  "<<rough1->ConvertTheta0GlassAir_to_AirGlass(ii)<<std::endl;
  //}
  //for(int ii=-90; ii<=90; ii+=10){
  //  for (int jj=-90; jj<=90; jj+=10){
  //    std::cerr<<ii<<"  "<<jj<<"  "<< rough1->InterpolatePowerValue(ii, jj)<<std::endl;
  //  }
  //}

  // screen testing
  //panel1->SetEdgeLength(10.);
  //panel1->SetCentralPoint(Position(.1,-.1,.1));
  //panel1->SetNormal(Vector(1,.1,.1));
  //Position startpos;
  //for (int ii=0; ii<fSCREEN_NUMPOINTS_EDGE*fSCREEN_NUMPOINTS_EDGE; ii++){
  //  startpos = panel1->GetNextPosition();
  //  std::cerr << startpos[0] << "  " << startpos[1] << "  " << startpos[2] << std::endl;
  //}

  // declare instance of trigger class.
  // this constructor reads from files with info to parameterize the
  // signal efficiency curves.
  // and reads in scalar vs. threshold scans
  
  time_t raw_start_time = time(NULL);
  struct tm * start_time = localtime(&raw_start_time);
  
  cout << "Date and time at start of run are: " << asctime (start_time) << "\n";
  
  if (settings1->FORSECKEL)
    SetupViewangles(sig1);// set up viewing angles for testing against jaime's parameterizations
  
  // for attenuation of neutrino in atmosphere
  // only important for black hole studies
  double col1[900];
  GetAir(col1);
  double myair;//air column density, kg/m^2
  
  // zeroing global variables.
  Tools::Zero(sum_frac, 3);
  Tools::Zero(sum_frac_db, 3);
  Tools::Zero(anita1->NRX_PHI, Anita::NLAYERS_MAX);
  for (int i=0;i<Anita::NLAYERS_MAX;i++) {
    Tools::Zero(anita1->PHI_EACHLAYER[i], Anita::NPHI_MAX);
  }
  Tools::Zero(anita1->PHI_OFFSET, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->THETA_ZENITH, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->LAYER_VPOSITION, Anita::NLAYERS_MAX);
  Tools::Zero(anita1->RRX, Anita::NLAYERS_MAX);
  
  
  //added djg ////////////////////////////////////////////////////////
  al_voltages_direct<<"antenna #"<<"   "<<"volts chan 1"<<"   "<<"volts chan 2"<<"    "<<"volts chan 3"<<"    "<<"volts chan 4"<<"    "<<"noise chan 1"<<"    "<<"noise chan 2"<<"    "<<"noise chan 3"<<"   "<<"noise chan 4"<<"  "<<"weight"<<endl;
  ////////////////////////////////////////////////////////////////////
  
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
  double vmmhz_lowfreq=0.; // V/m/MHz after 1/r,  attenuation at the lowest freq.
  double bestcase_atten=0;// attenuation factor,  best case
  double vmmhz1m_fresneledonce=0; // above,  after fresnel factor applied for ice-air interface
  double vmmhz1m_fresneledtwice=0; // above,  after fresnel factor applied for firn
  double vmmhz[Anita::NFREQ];                        //  V/m/MHz at balloon (after all steps)
  // given the angle you are off the Cerenkov cone,  the fraction of the observed e field that comes from the em shower
  double vmmhz_em[Anita::NFREQ];
  
  double vmmhz_temp=1.;
  double vmmhz_min_thatpasses=1000;
  double vmmhz_min=0;   // minimum of the above array
  double vmmhz_max=0;                        // maximum of the above array
  double deltheta_em[Anita::NFREQ], deltheta_had[Anita::NFREQ];     // for ch angular distribution
  double deltheta_em_max, deltheta_had_max;     // maximum value of above array angular distribution
  double deltheta_em_mid2, deltheta_had_mid2;     // widths of cones for the mid2 band
  
  // shower properties
  double emfrac, hadfrac, sumfrac;               // em and had fractions
  int n_interactions=1;            // count number of interactions for this event,  including secondaries.
  double emfrac_db, hadfrac_db;
  int nuflavorint2=interaction1->nuflavorint;
  double costheta_nutraject2=interaction1->costheta_nutraject;
  double phi_nutraject2=interaction1->phi_nutraject;
  double altitude_int2=interaction1->altitude_int;
  int currentint2=interaction1->currentint;
  double d12=interaction1->d1;
  double d22=interaction1->d2;
  double dtryingdirection2=interaction1->dtryingdirection;
  double logchord2=interaction1->logchord;
  double r_fromballoon2=interaction1->r_fromballoon[0];
  double chord_kgm2_bestcase2=interaction1->chord_kgm2_bestcase;
  double chord_kgm2_ice2=interaction1->chord_kgm2_ice;
  double weight_bestcase2=interaction1->weight_bestcase;
  double r_exit2bn2=interaction1->r_exit2bn;
  double r_exit2bn_measured2=interaction1->r_exit2bn_measured;

  string taudecay;                   // tau decay type: e, m, h
  
  double elast_y=0;                   // inelasticity
  double volts_rx_0=0;              // voltage on an individual antenna,  lc polarization
  double volts_rx_1=0;              // voltage on an individual antenna,  rc polarization
  double volts_rx_max=0; // max voltage seen on an antenna - just for debugging purposes
  double volts_rx_ave=0; // ave voltage seen on an antenna,  among hit antennas
  double volts_rx_sum=0; // ave voltage seen on an antenna,  among hit antennas
  
  double volts_rx_max_highband; // max voltage seen on an antenna - just for debugging purposes
  double volts_rx_max_lowband; // max voltage seen on an antenna - just for debugging purposes
  double volts_rx_rfcm_lab_e_all[48][512];
  double volts_rx_rfcm_lab_h_all[48][512];
  
 // variable declarations for functions GetEcompHcompEvector and GetEcompHcompkvector - oindree
  double e_component=0; // E comp along polarization
  double h_component=0; // H comp along polarization
  double n_component=0; // normal comp along polarization

  double e_component_kvector=0; // component of e-field along the rx e-plane
  double h_component_kvector=0; // component of the e-field along the rx h-plane
  double n_component_kvector=0; // component of the e-field along the normal 

 
  Vector n_eplane = const_z; 
  Vector n_hplane = -const_y;
  Vector n_normal = const_x;

  double chengji = 0;
  Vector ant_normal; //Vector normal to the face of the antenna
  
  double hitangle_e, hitangle_h;       // angle the ray hits the antenna wrt e-plane, h-plane
  double hitangle_e_all[Anita::NANTENNAS_MAX];         // hit angles rel. to e plane stored for each antenna
  double hitangle_h_all[Anita::NANTENNAS_MAX];         // hit angles rel. to h plane stored for each antenna
  
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
  double mytheta; //!< alpha minus 90 degrees
  double cosbeta0; // angle between nu momentum and surface normal at interaction point.
  double mybeta; //!< beta minus 90 degrees
  double nuexitlength=0; // distance from interaction to where neutrino would leave
  double nuexitice=0; 
  double nuentrancelength=0; // for taus
  double taulength=0;  // distance tau travels in ice before decaying
  double icethickness=0; // for taus
  double theta_pol_measured; // theta of the polarization as measured at the payload (for now we don't correct for the 10 degree cant)
  
  double ptaui=0;
  double ptauf =0;
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
 
  Vector n_nutraject_ontheground; //direction of the neutrino from the person standing on the ground just below the balloon.
  Vector n_pol; // direction of polarization
  Vector n_pol_eachboresight[Anita::NLAYERS_MAX][Anita::NPHI_MAX]; // direction of polarization of signal seen at each antenna
  Vector n_pol_db; // same,  double bangs
  
  int l3trig[Anita::NPOL];  // 16 bit number which says which phi sectors pass L3 V-POL
    // For each trigger layer,  which "clumps" pass L2.  16 bit,  16 bit and 8 bit for layers 1 & 2 and nadirs  
  int l2trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];
   //For each trigger layer,  which antennas pass L1.  16 bit,  16 bit and 8 bit and layers 1,  2 and nadirs  
  int l1trig[Anita::NPOL][Anita::NTRIGGERLAYERS_MAX];
  
  // these are declared here so that they can be stuck into trees
  int loctrig[Anita::NPOL][Anita::NLAYERS_MAX][Anita::NPHI_MAX]; //counting how many pass trigger requirement  
  
  int loctrig_nadironly[Anita::NPOL][Anita::NPHI_MAX]; //counting how many pass trigger requirement  

 
  UShort_t phiTrigMask;
  UShort_t phiTrigMaskH;
  UShort_t l1TrigMask;
  UShort_t l1TrigMaskH;
  
  int nchannels_triggered = 0; // total number of channels triggered
  int nchannels_perrx_triggered[48]; // total number of channels triggered
  
  // zeroing
  interaction1->dnutries=0;
  eventsfound=0.; // sums weights for events that pass
  
  Tools::Zero(count1->npass, 2); // sums events that pass,  without weights
  Tools::Zero(sum, 3);
  eventsfound_db=0;
  eventsfound_nfb=0;
  
  Tools::Zero(eventsfound_binned, NBINS);
  Tools::Zero(eventsfound_binned_e, NBINS);
  Tools::Zero(eventsfound_binned_mu, NBINS);
  Tools::Zero(eventsfound_binned_tau, NBINS);
  Tools::Zero(eventsfound_nfb_binned, NBINS);

  //we pick both the interaction point and its corresponding mirror point
  
  //for drawing the events on the events map
  TH2F *ref_int_coord=new TH2F("ref_int_coord",  "",  600,  -3000,  3000,  500,  -2500,  2500);
  ref_int_coord->SetMarkerSize(0.7);;
  ref_int_coord->SetMarkerStyle(30);
  ref_int_coord->SetMarkerColor(kBlue);

  TH2F *dir_int_coord=new TH2F("dir_int_coord",  "",  600,  -3000,  3000,  500,  -2500,  2500);
  dir_int_coord->SetMarkerSize(0.7);
  dir_int_coord->SetMarkerStyle(30);
  
  TH1F *prob_eachphi_bn=new TH1F("prob_eachphi_bn",  "prob_eachphi_bn",  100,  0.,  6.3);
  TH1F *prob_eachilon_bn=new TH1F("prob_eachilon_bn",  "prob_eachilon_bn",  180,  0.,  180.);
  TH2F *h6=new TH2F("theta_vs_hitangle_h",  "theta_vs_hitangle_h",  100,  -3.14,  3.14,  100,  -1.1,  1.1);
  TH1F *h10=new TH1F("hitangle_e",  "hitangle_e",  20,  -1.6,  1.6);
  TH1F *hy=new TH1F("hy",  "hy",  100,  0.,  1.);
  TH1F *fraction_sec_muons = new TH1F("fraction_sec_muons",  "fraction_sec_muons",  100,  0.,  1.);
  TH1F *fraction_sec_taus=new TH1F("fraction_sec_taus",  "fraction_sec_taus",  100,  0.,  1.);
  TH1F *n_sec_muons= new TH1F("n_sec_muons",  "n_sec_muons",  100,  0.,  10.);
  TH1F *n_sec_taus= new TH1F("n_sec_taus",  "n_sec_taus",  100,  0.,  10.);
  TH1F *sampleweights=new TH1F("sampleweights",  "sampleweights",  100,  -5.,  0.); // we sample weights for early events and save
  //them in this histogram,  
  
  //to determine where the cut should be.
  stemp=settings1->outputdir+"/icefinal"+run_num+".root";
  TFile *hfile = new TFile(stemp.c_str(), "RECREATE", "ice");
  TTree *tree2 = new TTree("h2000", "h2000"); // tree2 filled for each event that is beyond the horizon.
  
  tree2->Branch("inu", &inu, "inu/I");
  tree2->Branch("horizcoord", &horizcoord, "horizcoord/D");
  tree2->Branch("vertcoord", &vertcoord, "vertcoord/D");
  tree2->Branch("scalefactor_distance", &scalefactor_distance, "scalefactor_distance/D");
  tree2->Branch("scalefactor_attenuation", &scalefactor_attenuation, "scalefactor_attenuation/D");
  
  TTree *tree3 = new TTree("h3000", "h3000"); // tree3 if signal is detectable.
  tree3->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  tree3->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  tree3->Branch("theta_threshold_deg", &theta_threshold_deg, "theta_threshold_deg/D");
  tree3->Branch("nsigma_em_threshold", &nsigma_em_threshold, "nsigma_em_threshold/D");
  tree3->Branch("nsigma_had_threshold", &nsigma_had_threshold, "nsigma_had_threshold/D");
  tree3->Branch("horizcoord", &horizcoord, "horizcoord/D");
  tree3->Branch("vertcoord", &vertcoord, "vertcoord/D");
  tree3->Branch("vmmhz_max", &vmmhz_max, "vmmhz_max/D");
  tree3->Branch("vmmhz_min", &vmmhz_min, "vmmhz_min/D");
  tree3->Branch("dviewangle_deg", &dviewangle_deg, "dviewangle_deg/D");
  tree3->Branch("viewangle_deg", &viewangle_deg, "viewangle_deg/D");
  tree3->Branch("changle_deg", &changle_deg, "changle_deg/D");
  tree3->Branch("cosviewangle", &cosviewangle, "cosviewangle/D");
  tree3->Branch("emfrac", &emfrac, "emfrac/D");
  tree3->Branch("hadfrac", &hadfrac, "hadfrac/D");
      
  TTree *tree5 = new TTree("h5000", "h5000"); // tree5 filled for each nutau.
  tree5->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  tree5->Branch("inu", &inu, "inu/I");
  tree5->Branch("nuexitlength", &nuexitlength, "nuexitlength/D");
  tree5->Branch("nuexitice",  &nuexitice,  "nuexitice");
  tree5->Branch("vmmhz_max", &vmmhz_max, "vmmhz_max");
  tree5->Branch("maxtaper", &maxtaper, "maxtaper");
  tree5->Branch("inu", &inu, "inu/I");
  tree5->Branch("whichray", &whichray, "whichray/I");
  tree5->Branch("pnu", &pnu, "pnu/D");
  tree5->Branch("costhetanu", &costhetanu, "costhetanu/D");
  tree5->Branch("viewangle", &viewangle, "viewangle/D");
  tree5->Branch("offaxis", &offaxis, "offaxis/D");
  tree5->Branch("nsigma_offaxis", &nsigma_offaxis, "nsigma_offaxis/D");
  tree5->Branch("hadfrac", &hadfrac, "hadfrac/D");
  tree5->Branch("emfrac", &emfrac, "emfrac/D");
  tree5->Branch("sumfrac", &sumfrac, "sumfrac/D");
  tree5->Branch("horizcoord", &horizcoord, "horizcoord/D");
  tree5->Branch("vertcoord", &vertcoord, "vertcoord/D");
  tree5->Branch("weight1", &weight1, "weight1/D");
  tree5->Branch("nearthlayers", &nearthlayers, "nearthlayers/D");
  tree5->Branch("logchord", &logchord2, "interaction1->logchord/D");
  tree5->Branch("diff_3tries", &diff_3tries, "diff_3tries/D");
  tree5->Branch("fresnel2", &fresnel2, "fresnel2/D");
  tree5->Branch("costheta_inc", &costheta_inc, "costheta_inc/D");
  tree5->Branch("costheta_exit", &costheta_exit, "costheta_exit/D");
  tree5->Branch("deltheta_em", &deltheta_em[0], "deltheta_em/D");
  tree5->Branch("deltheta_had", &deltheta_had[0], "deltheta_had/F");
  tree5->Branch("r_fromballoon", &interaction1->r_fromballoon[0], "r_fromballoon/F");
  tree5->Branch("theta_in", &theta_in, "theta_in/D");
  tree5->Branch("lat_in", &lat_in, "lat_in/D");
  
    
  TTree *tree6 = new TTree("h6000", "h6000"); // tree6 filled for neutrinos that enter S of 60 deg S latitude.
  tree6->Branch("volts_rx_0", &volts_rx_0, "volts_rx_0/D");
  tree6->Branch("volts_rx_1", &volts_rx_1, "volts_rx_1/D");
  tree6->Branch("horizcoord", &horizcoord, "horizcoord/D");
  tree6->Branch("vertcoord", &vertcoord, "vertcoord/D");
  tree6->Branch("theta_in", &theta_in, "theta_in/D");
  tree6->Branch("chord_kgm2_bestcase", &chord_kgm2_bestcase2, "chord_kgm2_bestcase/D");
  tree6->Branch("chord_kgm2_ice", &interaction1->chord_kgm2_ice, "chord_kgm2_ice/D");
  tree6->Branch("costheta_nutraject", &interaction1->costheta_nutraject, "costheta_nutraject/D");
  tree6->Branch("weight1", &weight1, "weight1/D");
  tree6->Branch("weight_bestcase", &weight_bestcase2, "weight_bestcase/D");
  tree6->Branch("whichray", &whichray, "whichray/I");
  tree6->Branch("mybeta", &mybeta, "mybeta/D");
  tree6->Branch("longitude", &longitude_this, "longitude/D");
  
  
  TTree *tree6b = new TTree("h6001", "h6001"); // tree6b filled for the closest antenna to the interaction
  tree6b->Branch("bwslice_vnoise", bwslice_vnoise_thislayer, "bwslice_vnoise_thislayer[4]/D");
  
  TTree *tree7 = new TTree("h7000", "h7000"); // tree6 filled just after flavor is set
  tree7->Branch("emfrac", &emfrac, "emfrac/D");
  tree7->Branch("hadfrac", &hadfrac, "hadfrac/D");
  tree7->Branch("current", &interaction1->currentint, "currentint/I");
  tree7->Branch("nuflavor", &interaction1->nuflavorint, "nuflavorint/I");
  tree7->Branch("sumfrac", &sumfrac, "sumfrac/D");
  tree7->Branch("slopeyangle", &slopeyangle, "slopeyangle/D");
  
  TTree *jaimetree=new TTree("jaimetree", "jaimetree"); // signal as it is produced at the interaction
  jaimetree->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  jaimetree->Branch("emfrac", &emfrac, "emfrac/D");
  jaimetree->Branch("hadfrac", &hadfrac, "hadfrac/D");
  jaimetree->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  jaimetree->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  jaimetree->Branch("sumfrac", &sumfrac, "sumfrac/D");
  jaimetree->Branch("vmmhz1m_visible", &vmmhz1m_visible, "vmmhz1m_visible/D");
  
  TTree *viewangletree=new TTree("viewangletree", "viewangletree"); // signal as it is produced at the interaction
  viewangletree->Branch("dviewangle_deg", &dviewangle_deg, "dviewangle_deg/D");
  viewangletree->Branch("emfrac", &emfrac, "emfrac/D");
  viewangletree->Branch("hadfrac", &hadfrac, "hadfrac/D");
  viewangletree->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  viewangletree->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  viewangletree->Branch("theta_threshold_deg", &theta_threshold_deg, "theta_threshold_deg/D");
  viewangletree->Branch("dnutries", &interaction1->dnutries, "dnutries/D");
  viewangletree->Branch("viewangle", &viewangle, "viewangle/D");
  viewangletree->Branch("chord", &interaction1->chord, "dnutries/D");
  
  TTree *neutrino_positiontree=new TTree("neutrino_positiontree", "neutrino_positiontree");
  neutrino_positiontree->Branch("nnu", &interaction1->nnu, "nnu[3]/D");
  neutrino_positiontree->Branch("dtryingdirection", &interaction1->dtryingdirection, "dtryingdirection/D");
  neutrino_positiontree->Branch("bn1->dtryingposition", &bn1->dtryingposition, "bn1->dtryingposition/D");
  
  //Filled just after Getchord,  where we find the neutrino's path through the Earth
  TTree *nupathtree=new TTree("nupathtree", "nupathtree");
  nupathtree->Branch("total_kgm2", &total_kgm2, "total_kgm2/D");
  nupathtree->Branch("chord", &interaction1->chord, "chord/D");
  nupathtree->Branch("crust_entered", &crust_entered, "crust_entered/I");
  nupathtree->Branch("mantle_entered", &mantle_entered, "mantle_entered/I");
  nupathtree->Branch("core_entered", &core_entered, "core_entered/I");
  nupathtree->Branch("mybeta", &mybeta, "mybeta/D");
  nupathtree->Branch("costheta_nutraject", &interaction1->costheta_nutraject, "costheta_nutraject/D");
  
  
  
  TTree *finaltree = new TTree("passing_events", "passing_events"); // finaltree filled for all events that pass
  finaltree->Branch("inu", &inu, "inu/I");
  finaltree->Branch("vmmhz_min", &vmmhz_min, "vmmhz_min/D");
  finaltree->Branch("vmmhz_max", &vmmhz_max, "vmmhz_max/D");
  finaltree->Branch("horizcoord", &horizcoord, "horizcoord/D");
  finaltree->Branch("vertcoord", &vertcoord, "vertcoord/D");
  finaltree->Branch("horizcoord_bn", &bn1->horizcoord_bn, "horizcoord_bn/D");
  finaltree->Branch("vertcoord_bn", &bn1->vertcoord_bn, "vertcoord_bn/D");
  finaltree->Branch("r_bn", &r_bn_array, "r_bn_array[3]/D");
  finaltree->Branch("n_bn", &n_bn_array, "n_bn_array[3]/D");
  finaltree->Branch("longitude_bn", &longitude_this, "longitude_bn/D");
  finaltree->Branch("heading_bn", &heading_this, "heading_bn/D");
  finaltree->Branch("gps_offset", &gps_offset, "gps_offset/D");
  // this one is just weight due to earth absorption
  finaltree->Branch("weight1", &weight1, "weight1/D");
  // this is the total weight - the one you want to use!
  finaltree->Branch("weight", &weight, "weight/D");
  finaltree->Branch("logweight", &logweight, "logweight/D");
  finaltree->Branch("posnu", &posnu_array, "posnu_array[3]/D");
  finaltree->Branch("costheta_nutraject", &costheta_nutraject2, "costheta_nutraject/D");
  finaltree->Branch("chord_kgm2_ice",  &chord_kgm2_ice2, "chord_kgm2_ice/D");
  finaltree->Branch("phi_nutraject", &phi_nutraject2, "phi_nutraject/D");
  finaltree->Branch("altitude_int", &altitude_int2, "altitude_int/D");
  finaltree->Branch("nnu", &nnu_array, "nnu_array[3]/D");
  finaltree->Branch("n_exit2bn", &n_exit2bn_array, "n_exit2bn_array[5][3]/D");
  finaltree->Branch("n_exit_phi", &n_exit_phi, "n_exit_phi/D");
  finaltree->Branch("rfexit", &rfexit_array, "rfexit_array[5][3]/D");
  
  finaltree->Branch("pnu", &pnu, "pnu/D");
  finaltree->Branch("elast_y", &elast_y, "elast_y/D");
  finaltree->Branch("emfrac", &emfrac, "emfrac/D");
  finaltree->Branch("hadfrac", &hadfrac, "hadfrac/D");
  finaltree->Branch("sumfrac", &sumfrac, "sumfrac/D");
  finaltree->Branch("nuflavor", &nuflavorint2, "nuflavorint/I");//1=electron,  2=muon,  3=tau
  finaltree->Branch("current", &currentint2, "currentint/I");//0=charged current,  1=neutral current
  finaltree->Branch("logchord",  &logchord2,  "logchord/D");
  finaltree->Branch("nuexitice",  &nuexitice,  "nuexitice/D");
  finaltree->Branch("weight_bestcase",  &weight_bestcase2,  "weight_bestcase/D");
  finaltree->Branch("chord_kgm2_bestcase",  &chord_kgm2_bestcase2,  "chord_kgm2_bestcase/D");
  finaltree->Branch("dtryingdirection",  &dtryingdirection2,  "dtryingdirection/D");
  finaltree->Branch("l3trig", &l3trig, "l3trig[2]/I");
  finaltree->Branch("l2trig", &l2trig, "l2trig[2][3]/I");
  finaltree->Branch("l1trig", &l1trig, "l1trig[2][3]/I");
  finaltree->Branch("phiTrigMask", &phiTrigMask, "phiTrigMask/s");
  finaltree->Branch("phiTrigMaskH", &phiTrigMaskH, "phiTrigMaskH/s");
  finaltree->Branch("l1TrigMask", &l1TrigMask, "l1TrigMask/s");
  finaltree->Branch("l1TrigMaskH", &l1TrigMaskH, "l1TrigMaskH/s");
  finaltree->Branch("max_antenna0", &max_antenna0, "max_antenna0/I");
  finaltree->Branch("max_antenna1", &max_antenna1, "max_antenna1/I");
  finaltree->Branch("max_antenna2", &max_antenna2, "max_antenna2/I");
  
  finaltree->Branch("viewangle", &viewangle, "viewangle/D");
  finaltree->Branch("offaxis", &offaxis, "offaxis/D");
  finaltree->Branch("rx0_signal_eachband", &rx0_signal_eachband, "rx0_signal_eachband[2][5]/D");
  finaltree->Branch("rx0_threshold_eachband", &rx0_threshold_eachband, "rx0_threshold_eachband[2][5]/D");
  finaltree->Branch("rx0_noise_eachband", &rx0_noise_eachband, "rx0_noise_eachband[2][5]/D");
  finaltree->Branch("rx0_passes_eachband", &rx0_passes_eachband, "rx0_passes_eachband[2][5]/I");
  finaltree->Branch("e_component", &e_component, "e_component/D");
  finaltree->Branch("h_component", &h_component, "h_component/D");
  finaltree->Branch("dist_int_bn_2d", &dist_int_bn_2d, "dist_int_bn_2d/D");
  finaltree->Branch("d1", &d12, "d1/D");
  
  finaltree->Branch("cosalpha", &cosalpha, "cosalpha/D");
  finaltree->Branch("mytheta", &mytheta, "mytheta/D");
  finaltree->Branch("cosbeta0", &cosbeta0, "cosbeta0/D");
  finaltree->Branch("mybeta", &mybeta, "mybeta/D");
  finaltree->Branch("d1", &d12, "d1/D");
  finaltree->Branch("d2", &d22, "d2/D");
  
  //Begin block added by Stephen for verification plots
  finaltree->Branch("fresnel1", &fresnel1, "fresnel1/D");
  finaltree->Branch("fresnel2", &fresnel2, "fresnel2/D");
  finaltree->Branch("mag1", &mag1, "mag1/D");
  finaltree->Branch("mag2", &mag2, "mag2/D");
  finaltree->Branch("t_coeff_pokey", &t_coeff_pokey, "t_coeff_pokey/D");
  finaltree->Branch("t_coeff_slappy", &t_coeff_slappy, "t_coeff_slappy/D");
  finaltree->Branch("exponent", &settings1->EXPONENT, "EXPONENT/D");
  
  finaltree->Branch("hitangle_e_all", &hitangle_e_all, "hitangle_e_all[48]/D");
  finaltree->Branch("hitangle_h_all", &hitangle_h_all, "hitangle_h_all[48]/D");
  
  finaltree->Branch("e_comp_max1", &e_comp_max1, "e_comp_max1/D");
  finaltree->Branch("h_comp_max1", &h_comp_max1, "h_comp_max1/D");
  finaltree->Branch("e_comp_max2", &e_comp_max2, "e_comp_max2/D");
  finaltree->Branch("h_comp_max2", &h_comp_max2, "h_comp_max2/D");
  finaltree->Branch("e_comp_max3", &e_comp_max3, "e_comp_max3/D");
  finaltree->Branch("h_comp_max3", &h_comp_max3, "h_comp_max3/D");
  finaltree->Branch("max_antenna_volts0", &max_antenna_volts0, "max_antenna_volts0/D");
  finaltree->Branch("max_antenna_volts0_em", &max_antenna_volts0_em, "max_antenna_volts0_em/D");
  finaltree->Branch("max_antenna_volts1", &max_antenna_volts1, "max_antenna_volts1/D");
  finaltree->Branch("max_antenna_volts2", &max_antenna_volts2, "max_antenna_volts2/D");
  finaltree->Branch("triggers", &nchannels_perrx_triggered, "nchannels_perrx_triggered[48]/I");
  finaltree->Branch("nchannels_triggered", &nchannels_triggered, "nchannels_triggered/I");
  finaltree->Branch("volts_rx_max", &volts_rx_max, "volts_rx_max/D");
  finaltree->Branch("volts_rx_ave", &volts_rx_ave, "volts_rx_ave/D");
  finaltree->Branch("volts_rx_sum", &volts_rx_sum, "volts_rx_sum/D");
  
  finaltree->Branch("volts_rx_max_highband", &volts_rx_max_highband, "volts_rx_max_highband/D");
  finaltree->Branch("volts_rx_max_lowband", &volts_rx_max_lowband, "volts_rx_max_lowband/D");
  finaltree->Branch("theta_pol_measured", &theta_pol_measured, "theta_pol_measured/D");
  finaltree->Branch("theta_rf_atbn", &theta_rf_atbn, "theta_rf_atbn/D");
  finaltree->Branch("theta_rf_atbn_measured", &theta_rf_atbn_measured, "theta_rf_atbn_measured/D");
  finaltree->Branch("voltage", &voltagearray, "voltagearray[48]/D");
  finaltree->Branch("nlayers", &settings1->NLAYERS, "settings1->NLAYERS/I");
  
  finaltree->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  finaltree->Branch("vmmhz_lowfreq", &vmmhz_lowfreq, "vmmhz_lowfreq/D");
  
  finaltree->Branch("deltheta_em_max", &deltheta_em_max, "deltheta_em_max/D");
  finaltree->Branch("deltheta_had_max", &deltheta_had_max, "deltheta_had_max/D");
  finaltree->Branch("r_enterice", &r_enterice_array, "r_enterice_array[3]/D");
  finaltree->Branch("n_exit2bn_db", &n_exit2bn_db_array, "n_exit2bn_db_array[5][3]/D");
  
  finaltree->Branch("rfexit_db", &rfexit_db_array, "rfexit_db_array[5][3]/D");
  finaltree->Branch("r_in", &r_in_array, "r_in_array[3]/D");
  finaltree->Branch("nsurf_rfexit", &nsurf_rfexit_array, "nsurf_rfexit_array[3]/D");
  finaltree->Branch("nsurf_rfexit_db", &nsurf_rfexit_db_array, "nsurf_rfexit_db_array[3]/D");
  finaltree->Branch("r_fromballoon", &r_fromballoon2, "r_fromballoon/D");
  finaltree->Branch("r_fromballoon_db", &interaction1->r_fromballoon_db, "r_fromballoon_db/D");
  
  finaltree->Branch("nuexitlength", &nuexitlength, "nuexitlength/D");
  finaltree->Branch("nuentrancelength", &nuentrancelength, "nuentrancelength/D");
  finaltree->Branch("taulength", &taulength, "taulength/D");
  finaltree->Branch("icethickness", &icethickness, "icethickness/D");
  finaltree->Branch("nrf_iceside", &nrf_iceside_array, "nrf_iceside_array[5][3]/D");
  finaltree->Branch("nrf_iceside_db", &nrf_iceside_db_array, "nrf_iceside_db_array[5][3]/D");
  finaltree->Branch("ant_normal0", &ant_max_normal0_array, "ant_max_normal0_array[3]/D");
  finaltree->Branch("ant_normal1", &ant_max_normal1_array, "ant_max_normal1_array[3]/D");
  finaltree->Branch("ant_normal2", &ant_max_normal2_array, "ant_max_normal2_array[3]/D");
  finaltree->Branch("vmmhz1m_visible", &vmmhz1m_visible, "vmmhz1m_visible/D");
  finaltree->Branch("freq_bins", &freq_bins, "freq_bins/I");
  finaltree->Branch("vmmhz", &vmmhz, "vmmhz[freq_bins]/D");
  
  finaltree->Branch("dist_int_bn_2d_chord", &dist_int_bn_2d_chord, "dist_int_bn_2d_chord/D");
  
  finaltree->Branch("dviewangle_deg", &dviewangle_deg, "dviewangle_deg/D");
  finaltree->Branch("theta_threshold_deg", &theta_threshold_deg, "theta_threshold_deg/D");
  finaltree->Branch("total_kgm2", &total_kgm2, "total_kgm2/D");
  finaltree->Branch("chord", &interaction1->chord, "chord/D");
  finaltree->Branch("crust_entered", &crust_entered, "crust_entered/I");
  finaltree->Branch("mantle_entered", &mantle_entered, "mantle_entered/I");
  finaltree->Branch("core_entered", &core_entered, "core_entered/I");
  finaltree->Branch("n_pol", &n_pol_array, "n_pol_array[3]/D");
  finaltree->Branch("vmmhz_min_thatpasses", &vmmhz_min_thatpasses, "vmmhz_min_thatpasses/D");
  
  finaltree->Branch("pieceofkm2sr", &pieceofkm2sr, "pieceofkm2sr/D");
  //finaltree->Branch("volts_original", &volts_original, "volts_original[10][20][2]/D");
  finaltree->Branch("r_exit2bn", &r_exit2bn2, "r_exit2bn/D");
  finaltree->Branch("r_exit2bn_measured", &r_exit2bn_measured2, "r_exit2bn_measured/D");
  finaltree->Branch("scalefactor_attenuation", &scalefactor_attenuation, "scalefactor_attenuation/D");
  finaltree->Branch("anita1->PHI_OFFSET", &anita1->PHI_OFFSET, "anita1->PHI_OFFSET/D");
  finaltree->Branch("igps", &bn1->igps, "igyps/I");
  finaltree->Branch("volts_rx_rfcm_lab_e_all", &volts_rx_rfcm_lab_e_all, "volts_rx_rfcm_lab_e_all[48][512]/D");
  finaltree->Branch("volts_rx_rfcm_lab_h_all", &volts_rx_rfcm_lab_h_all, "volts_rx_rfcm_lab_h_all[48][512]/D");
  finaltree->Branch("ptaui", &ptaui, "ptaui/D");
  finaltree->Branch("ptauf", &ptauf, "ptauf/D");
  finaltree->Branch("sourceLon", &sourceLon, "sourceLon/D");
  finaltree->Branch("sourceLat", &sourceLat, "sourceLat/D");
  finaltree->Branch("sourceAlt", &sourceAlt, "sourceAlt/D");
  finaltree->Branch("sourceMag", &sourceMag, "sourceMag/D");
  
  TTree *mytaus_tree = new TTree("mytaus", "mytaus");
  mytaus_tree->Branch("taus",  &TauPtr);
  
  double rms_rfcm_e;
  double rms_rfcm_h;
  double rms_lab_e;
  double rms_lab_h;
  
  double avgfreq_rfcm[Anita::NFREQ];
  double avgfreq_rfcm_lab[Anita::NFREQ];
  double freq[Anita::NFREQ];
  
  TTree *summarytree = new TTree("summarytree", "summarytree"); // finaltree filled for all events that pass
  summarytree->Branch("NNU", &NNU, "NNU/I");
  summarytree->Branch("EXPONENT", &settings1->EXPONENT, "EXPONENT/D");
  summarytree->Branch("eventsfound_beforetrigger", &eventsfound_beforetrigger, "eventsfound_beforetrigger/D");
  summarytree->Branch("rms_rfcm_e", &rms_rfcm_e, "rms_rfcm_e/D");
  summarytree->Branch("rms_rfcm_h", &rms_rfcm_h, "rms_rfcm_h/D");
  summarytree->Branch("rms_lab_e", &rms_lab_e, "rms_lab_e/D");
  summarytree->Branch("rms_lab_h", &rms_lab_h, "rms_lab_h/D");
  summarytree->Branch("avgfreq_rfcm", &avgfreq_rfcm, "avgfreq_rfcm[128]/D");
  summarytree->Branch("avgfreq_rfcm_lab", &avgfreq_rfcm_lab, "avgfreq_rfcm_lab[128]/D");
  summarytree->Branch("freq", &freq, "freq[128]/D");
  
  TTree *banana_tree = new TTree("banana_tree", "banana_tree");  //To record banana plot info - Stephen
  banana_tree->Branch("r_bn", &bn1->r_bn, "r_bn[3]/D");
  
  TTree *ytree = new TTree("ytree", "ytree"); //To record y distributions
  ytree->Branch("elast_y", &elast_y, "elast_y/D");
  
  double icethck;
  double elev;
  double lon_ground;
  double lat_ground;
  double lon_ice;
  double lat_ice;
  double h20_depth;
  double lon_water;
  double lat_water;
  
  TTree *icetree = new TTree("icetree", "icetree");
  icetree->Branch("icethck", &icethck, "icethck/D");
  icetree->Branch("lon_ice", &lon_ice, "lon_ice/D");
  icetree->Branch("lat_ice", &lat_ice, "lat_ice/D");
  icetree->Branch("lon_water", &lon_water, "lon_water/D");
  icetree->Branch("lat_water", &lat_water, "lat_water/D");
  icetree->Branch("h20_depth", &h20_depth, "h20_depth/D");
  
  TTree *groundtree = new TTree("groundtree", "groundtree");
  groundtree->Branch("elev", &elev, "elev/D");
  groundtree->Branch("lon_ground", &lon_ground, "lon_ground/D");
  groundtree->Branch("lat_ground", &lat_ground, "lat_ground/D");
  
  //End block added by Stephen
    
  TTree *tree11 = new TTree("h11000", "h11000"); // tree11
  tree11->Branch("loctrig00", &loctrig[0][0], "loctrig0/D");
  tree11->Branch("loctrig10", &loctrig[1][0], "loctrig0/D");
  tree11->Branch("loctrig20", &loctrig[2][0], "loctrig0/D");
  tree11->Branch("loctrig_nadironly0", &loctrig_nadironly[0], "loctrig_nadironly0/D");
  tree11->Branch("loctrig01", &loctrig[0][1], "loctrig1/D");
  tree11->Branch("loctrig11", &loctrig[1][1], "loctrig1/D");
  tree11->Branch("loctrig21", &loctrig[2][1], "loctrig1/D");
  tree11->Branch("loctrig_nadironly1", &loctrig_nadironly[1], "loctrig0/D");
    
  TTree *tree16 = new TTree("h16000", "h16000");
  tree16->Branch("pnu", &pnu, "pnu/D");
  tree16->Branch("ptau", &ptau, "ptau/D");
  tree16->Branch("taulength", &taulength, "taulength/D");
  tree16->Branch("weight1", &weight1, "weight1/D");
  tree16->Branch("emfrac", &emfrac, "emfrac/D");
  tree16->Branch("hadfrac", &hadfrac, "hadfrac/D");
  tree16->Branch("nuentrancelength", &nuentrancelength, "nuentrancelength/D");
  
  int pdgcode;
  
  TTree *tree18 = new TTree("h18000", "h18000");
  tree18->Branch("emfrac",  &emfrac,  "emfrac/D");
  tree18->Branch("hadfrac", &hadfrac, "hadfrac/D");
  tree18->Branch("pdgcode", &pdgcode, "pdgcode/I");
  
    
  TH1D *h1mybeta = new TH1D("betaforall", "betaforall(deg)", 180, -15, 15);
  TH1D *h1mytheta= new TH1D("mytheta", "mytheta(deg)", 180, -90, 90);//90-incidentangle when neutrinos enter the Earth.
  TH1F *hundogaintoheight_e=new TH1F("undogaintoheight_e", "undogaintoheight_e", 100, 0., 1.);
  TH1F *hundogaintoheight_h=new TH1F("undogaintoheight_h", "undogaintoheight_h", 100, 0., 1.);
  TH1F *rec_diff=new TH1F("rec_diff", "rec_diff", 100, -1., 1.);
  TH1F *recsum_diff=new TH1F("recsum_diff", "recsum_diff", 100, -1., 1.);
  TH1F *rec_diff0=new TH1F("rec_diff0", "rec_diff0", 100, -1., 1.);
  TH1F *rec_diff1=new TH1F("rec_diff1", "rec_diff1", 100, -1., 1.);
  TH1F *rec_diff2=new TH1F("rec_diff2", "rec_diff2", 100, -1., 1.);
  TH1F *rec_diff3=new TH1F("rec_diff3", "rec_diff3", 100, -1., 1.);
    
  TTree *vmmhz_tree = new TTree("vmmhz_tree", "vmmhz_tree"); //To record frequency spread at point where it is first filled
  vmmhz_tree->Branch("freq_bins", &freq_bins, "freq_bins/I");
  vmmhz_tree->Branch("vmmhz", &vmmhz, "vmmhz[freq_bins]/D");
  
  TTree *tree1 = new TTree("h1000", "h1000"); // tree1 filled for each neutrino
  tree1->Branch("inu", &inu, "inu/I");
  tree1->Branch("diffexit", &diffexit, "diffexit/D");
  tree1->Branch("diffrefr", &diffrefr, "diffrefr/D");
  tree1->Branch("horizcoord", &horizcoord, "horizcoord/D");
  tree1->Branch("vertcoord", &vertcoord, "vertcoord/D");
  tree1->Branch("costhetanu", &costhetanu, "costhetanu/D");
  tree1->Branch("vmmhz1m_max", &vmmhz1m_max, "vmmhz1m_max/D");
  tree1->Branch("volume_thishorizon", &volume_thishorizon, "volume_thishorizon/D");
  tree1->Branch("realtime", &realtime_this, "realtime/D");
  tree1->Branch("longitude", &longitude_this, "longitude/D");
  tree1->Branch("latitude", &latitude_this, "latitude/D");
  tree1->Branch("MAXHORIZON", &bn1->MAXHORIZON, "MAXHORIZON/D");
  tree1->Branch("igps", &bn1->igps, "igps/I");
  tree1->Branch("passes_thisevent", &passes_thisevent, "passes_thisevent/I");
  tree1->Branch("igps", &bn1->igps, "igps/I");
  tree1->Branch("weight", &weight, "weight/D");
  tree1->Branch("r_exit2bn", &interaction1->r_exit2bn, "r_exit2bn/D");
  tree1->Branch("bn1->igps", &bn1->igps, "bn1->igps/I");
    
  // set up balloontree 

  TTree *balloontree = new TTree("balloon", "balloon"); //filled for all events
  balloontree->Branch("heading", &bn1->heading, "heading/D");
  balloontree->Branch("pitch", &bn1->pitch, "pitch/D");
  balloontree->Branch("roll", &bn1->roll, "roll/D"); 
  balloontree->Branch("realTime_flightdata", &bn1->realTime_flightdata, "realTime_flightdata/I"); 
  balloontree->Branch("latitude", &bn1->latitude, "latitude/D");
  balloontree->Branch("longitude", &bn1->longitude, "longitude/D");
  balloontree->Branch("altitude", &bn1->altitude, "altitude/D"); 
  balloontree->Branch("horizcoord_bn", &bn1->horizcoord_bn, "horizcoord_bn/D"); 
  balloontree->Branch("vertcoord_bn", &bn1->vertcoord_bn, "vertcoord_bn/D"); 

  
  // these variables are for energy reconstruction studies
  double undogaintoheight_e=0;
  double undogaintoheight_h=0;
  
  double undogaintoheight_e_array[4];
  double undogaintoheight_h_array[4];
  
  double nbins_array[4];
    
  double rec_efield=0;
  double true_efield=0;
  
  double rec_efield_array[4];
  double true_efield_array[4];
  // end energy reconstruction variables
  
#ifdef ANITA_UTIL_EXISTS

  string outputAnitaFile =settings1->outputdir+"/SimulatedAnitaEventFile"+run_num+".root";
  TFile *anitafileEvent = new TFile(outputAnitaFile.c_str(), "RECREATE");

  TTree *eventTree = new TTree("eventTree", "eventTree");
  eventTree->Branch("event",             &realEvPtr           );
  eventTree->Branch("run",               &run_no,   "run/I"   );
  eventTree->Branch("weight",            &weight,   "weight/D");

  outputAnitaFile =settings1->outputdir+"/SimulatedAnitaHeadFile"+run_num+".root";
  TFile *anitafileHead = new TFile(outputAnitaFile.c_str(), "RECREATE");
  
  TTree *headTree = new TTree("headTree", "headTree");
  headTree->Branch("header",  &rawHeaderPtr           );
  headTree->Branch("weight",  &weight,      "weight/D");

  outputAnitaFile =settings1->outputdir+"/SimulatedAnitaGpsFile"+run_num+".root";
  TFile *anitafileGps = new TFile(outputAnitaFile.c_str(), "RECREATE");

  UInt_t eventNumber;
  TTree *adu5PatTree = new TTree("adu5PatTree", "adu5PatTree");
  adu5PatTree->Branch("pat",          &Adu5PatPtr                   );
  adu5PatTree->Branch("eventNumber",  &eventNumber,  "eventNumber/I");
  adu5PatTree->Branch("weight",       &weight,       "weight/D"     );

  AnitaGeomTool *AnitaGeom1 = AnitaGeomTool::Instance();
  
#ifdef ANITA3_EVENTREADER
  outputAnitaFile =settings1->outputdir+"/SimulatedAnitaTruthFile"+run_num+".root";
  TFile *anitafileTruth = new TFile(outputAnitaFile.c_str(), "RECREATE");

  TTree *truthAnitaTree = new TTree("truthAnitaTree", "Truth Anita Tree");
  truthAnitaTree->Branch("truth",     &truthEvPtr                   );
#endif
  
#endif
  
  
  
  //end ROOT variable definitions
  ///////////////////////////////////////////////////////////////////////
  
  if (bn1->WHICHPATH==3) {
    settings1->CONSTANTCRUST = 1;  //Set ice depths and elevations all the same
    settings1->CONSTANTICETHICKNESS = 1;
    settings1->FIXEDELEVATION = 1;
  }
  
  IceModel *antarctica = new IceModel(settings1->ICE_MODEL + settings1->NOFZ*10, settings1->CONSTANTICETHICKNESS * 1000 + settings1->CONSTANTCRUST * 100 + settings1->FIXEDELEVATION * 10 + 0, settings1->WEIGHTABSORPTION);
  cout << "area of the earth's surface covered by antarctic ice is " << antarctica->ice_area << "\n";
  
  for (int i=0;i<antarctica->nRows_ice;i++) {
    for (int j=0;j<antarctica->nCols_ice;j++) {
      antarctica->IceENtoLonLat(j, i, lon_ice, lat_ice);
      icethck=antarctica->IceThickness(lon_ice, lat_ice);
      h20_depth=antarctica->water_depth[j][i];
      icetree->Fill();
    }
  }
  for (int i=0;i<antarctica->nRows_ground;i++) {
    for (int j=0;j<antarctica->nCols_ground;j++) {
      antarctica->GroundENtoLonLat(j, i, lon_ground, lat_ground);
      elev=antarctica->SurfaceAboveGeoid(lon_ground, lat_ground);
      groundtree->Fill();
    }
  }
  
  // fills arrays according to antenna specs
  anita1->GetBeamWidths(settings1); // this is used if GAINS set to 0
  // Antenna measured gain vs. frequency
  anita1->ReadGains(); // this is used if GAINS set to 1
  anita1->Set_gain_angle(settings1, sig1->NMEDIUM_RECEIVER);
  if(settings1->WHICH == 2 || settings1->WHICH == 6) anita1->SetDiffraction(); // for the upper ring
  
  TCanvas *cgains=new TCanvas("cgains", "cgains", 880, 800);
  TGraph *ggains=new TGraph(anita1->NPOINTS_GAIN, anita1->frequency_forgain_measured, anita1->vvGaintoHeight);
  ggains->Draw("al");
  stemp=settings1->outputdir+"/gains.eps";
  cgains->Print((TString)stemp);  
  
  // sets position of balloon and related quantities
  // these are all passed as pointers
  // theta,  phi,  altitude of balloon
  // position of balloon,  altitude and position of surface of earth (relative to the center of the earth) under balloon
  bn1->SetDefaultBalloonPosition(antarctica);
  
  anita1->SetNoise(settings1, bn1, antarctica);
  //pathtree->Fill(); //Added by Stephen for verification of path
  
  // find the maximum distance the interaction could be from the balloon and still be within the horizon.
  antarctica->GetMAXHORIZON(bn1);
  
  // calculate the volume of ice seen by the balloon for all balloon positions
  antarctica->CreateHorizons(settings1, bn1, bn1->theta_bn, bn1->phi_bn, bn1->altitude_bn, foutput);
  cout << "Done with CreateHorizons.\n";
  
  // get energy at which LPM effect turns on.
  //elpm=sig1->GetLPM();
  // elpm=sig1->GetELPM();
  
  // sets neutrino energy
  //  if(EXPONENT>15&&EXPONENT<25)//if EXPONENT is set to be a standard energy within ANITA's energy sensitivity
  if ( spectra1->IsMonoenergetic() ){
    pnu=pow(10., settings1->EXPONENT);
    primary1->GetSigma(pnu, sigma, len_int_kgm2, settings1, 0, 0);    // get cross section and interaction length.
    cout << "pnu,  sigma,  len_int_kgm2 are " << pnu << " " << sigma << " " << len_int_kgm2 << "\n";
  }
  
  if (settings1->WRITEPOSFILE==1) {
    nu_out << "Neutrinos with energy " << pnu << "\n\n"; //Write header to file of neutrino positions
    nu_out << "nu #,  position of nu interaction,  depth of int.,  Direction of nu momentum,  Position of balloon,  nu flavour,  current type,  elasticity\n\n\n\n";
  }
  
  if (bn1->WHICHPATH==3) { //Force certain parameters if we're doing a banana plot
    NNU = settings1->horizontal_banana_points*settings1->vertical_banana_points; //Total number of points to look at
    anita1->maxthreshold = Interaction::banana_sigma;
    settings1->SIGNAL_FLUCT = Interaction::banana_signal_fluct;
    settings1->CONSTANTCRUST=1;  //Set ice depths and elevations all the same
    settings1->CONSTANTICETHICKNESS=1;
    settings1->FIXEDELEVATION=1;
    pnu=Interaction::pnu_banana;
    int_banana->surface_over_banana_nu=antarctica->Surface(int_banana->nu_banana);
    int_banana->nu_banana_surface=(int_banana->surface_over_banana_nu) * (int_banana->nu_banana);
  } //Done setting parameters for banana plot

  cout << "reminder that I took out ChangeCoord.\n";  
  
  // builds payload based on read inputs
  anita1->GetPayload(settings1,  bn1);

  if (settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME==5) {
    Vector plusz(0., 0., 1.);
    bn1->PickBalloonPosition(plusz, antarctica, settings1, anita1);
    anita1->calculate_all_offsets();
    double angle_theta=16.;
    double angle_phi=0.;
    Vector x = Vector(cos(angle_theta * RADDEG) * cos((angle_phi+11.25) * RADDEG),
                      cos(angle_theta * RADDEG) * sin((angle_phi+11.25) * RADDEG),
                      sin(angle_theta * RADDEG));  
    anita1->GetArrivalTimes(x,bn1,settings1);
    cout << "end of getarrivaltimes\n";
  }

  // comment this out for running on auger cluster
  // get positions of the anita payload during the slac test
  if (settings1->SLAC)
    bn1->GetSlacPositions(anita1);
  
  for (int j=0;j<settings1->NLAYERS;j++) { // loop over the layers of the payload
    // noise depends on cant angle of antenna
    
  } //for loop over antenna layers
  
  if (settings1->WHICHRAYS==1) {
    settings1->MINRAY=0;
    settings1->MAXRAY=0;
  } //if (settings1->WHICHRAYS==1)
  if (settings1->WHICHRAYS==2) {
    settings1->MINRAY=0;
    settings1->MAXRAY=1;
  } //if (settings1->WHICHRAYS==2)
  if (settings1->WHICHRAYS==3) {
    settings1->MINRAY=1;
    settings1->MAXRAY=1;
  } //if (settings1->WHICHRAYS==3)
  
  time_t raw_loop_start_time = time(NULL);
  cout<<"Starting loop over events.  Time required for setup is "<<(int)((raw_loop_start_time - raw_start_time)/60)<<":"<< ((raw_loop_start_time - raw_start_time)%60)<<endl;
  
  // Spectra *spectra1 = new Spectra((int)EXPONENT);
  //   TGraph *gE2=new TGraph(10, spectra1->energy, spectra1->E2dNdEdAdt);
  //   TCanvas *cE2=new TCanvas("cE2", "cE2", 880, 800);
  //cout<<"assign Spectra"<<endl;
  
  TCanvas *ctest1 = new TCanvas("ctest1", "", 880, 800);
  
  spectra1->GetGEdNdEdAdt()->Draw("al");
  spectra1->GetGEdNdEdAdt()->GetHistogram()->SetMinimum(-0.2*spectra1->Getmaxflux());
  spectra1->GetSEdNdEdAdt()->SetLineColor(2);
  spectra1->GetSEdNdEdAdt()->Draw("l same");
    
  stemp=settings1->outputdir+"/GetG_test1.pdf";
  ctest1->Print((TString)stemp);
  
  //   gE2->Draw("al");
  //   //  spectra1->gspectrum[(int)EXPONENT]->Draw("al");
  //   stemp=settings1->outputdir+"cE2.pdf";
  //   cE2->Print((TString)stemp);
  
  
  
  // for averaging balloon altitude and distance from center of earth
  // for comparing with Peter
  double average_altitude=0.;
  double average_rbn=0.;

  //  TRandom r(0); // use a seed generated using machine clock (different every second)
  TRandom r(settings1->SEED); // use seed set as input

  signal(SIGINT,  interrupt_signal_handler);     // This function call allows icemc to gracefully abort and write files as usual rather than stopping abruptly.

  // Setting gps offset for plotting direction wrt north
  if (settings1->WHICH==7){
    gps_offset=atan2(-0.7042,0.71)*DEGRAD;
  } else if(settings1->WHICH==8){
    gps_offset=atan2(-0.7085,0.7056)*DEGRAD;
  } else if (settings1->WHICH==9 || settings1->WHICH==10){
    gps_offset=45;
  } else gps_offset=0;

  // begin looping over NNU neutrinos doing the things
  for (inu = 0; inu < NNU; inu++) {
    
    if (NNU >= 100) {
      if (inu % (NNU / 100) == 0)
        cout << inu << " neutrinos. " << (double(inu)/double(NNU)) * 100 << "% complete.\n";
    }
    else
      cout << inu << " neutrinos.  " << (double(inu) / double(NNU)) * 100 << "% complete.\n";
    


    for (whichray = settings1->MINRAY; whichray <= settings1->MAXRAY; whichray++) {
      anita1->passglobtrig[0]=0;
      anita1->passglobtrig[1]=0;
      passes_thisevent=0;
      unmasked_thisevent=1;
      vmmhz_min_thatpasses=1000; // initializing.  want to find the minumum voltage that passes a
      if ( spectra1->IsSpectrum() ){//if using energy spectrum
        //pnu=spectra1->GetNuEnergy();
        cout<<"using spectrum \n";
        pnu=spectra1->GetCDFEnergy();
        cout<<"pnu is "<<pnu<<"\n";
        ierr=primary1->GetSigma(pnu, sigma, len_int_kgm2, settings1, 0, 0);	// given neutrino momentum,  cross section and interaction length of neutrino.
        // ierr=0 if the energy is too low for the parameterization
        // ierr=1 otherwise
        len_int=1.0/(sigma*sig1->RHOH20*(1./M_NUCL)*1000); // in km (why interaction length in water?) //EH
      }// end IsSpectrum
      n_interactions=1;
      count_pass=0;
      passestrigger=0;
      chanceinhell2=0;
      sec1->secondbang=false;
      count_total++;
      // initializing the voltage seen by each polarization of each antenna
      interaction1->dtryingdirection=0;
      bn1->dtryingposition=0;
      for (int i=0; i<Anita::NFREQ;i++) {
        vmmhz[i] = 0.; // the full signal with all factors accounted for (1/r,  atten. etc.)
        vmmhz_em[i]=0.; // for keeping track of just the em component of the shower
      } //Zero the vmmhz array - helpful for banana plots,  shouldn't affect anything else - Stephen

      // Picks the balloon position and at the same time sets the masks and thresholds
      bn1->PickBalloonPosition(antarctica,  settings1,  inu,  anita1,  r.Rndm());
      // also sets phiTrigMask
      phiTrigMask=anita1->phiTrigMask;
      phiTrigMaskH=anita1->phiTrigMaskH;
      l1TrigMask=anita1->l1TrigMask;
      l1TrigMaskH=anita1->l1TrigMaskH;

      // find average balloon altitude and distance from center of earth for
      // making comparisons with Peter
      average_altitude+=bn1->altitude_bn/(double)NNU;
      average_rbn+=bn1->r_bn.Mag()/(double)NNU;
      
      realtime_this=bn1->realTime_flightdata;
      longitude_this=bn1->longitude;
      latitude_this=bn1->latitude;
      altitude_this=bn1->altitude;
      heading_this=bn1->heading;
      
      if (settings1->HIST && !settings1->ONLYFINAL 
        && prob_eachphi_bn->GetEntries() < settings1->HIST_MAX_ENTRIES) {
        prob_eachphi_bn->Fill(bn1->phi_bn);
        prob_eachilon_bn->Fill(bn1->r_bn.Lon());
      }
      
      if (bn1->WHICHPATH==3) { // for banana plot
        //Set observation location
        bn1->setObservationLocation(int_banana, inu, antarctica, settings1);
      } //End else if (WHICHPATH==3) : Banana plot locations
      balloontree->Fill();
      
      // pick random point in ice.
      // also get initial guess shower exit position
      // unit vector from shower exit to balloon.
      // get both positions of interaction point and its mirror point.
      //-------------------------------------------------------
      beyondhorizon = 0;
      
      if (interaction1)
        delete interaction1;
      interaction1 = new Interaction("nu",  primary1,  settings1,  whichray,  count1);

      if(taus1)
        delete taus1;
      taus1 = new Taumodel();
      //pnu=pow(10., settings1->EXPONENT); // pnu already obtained above
 
      int taumodes = settings1->taumodes;
      tauweighttrigger=0;
      interaction1->weight_nu=0;
      interaction1->weight_nu_prob=0;
      taus1->weight_tau_prob=0;
      
      if (taumodes==1 && interaction1->nuflavor=="nutau" && interaction1->current=="cc")
        tautrigger=1;//!< tau trigger sets the chance to create tau particle
      else
        tautrigger=0;
      
      bn1->PickDownwardInteractionPoint(interaction1,  anita1,  settings1,  antarctica,  ray1,  beyondhorizon);
           
      if (interaction1->noway)
        continue;
      count1->noway[whichray]++;
      
      if (interaction1->wheredoesitleave_err)
        continue;
      count1->wheredoesitleave_err[whichray]++;
      
      if (interaction1->neverseesice)
        continue;
      count1->neverseesice[whichray]++;
      
      if (interaction1->wheredoesitenterice_err)
        continue;
      count1->wheredoesitenterice_err[whichray]++;
      
      if (interaction1->toohigh)
        continue;
      count1->toohigh[whichray]++;
      
      if (interaction1->toolow)
        continue;
      count1->toolow[whichray]++;
      
      
      
      if (bn1->WHICHPATH==3)
        interaction1=int_banana;
      //we pick both the interaction point (interaction1->posnu) and its corresponding mirror point (interaction1->posnu_down)
      
      
      if (!interaction1->iceinteraction)
        continue;
      count1->iceinteraction[whichray]++;
      
      if (beyondhorizon) {
        //delete ray1;
        //delete interaction1;
        continue;
      }
      count1->inhorizon[whichray]++;
      
      
      // cerenkov angle depends on depth because index of refraction depends on depth.
       
      if (settings1->FIRN) {
        sig1->SetNDepth(antarctica->GetN(interaction1->altitude_int));
        //      changle = acos(1/N_DEPTH);
        changle_deg=sig1->changle*DEGRAD;
      }
      
      if (settings1->FORSECKEL==1)
        sig1->SetChangle(acos(1/sig1->NICE));
      
      //      cout << "altitude_int,  n_depth are " << altitude_this << " " << latitude_this << " " << longitude_this << " " << sig1->N_DEPTH << "\n";
      
      // x and y components of interaction in km.
      
      horizcoord=interaction1->posnu[0]/1000;
      vertcoord=interaction1->posnu[1]/1000;
      
            
      ray1->GetSurfaceNormal(settings1, antarctica, interaction1->posnu, slopeyangle, 0);

      // *** warning **** for Snell's law,  I call the ray on the air-side
      // the incident angle and the ice-side ray the refracted
      
      // ray's angle of incidence (in the air) onto ice
      costheta_inc=ray1->n_exit2bn[0]*ray1->nsurf_rfexit;    // just for plotting

      // just for plotting     
      costheta_exit=cos(ray1->rfexit[0].Theta()); // just for plotting
      
      if (!ray1->TraceRay(settings1, anita1, 1, sig1->N_DEPTH)) {
        //delete ray1;
        //delete interaction1;
        continue;
      }
      
      //       // use snell's law to get the first guess at the
      //       // direction of the rf as it leaves ice surface.
      //       // 0th guess was simply radially outward from interaction position
      //       // this now takes into account balloon position and surface normal.
      
      // fills ray1->n_exit2bn[1]
      ray1->GetRFExit(settings1, anita1, whichray, interaction1->posnu, interaction1->posnu_down, bn1->r_bn, bn1->r_boresights, 1, antarctica);
            
      ray1->GetSurfaceNormal(settings1, antarctica, interaction1->posnu, slopeyangle, 1);
      
      if (!ray1->TraceRay(settings1, anita1, 2, sig1->N_DEPTH)) {; // trace ray,  2nd iteration.
        //delete ray1;
        //delete interaction1;
        continue;
      }
      
      // fills ray1->n_exit2bn[2] ?
      ray1->GetRFExit(settings1, anita1, whichray, interaction1->posnu, interaction1->posnu_down, bn1->r_bn, bn1->r_boresights, 2, antarctica);
      
      ray1->GetSurfaceNormal(settings1, antarctica, interaction1->posnu, slopeyangle, 2);
            
      if (bn1->WHICHPATH==4)  // if this is for comparison with Peter,  print angles of incidence
        ray1->PrintAnglesofIncidence();
      
      // intermediate counter
      count1->nraypointsup1[whichray]++;
       
      
      sec1->GetTauDecay(interaction1->nuflavor, interaction1->current, taudecay,  emfrac_db,  hadfrac_db);
      
      // pick elasticity
        
      // EH check
      //cout<<"pnu right before it first used: "<<pnu<<endl;

      elast_y=primary1->pickY(settings1, pnu, 0, 0);
      if (settings1->CONSTANTY==1) { // if we ask to make y a constant=0.2
        elast_y=0.2;
        interaction1->nuflavor="nue";
        interaction1->current="cc";
      }
      if (bn1->WHICHPATH==3)
        elast_y = Interaction::banana_y;
      if (bn1->WHICHPATH==4)
        elast_y=1.;
      if (settings1->FORSECKEL==1) {
        if (settings1->SHOWERTYPE==0) // all hadronic shower
          elast_y=1.;
        if (settings1->SHOWERTYPE==1) // all em shower
          elast_y=0.;
      } //if (settings1->FORSECKEL)
      
      if (ytree->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
        ytree->Fill();
      
      
      //TAU STUFF. Pick whether it will stay as a neutrino or create tau
      if(tautrigger==1){
        if ((settings1->ROUGHNESS==0 || !settings1->UNBIASED_SELECTION) && !settings1->SLAC ) {
          err=GetDirection(settings1, interaction1, ray1->nrf_iceside[4], deltheta_em_max, deltheta_had_max, emfrac, hadfrac, vmmhz1m_max*bestcase_atten, interaction1->r_fromballoon[whichray], ray1, sig1, interaction1->posnu, anita1, bn1, interaction1->nnu, costhetanu, theta_threshold);
          //cout<<"UNBIASED_SELECTION IS "<<settings1->UNBIASED_SELECTION<<"\n";
        }
        else if (settings1->SLAC) {
          Vector xaxis(1., 0., 0.);
          //nnu=(rfexit[0].Unit()).Rotate(-10.*RADDEG, interaction1->posnu.Cross(zaxis));
          interaction1->nnu = xaxis.RotateY(bn1->theta_bn-settings1->SLAC_HORIZDIST/EarthModel::R_EARTH);  //direction of neutrino- for slac,  that's the direction of the beam
          interaction1->nnu = interaction1->nnu.RotateZ(bn1->phi_bn);
          costhetanu=cos(interaction1->nnu.Theta());
          theta_threshold=1.; // this is a bogus theta_threshold but it is only used for plotting anyway
          if (settings1->BORESIGHTS) {
            fslac_viewangles << bn1->sslacpositions[bn1->islacposition] << "\n";
            for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
              for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
                viewangle_eachboresight[ilayer][ifold]=acos(interaction1->nnu.Dot(ray1->nrf_iceside_eachboresight[4][ilayer][ifold]));
                fslac_viewangles << ilayer << "\t" << ifold << "\t" << (viewangle_eachboresight[ilayer][ifold]-sig1->changle)*DEGRAD << "\n";
              }//end for ifold
            }//end for ilayer
          }//end if boresights
          err=1; // everything is a-okay
        }// end else if slac
	
        if(err==0)
          continue;//bad stuff has happened.

        interaction1->r_in = antarctica->WhereDoesItEnter(interaction1->posnu, interaction1->nnu);
        
        taus1->GetTauWeight(primary1,  settings1,  antarctica,  interaction1,  pnu,  1,  ptauf, crust_entered);


        antarctica->Getchord(settings1, len_int_kgm2, interaction1->r_in, interaction1->r_enterice, interaction1->nuexitice, interaction1->posnu, inu, interaction1->chord, interaction1->weight_nu_prob, interaction1->weight_nu, nearthlayers, myair, total_kgm2, crust_entered,  mantle_entered, core_entered);
	
        nutauweight = interaction1->weight_nu_prob;
        tauweight = taus1->weight_tau_prob;

        nutauweightsum +=nutauweight;
        tauweightsum +=tauweight;
        double xrndm=gRandom->Rndm();

        if(xrndm <=taus1->weight_tau_prob/(taus1->weight_tau_prob+interaction1->weight_nu_prob)){
          pnu=ptauf;//set the energy we are looking at to the final energy of the tau. cuts out alot of if-else statements
          tauweighttrigger=1;//From now on,  signifies a tau particle  
        }
        if(TauPtr)
          delete TauPtr;

        TauPtr = new Taumodel();
        TauPtr->inu = inu;
        TauPtr->ptauf = ptauf;
        TauPtr->weight_nu_prob = interaction1->weight_nu_prob;
        TauPtr->weight_tau_prob = taus1->weight_tau_prob;
        mytaus_tree->Fill();

        //delete TauPtr;
      }//end tautrigger ==1
      /////////////////////// end of tau stuff
      
     
      
      // get fraction of shower that is electromagnetic.
      // pi^0's are counted as hadronic.
      sec1->GetEMFrac(settings1, interaction1->nuflavor, interaction1->current, taudecay, elast_y, hy, pnu, inu,emfrac, hadfrac, n_interactions, tauweighttrigger);

      if (emfrac+hadfrac>1.000001) {
        cout << "Warning:  " << inu << " " << emfrac+hadfrac << "\n";
      }
      
      // for plotting
      sumfrac=emfrac+hadfrac;
      //cout << "tree7 check" <<interaction1->nuflavorint << endl;
      if (tree7->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
        tree7->Fill();
      
      // EH check
      //cout<<"vmmhz1m_max : "<<vmmhz1m_max<<endl; // why vmmhz1m_max is always 1.0909? any meaning?

      vmmhz1m_visible = (emfrac+hadfrac)*vmmhz1m_max; //Stephen - Record actual V/m/Mhz for display
      
      // plots for debugging.
      
      if (interaction1->nuflavor=="numu" && bn1->WHICHPATH != 3 && !settings1->ONLYFINAL && settings1->HIST==1 && fraction_sec_muons->GetEntries()<settings1->HIST_MAX_ENTRIES) {
        fraction_sec_muons->Fill(emfrac+hadfrac, weight);
        n_sec_muons->Fill((double)n_interactions);
      }
      
      if (interaction1->nuflavor=="nutau" && bn1->WHICHPATH != 3 && !settings1->ONLYFINAL && settings1->HIST==1 && fraction_sec_taus->GetEntries()<settings1->HIST_MAX_ENTRIES) {
        fraction_sec_taus->Fill(emfrac+hadfrac, weight);
        n_sec_taus->Fill((double)n_interactions);
      }
      
      // for double bangs
      if(sec1->secondbang && sec1->interestedintaus) {
        ptau=(1-elast_y)*pnu;
        emfrac=emfrac_db;
        hadfrac=hadfrac_db;
      }
      
      // Find the highest possible electric field emitted at this energy
      // (this corresponds to the electric field at the highest frequency
      // detectable by the antennas.
      // Also find the maximum width of Cerenkov cone (which is at lowest frequency)
      // These are used to find the maximum angular deviation from Cerenkov
      // cone where signal is still detectable.
      if(sec1->secondbang && sec1->interestedintaus) {
        vmmhz1m_max=sig1->GetVmMHz1m(ptau, anita1->FREQ_HIGH);
        sig1->GetSpread(ptau, emfrac, hadfrac, anita1->FREQ_LOW, deltheta_em_max, deltheta_had_max);
        //EH check
        //cout<<"tau case ";
      } //if (secondbang && interestedintaus)
      else {// get peak signal at highest edge of frequency band because that is where it is highest
        vmmhz1m_max=sig1->GetVmMHz1m(pnu, anita1->FREQ_HIGH);
        sig1->GetSpread(pnu, emfrac, hadfrac, anita1->FREQ_LOW, 
        deltheta_em_max, deltheta_had_max);
        //cout << "pnu,  emfrac,  hadfrac,  deltheta_em_max,  deltheta_had_max are " << pnu << " " << emfrac << " " << hadfrac << " " << deltheta_em_max*DEGRAD << " " << deltheta_had_max*DEGRAD << "\n";
      } //end else (not secondbang or not interested in taus)

      // EH check
      //cout<<"vmmhz1m_max after GetVmMHz1m : "<<vmmhz1m_max<<endl;
      
      if (jaimetree->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
        jaimetree->Fill();
      
      //  Using highest possible signal and minimum noise, 
      // and accounting for distance from balloon, 
      // and best case attenutation, 
      // reject if signal is undetectable.
      //cout << "altitude_int is " << interaction1->altitude_int << "\n";
      
      if (whichray==0) // if we are looking at direct rays
        bestcase_atten=exp(interaction1->altitude_int/MAX_ATTENLENGTH); // the attenuation is obtained from the altitude of the interaction (shortest path is if the signal went straight up)
      if (whichray==1) // if we are looking at reflected rays
        bestcase_atten=exp(interaction1->altitude_int_mirror/MAX_ATTENLENGTH);//use the real path which seems from the mirror point.
      
      //      cout << "this is being tested. " << anita1->VNOISE[0]/10.*anita1->maxthreshold/((hadfrac+emfrac)*vmmhz1m_max*bestcase_atten/interaction1->r_fromballoon[whichray]*heff_max*anita1->bwmin/1.E6) << " against " << settings1->CHANCEINHELL_FACTOR << " " << settings1->SKIPCUTS << "\n";
      //cout << "hadfrac+emfrac,  vmmhz1m_max are " << hadfrac << "\t" << emfrac << "\t" << vmmhz1m_max << "\n";


      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/((hadfrac+emfrac)*vmmhz1m_max*bestcase_atten/interaction1->r_fromballoon[whichray]*heff_max*anita1->bwmin/1.E6)>settings1->CHANCEINHELL_FACTOR && !settings1->SKIPCUTS) {
      //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/((hadfrac+emfrac)*vmmhz1m_max*bestcase_atten/r_fromballoon*heff_max*bw/1.E6)>1. && !settings1->SKIPCUTS) {
        //delete ray1;
        //delete interaction1;
        continue; // by comparing highest possible signal to the lowest possible noise,  reject if there is just no way we could detect this event.
        // vmmhz1m_max=signal at highest frequency
        // bestcase_atten=best case attenuation
        // r_fromballoon=distance from interaction to balloon
        //heff_max=maximum effective height over the frequency range
      }
      
      
      // intermediate counter
      count1->nnottoosmall[whichray]++;
      
      // pick neutrino direction.
      // This function picks a neutrino direction such that its cerenkov cone
      // is close enough to the balloon line of sight that you have a chance in hell of seeing the signal.
      // It returns zero if emfrac and hadfrac are both zero (should never be true)
      // It returns 2 if the signal is undetectable. // It does not!
      
      if (whichray==1) {
        chengji=ray1->nrf_iceside[4]*ray1->nrf_iceside[0];//the projection of nrf_iceside[2] on the direction of radius direction of interaction1->posnu
        //original nrf_iceside[4] is the upgoing direction of signals after being reflected.
        //now I get the corresponding downward direction of real signals in my case.
        //The two vectors are symmetric to the tangent plane of the Earth at interaction point
        ray1->nrf_iceside[4] = ray1->nrf_iceside[4] - 2*chengji*ray1->nrf_iceside[0];
      } //if whichray==1
      
      if(tautrigger==0){//did this for cc- taus already,  do again for all other particles
        if ((settings1->ROUGHNESS==0 || !settings1->UNBIASED_SELECTION) && !settings1->SLAC )
          err=GetDirection(settings1, interaction1, ray1->nrf_iceside[4], deltheta_em_max, deltheta_had_max, emfrac, hadfrac, vmmhz1m_max*bestcase_atten, interaction1->r_fromballoon[whichray], ray1, sig1, interaction1->posnu, anita1, bn1, interaction1->nnu, costhetanu, theta_threshold);
        else if (settings1->SLAC) {
          Vector xaxis(1., 0., 0.);
          //nnu=(rfexit[0].Unit()).Rotate(-10.*RADDEG, interaction1->posnu.Cross(zaxis));
          interaction1->nnu = xaxis.RotateY(bn1->theta_bn-settings1->SLAC_HORIZDIST/EarthModel::R_EARTH);  //direction of neutrino- for slac,  that's the direction of the beam
          interaction1->nnu = interaction1->nnu.RotateZ(bn1->phi_bn);
          costhetanu=cos(interaction1->nnu.Theta());
          theta_threshold=1.; // this is a bogus theta_threshold but it is only used for plotting anyway
          if (settings1->BORESIGHTS) {
            fslac_viewangles << bn1->sslacpositions[bn1->islacposition] << "\n";
            for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
              for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
                viewangle_eachboresight[ilayer][ifold]=acos(interaction1->nnu.Dot(ray1->nrf_iceside_eachboresight[4][ilayer][ifold]));
                fslac_viewangles << ilayer << "\t" << ifold << "\t" << (viewangle_eachboresight[ilayer][ifold]-sig1->changle)*DEGRAD << "\n";
              }//end ifold
            }//end ilayer
          }//end boresight
          err=1; // everything is a-okay
        }//end else if slac
      }//end tau trigger ==0
      
      //if (settings1->ROUGHNESS==1)
      //  rough1->GetExitforRoughness(settings1, antarctica, emfrac, hadfrac, deltheta_em_mid2, deltheta_had_mid2, ray1, interaction1->nnu, bn1->r_bn, interaction1->posnu); // inputs
      //ray1 is modified in function

        // gets angle between ray and neutrino direction
      viewangle = GetViewAngle(ray1->nrf_iceside[4], interaction1->nnu);
      if(viewangle>1.57 && !settings1->SKIPCUTS) { //discard the event if viewangle is greater than 90 degrees
        //delete interaction1;
        continue;
      }
      
      count1->nviewangle_lt_90[whichray]++; // add to counter
      
      if (!Ray::WhereDoesItLeave(interaction1->posnu, interaction1->nnu, antarctica, interaction1->nuexit))
        continue; // doesn't give a real value from quadratic formula
      
      GetBalloonLocation(interaction1, ray1, bn1, antarctica);
      
      nuexitlength=interaction1->posnu.Distance(interaction1->nuexit);
      // probability a tau would decay within this length at this
      // energy
      nuexitice=interaction1->posnu.Distance(interaction1->nuexitice);
      theta_threshold_deg=theta_threshold*DEGRAD;
      
      // neutrino direction in frame where balloon is up,  0=east, 1=north, 2=up
      n_nutraject_ontheground = Vector(bn1->n_east*interaction1->nnu,  bn1->n_north*interaction1->nnu,  bn1->n_bn*interaction1->nnu);
      
      cosviewangle=cos(viewangle); // cosine angle
      viewangle_deg=viewangle*DEGRAD; // same angle but in degrees
      dviewangle_deg=(sig1->changle-viewangle)*DEGRAD; // deviation from cerenkov angle
      
      if (viewangletree->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
        viewangletree->Fill(); // fills variables related to viewing angle
      
      if (neutrino_positiontree->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
        neutrino_positiontree->Fill(); // fills variables related to neutrino position

      if (whichray==1) {
        //return it to the upgoing direction that is after being reflected
        ray1->nrf_iceside[4] = ray1->nrf_iceside[4] + 2*chengji*ray1->nrf_iceside[0];
      }
      
      if (err==0) {
        count1->nbadfracs[whichray]++;
        cout<<"err==0,  so leaving.\n";
        //delete ray1;
        //delete interaction1;
        continue;
      }
      
      count1->ngoodfracs[whichray]++;

      // these variables are just for plotting
      nsigma_em_threshold=theta_threshold/deltheta_em_max;
      nsigma_had_threshold=theta_threshold/deltheta_had_max;
      
      // For each neutrino,  multiply the number of tries
      // necessary to generate the appropriate direction, 
      // and the number of tries necessary to generate
      // an appropriate position,  and assume they are
      // independent.
      interaction1->dnutries=interaction1->dtryingdirection*bn1->dtryingposition;
      
      // for plotting aperture per ring radius from balloon
      index_distance=(int)(bn1->r_bn.SurfaceDistance(interaction1->posnu, bn1->surface_under_balloon) / (700000./(double)NBINS_DISTANCE));
      
      // where the neutrino enters the earth
      if (tautrigger==0){//did for cc-taus already,  do for all other particles
        interaction1->r_in = antarctica->WhereDoesItEnter(interaction1->posnu, interaction1->nnu);
       }
      
      // total chord
      double chord_kgm2_test=interaction1->posnu.Distance(interaction1->r_in)*sig1->RHOMEDIUM;
      
      double weight_test=0;  // weight if the whole chord from interaction to earth entrance is ice.
      // take best case scenario chord length and find corresponding weight
      
      IsAbsorbed(chord_kgm2_test, len_int_kgm2, weight_test);
      // if the probably the neutrino gets absorbed is almost 1,  throw it out.

      if (bn1->WHICHPATH!=4 && settings1->FORSECKEL!=1 && !settings1->SKIPCUTS) {
        if (weight_test<CUTONWEIGHTS) {
          //delete ray1;
          //delete interaction1;
	  continue;
        }
      }

      count_chanceofsurviving++;
      
      
      // if (tree8->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1 && bn1->WHICHPATH != 3)
      //   tree8->Fill();
      
      // theta of nu entrance point,  in earth frame
      // and latitude
      theta_in=interaction1->r_in.Theta();
      lat_in=-90+theta_in*DEGRAD;
      
      // find quantities relevent for studying impact of atmosphere
      // for black hole studies
      // costheta and mytheta: theta of neutrino wrt surface normal where neutrino enters earth
      // cosbeta0, mybeta: theta of neutrino wrt surface normal for a person standing above the interaction point
      myair=GetThisAirColumn(settings1,  interaction1->r_in, interaction1->nnu, interaction1->posnu, col1, cosalpha, mytheta, cosbeta0, mybeta);
      
      // where the neutrino enters the ice
      // reject if it enters beyond the borders of the continent.
      // step size is 1/10 of interaction length
      if (!settings1->FORSECKEL && !settings1->UNBIASED_SELECTION) {
        if (!antarctica->WhereDoesItEnterIce(interaction1->posnu, interaction1->nnu, len_int_kgm2/sig1->RHOMEDIUM/10., interaction1->r_enterice)) {
          //r_enterice.Print();
          if (antarctica->OutsideAntarctica(interaction1->r_enterice)) {
            cout<<"Warning!  Neutrino enters beyond continent,  program is rejecting neutrino! inu = "<<inu<<endl;
            if (bn1->WHICHPATH==3)
              cout<<"Warning!  Neutrino enters beyond continent,  program is rejecting neutrino!"<<endl;

            //delete ray1;
            //delete interaction1;
            continue;
          }// end outside antarctica
        }// end wheredoesitenterice
      }// end if !settings forseckel && unbiased
      
      // intermediate counter
      count1->nentersice[whichray]++;
      
      // d1=earth entrance to rock-ice interface
      // d2=rock-ice interface to position of neutrino interaction
      interaction1->d1=interaction1->r_enterice.Distance(interaction1->r_in);
      interaction1->d2=interaction1->r_enterice.Distance(interaction1->posnu);
      

      // get a lower limit on the chord that the neutrino traverses, 
      // so that later we can see if the signal is detectable in
      // the best case scenario.
      
      if(sec1->secondbang && sec1->interestedintaus) {
        cout << "Need to bring back GetFirstBang before you can simulate taus.\n";
        cout << "I removed it because it required EarthModel and I wanted Secondaries to be a stand-alone class to use in the embedded simulation.\n";
        //sec1->GetFirstBang(interaction1->r_in, nnu, interaction1->posnu2, len_int_kgm2, d1, nuentrancelength);
        //taulength=interaction1->posnu.Distance(interaction1->posnu2);
        icethickness=interaction1->r_enterice.Distance(interaction1->nuexit);
        interaction1->chord_kgm2_bestcase=nuentrancelength*Tools::dMin(densities, 3);
      }
      else {
        // finds minimum chord (in kg/m^2) traversed by neutrino
        // only keeping events with weight > 10^-3
        // periodically need to make sure this is still valid
        // chord_kgm2_bestcase=(d1+d2)*sig1->RHOMEDIUM;
        interaction1->chord_kgm2_bestcase=(interaction1->d1+interaction1->d2)*Tools::dMin(densities, 3);
      }

      // chord just through ice.
      interaction1->chord_kgm2_ice=interaction1->d2*sig1->RHOMEDIUM;
      
      
      if (tree6->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
        tree6->Fill();
      
      // take best case scenario chord length and find corresponding weight
      IsAbsorbed(interaction1->chord_kgm2_bestcase, len_int_kgm2, interaction1->weight_bestcase);
      
      // if the probability that the neutrino gets absorbed is almost 1,  throw it out.
      if (bn1->WHICHPATH!=4 && interaction1->weight_bestcase<CUTONWEIGHTS && !settings1->SKIPCUTS && !settings1->FORSECKEL && !settings1->ROUGHNESS) {
        if (bn1->WHICHPATH==3)
          cout<<"Neutrino is getting absorbed and thrown out!"<<endl;

        continue;
      }
      
      //intermediate counter
      count1->nabsorbed[whichray]++;
      
      // intermediate counter
      count1->nraywithincontinent1[whichray]++;
      
      // now we have second guess for rf exit point,  which is
      // pretty close to the right answer.
      // now get our best case attenuation again, 
      // and see if we can reject the event.
      if (settings1->ROUGHNESS==0) {
        if (whichray==0)
          bestcase_atten=exp(-1*ray1->rfexit[1].Distance(interaction1->posnu)/MAX_ATTENLENGTH);
        if (whichray==1)
          bestcase_atten=exp(-1*ray1->rfexit[1].Distance(interaction1->posnu_down)/MAX_ATTENLENGTH);//use the real distance
      }// end if setting roughness
      else {
        bestcase_atten=exp(-1*(ray1->rfexit[1].Mag()-interaction1->posnu.Mag())/MAX_ATTENLENGTH);
      }//end else

      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/((hadfrac+emfrac)*vmmhz1m_max*bestcase_atten/interaction1->r_fromballoon[whichray]*heff_max*anita1->bwmin/1.E6)>settings1->CHANCEINHELL_FACTOR && !settings1->SKIPCUTS && !settings1->FORSECKEL) {
      //if (CHANCEINHELL_FACTOR*Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/((hadfrac+emfrac)*vmmhz1m_max*bestcase_atten/r_fromballoon*heff_max*bw/1.E6)>1. && !settings1->SKIPCUTS && !settings1->FORSECKEL) {
        if (bn1->WHICHPATH==3)
          cout<<"Event rejected.  Check."<<endl;
;
        continue;
      }
      
      count_chanceinhell0++;
      
      // intermediate counting
      count1->nraypointsup2[whichray]++;
      
      double nbelowsurface;
      // reject if it is totally internally reflected at the surface
      // !!  AND NOT CONSIDERING ROUGHNESS
      if (!settings1->ROUGHNESS) {
        if (settings1->FIRN)
          nbelowsurface=NFIRN;
        else
          nbelowsurface=sig1->NICE;
        // this is purely a sanity check.
        // if everything is working,  events should pass with 100% efficiency
        if (TIR(ray1->nsurf_rfexit, ray1->nrf_iceside[3], nbelowsurface, sig1->N_AIR)) {
          //delete interaction1;
          continue;
        }
      }// end !settings roughness
      
      count1->nnottir[whichray]++;
      
      // this sets n_exit2bn[2] to the ray from the exit point to the balloon, 
      // last iteration.  Now we're ready to do some calculations!!!!
      ray1->GetRFExit(settings1, anita1, whichray, interaction1->posnu, interaction1->posnu_down, bn1->r_bn, bn1->r_boresights, 2, antarctica);
      
      count1->nraywithincontinent2[whichray]++;
      
      // for plotting- cos(theta) of neutrino direction standing on earth
      // below balloon.
      interaction1->costheta_nutraject=(interaction1->nnu*bn1->r_bn)/sqrt(bn1->r_bn*bn1->r_bn);
      
      theta_rf_atbn = ray1->n_exit2bn[2].Angle(bn1->r_bn); // polar angle of the rf signal as seen at the balloon.
      // measured theta of the rf,  which is actual smeared by SIGMA_THETA,  whose default is 0.5 degrees.
      theta_rf_atbn_measured = theta_rf_atbn+gRandom->Gaus()*anita1->SIGMA_THETA;
      interaction1->r_exit2bn=bn1->r_bn.Distance(ray1->rfexit[2]);
      interaction1->r_exit2bn_measured=bn1->altitude_bn/cos(theta_rf_atbn_measured);
      
      if((settings1->WHICH == 2 || settings1->WHICH == 6) && theta_rf_atbn < 0.3790091) {
        //delete interaction1;
        continue; // the deck will mess up the arrival times in the top ring
      }
      // reject if the rf leaves the ice where there is water,  for example.
      //cout << "rfexit[2] is ";ray1->rfexit[2].Print();
      if (!antarctica->AcceptableRfexit(ray1->nsurf_rfexit, ray1->rfexit[2], ray1->n_exit2bn[2])){
        if (bn1->WHICHPATH==3)
          cout<<"Should look at this. Not expecting to be here."<<endl;
        continue;
      }//end if acceptableRFexit
      
      
      // intermediate counting
      count1->nacceptablerf[whichray]++;
      
      // difference between exit points of 2nd and 3rd iterations.
      diff_3tries=ray1->rfexit[1].Distance(ray1->rfexit[2]);
      
      if (tree5->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1  && bn1->WHICHPATH != 3)
        tree5->Fill();
      
      // reject if 2nd and 3rd tries
      // don't converge within 10m.
      if (settings1->ROUGHNESS==0) {
        if (diff_3tries>10) {
          //delete ray1;
          //delete interaction1;
          continue;
        }
      }
      count1->nconverges[whichray]++;
            
      // Get Polarization vector.  See Jackson,  Cherenkov section.
      n_pol = GetPolarization(interaction1->nnu, ray1->nrf_iceside[4]);
      
      if (settings1->BORESIGHTS) {
        for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
          for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
            n_pol_eachboresight[ilayer][ifold]=GetPolarization(interaction1->nnu, ray1->nrf_iceside_eachboresight[4][ilayer][ifold]);
          } // end looping over antennas in phi
        } // end looping over layers
      } // if we are calculating for all boresights
      /*
      ofstream roughout; // length of chord in air vs. theta (deg)
      stemp = settings1->outputdir+"/scan_Npts_Ledge.dat";
      roughout.open(stemp);
      ofstream props; // length of chord in air vs. theta (deg)
      stemp = settings1->outputdir+"/eventproperties.dat";
      props.open(stemp);
*/
      if(!settings1->ROUGHNESS) {  // IF NO ROUGHNESS THEN DO THIS (FOR CONSISTENCY CHANGE NOTHING BELOW HERE IN THE if !rough)
/*
        // output some general event parameters
        props<<ray1->rfexit[2].Lon()<<"  "<<-90+ray1->rfexit[2].Lat()<<std::endl;
        props<<bn1->r_bn.Lon()<<"  "<<-90+bn1->r_bn.Lat()<<std::endl;
        props<<interaction1->posnu.Lon()<<"  "<<-90+interaction1->posnu.Lat()<<std::endl;
*/
        if (settings1->FIRN) {
          // now rotate that polarization vector according to ray paths in firn and air.
          // fresnel factor at ice-firn interface
          GetFresnel(rough1, settings1->ROUGHNESS, ray1->nsurf_rfexit, ray1->nrf_iceside[3], n_pol, ray1->nrf_iceside[4], vmmhz1m_max, emfrac, hadfrac, deltheta_em_max, deltheta_had_max, t_coeff_pokey, t_coeff_slappy, fresnel1, mag1);
          if (bn1->WHICHPATH==4)
            cout << "Lentenin factor is " << 1./mag1 << "\n";
  
          //The gradual transition in the firn means that there is no fresnel factor,  only magnification
          // and the magnification factor is upside down compared to what it is
          // for the firn-air interface
          vmmhz1m_fresneledonce = vmmhz1m_max/mag1;

          //  get fresnel factor at firn-air interface
          GetFresnel(rough1, settings1->ROUGHNESS, ray1->nsurf_rfexit, ray1->n_exit2bn[2], n_pol, ray1->nrf_iceside[3], vmmhz1m_fresneledonce, emfrac, hadfrac, deltheta_em_max, deltheta_had_max, t_coeff_pokey, t_coeff_slappy, fresnel2, mag2);
          // use both fresnel and magnification factors at firn-air interface.  Notice that magnification factor is
          //upside-down compared to what it is in the firn.
          vmmhz1m_fresneledtwice=vmmhz1m_fresneledonce*fresnel2*mag2;

          if (settings1->BORESIGHTS) { 
            for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
              for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
          	GetFresnel(rough1, settings1->ROUGHNESS, ray1->nsurf_rfexit, ray1->n_exit2bn_eachboresight[2][ilayer][ifold], 
          		   n_pol_eachboresight[ilayer][ifold], ray1->nrf_iceside_eachboresight[3][ilayer][ifold], 
          		   vmmhz1m_max, emfrac, hadfrac, deltheta_em_max, deltheta_had_max, t_coeff_pokey, t_coeff_slappy, 
          		   fresnel1_eachboresight[ilayer][ifold], mag1_eachboresight[ilayer][ifold]);
          	//		std::cout << fresnel1_eachboresight[ilayer][ifold] << std::endl;
              } // end looping over phi sectors
            } // end looping over layers
          } // end if we are calculating for all boresights

          
          if (bn1->WHICHPATH==4)
            cout << "firn-air interface:  fresnel2,  mag2 are " << fresnel2 << " " << mag2 << "\n";
          
        }//end if firn
        else {
          sig1->GetSpread(pnu, emfrac, hadfrac, (anita1->bwslice_min[2]+anita1->bwslice_max[2])/2., deltheta_em_mid2, deltheta_had_mid2);
          
          GetFresnel(rough1, settings1->ROUGHNESS, ray1->nsurf_rfexit, ray1->n_exit2bn[2], n_pol, ray1->nrf_iceside[4], vmmhz1m_max, emfrac, hadfrac, deltheta_em_mid2, deltheta_had_mid2, t_coeff_pokey, t_coeff_slappy,  fresnel1, mag1);
          
          vmmhz1m_fresneledtwice = vmmhz1m_max*fresnel1*mag1;  //  only the ice-air interface
          
          if (settings1->BORESIGHTS) {
            for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) { // loop over layers on the payload
              for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
                GetFresnel(rough1, settings1->ROUGHNESS, ray1->nsurf_rfexit, ray1->n_exit2bn_eachboresight[2][ilayer][ifold], n_pol_eachboresight[ilayer][ifold], ray1->nrf_iceside_eachboresight[4][ilayer][ifold], vmmhz1m_max, emfrac, hadfrac, deltheta_em_max, deltheta_had_max, t_coeff_pokey, t_coeff_slappy, fresnel1_eachboresight[ilayer][ifold], mag1_eachboresight[ilayer][ifold]);
              } // end looping over phi sectors
            } // end looping over layers
          } // end if we are calculating for all boresights
        }//end else firn
      }
      // OTHERWISE THERE IS ROUGHNESS SO DO MAGIC
      else{
        //(vector) ray1->nsurf_rfexit:  surface normal at RFexit position
        //(pos)        ->rfexit[2]:     final iterated position of RF exit
        //(vector)     ->n_exit2bn[2]:  vector from RF exit position TO balloon
        //(pos)    bn1->r_bn:           position of balloon
        //(vector) n_pol:               polarization vector
        //(pos)    posnu:               position of neutrino interaction

        //IMPORTANT NOTE: WE ARE IGNORING THE FIRN UNTIL FURTHER NOTICE FOR SIMPLICITY

        /////////////
        //  THESE ARE FOR DETERMINING SOME SPECULAR PARAMETERS (POL, MAGNITUDE, ..)
        //
        Vector pol_specular = GetPolarization(interaction1->nnu, ray1->nrf_iceside[4]);;

        double pathlength_specular = interaction1->posnu.Distance(ray1->rfexit[2]) + ray1->rfexit[2].Distance(bn1->r_bn);
        Vector vec_bln_vert = antarctica->GetSurfaceNormal( bn1->r_bn ).Unit();

        sig1->GetSpread(pnu, emfrac, hadfrac, (anita1->bwslice_min[2]+anita1->bwslice_max[2])/2., deltheta_em_mid2, deltheta_had_mid2);
        GetFresnel(rough1, 0, ray1->nsurf_rfexit, ray1->n_exit2bn[2], pol_specular, ray1->nrf_iceside[4], vmmhz1m_max, emfrac, hadfrac, deltheta_em_mid2, deltheta_had_mid2, t_coeff_pokey, t_coeff_slappy,  fresnel1, mag1);
        vmmhz1m_fresneledtwice = vmmhz1m_max*fresnel1*mag1;  //  only the ice-air interface
        vmmhz1m_fresneledtwice /= ( interaction1->posnu.Distance(ray1->rfexit[2]) + ray1->rfexit[2].Distance(bn1->r_bn) );
        Attenuate(antarctica, settings1, vmmhz1m_fresneledtwice,  interaction1->posnu.Distance(ray1->rfexit[2]),  interaction1->posnu);

        /////////////

        //double minimumScreenLength = 8.; // screen minimum size limit, default 8.
        double subscreenFractionLimit = 1e-4; // default 1e-4
        int maximumSubscreenGeneration = 2;  // value is inclusive; default 3
        int subscreenDivisions = 4;
        int num_validscreenpoints = 0;
        Position pos_basescreen_centralpos;
        Position pos_current;
        Vector pos_current_localnormal; // local normal at point on ground below screen current point
        Vector vec_pos_current_to_balloon;

        Position pos_projectedImpactPoint;
        Vector vec_specularnormal;      //normalized, normal vector at specular point
        Vector vec_localnormal;         //normalized, normal vector at projected ground point
        Vector vec_nnu_to_impactPoint;  //normalized

        double viewangle_local;
        double theta_local;             //angle between local surface normal and vector to balloon [radians]
        double theta_0_local;                 //angle between 'inverted' local surface normal and incident direction [radians]

        double element_sa;
        double interpolatedPower;       // \muW for these angles

        Vector npol_local_inc, npol_local_trans;
        Vector temp_a;

        double Emag_local, Emag_local_notaper;
        Vector Efield_local;
        Vector Efield_screentotal = Vector(0,0,0);

        std::vector<double> basescrn_Emags;             // container for base screen Efield
        std::vector<Position> basescrn_pos;             // container for base screen points, will sample to pass along to seedpositions below
        std::vector<double> basescrn_length;            // container for base screen edge lengths
        std::vector<Position> seedpositions;            // container for seed positions for additional passes
        std::vector<double> seedEdgeLengths;            // container for edge lengths for additional passes
        std::vector<int> seedGeneration;                // container for tracking generations for additional passes

        double pathlength_local;        // set for each screen point
        double vmmhz_local;   //  V/m/MHz at balloon (after all steps)
        double tempvar;

        //iterate points on the screen, get their position and project back to find ground impact
        //calculate incident and transmitted angles, look up power fraction and loss correction factor
        //add to running total, and output observables for checks

        //reset the counter and set the BASE screen properties based on current event/balloon geometry
        panel1->ResetPositionIndex();
        panel1->SetNsamples(20);  //default 20
        panel1->SetEdgeLength( settings1->SCREENEDGELENGTH );
        panel1->SetCentralPoint( ray1->rfexit[2] + 0.8*bn1->r_bn.Distance(ray1->rfexit[2])*ray1->n_exit2bn[2].Unit() ); //move towards balloon by 80% of distance to balloon
        panel1->SetNormal( ray1->n_exit2bn[2].Unit() );
        panel1->SetCosineProjectionFactor( ray1->n_exit2bn[2].Unit().Dot(antarctica->GetSurfaceNormal(ray1->rfexit[2]).Unit()) );
        pos_current_localnormal = antarctica->GetSurfaceNormal( panel1->GetCentralPoint() ).Unit();
        temp_a = pos_current_localnormal - ray1->n_exit2bn[2].Dot(pos_current_localnormal) * ray1->n_exit2bn[2].Unit();
        panel1->SetUnitX( ray1->n_exit2bn[2].Cross(temp_a).Unit() );  // set screen y-hat as up 'perp.' to ground in the screen plane, define x-hat as perp.
        panel1->SetUnitY( -1.*ray1->n_exit2bn[2].Cross(ray1->n_exit2bn[2].Cross(temp_a).Unit()) );
        pos_basescreen_centralpos = panel1->GetCentralPoint();

        // output some general event parameters
        temp_a = interaction1->nnu - pos_current_localnormal.Dot(interaction1->nnu) * pos_current_localnormal;
        Position tpos = interaction1->posnu + 500.*temp_a.Unit();
        //std::cerr<<"[ ";
        //std::cerr<<"["<<ray1->rfexit[2].Lon()<<",  "<<-90+ray1->rfexit[2].Lat()<<"],  ";
        //std::cerr<<"["<<bn1->r_bn.Lon()<<",  "<<-90+bn1->r_bn.Lat()<<"],  ";
        //std::cerr<<"["<<interaction1->posnu.Lon()<<",  "<<-90+interaction1->posnu.Lat()<<"],  ";
        //std::cerr<<"["<<tpos.Lon()<<",  "<<-90+tpos.Lat()<<"],  ";
        //std::cerr<<"["<<pol_specular.Dot(antarctica->GetSurfaceNormal(bn1->r_bn)) <<",  "<<(pol_specular - pol_specular.Dot(antarctica->GetSurfaceNormal(bn1->r_bn)) * antarctica->GetSurfaceNormal(bn1->r_bn)).Mag()<<",  "<< vmmhz1m_fresneledtwice <<"]  ";
        //std::cerr<<"],"<<std::endl;
/*        props<<ray1->rfexit[2].Lon()<<"  "<<-90+ray1->rfexit[2].Lat()<<std::endl;
        props<<bn1->r_bn.Lon()<<"  "<<-90+bn1->r_bn.Lat()<<std::endl;
        props<<interaction1->posnu.Lon()<<"  "<<-90+interaction1->posnu.Lat()<<std::endl;
        props<<tpos.Lon()<<"  "<<-90+tpos.Lat()<<std::endl;
        props<<pol_specular.Dot(antarctica->GetSurfaceNormal(bn1->r_bn)) <<"  "<<(pol_specular - pol_specular.Dot(antarctica->GetSurfaceNormal(bn1->r_bn)) * antarctica->GetSurfaceNormal(bn1->r_bn)).Mag()<<"  "<< vmmhz1m_fresneledtwice <<std::endl;*/


        //#########
        // First treat the base 5x5 screen to establish seed points for higher-res scans (we know that there will be sufficient variation across the base screen)
        //std::cerr<<"Base screen..."<<std::endl;
        double summedE=0.;
        for (int ii=0; ii<panel1->GetNsamples()*panel1->GetNsamples(); ii++){
          Emag_local = 0.;
          pos_current = panel1->GetNextPosition(ii);        // this gets the new screen position
          pos_projectedImpactPoint = Position(1,1,1);     // placeholder, is set below in WhereDoesItEnterIce()
          vec_pos_current_to_balloon = Vector( bn1->r_bn[0] - pos_current[0], bn1->r_bn[1] - pos_current[1], bn1->r_bn[2] - pos_current[2] );

          //Determine ground impact position where the projected ray enters the ice
          // reject if it enters beyond the borders of the continent.
          // step size is 10 meters
          if (!antarctica->WhereDoesItExitIce(pos_current, panel1->GetNormal(), 10., pos_projectedImpactPoint)){
            //std::cerr<<"Warning!  Projected ground impact position of screen point does NOT enter ice. Skipping this screen point."<<std::endl;
            continue;
          }
          if (antarctica->OutsideAntarctica(pos_projectedImpactPoint)) {
              //std::cerr<<"Warning!  Projected ground impact position of screen point is off-continent. Skipping this screen point."<<std::endl;
            continue;
          }// end outside antarctica

          // local angles of transmission and incidence in their respective planes
          vec_localnormal = antarctica->GetSurfaceNormal(pos_projectedImpactPoint).Unit();
          vec_nnu_to_impactPoint =  Vector( pos_projectedImpactPoint[0]-interaction1->posnu[0], pos_projectedImpactPoint[1]-interaction1->posnu[1], pos_projectedImpactPoint[2]-interaction1->posnu[2] ).Unit();
          theta_local = vec_localnormal.Angle( (const Vector)vec_pos_current_to_balloon ); //[rad]
          theta_0_local = vec_localnormal.Angle(vec_nnu_to_impactPoint); //[rad]

          /////
          // Field Magnitude
          element_sa = ( (panel1->GetEdgeLength()/panel1->GetNsamples()) / pos_projectedImpactPoint.Distance(pos_current) ) * 180./PI;
          interpolatedPower = rough1->InterpolatePowerValue(theta_0_local*180./PI, theta_local*180./PI);
          Emag_local = vmmhz1m_max
                        * sqrt(interpolatedPower / rough1->GetLaserPower() / rough1->GetLossCorrectionFactor(theta_0_local));
          // account for 1/r for 1)interaction point to impact point and 2)impact point to balloon, and attenuation in ice
          pathlength_local = interaction1->posnu.Distance(pos_projectedImpactPoint) + pos_projectedImpactPoint.Distance(bn1->r_bn);
          Emag_local /= pathlength_local ;
          Attenuate(antarctica, settings1, Emag_local,  interaction1->posnu.Distance(pos_projectedImpactPoint),  interaction1->posnu);

          summedE += Emag_local;
          basescrn_Emags.push_back(Emag_local);
          basescrn_pos.push_back(pos_current);
          basescrn_length.push_back(panel1->GetEdgeLength() / (float)panel1->GetNsamples());
        }

        //#########
        // Second, now select those points in the base screen that contribute most of the signal strength
        ofstream basescreenout; // length of chord in air vs. theta (deg)
        stemp = settings1->outputdir+"/basescreen.dat";
        basescreenout.open(stemp);
        for (int ii=0; ii< panel1->GetNsamples()*panel1->GetNsamples(); ii++){
          pos_current = panel1->GetNextPosition(ii);
          temp_a = Vector( pos_current[0] - pos_basescreen_centralpos[0], pos_current[1] - pos_basescreen_centralpos[1], pos_current[2] - pos_basescreen_centralpos[2] );
          basescreenout<<basescrn_Emags[ii]/summedE<<"  "
                  <<temp_a.Dot(panel1->GetUnitX())<<"  "
                  <<temp_a.Dot(panel1->GetUnitY())<<std::endl;
          //if (basescrn_Emags[ii]/summedE < 1.e-50){
          //  continue;
          //}
          seedpositions.push_back(basescrn_pos[ii]);
          seedEdgeLengths.push_back(basescrn_length[ii]);
          seedGeneration.push_back(1);
        }
        basescreenout.close();
        //std::cerr<<seedpositions.size()<<" good seeds."<<std::endl;

        //#########
        // Third, now loop over the seed positions, resetting the screen parameters for each seed position.
        //if the min(Efield) < fraction * max(Efield), then add these points to the seedposition vector to sample further
        for (int ii=0; ii< seedpositions.size(); ii++){
          //std::cerr<<"  Screen seed "<<ii;

          panel1->ResetPositionIndex();
          panel1->SetNsamples(subscreenDivisions); //default 5
          panel1->SetEdgeLength( seedEdgeLengths[ii] );
          panel1->SetCentralPoint( seedpositions[ii] );
          // the Nsamples, normal, cosine projection factor, and unit vectors are the same as for the base screen

          double minE=1.e22;
          double maxE = -1.;

          //these are per-screen and will be overwritten for each seed screen, use to store values before deciding if screen was 'valid'
          std::vector<Position> seedscreens_pos;
          std::vector<double> seedscreens_vmmhzlocal;  //use this to save Emag_local aka vmmhz[0] ONLY
          std::vector<Vector> seedscreens_2bln;
          std::vector<Vector> seedscreens_pols;
          std::vector<double> seedscreens_phasedelay;
          std::vector<Position> seedscreens_impactpt;
          std::vector<double> seedscreens_viewangle;

          for (int jj=0; jj< panel1->GetNsamples()*panel1->GetNsamples(); jj++){
            Emag_local = 0.;
            pos_current = panel1->GetNextPosition(jj);        // this gets the new screen position
            pos_projectedImpactPoint = Position(1,1,1);     // placeholder, is set below in WhereDoesItEnterIce()
            vec_pos_current_to_balloon = Vector( bn1->r_bn[0] - pos_current[0], bn1->r_bn[1] - pos_current[1], bn1->r_bn[2] - pos_current[2] );

            //Determine ground impact position where the projected ray enters the ice
            // reject if it enters beyond the borders of the continent.
            // step size is 10 meters
            if (!antarctica->WhereDoesItExitIce(pos_current, panel1->GetNormal(), 10., pos_projectedImpactPoint)){
              //std::cerr<<"Warning!  Projected ground impact position of screen point does NOT enter ice. Skipping this screen point."<<std::endl;
              continue;
            }
            if (antarctica->OutsideAntarctica(pos_projectedImpactPoint)) {
                //std::cerr<<"Warning!  Projected ground impact position of screen point is off-continent. Skipping this screen point."<<std::endl;
              continue;
            }// end outside antarctica

            // local angles of transmission and incidence in their respective planes
            vec_localnormal = antarctica->GetSurfaceNormal(pos_projectedImpactPoint).Unit();
            vec_nnu_to_impactPoint =  Vector( pos_projectedImpactPoint[0]-interaction1->posnu[0], pos_projectedImpactPoint[1]-interaction1->posnu[1], pos_projectedImpactPoint[2]-interaction1->posnu[2] ).Unit();
            theta_local = vec_localnormal.Angle( (const Vector)vec_pos_current_to_balloon ); //[rad]
            theta_0_local = vec_localnormal.Angle(vec_nnu_to_impactPoint); //[rad]

            viewangle_local = GetViewAngle(vec_nnu_to_impactPoint, interaction1->nnu);

            /////
            // Field Magnitude
            element_sa = ( panel1->GetEdgeLength() / pos_projectedImpactPoint.Distance(pos_current) ) * 180./PI;
            interpolatedPower = rough1->InterpolatePowerValue(theta_0_local*180./PI, theta_local*180./PI);
            Emag_local = vmmhz1m_max
                          * sqrt(interpolatedPower*element_sa / rough1->GetLaserPower() / rough1->GetLossCorrectionFactor(theta_0_local));
            // account for 1/r for 1)interaction point to impact point and 2)impact point to balloon, and attenuation in ice
            pathlength_local = interaction1->posnu.Distance(pos_projectedImpactPoint) + pos_projectedImpactPoint.Distance(bn1->r_bn);
            Emag_local /= pathlength_local ;
            Attenuate(antarctica, settings1, Emag_local,  interaction1->posnu.Distance(pos_projectedImpactPoint),  interaction1->posnu);
            Emag_local_notaper = Emag_local;
            //check the vmmhz[0] component ONLY here
            deltheta_em[0]=deltheta_em_max*anita1->FREQ_LOW/anita1->freq[0];
            deltheta_had[0]=deltheta_had_max*anita1->FREQ_LOW/anita1->freq[0];
            sig1->TaperVmMHz(viewangle_local, deltheta_em[0], deltheta_had[0], emfrac, hadfrac, Emag_local, vmmhz_em[0]);// this applies the angular dependence.
            //check if Taper would remove this
            if(Emag_local==0){
              //std::cerr<<"VmMHz == 0"<<std::endl;
              //seedpositions.push_back(pos_current);
              //seedEdgeLengths.push_back(panel1->GetEdgeLength()/panel1->GetNsamples());
              //std::cerr<<" -------"<<std::endl;
              continue;
            }

            /////
            // Transmitted Polarization
            npol_local_inc = GetPolarization(interaction1->nnu, vec_nnu_to_impactPoint).Unit();
            // calculate polarization transformation matrix coefficients (from Bahar 1995, page 535)
            // treat a_y = normal at specular exit point
            vec_specularnormal = antarctica->GetSurfaceNormal(ray1->rfexit[2]).Unit();
            // transmitted polarization needs to be perpendicular to to-balloon vector, and the horizontal component is 'set', so need to find appropriate vector for the vertical component to ensure perpendicularity
            npol_local_trans = panel1->CalculateTransmittedPolarization(interaction1->nnu, vec_specularnormal, vec_localnormal, vec_pos_current_to_balloon, vec_nnu_to_impactPoint, npol_local_inc);
            if(npol_local_trans[0]!=npol_local_trans[0]){
              //std::cerr<<"Bad transmitted polarization"<<std::endl;
              //seedpositions.push_back(pos_current);
              //seedEdgeLengths.push_back(panel1->GetEdgeLength()/panel1->GetNsamples());
              continue;   // skip if transmitted polarization is undefined
            }

            minE = Tools::dMin(minE, Emag_local_notaper);
            maxE = Tools::dMax(maxE, Emag_local_notaper);
            seedscreens_pos.push_back(pos_current);
            seedscreens_vmmhzlocal.push_back( Emag_local_notaper );
            //for (int kk=0; kk<Anita::NFREQ; kk++){
            //  seedscreens_vmmhzlocal.push_back( vmmhz_local[kk] );
            //}
            seedscreens_2bln.push_back(vec_pos_current_to_balloon);
            seedscreens_pols.push_back(npol_local_trans);
            for (int k=0;k<Anita::NFREQ;k++) {
              seedscreens_phasedelay.push_back( TWOPI*(pathlength_specular-pathlength_local)*anita1->freq[k]/CLIGHT * 180./PI );
            }
            seedscreens_impactpt.push_back(pos_projectedImpactPoint);
            seedscreens_viewangle.push_back(viewangle_local);
          }// end for jj for this seed screen

          // simply drop this screen if the maximum vmmhz[0] is below some threshold
          if (maxE < 1.e-50){
            continue;
          }

          //if( (minE >= subscreenFractionLimit*maxE) || (minimumScreenLength > panel1->GetEdgeLength()/panel1->GetNsamples()) ){ // store these and move on to the next seed screen
          if( (minE >= subscreenFractionLimit*maxE) || (seedGeneration[ii] == maximumSubscreenGeneration) ){ // store these and move on to the next seed screen
            //std::cerr<<ii<<"  "<<seedEdgeLengths[ii]<<"  saved."<<std::endl;
            double pol_bln_vert, pol_bln_horiz;
            for (int jj=0; jj<seedscreens_pos.size(); jj++){
              // increment the valid point counter so we can track the size of the screen's vmmhz_freq vector
              num_validscreenpoints++;

              //add the contribution to the running total
              panel1->AddVmmhz0(seedscreens_vmmhzlocal[jj]);
              panel1->AddVec2bln(seedscreens_2bln[jj]);
              panel1->AddPol(seedscreens_pols[jj]);
              for (int k=0;k<Anita::NFREQ;k++) {
                panel1->AddDelay( seedscreens_phasedelay[jj*Anita::NFREQ + k] );
              }
              panel1->AddImpactPt(seedscreens_impactpt[jj]);
              panel1->AddViewangle(seedscreens_viewangle[jj]);
              panel1->AddWeight( panel1->GetEdgeLength() / panel1->GetNsamples() );

              pol_bln_vert = npol_local_trans.Dot(vec_bln_vert);
              pol_bln_horiz = (npol_local_trans - pol_bln_vert*vec_bln_vert).Mag();
              // use temp_a to determine position of this screen point relative to the base screen's central position
              temp_a = Vector( seedscreens_pos[jj][0] - pos_basescreen_centralpos[0], seedscreens_pos[jj][1] - pos_basescreen_centralpos[1], seedscreens_pos[jj][2] - pos_basescreen_centralpos[2] );


/*              roughout<<inu<<"  "
                      <<seedscreens_vmmhzlocal[jj]<<"  "
                      <<seedscreens_impactpt[jj].Lon()<<"  "
                      <<-90+seedscreens_impactpt[jj].Lat()<<"  "
                      <<seedscreens_phasedelay[jj*Anita::NFREQ]*PI/180.<<"  "
                      <<panel1->GetEdgeLength()<<"  "
                      <<pol_bln_vert<<"  "
                      <<pol_bln_horiz<<"  "
                      <<temp_a.Dot(panel1->GetUnitX())<<"  "
                      <<temp_a.Dot(panel1->GetUnitY())<<"  "
                      <<std::endl;*/

            }// end for jj<seedscreens_pos
          }// end else
          else { //if true, reject these points and pass the positions back into seedpositions so we sample at higher resolutions (smaller edge length)
            //std::cerr<<ii<<"  "<<seedEdgeLengths[ii]<<"  <-- ! redo."<<std::endl;
            for (int jj=0; jj<seedscreens_pos.size(); jj++){
              seedpositions.push_back(seedscreens_pos[jj]);
              seedEdgeLengths.push_back(panel1->GetEdgeLength()/panel1->GetNsamples());
              seedGeneration.push_back(seedGeneration[ii]+1);
              /*roughout<<inu<<"  "
                      <<0<<"  "
                      <<seedscreens_impactpt[jj].Lon()<<"  "
                      <<-90+seedscreens_impactpt[jj].Lat()<<"  "
                      <<0<<"  "
                      <<panel1->GetEdgeLength()<<"  "
                      <<std::endl;*/
            }
          }
        }// end for ii loop over seedpositions

        panel1->SetNvalidPoints(num_validscreenpoints);

        //now construct the Screen's vmmhz array for all points, so it gets passed to the trigger object later to make the waveforms
        // here we get the array vmmhz by taking vmmhz1m_max (signal at lowest frequency bin) and vmmhz_max (signal at lowest frequency after applying 1/r factor and attenuation factor) and making an array across frequency bins by putting in frequency dependence.
        double validScreenSummedArea = 0.;
        double vmmhz_local_array[Anita::NFREQ];
        for (int jj=0; jj<panel1->GetNvalidPoints(); jj++){
          sig1->GetVmMHz(panel1->GetVmmhz0(jj), vmmhz1m_max, pnu, anita1->freq, anita1->NOTCH_MIN, anita1->NOTCH_MAX, vmmhz_local_array, Anita::NFREQ);
          // apply the off-angle tapering
          for (int k=0;k<Anita::NFREQ;k++) {
            deltheta_em[k]=deltheta_em_max*anita1->FREQ_LOW/anita1->freq[k];
            deltheta_had[k]=deltheta_had_max*anita1->FREQ_LOW/anita1->freq[k];
            sig1->TaperVmMHz(panel1->GetViewangle(jj), deltheta_em[k], deltheta_had[k], emfrac, hadfrac, vmmhz_local_array[k], vmmhz_em[k]);// this applies the angular dependence.
            panel1->AddVmmhz_freq(vmmhz_local_array[k]);
            //vmmhz[k] += vmmhz_local_array[k];
          }
          validScreenSummedArea += panel1->GetWeight(jj);
          Efield_local = panel1->GetVmmhz_freq(jj*Anita::NFREQ) * panel1->GetPol(jj);
          Efield_screentotal = Efield_screentotal + Efield_local;
        }//end jj over panel Nvalid points
        panel1->SetWeightNorm(validScreenSummedArea);
        vmmhz_max = Efield_screentotal.Mag();
        n_pol = Efield_screentotal.Unit();
      }//end else roughness

//      roughout.close();

      // reject if the event is undetectable.
      // THIS ONLY CHECKS IF ROUGHNESS == 0, WE WILL SKIP THIS IF THERE IS ROUGHNESS
      //if (vmmhz1m_fresneledtwice*heff_max*0.5*(bw/1.E6)<CHANCEINHELL_FACTOR*anita1->maxthreshold*Tools::dMin(VNOISE, settings1->NLAYERS) && !settings1->SKIPCUTS) {
      if (settings1->ROUGHNESS==0){
        if(settings1->CHANCEINHELL_FACTOR*vmmhz1m_fresneledtwice*heff_max*0.5*(anita1->bwmin/1.E6)<anita1->maxthreshold*anita1->VNOISE[0]/10.&& !settings1->SKIPCUTS) {
          if (bn1->WHICHPATH==3)
            cout<<"Event is undetectable.  Leaving loop."<<endl;

          continue;
        }
        count1->nchanceinhell_fresnel[whichray]++;
      } //end if CHANCEINHELL factor and SKIPCUTS
      
      //
      // for plotting
      diffexit=ray1->rfexit[0].Distance(ray1->rfexit[1]);
      diffnorm=acos(ray1->nsurf_rfexit[0]*ray1->nsurf_rfexit[1]);
      diffrefr=acos(ray1->nrf_iceside[4]*ray1->nrf_iceside[0]);
      
      // scale by 1/r once you've found the 3rd iteration exit point
      // ALREADY DEALT WITH IN CASE OF ROUGHNESS
      if (settings1->ROUGHNESS==0) {
        if (whichray==0)
          vmmhz_max=ScaleVmMHz(vmmhz1m_fresneledtwice, interaction1->posnu, bn1->r_bn);
        if (whichray==1)
          vmmhz_max=ScaleVmMHz(vmmhz1m_fresneledtwice, interaction1->posnu_down, bn1->r_bn);//use the mirror point
        // EH
        //cout<<"after applying ScaleVmMHz,  vmmhz_max: "<<vmmhz_max<<",  vmmhz1m_fresneledtwice: "<<vmmhz1m_fresneledtwice<<endl;
      }
      //else
      //  vmmhz_max=vmmhz1m_max;
      
      // reject if the event is undetectable.
      if (settings1->ROUGHNESS==0){
        if (settings1->CHANCEINHELL_FACTOR*vmmhz_max*heff_max*0.5*(anita1->bwmin/1.E6)<anita1->maxthreshold*anita1->VNOISE[0]/10. && !settings1->SKIPCUTS) {
          if (bn1->WHICHPATH==3)
            cout<<"Event is undetectable.  Leaving loop."<<endl;

          //delete ray1;
          //delete interaction1;
          continue;
        } //if
      }
      count1->nchanceinhell_1overr[whichray]++;
      
      
      // distance ray travels through ice.
      if (settings1->ROUGHNESS==0) {
        if (whichray==0) {
          rflength=interaction1->posnu.Distance(ray1->rfexit[2]);
        }
        if (whichray==1) {
          rflength=interaction1->posnu_down.Distance(ray1->rfexit[2]);//use the real distance that singals pass
        }
      }// if this is not for roughness
      //else
      //  rflength=ray1->rfexit[1].Mag()-interaction1->posnu.Mag();
      
      // take into account attenuation through the ice.
      // here,  vmmhz is input and output
      if (bn1->WHICHPATH==4)
        cout << "rflength is " << rflength << "\n";
      
      // applying ice attenuation factor
      if (settings1->ROUGHNESS==0) {
        if (whichray==0)
          Attenuate(antarctica, settings1, vmmhz_max,  rflength,  interaction1->posnu);
        if (whichray==1)
          Attenuate_down(antarctica, settings1, vmmhz_max,  ray1->rfexit[2],  interaction1->posnu, interaction1->posnu_down);
      }
      // roughness attenuation already dealt with

      // fill for just 1/10 of the events.
      if (tree2->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1 && bn1->WHICHPATH != 3)
        tree2->Fill();
      
      // intermediate counting
      count_dbexitsice++;
      
      // reject if the event is undetectable.
      //     if (vmmhz_max*heff_max*0.5*(bw/1.E6)<CHANCEINHELL_FACTOR*anita1->maxthreshold*Tools::dMin(VNOISE, settings1->NLAYERS) && !settings1->SKIPCUTS) {
      if (settings1->ROUGHNESS==0){
        if (settings1->CHANCEINHELL_FACTOR*vmmhz_max*heff_max*0.5*(anita1->bwmin/1.E6)<anita1->maxthreshold*anita1->VNOISE[0]/10. && !settings1->SKIPCUTS) {
          if (bn1->WHICHPATH==3)
            cout<<"Event is undetectable.  Leaving loop."<<endl;

          //delete ray1;
          //delete interaction1;
          continue;
        } //if
      }
      
      count1->nchanceinhell[whichray]++;
      
      // for plotting
      if (tree3->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1 && bn1->WHICHPATH != 3)
        tree3->Fill();
      
      // index for each antenna so you can use it to fill arrays
      count_rx=0;
      // keeps track of maximum voltage seen on either polarization of any antenna
      volts_rx_max=0;
      
      //std::cerr<<vmmhz_max<<std::endl;

      // Make a vector of V/m/MHz scaled by 1/r and attenuated.
      // Calculates Jaime's V/m/MHz at 1 m for each frequency
      // then multiplies by scale factor vmmhz_max/vmmhz1m_max
      // this will need to be improved once frequency-dependent
      // attenuation length is included.
      if (settings1->ROUGHNESS==0){
        if (settings1->FORSECKEL==1)
          sig1->SetNDepth(sig1->NICE); // for making array of signal vs. frequency,  viewangle
        
        sig1->GetVmMHz(vmmhz_max, vmmhz1m_max, pnu, anita1->freq, anita1->NOTCH_MIN, anita1->NOTCH_MAX, vmmhz, Anita::NFREQ); // here we get the array vmmhz by taking vmmhz1m_max (signal at lowest frequency bin) and
        //   vmmhz_max (signal at lowest frequency after applying 1/r factor and attenuation factor)
        // and making an array across frequency bins by putting in frequency dependence.
      }
      
      // For each frequency,  get the width of Cerenkov cone
      // and size of signal once position of viewing angle is taken into account
      
      // these variables are for energy reconstruction studies
      undogaintoheight_e=0;
      undogaintoheight_h=0;
      
      for (int k=0;k<4;k++) {
        undogaintoheight_e_array[k]=0.;
        undogaintoheight_h_array[k]=0.;
        nbins_array[k]=0;
        true_efield_array[k]=0.;
        rec_efield_array[k]=0.;
      }
      
      rec_efield=0;
      true_efield=0;
      
      
      if (settings1->ROUGHNESS==0){
        // don't loop over frequencies if the viewing angle is too far off
        double rtemp=Tools::dMin((viewangle-sig1->changle)/(deltheta_em_max), (viewangle-sig1->changle)/(deltheta_had_max));
        if (rtemp>Signal::VIEWANGLE_CUT && !settings1->SKIPCUTS) {
          //delete interaction1;
          continue;
        }
        count1->nviewanglecut[whichray]++;
        
        for (int k=0;k<Anita::NFREQ;k++) {
          //GetSpread(pnu, emfrac, hadfrac, freq[k], N_DEPTH, X0ICE, 
          //deltheta_em[k], deltheta_had[k]); // this gets the angular width of the em and had components of the cerenkov cone
          deltheta_em[k]=deltheta_em_max*anita1->FREQ_LOW/anita1->freq[k];
          deltheta_had[k]=deltheta_had_max*anita1->FREQ_LOW/anita1->freq[k];

          if (settings1->FORSECKEL==1) {// this is for making plots of the signal
            for (int iviewangle=0;iviewangle<NVIEWANGLE;iviewangle++) {// loop over viewing angles
              // remove the 1/r and attenuation factors that are contained in the ratio vmmhz1m_max/vmmhz_max
                    //
              vmmhz_temp=vmmhz[k]*vmmhz1m_max/vmmhz_max; // EH,  I guess temporally remove 1/r and attenuation factors,  not for general use?

              viewangle_temp=viewangles[iviewangle]; //grab the viewing angle from this array
              //apply the gaussian dependence away from the cerenkov angle.  vmmhz_temp is both an input and an output.
              //vmmhz_temp as an output is the signal with the angular dependence applied.
              sig1->TaperVmMHz(viewangle_temp, deltheta_em[k], deltheta_had[k], emfrac, hadfrac, vmmhz_temp, djunk);
              forseckel[iviewangle][k]=vmmhz_temp;// put this in an array which we will plot later.
            } //for (loop over viewing angles)
          } //if (settings1->FORSECKEL==1)

          //if (settings1->ROUGHNESS==0) {// we apply this normal angular dependence if we are not dealing with roughness here
          sig1->TaperVmMHz(viewangle, deltheta_em[k], deltheta_had[k], emfrac, hadfrac, vmmhz[k], vmmhz_em[k]);// this applies the angular dependence.
            // viewangle is which viewing angle we are at
            // deltheta_em is the width of the em component at this frequency
            // deltheta_had is the width of the had component at this frequency
            // emfrac is the em fraction of the shower
            // hadfrac is the hadronic fraction of the shower
            // vmmhz is the strength of the signal in V/m/MHz at this viewing angle
            // vmmhz_em is the strength of the em component

          // MS - moved this stuff to just before globaltrigger so it uses screen
          //Tools::Zero(sumsignal_aftertaper, 5);
          //if (bn1->WHICHPATH==4)
          //  IntegrateBands(anita1, k, panel1, anita1->freq, vmmhz1m_max/(vmmhz_max*1.E6), sumsignal_aftertaper);

          vmmhz_lowfreq=vmmhz[0]; // for plotting,  vmmhz at the lowest frequency
            // EH,  now here,  vmmhz array is the spectrum right infront of antennas (so not yet any detector properties are applied)
          //GetVmMHz(vmmhz1m_max, vmmhz1m_max, pnu, vmmhz); //Get vmmhz in frequency domain before the antenna response,  for modeling time dependent pulses

          // just want to see the maximum effect of viewing angle being off cerenkov cone
          // should be at highest frequency
          // just for plotting
          maxtaper=-1000;
          if (sig1->logscalefactor_taper>maxtaper)
            maxtaper=sig1->logscalefactor_taper;

          if (interaction1->nuflavor=="nue")        pdgcode = 12;
          else if (interaction1->nuflavor=="numu")  pdgcode = 14;
          else if (interaction1->nuflavor=="nutau") pdgcode = 16;

          if (settings1->HIST==1 && !settings1->ONLYFINAL && bn1->WHICHPATH != 3 && k==Anita::NFREQ/2 && tree18->GetEntries()<settings1->HIST_MAX_ENTRIES) {
            tree18->Fill();
          }

          if (bn1->WHICHPATH == 3)
            interaction1->banana_volts += vmmhz[k]*(settings1->BW/(double)Anita::NFREQ/1.E6);
        }//end for (int k=0;k<Anita::NFREQ;k++)

        if (bn1->WHICHPATH==3 && interaction1->banana_volts != 0 && settings1->HIST && banana_tree->GetEntries()<settings1->HIST_MAX_ENTRIES) {
          banana_tree->Fill();
          //delete ray1;
          //delete interaction1;
          continue;
        } //This is all the data needed for the banana plot - we now have the final value of vmmhz[]
        else if (bn1->WHICHPATH==3 && interaction1->banana_volts == 0) {
          //delete ray1;
          //delete interaction1;
          continue; //Exit the loop if there's no voltage here - no value at a point is the same as zero,  and this will save HD space
        }
        // reject if it is undetectable now that we have accounted for viewing angle
        
        //if (Tools::dMax(vmmhz, Anita::NFREQ)*heff_max*0.5*(bw/1.E6)<CHANCEINHELL_FACTOR*anita1->maxthreshold*Tools::dMin(VNOISE, settings1->NLAYERS) && !settings1->SKIPCUTS) {
        if (settings1->CHANCEINHELL_FACTOR*Tools::dMax(vmmhz, Anita::NFREQ)*heff_max*0.5*(anita1->bwmin/1.E6)<anita1->maxthreshold*anita1->VNOISE[0]/10. && !settings1->SKIPCUTS) {
          //delete interaction1;
          //delete ray1;
          continue;
        }
        
      }//end if roughness==0 before the Anita::NFREQ k loop, this isolates the TaperVmMHz()

      // just for plotting
      vmmhz_max=Tools::dMax(vmmhz, Anita::NFREQ);
      vmmhz_min=Tools::dMin(vmmhz, Anita::NFREQ);
      
      //std::cerr<<"++ "<<fSCREEN_NUMPOINTS_EDGE<<"   "<<panel1->GetEdgeLength()<<"   "<<vmmhz_min<<std::endl;
      
      // intermediate counting
      count1->nchanceinhell2[whichray]++;
      chanceinhell2=1;

      //if no-roughness case, add its parameters to the saved screen parameters
      if(!settings1->ROUGHNESS){
        panel1->SetNvalidPoints(1);
        for (int k=0;k<Anita::NFREQ;k++) {
          panel1->AddVmmhz_freq(vmmhz[k]);
          panel1->AddDelay( 0. );
        }
        panel1->AddVec2bln(ray1->n_exit2bn[2]);
        panel1->AddPol(n_pol);
        panel1->AddWeight( 1. );
        panel1->SetWeightNorm( 1. );

        for (int k=0;k<Anita::NFREQ;k++) {
          Tools::Zero(sumsignal_aftertaper, 5);
          if (bn1->WHICHPATH==4)
            IntegrateBands(anita1, k, panel1, anita1->freq, vmmhz1m_max/(vmmhz_max*1.E6), sumsignal_aftertaper);
        }
      }


      // make a global trigger object (but don't touch the electric fences)
      globaltrig1 = new GlobalTrigger(settings1, anita1);
      
      Tools::Zero(anita1->arrival_times, Anita::NLAYERS_MAX*Anita::NPHI_MAX);
      if(settings1->BORESIGHTS)
        anita1->GetArrivalTimesBoresights(ray1->n_exit2bn_eachboresight[2],bn1,settings1);
      else
        anita1->GetArrivalTimes(ray1->n_exit2bn[2],bn1,settings1);
      
      anita1->rx_minarrivaltime=Tools::WhichIsMin(anita1->arrival_times, settings1->NANTENNAS);

      //Zeroing
      for (int i=0;i<settings1->NANTENNAS;i++) {
        voltagearray[i]=0;
        discones_passing=0;
      } //Zero the trigger array
      
      max_antenna0=-1;
      max_antenna1=-1;
      max_antenna2=-1;
      max_antenna_volts0 =0;
      max_antenna_volts1 =0;
      max_antenna_volts2 =0;
      e_comp_max1 = 0;
      h_comp_max1 = 0;
      e_comp_max2 = 0;
      h_comp_max2 = 0;
      e_comp_max3 = 0;
      h_comp_max3 = 0;
      //End zeroing
      
        
      // start looping over antennnas.
      // ilayer loops through vertical layers
      
      if (settings1->SLAC)
        fslac_hitangles << bn1->sslacpositions[bn1->islacposition] << "\n";
      
      if (RANDOMISEPOL) {
        double rotateangle=gRandom->Gaus(RANDOMISEPOL*RADDEG);
        n_pol=n_pol.Rotate(rotateangle, ray1->n_exit2bn[2]);
      }
    
      if (bn1->WHICHPATH==4) {
        Tools::Zero(sumsignal, 5);
        for (int k=0;k<Anita::NFREQ;k++)
          IntegrateBands(anita1, k, panel1, anita1->freq, bn1->r_bn.Distance(interaction1->posnu)/1.E6, sumsignal);
      }//end if whichpath==4
      
      if (bn1->CENTER){
        //bn1->CenterPayload(settings1, anita1, ray1->n_exit2bn[2], n_pol, anita1->GetLayer(anita1->rx_minarrivaltime), anita1->GetIfold(anita1->rx_minarrivaltime));
        bn1->CenterPayload(hitangle_e);
      }

      if (ray1->MAKEVERTICAL) {
        //	cout << n_pol << " before and after "; 
        n_pol=bn1->n_bn;
        //	cout << n_pol << endl;
          // rotate n_exit2bn too
          // rotation axis n_bn crossed with n_exit2bn
        Vector rotationaxis=ray1->n_exit2bn[2].Cross(bn1->n_bn);
        double rotateangle=PI/2.-ray1->n_exit2bn[2].Dot(bn1->n_bn);
        ray1->n_exit2bn[2]=ray1->n_exit2bn[2].Rotate(rotateangle, rotationaxis);
        
        for (int ilayer=0; ilayer < settings1->NLAYERS; ilayer++) { // loop over layers on the payload
          // ifold loops over phi
          for (int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
            Vector rotationaxis2=ray1->n_exit2bn_eachboresight[2][ilayer][ifold].Cross(n_pol_eachboresight[ilayer][ifold]);
            double rotateangle2=PI/2.-ray1->n_exit2bn_eachboresight[2][ilayer][ifold].Dot(n_pol_eachboresight[ilayer][ifold]);
            ray1->n_exit2bn_eachboresight[2][ilayer][ifold].Rotate(rotateangle2, rotationaxis2);
          } // end loop over phi
        } // end loop over layers
      }//end if ray1->makevertical

      globaltrig1->volts_rx_rfcm_trigger.assign(16,  vector <vector <double> >(3,  vector <double>(0)));
      anita1->rms_rfcm_e_single_event = 0;

      for (int ilayer=0; ilayer < settings1->NLAYERS; ilayer++) { // loop over layers on the payload
        for (int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
          
          AntTrigger *anttrig1 = new AntTrigger();
          anttrig1->InitializeEachBand(anita1);

          bn1->GetAntennaOrientation(settings1,  anita1,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);


          // for this (hitangle_h_all[count_rx]=hitangle_h;) and histogram fill, use specular case
          //although the GetEcomp..() functions are called in ConvertInputWFtoAntennaWF() to calculate the actual waveforms
          if (!settings1->BORESIGHTS) {
            bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  ray1->n_exit2bn[2], e_component_kvector,  h_component_kvector,  n_component_kvector);
            bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  n_pol,  e_component,  h_component,  n_component);
          }
          else{ // i.e. if BORESIGHTS is true
            bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  ray1->n_exit2bn_eachboresight[2][ilayer][ifold],  e_component_kvector,  h_component_kvector,  n_component_kvector);
            bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  n_pol_eachboresight[ilayer][ifold], e_component,  h_component,  n_component);
            fslac_hitangles << ilayer << "\t" << ifold << "\t" << hitangle_e << "\t" << hitangle_h << "\t" << e_component_kvector << "\t" << h_component_kvector << "\t" << fresnel1_eachboresight[ilayer][ifold] << " " << mag1_eachboresight[ilayer][ifold] << "\n";
          }
          bn1->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);
          // store hitangles for plotting
          hitangle_h_all[count_rx]=hitangle_h;
          hitangle_e_all[count_rx]=hitangle_e;
          // for debugging
          if (h6->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
            h6->Fill(hitangle_h, ray1->n_exit2bn[2]*bn1->n_bn);


          anttrig1->ConvertInputWFtoAntennaWF(settings1, anita1, bn1, panel1, n_eplane,  n_hplane,  n_normal);

          // AntTrig::ImpulseResponse needs to be outside the ray loop
          anttrig1->ImpulseResponse(settings1, anita1, ilayer, ifold);
          
          anttrig1->TimeShiftAndSignalFluct(settings1, anita1, ilayer, ifold, volts_rx_rfcm_lab_e_all,  volts_rx_rfcm_lab_h_all);


          //+++++//+++++//+++++//+++++//+++++//+++++//+++++

          // THIS IS WHERE WE ACTUALLY CONSTRUCT THE WAVEFORMS THAT GET PASSED TO THE TRIGGER
          //anttrig1->Banding(settings1,  anita1, vmmhz);

          //temporary containers for each screen point to construct the summed waveforms
          double TMP_v_banding_rfcm_e_forfft[5][Anita::HALFNFOUR]={0.};
          double TMP_v_banding_rfcm_h_forfft[5][Anita::HALFNFOUR]={0.};

          Tools::Zero(sumsignal, 5);
          double integrate_energy_freq[5]={0.,0.,0.,0.,0.};

          for (int iband=0;iband<5;iband++) { // loop over bands
            
            anita1->iminbin[iband]=0.;
            anita1->imaxbin[iband]=anita1->NFOUR/2;

            for (int i=0;i<anita1->NFOUR/2;i++) {
              anttrig1->v_banding_rfcm_e_forfft[iband][i] = 0.;
              anttrig1->v_banding_rfcm_h_forfft[iband][i] = 0.;
            }
            
            for (int jpt=0; jpt<panel1->GetNvalidPoints(); jpt++){
              //get the orientation for this screen point
              bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal, panel1->GetVec2bln(jpt), e_component_kvector,  h_component_kvector,  n_component_kvector);
              bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  panel1->GetPol(jpt),  e_component,  h_component,  n_component);
              bn1->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);

              for (int i=0;i<Anita::NFREQ;i++) {
                anita1->vmmhz_banding[i]=panel1->GetVmmhz_freq(jpt*Anita::NFREQ + i);//; // now copy vmmhz to vmmhz_bak instead, which we now play with to get the time domain waveforms for each subband
                // remember vmmhz is V/m/MHz at the face of the antenna
              }
              // Don't we need to apply antenna gains here?
                
              // impose banding on the incident signal
              anita1->Banding(iband,anita1->freq,anita1->vmmhz_banding,Anita::NFREQ); // impose banding whatever the trigger scheme
                
              for (int i=0;i<Anita::NFREQ;i++) {
                anita1->vmmhz_banding_rfcm[i]=anita1->vmmhz_banding[i];
              }
                
              // for frequency-domain voltage-based trigger (triggerscheme==0)
              // we do not apply rfcm's
              // for other trigger types we do
              if (settings1->TRIGGERSCHEME==1 || settings1->TRIGGERSCHEME==2 || settings1->TRIGGERSCHEME == 3 || settings1->TRIGGERSCHEME == 4 || settings1->TRIGGERSCHEME == 5)
                anita1->RFCMs(1,1,anita1->vmmhz_banding_rfcm);
                
                
              if (settings1->TRIGGERSCHEME >=2) { // we need to prepar the signal for the diode integration
                for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {

                  anita1->vmmhz_banding_rfcm[ifreq]=anita1->vmmhz_banding_rfcm[ifreq]/sqrt(2.)/(anita1->TIMESTEP*1.E6);
                  // vmmhz was set to account for both negative and positive frequencies
                  // now it has units of volts/(meter*s) so below we copy it to vm_banding_rfcm_e,h
                    
                  // here let's treat it as just one side of a double sided function though
                  // divide by timestep so that we can use it for discrete fourier transform.
                  // as in numerical recipes, F = (Delta t) * H
                }
              }
              double integral=0.;
              for (int ifreq=0;ifreq<Anita::NFREQ;ifreq++) {
                anttrig1->vm_banding_rfcm_e[iband][ifreq]=anita1->vmmhz_banding_rfcm[ifreq]; // this is now Volts/(m*s) vs. frequency with banding and rfcm's applied
                anttrig1->vm_banding_rfcm_h[iband][ifreq]=anita1->vmmhz_banding_rfcm[ifreq];
                integral += anttrig1->vm_banding_rfcm_e[iband][ifreq]*anttrig1->vm_banding_rfcm_e[iband][ifreq];
              } // end loop over nfreq

              for (int k=0;k<Anita::NFREQ;k++) {
                if (anita1->freq[k]>=settings1->FREQ_LOW_SEAVEYS && anita1->freq[k]<=settings1->FREQ_HIGH_SEAVEYS){
                  // need to calculate lcp and rcp components after antenna voltages are recorded.
                    anttrig1->v_banding_rfcm_e[iband][k]=anttrig1->vm_banding_rfcm_e[iband][k];
                    anttrig1->v_banding_rfcm_h[iband][k]=anttrig1->vm_banding_rfcm_h[iband][k];

                    anita1->AntennaGain(settings1, hitangle_e, hitangle_h, e_component, h_component, k, anttrig1->v_banding_rfcm_e[iband][k], anttrig1->v_banding_rfcm_h[iband][k]);
                    if (bn1->WHICHPATH==4 && anita1->Match(ilayer, ifold, anita1->rx_minarrivaltime))
                      Integrate(anita1, iband, k, anttrig1->v_banding_rfcm_e[iband], anita1->freq, 1./1.E6, sumsignal[iband]);
                    anttrig1->addToChannelSums(settings1, anita1, iband, k);
                } // end if (seavey frequencies)
              } // end looping over frequencies.

              Tools::Zero(TMP_v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
              Tools::Zero(TMP_v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
            
              anita1->MakeArraysforFFT(anttrig1->v_banding_rfcm_e[iband], anttrig1->v_banding_rfcm_h[iband], TMP_v_banding_rfcm_e_forfft[iband], TMP_v_banding_rfcm_h_forfft[iband], 0., false);
            
              // for some reason I'm averaging over 10 neighboring bins
              // to get rid of the zero bins
              for (int i=0;i<anita1->NFOUR/4;i++) {
                anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i]  =0.;
                anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i+1]=0.;
                anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i]  =0.;
                anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i+1]=0.;

                int tempcount = 0;
                for (int k=i;k<i+10;k++) {
                  if (k<anita1->NFOUR/4) {
                    anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i]  += TMP_v_banding_rfcm_e_forfft[iband][2*k];
                    anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i+1]+= TMP_v_banding_rfcm_e_forfft[iband][2*k+1];
                    anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i]  += TMP_v_banding_rfcm_h_forfft[iband][2*k];
                    anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i+1]+= TMP_v_banding_rfcm_h_forfft[iband][2*k+1];
                    tempcount++;
                  }
                }

                TMP_v_banding_rfcm_e_forfft[iband][2*i]  = anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i]/tempcount;
                TMP_v_banding_rfcm_e_forfft[iband][2*i+1]= anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i+1]/tempcount;
                TMP_v_banding_rfcm_h_forfft[iband][2*i]  = anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i]/tempcount;
                TMP_v_banding_rfcm_h_forfft[iband][2*i+1]= anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i+1]/tempcount;

                anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i]  = TMP_v_banding_rfcm_e_forfft[iband][2*i];
                anttrig1->v_banding_rfcm_e_forfft_temp[iband][2*i+1]= TMP_v_banding_rfcm_e_forfft[iband][2*i+1];
                anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i]  = TMP_v_banding_rfcm_h_forfft[iband][2*i];
                anttrig1->v_banding_rfcm_h_forfft_temp[iband][2*i+1]= TMP_v_banding_rfcm_h_forfft[iband][2*i+1];
              }

              for (int ifour=0;ifour<Anita::NFOUR/4;ifour++) {
                TMP_v_banding_rfcm_e_forfft[iband][2*ifour]*=cos(panel1->GetDelay(jpt)*PI/180.);
                TMP_v_banding_rfcm_e_forfft[iband][2*ifour+1]*=sin(panel1->GetDelay(jpt)*PI/180.);
                TMP_v_banding_rfcm_h_forfft[iband][2*ifour]*=cos(panel1->GetDelay(jpt)*PI/180.);
                TMP_v_banding_rfcm_h_forfft[iband][2*ifour+1]*=sin(panel1->GetDelay(jpt)*PI/180.);  
              }

              Tools::realft(TMP_v_banding_rfcm_e_forfft[iband],-1,anita1->NFOUR/2);
              // now v_banding_rfcm_e_forfft is in the time domain
              // and now it is really in units of V

              Tools::realft(TMP_v_banding_rfcm_h_forfft[iband],-1,anita1->NFOUR/2);
              // now v_banding_rfcm_h_forfft is in the time domain
              // and now it is really in units of V
            
              // put it in normal time ording -T to T
              // instead of 0 to T, -T to 0
              Tools::NormalTimeOrdering(anita1->NFOUR/2,TMP_v_banding_rfcm_e_forfft[iband]);
              Tools::NormalTimeOrdering(anita1->NFOUR/2,TMP_v_banding_rfcm_h_forfft[iband]);

              if (settings1->ZEROSIGNAL) {
                Tools::Zero(TMP_v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
                Tools::Zero(TMP_v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
              }

              for (int i=0;i<anita1->NFOUR/4;i++) {
                integrate_energy_freq[iband] += TMP_v_banding_rfcm_e_forfft[iband][2*i]*TMP_v_banding_rfcm_e_forfft[iband][2*i]
                 + TMP_v_banding_rfcm_h_forfft[iband][2*i+1]*TMP_v_banding_rfcm_h_forfft[iband][2*i+1];
              }

              // now add the screen point's waveform to the total including the weighting and phase delays
              for (int i=0;i<anita1->NFOUR/2;i++) {
                anttrig1->v_banding_rfcm_e_forfft[iband][i] += TMP_v_banding_rfcm_e_forfft[iband][i] * panel1->GetWeight(jpt) * panel1->GetWeightNorm();
                anttrig1->v_banding_rfcm_h_forfft[iband][i] += TMP_v_banding_rfcm_h_forfft[iband][i] * panel1->GetWeight(jpt) * panel1->GetWeightNorm();
              }

            } // end loop over screen points
          } // end loop over bands
          
          // now hopefully we have converted the signal to time domain waveforms
          // for all the bands of the antenna and screen points


          // write the signal events to a tree
          for (int iband=0;iband<5;iband++) {
            for (int k=0;k<anita1->NFOUR/2;k++) {
              anita1->signal_vpol_inanita[iband][k]=anttrig1->v_banding_rfcm_e_forfft[iband][k];
            }
          }
          anita1->integral_vmmhz_foranita=anttrig1->integral_vmmhz;
        
          // Find the p2p value before adding noise
          for (int iband=0;iband<5;iband++) {
            anita1->peak_v_banding_rfcm_e[iband]=anttrig1->FindPeak(anttrig1->v_banding_rfcm_e_forfft[iband],anita1->NFOUR/2);
            anita1->peak_v_banding_rfcm_h[iband]=anttrig1->FindPeak(anttrig1->v_banding_rfcm_h_forfft[iband],anita1->NFOUR/2);
          }
/*
  ofstream props;
  std::string stemp;
  if(settings1->ROUGHNESS>0)
    stemp = settings1->outputdir+"/rough"+to_string(settings1->ROUGHSIZE)+"_"+to_string(ilayer)+"_"+to_string(ifold)+".dat";
  else
    stemp = settings1->outputdir+"/rough"+to_string(settings1->ROUGHNESS)+"_"+to_string(ilayer)+"_"+to_string(ifold)+".dat";
  props.open(stemp);
  for (int iband=0;iband<5;iband++) {
    for (int k=0;k<anita1->NFOUR/2;k++) {
      props << iband <<"  "<<k<<"  "<< anita1->signal_vpol_inanita[iband][k]<<"  "<<anttrig1->v_banding_rfcm_h_forfft[iband][k]<<std::endl;
    }
  }
  props.close();*/

          if (bn1->WHICHPATH==4 && ilayer==anita1->GetLayer(anita1->rx_minarrivaltime) && ifold==anita1->GetIfold(anita1->rx_minarrivaltime)) {
            for (int ibw=0;ibw<5;ibw++) {
              cout << "Just after Taper,  sumsignal is " << sumsignal_aftertaper[ibw] << "\n";
              cout << "Just after antennagain,  sumsignal is " << sumsignal[ibw] << "\n";
            }
          }

          // for energy reconstruction studies
          if (count_rx==anita1->rx_minarrivaltime) {
            undogaintoheight_e/=(double)Anita::NFREQ;
            undogaintoheight_h/=(double)Anita::NFREQ;
            for (int k=0;k<4;k++) {
              undogaintoheight_e_array[k]/=(double)nbins_array[k];
              undogaintoheight_h_array[k]/=(double)nbins_array[k];
            }
          }
          if (settings1->SCALEDOWNLCPRX1)
            globaltrig1->volts[0][ilayer][0]=globaltrig1->volts[0][ilayer][0]/sqrt(2.);
          if (settings1->RCPRX2ZERO)
            globaltrig1->volts[1][ilayer][1]=0.;
          if (settings1->LCPRX2ZERO)
            globaltrig1->volts[0][ilayer][1]=0.;
          if (settings1->SIGNAL_FLUCT) {
            if (settings1->WHICH==0) {
              globaltrig1->volts[ilayer][ifold][0]+=gRandom->Gaus(0., anita1->VNOISE_ANITALITE[ifold]);
              globaltrig1->volts[ilayer][ifold][1]+=gRandom->Gaus(0., anita1->VNOISE_ANITALITE[ifold]);
            } //else
          } //if adding noise
          if (count_rx==anita1->rx_minarrivaltime) {
            rec_efield=sqrt(pow(globaltrig1->volts_original[0][ilayer][ifold]/(undogaintoheight_e*0.5), 2)+pow(globaltrig1->volts_original[1][ilayer][ifold]/(undogaintoheight_h*0.5), 2));
            for (int ibw=0;ibw<4;ibw++) {
              rec_efield_array[ibw]=sqrt(pow(anttrig1->bwslice_volts_pole[ibw]/(undogaintoheight_e_array[ibw]*0.5), 2)+pow(anttrig1->bwslice_volts_polh[ibw]/(undogaintoheight_h_array[ibw]*0.5), 2));
              bwslice_vnoise_thislayer[ibw]=anita1->bwslice_vnoise[ilayer][ibw];// this is just for filling into a tree
            } // end loop over bandwidth slices
            tree6b->Fill();
          } // end if this is the closest antenna
          // for plotting
          if (ilayer==0 && ifold==0) {
            volts_rx_0=globaltrig1->volts[0][ilayer][ifold];
            if (settings1->SIGNAL_FLUCT)
              volts_rx_0+=gRandom->Gaus(0, anita1->VNOISE[ilayer]);
          } //if (first antenna,  top layer)

          // for debugging
          if (volts_rx_0>volts_rx_max) {
            volts_rx_max=volts_rx_0;
            volts_rx_max_highband=anttrig1->bwslice_volts_pol0[3];
            volts_rx_max_lowband=anttrig1->bwslice_volts_pol0[0];
            // theta of the polarization as measured at the antenna (approximately since we aren't correcting for the
            //cant of the antenna yet) =
            theta_pol_measured=atan(globaltrig1->volts_original[1][ilayer][ifold]/globaltrig1->volts_original[0][ilayer][ifold]);
          }
          // for plotting
          volts_rx_1=globaltrig1->volts[1][ilayer][ifold];


          //+++++//+++++//+++++//+++++//+++++//+++++//+++++

          anttrig1->WhichBandsPass(inu,settings1, anita1, globaltrig1, bn1, ilayer, ifold,  viewangle-sig1->changle, emfrac, hadfrac);









          if (Anita::GetAntennaNumber(ilayer, ifold)==anita1->rx_minarrivaltime) {
            for (int iband=0;iband<5;iband++) {
              for (int ipol=0;ipol<2;ipol++) {
                rx0_signal_eachband[ipol][iband]=anttrig1->signal_eachband[ipol][iband];
                rx0_threshold_eachband[ipol][iband]=anttrig1->threshold_eachband[ipol][iband];
                rx0_noise_eachband[ipol][iband]=anttrig1->noise_eachband[ipol][iband];
                rx0_passes_eachband[ipol][iband]=anttrig1->passes_eachband[ipol][iband];
              }
            }
          }

          //For verification plots: find antenna with max signal - added by Stephen
          if (ilayer == 0 && globaltrig1->volts[0][ilayer][ifold] > max_antenna_volts0) {
            max_antenna0 = count_rx;
            max_antenna_volts0 = globaltrig1->volts[0][ilayer][ifold];
            max_antenna_volts0_em=globaltrig1->volts_em[0][ilayer][ifold];
            ant_max_normal0 = ant_normal;
            e_comp_max1 = e_component;
            h_comp_max1 = h_component;
          }
          else if (ilayer == 0 && globaltrig1->volts[0][ilayer][ifold] == max_antenna_volts0 && globaltrig1->volts[0][ilayer][ifold] != 0){
            cout<<"Equal voltage on two antennas!  Event : "<<inu<<endl;
          }
          else if (ilayer == 1 && globaltrig1->volts[0][ilayer][ifold] > max_antenna_volts1) {
            max_antenna1 = count_rx;
            max_antenna_volts1 = globaltrig1->volts[0][ilayer][ifold];
            ant_max_normal1 = ant_normal;
            e_comp_max2 = e_component;
            h_comp_max2 = h_component;
          }
          else if (ilayer == 1 && globaltrig1->volts[0][ilayer][ifold] == max_antenna_volts1 && globaltrig1->volts[0][ilayer][ifold] != 0){
            cout<<"Equal voltage on two antennas!  Event : "<<inu<<endl;
          }
          else if (ilayer == 2 && globaltrig1->volts[0][ilayer][ifold] > max_antenna_volts2) {
            max_antenna2 = count_rx;
            max_antenna_volts2 = globaltrig1->volts[0][ilayer][ifold];
            ant_max_normal2 = ant_normal;
            e_comp_max3 = e_component;
            h_comp_max3 = h_component;
          }
          else if (ilayer == 2 && globaltrig1->volts[0][ilayer][ifold] == max_antenna_volts2 && globaltrig1->volts[0][ilayer][ifold] != 0){
            cout<<"Equal voltage on two antennas!  Event : "<<inu<<endl;
          }
          voltagearray[count_rx] = globaltrig1->volts[0][ilayer][ifold];
          //End verification plot block
          
          count_rx++; // counting antennas that we loop through,  for indexing
          
          
          if (settings1->TRIGTYPE==0 && ifold==1 && count_pass>=settings1->NFOLD) { //added djg --line below fills "direct" voltage output file
            al_voltages_direct<<"0 0 0"<<"   "<<"    "<<globaltrig1->volts_original[1][0][0]<<"    "<<(globaltrig1->volts_original[0][0][0]/sqrt(2.))<<"     "<<globaltrig1->volts_original[1][0][1]<<"     "<<globaltrig1->volts_original[0][0][1]<<"      "<<anita1->VNOISE[0]<<"     "<<anita1->VNOISE[0]<<"     "<<anita1->VNOISE[0]<<"     "<<anita1->VNOISE[0]<<"  "<<weight<<endl;
          }
          delete anttrig1;
        } //loop through the phi-fold antennas
      }  //loop through the layers of antennas


      anita1->rms_rfcm_e_single_event = sqrt(anita1->rms_rfcm_e_single_event / (anita1->HALFNFOUR * settings1->NANTENNAS));

      if (settings1->DISCONES==1) {
        // loop through discones
        for (int idiscone=0;NDISCONES;idiscone++) {
          AntTrigger *anttrig1=new AntTrigger();
          volts_discone=0.;
          polarfactor_discone=n_pol.Dot(bn1->n_bn); // beam pattern
          for (int k=0;k<Anita::NFREQ;k++) {
            if (anita1->freq[k]>=FREQ_LOW_DISCONES && anita1->freq[k]<=FREQ_HIGH_DISCONES) {
              thislambda=CLIGHT/sig1->N_AIR/anita1->freq[k];
              heff_discone= thislambda*sqrt(2*Zr*gain_dipole/Z0/4/PI*sig1->N_AIR);   // effective height of dipole,  using formula from Ped's note
              
              volts_discone+=vmmhz[k]*0.5*heff_discone*((settings1->BW/1E6)/(double)Anita::NFREQ)*polarfactor_discone;
            }
          }// end for k loop
          
          vnoise_discone=anita1->VNOISE[0]*sqrt(BW_DISCONES/settings1->BW_SEAVEYS);

          if (settings1->SIGNAL_FLUCT) {
            volts_discone+=gRandom->Gaus(0., vnoise_discone); // here I'm using the noise seen by an antenna pointed with a 10 degree cant.  Should be different for a discone but we'll change it later.
          }
          
          if (fabs(volts_discone)/vnoise_discone>anita1->maxthreshold)
            discones_passing++;
          
          delete anttrig1;
        } // end looping through discones
      } //end if settings discones==1
  
      for (int irx=0;irx<settings1->NANTENNAS;irx++) {
        nchannels_perrx_triggered[irx]=globaltrig1->nchannels_perrx_triggered[irx];
      }
    
      nchannels_triggered=Tools::iSum(globaltrig1->nchannels_perrx_triggered, settings1->NANTENNAS); // find total number of antennas that were triggered.
      volts_rx_ave=GetAverageVoltageFromAntennasHit(settings1, globaltrig1->nchannels_perrx_triggered, voltagearray, volts_rx_sum);
    
      // if it passes the trigger,  then go ahead and
      // calculate the chord length,  etc.
      // intermediate counter
      if(sec1->secondbang && sec1->interestedintaus)
        count_asktrigger_nfb++;  // just for taus
      else
        count_asktrigger++;
      dist_int_bn_2d_chord = ray1->rfexit[0].Distance(bn1->r_bn_shadow)/1000; // for sensitivity vs. distance plot
      //distance across the surface - Stephen
      dist_int_bn_2d = ray1->rfexit[0].SurfaceDistance(bn1->r_bn_shadow, bn1->surface_under_balloon) / 1000;
      //---------------------
      //just added this temporarily - will make it run slower
      //this gets the weight due to stopping in earth
      //returns 0 if chord<1m
      
    
      if (!antarctica->Getchord(settings1, len_int_kgm2, interaction1->r_in, interaction1->r_enterice, interaction1->nuexitice, interaction1->posnu, inu, interaction1->chord, interaction1->weight_nu_prob, interaction1->weight_nu, nearthlayers, myair, total_kgm2, crust_entered,  mantle_entered, core_entered)){
        interaction1->weight_nu_prob = -1.;
      }
      
      if(tauweighttrigger==1)
        weight1=interaction1->weight_nu_prob + taus1->weight_tau_prob;
      else
        weight1=interaction1->weight_nu_prob;
      
    
      weight = weight1 / interaction1->dnutries * settings1->SIGMA_FACTOR;  // total weight is the earth absorption factor
      // divided by the factor accounting for the fact that we only chose our interaction point within the horizon of the balloon
      // then multiply by the cross section multiplier,  to account for the fact that we get more interactions when the cross section is higher
      if (weight<CUTONWEIGHTS) {
        delete globaltrig1;
        continue;
      }

      eventsfound_beforetrigger+=weight;


      //////////////////////////////////////
      //       EVALUATE GLOBAL TRIGGER    //
      //          FOR VPOL AND HPOL       //
      //////////////////////////////////////
      

      int thispasses[Anita::NPOL]={0,0};
      //cout << "Event number " << inu << "\n";

      globaltrig1->PassesTrigger(settings1, anita1, discones_passing, 2, l3trig, l2trig, l1trig, settings1->antennaclump, loctrig, loctrig_nadironly, inu, 
				 thispasses);
      
      for (int i=0;i<2;i++) {
	for (int j=0;j<16;j++) {
	  for (int k=0;k<anita1->HALFNFOUR;k++) {
	    count1->nl1triggers[i][whichray]+=anita1->l1trig_anita3and4_inanita[i][j][k];
	  }
	}
      }
      ///////////////////////////////////////
      //       Require that it passes      //
      //            global trigger         //
      ///////////////////////////////////////
      // for Anita-lite,  Anita Hill, just L1 requirement on 2 antennas. This option is currently disabled
      // Save events that generate an RF trigger or that are part of the min bias sample
      // Minimum bias sample: save all events that we could see at the payload
      // Independentely from the fact that they generated an RF trigger

      if ( (thispasses[0]==1 && anita1->pol_allowed[0]==1)  
	   || (thispasses[1]==1 && anita1->pol_allowed[1]==1)
	   || (settings1->TRIGTYPE==0 && count_pass>=settings1->NFOLD)
	   || (settings1->MINBIAS==1)) {


        if (bn1->WHICHPATH==4)
          cout << "This event passes.\n";

        anita1->passglobtrig[0]=thispasses[0];
        anita1->passglobtrig[1]=thispasses[1];


        //calculate the phi angle wrt +x axis of the ray from exit to balloon
        n_exit_phi = Tools::AbbyPhiCalc(ray1->n_exit2bn[2][0], ray1->n_exit2bn[2][1]);

        // keep track of events passing trigger
        count1->npassestrigger[whichray]++;
        // tags this event as passing
        passestrigger=1;
	//	cout << "This event passes.\n";
  
        // for plotting
        if (tree11->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
          tree11->Fill();
        
        // for taus
        if(sec1->secondbang && sec1->interestedintaus)
          count_passestrigger_nfb++;
        
        crust_entered=0; //These are switches that let us tell how far a given neutrino penetrated.  Clear them before entering Getchord.
        mantle_entered=0;
        core_entered=0;

        // this gets the weight due to stopping in earth
        // returns 0 if chord<1m
        //cout << "inu is " << inu << "\n";
        if (tautrigger==1 || antarctica->Getchord(settings1, len_int_kgm2, interaction1->r_in, interaction1->r_enterice, interaction1->nuexitice, interaction1->posnu, inu, interaction1->chord, interaction1->weight_nu_prob, interaction1->weight_nu, nearthlayers, myair, total_kgm2, crust_entered, mantle_entered, core_entered)) {
          //cout << "passes chord.\n";
          if (nupathtree->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
            nupathtree->Fill();

          // counts how many have a good chord length
          count_chordgoodlength++;

          // divide phase space factor into weight1
          if(tauweighttrigger==1){
            weight_prob=interaction1->weight_nu_prob + taus1->weight_tau_prob;
          }
          else
            weight_prob=interaction1->weight_nu_prob;
          
          weight1=interaction1->weight_nu;
          weight=weight1/interaction1->dnutries*settings1->SIGMA_FACTOR;
          weight_prob=weight_prob/interaction1->dnutries*settings1->SIGMA_FACTOR;

          pieceofkm2sr=weight*antarctica->volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr/(double)NNU/len_int;
          if (h10->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST)
            h10->Fill(hitangle_e_all[0], weight);


////////////////////////////////////////////////////////////////////////
// this tau code hasn't been touched in a while so comment it out
//          // (correct me if I'm wrong) -AC
//          // //for taus
//          if (nuflavor=="nutau" && current=="cc") {
//            // //If it was a tay neutrino where the initial interaction point was the neutrino interaction,  look for the tau decay.
//            if(!secondbang && interestedintaus) {
//              if(IsDoubleBang(ray1, pnu, elast_y, elpm, hadfrac_db, emfrac_db, taudecay, theta_bn, phi_bn, bn1->phi_spin, interaction1->posnu, interaction1->posnu2, interaction1->nnu, deltheta_em_db, deltheta_had_db, vmmhz_db, volts_db, gain, freq, flare, n_pol_db, nuexit, ray1->rfexit_db, bn1->r_bn, n_exit2bn_db, nrf_iceside_db, altitude_int, ray1->nsurf_rfexit_db, ray1->nrf_iceside, ptau, deltheta_em_max_db, deltheta_had_max_db, vmmhz_max_db, vmmhz1m_max_db, volts_rx_db, volts_rx_max_db, viewangle_db, nuexitlength, r_fromballoon_db, ddistance_db, theta_db, phi_db, costheta_inc_db, costheta_exit_db, dtemp, diff_angle_db, e_component_db, h_component_db, sumfrac_db, hitangle_e_db, hitangle_h_db, sig1)) {
//                  if (tree15->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && HIST==1)
//                  tree15->Fill();
//
//                count_passestrigger_db++;
//                eventsfound_db+=weight_doublebang*weight1;
//              }
//              if (tree5->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && HIST==1) {
//                tree5->Fill();
//                tree14->Fill();
//              }
//            }
//          } // end if taus
//          //
//          if(secondbang && interestedintaus) {
//            if(interaction1->posnu2.Mag()<(antarctica->RockSurface(interaction1->posnu2)))
//              count_interacts_in_earth++;
//            
//            weight_doublebang=NFBWeight(ptau, taulength);
//            weight1=2*(1-weight1)*weight_doublebang*icethickness/(TAUDECAY_TIME*CLIGHT*ptau/MTAU);
//          }
//          else{
//            if((nuflavor=="nutau" && current=="cc") && taudecay!="m" && interestedintaus)
//              weight1*=2;
//          }
////////////////////////////////////////////////////////////////////////





          // log of weight and chord for plotting
          logweight=log10(weight);
          interaction1->logchord=log10(interaction1->chord);

          // if neutrino travels more than one meter in ice
          if (interaction1->d2>1) {
            // intermediate counter
            count_d2goodlength++;

            // for taus
            // add to tally of neutrinos found,  weighted.
            if(sec1->secondbang && sec1->interestedintaus) {
              eventsfound_nfb+=weight;
              index_weights=(int)(((logweight-MIN_LOGWEIGHT)/(MAX_LOGWEIGHT-MIN_LOGWEIGHT))*(double)NBINS);
              eventsfound_nfb_binned[index_weights]++;
              if (tree16->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
                tree16->Fill();
            }//end if secondbang & interestedintaus
            else {
              allcuts[whichray]++;
              allcuts_weighted[whichray]+=weight;
	      if (thispasses[0] && thispasses[1]) {
		allcuts_weighted_polarization[2]+=weight;
	      } else if (thispasses[0]){
		allcuts_weighted_polarization[0]+=weight;
	      } else if (thispasses[1]){
		allcuts_weighted_polarization[1]+=weight;
	      }
	      anita1->weight_inanita=weight;	      

              if (h1mybeta->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1)
                h1mybeta -> Fill(mybeta, weight); //get the angle distribution of mybeta
              
              eventsfound+=weight; // counting events that pass,  weighted.
              eventsfound_prob+=weight_prob; // counting events that pass,  probabilities.
              if (cosalpha>0)
                eventsfound_belowhorizon+=weight;
              
              count1->npass[whichray]++;	// counting events that pass,  unweighted.
              // for calculating errors on sensitivity
              // need to find how many events as a function of weight
              // here,  we find how to index weight
              if (logweight<MIN_LOGWEIGHT)  // underflows,  set to 0th bin
                index_weights=0;
              else if (logweight>MAX_LOGWEIGHT) // overflows,  set to last bin
                index_weights=NBINS-1;
              else // which index weight corresponds to.
                index_weights=(int)(((logweight-MIN_LOGWEIGHT)/(MAX_LOGWEIGHT-MIN_LOGWEIGHT))*(double)NBINS);

              // count number of events that pass,  binned in weight
              if (index_weights<NBINS)
                eventsfound_binned[index_weights]++;
              
              // number of events in a ring at distance from balloon
              if (index_distance<NBINS_DISTANCE)
                eventsfound_binned_distance[index_distance]+= weight;
              
              // same,  now binned in weight,  for calculating errors
              if (index_distance<NBINS_DISTANCE && index_weights<NBINS)
                eventsfound_binned_distance_forerror[index_distance][index_weights]++;
              // for debugging
              if (logweight>-3)
                eventsfound_weightgt01+=weight;
              
              // how many events just pass through crust,  for same purpose.
              if (nearthlayers==1)
                eventsfound_crust+=weight;
              
              if (h1mybeta->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && settings1->HIST==1) {
                h1mybeta -> Fill(mybeta, weight);
                h1mytheta -> Fill(mytheta, weight);//fill mytheta
              }
            }//end else secondbang & interestedintaus

            //for plotting events distribution map only
            if(weight>0.0001){
              double int_lon, int_lat;
              //here are the longitude and altitude which Amy defined
              int_lon = interaction1->posnu.Lon(); // what latitude,  longitude does interaction occur at
              int_lat = interaction1->posnu.Lat();
              //cout<<"int_lon="<<int_lon<<"  int_lat="<<int_lat<<endl;
              int event_e_coord=0, event_n_coord=0;
              float event_e, event_n;
              antarctica->IceLonLattoEN(int_lon, int_lat, event_e_coord, event_n_coord);
              event_e=float(antarctica->xLowerLeft_ice+event_e_coord*antarctica->cellSize)/1000.;
              event_n=float(-1*(antarctica->yLowerLeft_ice+(antarctica->cellSize*event_n_coord)))/1000.;
              if(whichray==0)//direct
                dir_int_coord->Fill(event_e, event_n);
              if(whichray==1)
                ref_int_coord->Fill(event_e, event_n);
            }

            // just for plotting.
            offaxis=(double)fabs(viewangle-sig1->changle);
            nsigma_offaxis=offaxis/deltheta_had_max;
            	      
            hundogaintoheight_e->Fill(undogaintoheight_e, weight);
            hundogaintoheight_h->Fill(undogaintoheight_h, weight);
            rec_diff->Fill((rec_efield-true_efield)/true_efield, weight);
            rec_diff0->Fill((rec_efield_array[0]-true_efield_array[0])/true_efield_array[0], weight);
            rec_diff1->Fill((rec_efield_array[1]-true_efield_array[1])/true_efield_array[1], weight);
            rec_diff2->Fill((rec_efield_array[2]-true_efield_array[2])/true_efield_array[2], weight);
            rec_diff3->Fill((rec_efield_array[3]-true_efield_array[3])/true_efield_array[3], weight);
            recsum_diff->Fill((rec_efield_array[0]+rec_efield_array[1]+rec_efield_array[2]+rec_efield_array[3]-true_efield)/true_efield, weight);


	    sourceLon = ray1->rfexit[2].Lon() - 180;
	    sourceLat = ray1->rfexit[2].Lat() - 90;
	    sourceAlt = antarctica->SurfaceAboveGeoid(sourceLon+180, sourceLat+90);
	      
	    //Now put data in Vectors and Positions into arrays for output to the ROOT file.
            if (settings1->HIST && finaltree->GetEntries()<settings1->HIST_MAX_ENTRIES) {
              for (int i=0;i<3;i++) {
                nnu_array[i] = interaction1->nnu[i];
                r_in_array[i] = interaction1->r_in[i];
                r_bn_array[i] = bn1->r_bn[i];
                n_bn_array[i] = bn1->n_bn[i];
                posnu_array[i] = interaction1->posnu[i];
                ant_max_normal0_array[i] = ant_max_normal0[i];
                ant_max_normal1_array[i] = ant_max_normal1[i];
                ant_max_normal2_array[i] = ant_max_normal2[i];
                n_pol_array[i] = n_pol[i];
                r_enterice_array[i] = interaction1->r_enterice[i];
                nsurf_rfexit_array[i] = ray1->nsurf_rfexit[i];
                nsurf_rfexit_db_array[i] = ray1->nsurf_rfexit_db[i];
              } //end for (fill arrays)
              for (int j=0;j<5;j++) {
                for (int i=0;i<3;i++) {
                  nrf_iceside_array[j][i] = ray1->nrf_iceside[j][i];
                  nrf_iceside_db_array[j][i] = nrf_iceside_db[j][i];
                  n_exit2bn_array[j][i] = ray1->n_exit2bn[j][i];
                  n_exit2bn_db_array[j][i] = n_exit2bn_db[j][i];
                  rfexit_array[j][i] = ray1->rfexit[j][i];
                  rfexit_db_array[j][i] = ray1->rfexit_db[j][i];
                } //end for
              } //end for
              if (vmmhz_tree->GetEntries()<20) {
                vmmhz_tree->Fill();
              }

              nuflavorint2 = interaction1->nuflavorint;
              costheta_nutraject2=interaction1->costheta_nutraject;
              phi_nutraject2=interaction1->phi_nutraject;
              altitude_int2=interaction1->altitude_int;
              currentint2=interaction1->currentint;
              d12=interaction1->d1;
              d22=interaction1->d2;
              dtryingdirection2=interaction1->dtryingdirection;
              logchord2=interaction1->logchord;
              r_fromballoon2=interaction1->r_fromballoon[0];
              chord_kgm2_bestcase2=interaction1->chord_kgm2_bestcase;
              chord_kgm2_ice2=interaction1->chord_kgm2_ice;
              weight_bestcase2=interaction1->weight_bestcase;
              r_exit2bn2=interaction1->r_exit2bn;
              r_exit2bn_measured2=interaction1->r_exit2bn_measured;

	      sourceMag = ray1->rfexit[2].Mag();
	      
              finaltree->Fill();
              count1->IncrementWeights_r_in(interaction1->r_in, weight);
            } //end if HIST & HISTMAXENTRIES

#ifdef ANITA_UTIL_EXISTS
              realEvPtr 	= new UsefulAnitaEvent();
              rawHeaderPtr	= new RawAnitaHeader();
              Adu5PatPtr	= new Adu5Pat();

              Adu5PatPtr->latitude= bn1->latitude;
              Adu5PatPtr->longitude=bn1->longitude;
              Adu5PatPtr->altitude=bn1->altitude;
              Adu5PatPtr->realTime=bn1->realTime_flightdata;
              Adu5PatPtr->heading = bn1->heading;
              Adu5PatPtr->pitch = bn1->pitch;
              Adu5PatPtr->roll = bn1->roll;
              Adu5PatPtr->run = run_no;

              int fNumPoints = 260;
              for (int i = 0; i < 96; i++){
                for (int j = 0; j < 260; j++){
                  realEvPtr->fVolts[i][j] = 0.;
                  realEvPtr->fTimes[i][j] = 0.;
                }
              }
            
              for (int iant = 0; iant < settings1->NANTENNAS; iant++){
                //int IceMCAnt = GetIceMCAntfromUsefulEventAnt(anita1,  AnitaGeom1,  iant);
                int IceMCAnt = GetIceMCAntfromUsefulEventAnt(settings1,  iant);
                int UsefulChanIndexH = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kHorizontal);
                int UsefulChanIndexV = AnitaGeom1->getChanIndexFromAntPol(iant,  AnitaPol::kVertical);
                realEvPtr->fNumPoints[UsefulChanIndexV] = fNumPoints;
                realEvPtr->fNumPoints[UsefulChanIndexH] = fNumPoints;

                realEvPtr->chanId[UsefulChanIndexV] = UsefulChanIndexV;
                realEvPtr->chanId[UsefulChanIndexH] = UsefulChanIndexH;

		for (int j = 0; j < fNumPoints; j++) {
		  // convert seconds to nanoseconds
		  realEvPtr->fTimes[UsefulChanIndexV][j] = j * anita1->TIMESTEP * 1.0E9;
		  realEvPtr->fTimes[UsefulChanIndexH][j] = j * anita1->TIMESTEP * 1.0E9;
		  // convert volts to millivolts
		  realEvPtr->fVolts[UsefulChanIndexH][j] =  volts_rx_rfcm_lab_h_all[IceMCAnt][j+128]*1000;
                  realEvPtr->fCapacitorNum[UsefulChanIndexH][j] = 0;

                  realEvPtr->fVolts[UsefulChanIndexV][j] =  volts_rx_rfcm_lab_e_all[IceMCAnt][j+128]*1000;
                  realEvPtr->fCapacitorNum[UsefulChanIndexV][j] = 0;
                }//end int j
              }// end int iant

              realEvPtr->eventNumber = inu;
       
              rawHeaderPtr->eventNumber = inu;
              rawHeaderPtr->surfSlipFlag = 0;
              rawHeaderPtr->errorFlag = 0;
	      if (settings1->MINBIAS==1) rawHeaderPtr->trigType = 8; // soft-trigger
	      else                       rawHeaderPtr->trigType = 1; // RF trigger
              rawHeaderPtr->run = run_no;

	      // put the vpol only as a placeholder - these are only used in Anita-2 anyway
              rawHeaderPtr->upperL1TrigPattern = l1trig[0][0];
              rawHeaderPtr->lowerL1TrigPattern = l1trig[0][1];
              rawHeaderPtr->nadirL1TrigPattern = l1trig[0][2];

              rawHeaderPtr->upperL2TrigPattern = l2trig[0][0];
              rawHeaderPtr->lowerL2TrigPattern = l2trig[0][1];
              rawHeaderPtr->nadirL2TrigPattern = l2trig[0][2];

              rawHeaderPtr->l3TrigPattern = (short) l3trig[0];
	      rawHeaderPtr->phiTrigMask  = (short) phiTrigMask;
	      if (settings1->WHICH==9 || settings1->WHICH==10) { 
                rawHeaderPtr->l3TrigPatternH = (short) l3trig[1];
                rawHeaderPtr->l1TrigMask   = (short) l1TrigMask;
                rawHeaderPtr->l1TrigMaskH  = (short) l1TrigMaskH;
                rawHeaderPtr->phiTrigMaskH = (short) phiTrigMaskH;
              }

              rawHeaderPtr->calibStatus = 15;
              rawHeaderPtr->realTime = bn1->realTime_flightdata;
              Adu5PatPtr->latitude= bn1->latitude;
              Adu5PatPtr->longitude=bn1->longitude;
              Adu5PatPtr->altitude=bn1->altitude;
              Adu5PatPtr->realTime=bn1->realTime_flightdata;
              Adu5PatPtr->heading = bn1->heading;
              Adu5PatPtr->pitch = bn1->pitch;
              Adu5PatPtr->roll = bn1->roll;
              Adu5PatPtr->run = run_no;
              eventNumber = inu;

#ifdef ANITA3_EVENTREADER
	      truthEvPtr        = new TruthAnitaEvent();
	      truthEvPtr->eventNumber      = inu;
	      truthEvPtr->realTime         = bn1->realTime_flightdata;
	      truthEvPtr->run              = run_no;
	      truthEvPtr->nuMom            = pnu;
	      truthEvPtr->nu_pdg           = pdgcode;
 	      truthEvPtr->e_component      = e_component;
	      truthEvPtr->h_component      = h_component;
	      truthEvPtr->n_component      = n_component;
 	      truthEvPtr->e_component_k    = e_component_kvector;
	      truthEvPtr->h_component_k    = h_component_kvector;
	      truthEvPtr->n_component_k    = n_component_kvector;
	      truthEvPtr->sourceLon        = sourceLon;    
	      truthEvPtr->sourceLat        = sourceLat;  
	      truthEvPtr->sourceAlt        = sourceAlt;
	      truthEvPtr->weight           = weight;
	      for (int i=0;i<3;i++){
		truthEvPtr->balloonPos[i]  = bn1->r_bn[i];
		truthEvPtr->balloonDir[i]  = bn1->n_bn[i];
		truthEvPtr->nuPos[i]       = interaction1->posnu[i];
		truthEvPtr->nuDir[i]       = interaction1->nnu[i];
	      }
	      for (int i=0;i<5;i++){
		for (int j=0;j<3;j++){
		  truthEvPtr->rfExitNor[i][j] = ray1->n_exit2bn[i][j];
		  truthEvPtr->rfExitPos[i][j] = ray1->rfexit[i][j];
		}
	      }
	      for (int i=0;i<48;i++){
		truthEvPtr->hitangle_e[i]  = hitangle_e_all[i];
		truthEvPtr->hitangle_h[i]  = hitangle_h_all[i];
	      }
	      for (int i=0;i<Anita::NFREQ;i++)
		truthEvPtr->vmmhz[i]       = vmmhz[i];

	      truthAnitaTree->Fill();
	      delete truthEvPtr;
#endif
      
              headTree->Fill();
              eventTree->Fill();
              adu5PatTree->Fill();

              delete realEvPtr;
              delete rawHeaderPtr;
              delete Adu5PatPtr;
#endif

            sum_weights+=weight;
            neutrinos_passing_all_cuts++;
            times_crust_entered_det+=crust_entered;  //Increment counter for neutrino numbers in each earth layer - passing neutrinos
            times_mantle_entered_det+=mantle_entered;
            times_core_entered_det+=core_entered;

            if (settings1->WRITEPOSFILE==1)
              WriteNeutrinoInfo(interaction1->posnu, interaction1->nnu, bn1->r_bn, interaction1->altitude_int, interaction1->nuflavor, interaction1->current, elast_y, nu_out);

            // sample first 1000 events that pass to see the distribution of weights
            if (settings1->HIST && !settings1->ONLYFINAL && sampleweights->GetEntries()<settings1->HIST_MAX_ENTRIES) {
              if (weight>1.E-6)
                sampleweights->Fill(log10(weight));
              else
                sampleweights->Fill(-6.);

              // on the 1000th one,  see how low you should make the cut so that you catch 99% of the events (weighted)
              if (sampleweights->GetEntries()==1000) {
                double sum_sampleintegral=0.;
                double sum_sample=0.;
                // first calculate total integral of all the weights
                for (int k=sampleweights->GetNbinsX();k>=1;k--) {
                  sum_sampleintegral+=sampleweights->GetBinContent(k)*pow(10., sampleweights->GetBinLowEdge(k));
                }
                // treat the underflow bin specially
                sum_sampleintegral+=sampleweights->GetBinContent(0)*pow(10., sampleweights->GetBinLowEdge(1));
                // now sum until you reach 99% of the integral.
                for (int k=sampleweights->GetNbinsX();k>=1;k--) {
                  sum_sample+=sampleweights->GetBinContent(k)*pow(10., sampleweights->GetBinLowEdge(k));
                  if (sum_sample>0.99*sum_sampleintegral) {
                    // reset the cut value.
                    CUTONWEIGHTS=pow(10., sampleweights->GetBinLowEdge(k));
                    cout << "CUTONWEIGHTS is " << CUTONWEIGHTS << "\n";
                    k=0;
                  }
                }
              }
            }//end if HIST & ONLYFINAL & sampleweights HISTMAXENTRIES
    
            // outputs to text file variables relevant to sky map.
            forbrian << interaction1->costheta_nutraject << " " << n_nutraject_ontheground.Phi() << " " << bn1->phi_bn << " " << logweight << "\n";
            // incrementing by flavor
            // also bin in weight for error calculation.
            if (interaction1->nuflavor=="nue") {
              sum[0]+=weight;
              eventsfound_binned_e[index_weights]++;
            } //if
            if (interaction1->nuflavor=="numu") {
              sum[1]+=weight;
              eventsfound_binned_mu[index_weights]++;
            } //if
            if(!sec1->secondbang || !sec1->interestedintaus) {
              if (interaction1->nuflavor=="nutau") {
                sum[2]+=weight;
                eventsfound_binned_tau[index_weights]++;
              } //if
            } //if

          } //end if interaction1->d2>1

        } //end if tautrigger || GetChord
        else {
          cout << "Chord is less than 1m.\n";
        } //end else GetChord


      //      cout <<"inu,  ston is " <<inu << "\t" << anita1->ston[4] << "\n";
      //      if (settings1->HIST==1 && !settings1->ONLYFINAL && anita1->tglob->GetEntries()<settings1->HIST_MAX_ENTRIES && anita1->ston[4]<-1.) {// all events
      if (settings1->HIST==1 && !settings1->ONLYFINAL && anita1->tglob->GetEntries()<settings1->HIST_MAX_ENTRIES) {// all events
	// cout << "Filling global trigger tree.  inu is " << inu << "\n";
        anita1->tglob->Fill();

      }
      //cout << "This event passes. " << inu << "\n";

      passes_thisevent=1; // flag this event as passing
      anita1->tdata->Fill();
      anita1->tgaryanderic->Fill();
      } // end if passing global trigger conditions
      else {
        passes_thisevent=0; // flag this event as not passing
        if (bn1->WHICHPATH==4)
          cout << "This event does not pass.\n";
      }// end else event does not pass trigger
/*
      props<<"weight_prob  "<< weight_prob<<std::endl;
      props<<"weight1  "<< weight1<<std::endl;
      props<<"weight_bestcase2  "<< weight_bestcase2<<std::endl;
      props<<"weight  "<< weight <<std::endl;
      props.close();*/
      ///////////////////////////////////////
      //  
      // WE GET HERE REGARDLESS OF WHETHER THE TRIGGER PASSES
      //
      /////////////


    
      delete globaltrig1;

      // keeping track of intermediate counters,  incrementing by weight1.
      // weight1 was not yet determined when integer counters were incremented.
      if (chanceinhell2)
        count_chanceinhell2_w += weight;

      if (passestrigger)
        count_passestrigger_w += weight;
      
      volume_thishorizon=antarctica->volume_inhorizon[bn1->Getibnposition()]/1.E9;
    
      if (settings1->HIST==1
        && !settings1->ONLYFINAL
        && tree1->GetEntries()<settings1->HIST_MAX_ENTRIES
        && bn1->WHICHPATH != 3){ // all events
        tree1->Fill();
      }//end if
      
    } // end for WHICHRAY 

    //looping over two types of rays - upgoing and downgoing.
    if (ABORT_EARLY){
      std::cout << "\n***********************************************************";
      std::cout << "\n* SIGINT received,  aborting loop over events early.";
      std::cout << "\n* Stopped after event " << inu << " instead of " << NNU;
      std::cout << "\n* Any output which relied on NNU should be corrected for.";
      std::cout << "\n***********************************************************\n";
      foutput << "\n***********************************************************";
      foutput << "\n* SIGINT received,  aborting loop over events early.";
      foutput << "\n* Stopped after event " << inu << " instead of " << NNU;
      foutput << "\n* Any output which relied on NNU should be corrected for.";
      foutput << "\n***********************************************************\n";
      break;
    }
  }//end NNU neutrino loop
  

  // Finished with individual neutrinos now ...


  gRandom=rsave;
  delete Rand3;


#ifdef ANITA_UTIL_EXISTS

  anitafileEvent->cd();
  eventTree->Write("eventTree");
  anitafileEvent->Close();
  delete anitafileEvent;
  
  anitafileHead->cd();
  headTree->Write("headTree");
  anitafileHead->Close();
  delete anitafileHead;
  
  anitafileGps->cd();
  adu5PatTree->Write("adu5PatTree");
  anitafileGps->Close();
  delete anitafileGps;

#ifdef ANITA3_EVENTREADER
  anitafileTruth->cd();
  truthAnitaTree->Write("truthAnitaTree");
  anitafileTruth->Close();
  delete anitafileTruth;
#endif  
  
#endif  

  
  
  cout << "about to close tsignals tree.\n";
  anita1->fsignals=anita1->tsignals->GetCurrentFile();
  anita1->fdata=anita1->tdata->GetCurrentFile();
  anita1->fdata=anita1->tgaryanderic->GetCurrentFile();
  anita1->fsignals->Write();
  anita1->fsignals->Close();
  
  anita1->fdata=anita1->tglob->GetCurrentFile();
  anita1->fdata->Write();
  anita1->fdata->Close();
  
  if (settings1->EVENTSMAP){
    //draw the S80-degree-latitude circle
    TH2F *lat80deg=new TH2F("lat80deg", "", 600, -3000, 3000, 500, -2500, 2500);
    lat80deg->SetMarkerColor(kRed);
    int tmp_e_coord=0, tmp_n_coord=0;
    float tmp_e, tmp_n=0;
    for(double lon=0;lon<360.;lon+=0.5){
        double lat=10.;
        antarctica->IceLonLattoEN(lon, lat, tmp_e_coord, tmp_n_coord);
        tmp_e=float(antarctica->xLowerLeft_ice+tmp_e_coord*antarctica->cellSize)/1000.;
        tmp_n=float(-1*(antarctica->yLowerLeft_ice+tmp_n_coord*antarctica->cellSize))/1000.;
        lat80deg->Fill(tmp_e, tmp_n);
      }//end for lon loop
  }// end if EVENTSMAP

  if (bn1->WHICHPATH==4) {// this is for comparing with Peter
    cout << "Earth radius at South Pole: " << antarctica->Geoid(0.) << "\n";
    Position posnu_temp=Position(180., 0., 1.); // points at the south pole
    cout << "Surface of ice at the South Pole: " << antarctica->Surface(posnu_temp) << "\n";
    cout << "Average balloon altitude is " << average_altitude << "\n";
    cout << "Average distance from earth center is " << average_rbn << "\n";
    cout << "Average height of balloon above ice surface is " << average_rbn-antarctica->Surface(bn1->r_bn) << "\n";
    cout << "theta_zenith are " << anita1->THETA_ZENITH[0]*DEGRAD << " " << anita1->THETA_ZENITH[1]*DEGRAD << " " << anita1->THETA_ZENITH[2]*DEGRAD << "\n";
    cout << "Index of refraction at this depth is " << sig1->N_DEPTH << "\n";
    cout << "Cerenkov angle is " << sig1->changle*DEGRAD << "\n";
    cout << "Nadir angle to surface exit point is " << DEGRAD*bn1->r_bn.Angle(ray1->n_exit2bn[2]) << "\n";
    cout << "Distance from rfexit to balloon is " << (bn1->r_bn+ -1*ray1->rfexit[2]).Mag() << "\n";
    cout << "Payload zenith angle at event source is " << DEGRAD*ray1->rfexit[2].Angle(ray1->n_exit2bn[2]) << "\n";
    cout << "Angle of incidence just below surface is " << DEGRAD*ray1->rfexit[2].Angle(ray1->nrf_iceside[4]) << "\n";
    cout << "Angle between neutrino and surface normal is " << DEGRAD*ray1->rfexit[0].Angle(interaction1->nnu)-90. << "\n";
    cout << "Angle of incidence below firn surface is " << DEGRAD*ray1->rfexit[2].Angle(ray1->nrf_iceside[3]) << "\n";
  }
  cout << "about to Summarize.\n";
  
  anita1->rms_rfcm_e = sqrt(anita1->rms_rfcm_e / (double)anita1->count_getnoisewaveforms)*1000.;
  anita1->rms_rfcm_h = sqrt(anita1->rms_rfcm_h / (double)anita1->count_getnoisewaveforms)*1000.;
  anita1->rms_lab_e = sqrt(anita1->rms_lab_e / (double)anita1->count_getnoisewaveforms)*1000.;
  anita1->rms_lab_h = sqrt(anita1->rms_lab_h / (double)anita1->count_getnoisewaveforms)*1000.;
  
  cout << "RMS noise in rfcm e-pol is " << anita1->rms_rfcm_e << " mV.\n";
  cout << "RMS noise in rfcm h-pol is " << anita1->rms_rfcm_h << " mV.\n";
  cout << "RMS noise in lab e-pol is " << anita1->rms_lab_e << "mV.\n";
  cout << "RMS noise in lab h-pol is " << anita1->rms_lab_h << "mV.\n";
  for (int i=0;i<Anita::NFREQ;i++) {
    anita1->avgfreq_rfcm[i]/=(double)anita1->count_getnoisewaveforms;
    anita1->avgfreq_rfcm_lab[i]/=(double)anita1->count_getnoisewaveforms;
  }
    
  rms_rfcm_e=anita1->rms_rfcm_e;
  rms_rfcm_h=anita1->rms_rfcm_h;
  rms_lab_e=anita1->rms_lab_e;
  rms_lab_h=anita1->rms_lab_h;
  for (int i=0;i<Anita::NFREQ;i++) {
    avgfreq_rfcm[i]=anita1->avgfreq_rfcm[i];
    avgfreq_rfcm_lab[i]=anita1->avgfreq_rfcm_lab[i];
    freq[i]=anita1->freq[i];
  }
  cout << "Filling summarytree.  rms_rfcm_e is " << rms_rfcm_e << "\n";
  summarytree->Fill();
  

  // maks the output file
  Summarize(settings1, anita1, count1, spectra1, sig1, primary1, pnu, eventsfound, eventsfound_db, eventsfound_nfb, sigma, sum, antarctica->volume, antarctica->ice_area, km3sr, km3sr_e, km3sr_mu, km3sr_tau, foutput, distanceout);

  veff_out << settings1->EXPONENT << "\t" << km3sr << "\t" << km3sr_e << "\t" << km3sr_mu << "\t" << km3sr_tau << "\t" << settings1->SIGMA_FACTOR << endl;//this is for my convenience
  
  // for each neutrino flavor,  fraction each contributes to sensitivity.
  sum_frac[0]=sum[0]/eventsfound;
  sum_frac[1]=sum[1]/eventsfound;
  sum_frac[2]=sum[2]/eventsfound;
  
  // for taus.
  sum_frac_db[0]=sum[0]/(eventsfound+eventsfound_db+eventsfound_nfb);
  sum_frac_db[1]=sum[1]/(eventsfound+eventsfound_db+eventsfound_nfb);
  sum_frac_db[2]=(sum[2]+eventsfound_db+eventsfound_nfb)/(eventsfound+eventsfound_db+eventsfound_nfb);
  //if (tree17->GetEntries()<settings1->HIST_MAX_ENTRIES && !settings1->ONLYFINAL && HIST==1)
    //tree17->Fill();


  cout << "closing file.\n";
  CloseTFile(hfile);

  time_t raw_end_time = time(NULL);
  struct tm * end_time = localtime(&raw_end_time);  
  cout << "Date and time at end of run are: " << asctime (end_time) << "\n";
  cout<<"Total time elapsed is "<<(int)((raw_end_time - raw_start_time)/60)<<":"<< ((raw_end_time - raw_start_time)%60)<<endl;
  
  foutput << "\nTotal time elapsed in run is " <<(int)((raw_end_time - raw_start_time)/60)<<":"<< ((raw_end_time - raw_start_time)%60)<<endl;
  std::cerr<<"CUT  "<< CUTONWEIGHTS<<std::endl;
  delete anita1;
  return 0;

} //END MAIN PROGRAM









//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//
//  Auxiliary functions
//
//
//


void IntegrateBands(Anita *anita1, int k, Screen *panel1, double *freq, double scalefactor, double *sumsignal) {
  for (int j=0;j<5;j++) {
  // if this frequency is in this bandwidth slice
    for (int jpt=0; jpt<panel1->GetNvalidPoints(); jpt++){
      if (anita1->bwslice_min[j]<=freq[k] && anita1->bwslice_max[j]>freq[k])
        sumsignal[j]+=panel1->GetVmmhz_freq(jpt*Anita::NFREQ + k)*(freq[k+1]-freq[k])*scalefactor;
    }
  }
}
//end IntegrateBands()


void Integrate(Anita *anita1, int j, int k, double *vmmhz, double *freq, double scalefactor, double sumsignal) {
  // if this frequency is in this bandwidth slice
  if (anita1->bwslice_min[j]<=freq[k] && anita1->bwslice_max[j]>freq[k])
    sumsignal+=vmmhz[k]*(freq[k+1]-freq[k])*scalefactor;
}
//end Integrate()


void WriteNeutrinoInfo(Position &posnu,  Vector &nnu,  Position &r_bn,  double altitude_int,  string nuflavor,  string current,  double elast_y,  ofstream &nu_out) {
  nu_out << "\n" << inu << "\t" << posnu[0] << " " << posnu[1] << " " << posnu[2] << "\t" << altitude_int << "\t" << nnu[0] << " " << nnu[1] << " " << nnu[2] << "\t" << r_bn[0] << " " << r_bn[1] << " " << r_bn[2] << "\t" << nuflavor << "\t" << current << "\t" << elast_y << "\n\n";
}
//end WriteNeutrinoInfo()


//void Summarize(Settings *settings1, Spectra *spectra1, Signal *sig1, Counting *count1, double pnu, double eventsfound, double eventsfound_db, double eventsfound_nfb, double sigma, double* sum, double volume, 
//	       double& km3sr, double& km3sr_e, double& km3sr_mu, double& km3sr_tau, ofstream &foutput, ofstream &distanceout, ofstream &foutseckel) {
void Summarize(Settings *settings1,  Anita* anita1,  Counting *count1, Spectra *spectra1, Signal *sig1, Primaries *primary1, double pnu, double eventsfound, double eventsfound_db, double eventsfound_nfb, double sigma, double* sum, double volume, double ice_area, double& km3sr, double& km3sr_e, double& km3sr_mu, double& km3sr_tau, ofstream &foutput, ofstream &distanceout) {
  double rate_v_thresh[NTHRESHOLDS];
  double errorup_v_thresh[NTHRESHOLDS];
  double errordown_v_thresh[NTHRESHOLDS];
  double rate_h_thresh[NTHRESHOLDS];
  double errorup_h_thresh[NTHRESHOLDS];
  double errordown_h_thresh[NTHRESHOLDS];
  double zeroes[NTHRESHOLDS];
  Tools::Zero(zeroes, NTHRESHOLDS);

  // plot result of threshold scan
  for (int i=0;i<NTHRESHOLDS;i++) {
    rate_v_thresh[i]=npass_v_thresh[i]/denom_v_thresh[i];
    // cout << "i,  npass_v_thresh are " << i << "\t" << npass_v_thresh[i] << " " << denom_v_thresh[i] << " " << rate_v_thresh[i] << "\n";
    if (npass_v_thresh[i]<=20) {
      errorup_v_thresh[i]=poissonerror_plus[(int)npass_v_thresh[i]]/denom_v_thresh[i];
      errordown_v_thresh[i]=poissonerror_minus[(int)npass_v_thresh[i]]/denom_v_thresh[i];
      // cout << errorup_v_thresh[i] << " " << errordown_v_thresh[i] << endl;
    }//end if
    if (settings1->WHICH==9 || settings1->WHICH==10){ // Anita-3
      rate_h_thresh[i]=npass_h_thresh[i]/denom_h_thresh[i];
      if (npass_h_thresh[i]<=20) {
        errorup_h_thresh[i]=poissonerror_plus[(int)npass_h_thresh[i]]/denom_h_thresh[i];
        errordown_h_thresh[i]=poissonerror_minus[(int)npass_h_thresh[i]]/denom_h_thresh[i];
      }//end if
    }//end if WHICH==9

  }//end for NTHRESHOLDS

  string stemp=settings1->outputdir+"/thresholds.root";
  TFile *fthresholds=new TFile(stemp.c_str(), "RECREATE");
  TCanvas *cthresh=new TCanvas("cthresh", "cthresh", 880, 800);
  cthresh->SetLogy();

  TGraph *gnpass=new TGraph(NTHRESHOLDS, thresholds, npass_v_thresh);
  gnpass->SetName("npass");
  TGraph *gdenom=new TGraph(NTHRESHOLDS, thresholds, denom_v_thresh);
  gdenom->SetName("denom");
 
  TGraphAsymmErrors *g=new TGraphAsymmErrors(NTHRESHOLDS, thresholds, rate_v_thresh, zeroes, zeroes, errorup_v_thresh, errordown_v_thresh);
  g->SetName("rate");

  g->SetLineWidth(2);
  g->SetMarkerStyle(21);
  g->Draw("ape");

  stemp = settings1->outputdir+"/thresholds.eps";
  cthresh->Print((TString)stemp);
  g->Write();
  gdenom->Write();
  gnpass->Write();


  if (settings1->WHICH==9 || settings1->WHICH==10){ // Anita-3 or Anita-4
    TGraph *gnpassH=new TGraph(NTHRESHOLDS, thresholds, npass_h_thresh);
    gnpassH->SetName("npassH");
    TGraph *gdenomH=new TGraph(NTHRESHOLDS, thresholds, denom_h_thresh);
    gdenomH->SetName("denomH");
    TGraphAsymmErrors *gH=new TGraphAsymmErrors(NTHRESHOLDS, thresholds, rate_h_thresh, zeroes, zeroes, errorup_h_thresh, errordown_h_thresh);
    gH->SetName("rate");
    gH->SetLineWidth(2);
    gH->SetMarkerStyle(21);
    gH->Draw("ape");
    stemp = settings1->outputdir+"/thresholds_HPOL.eps";
    cthresh->Print((TString)stemp);
    
    gH->Write();
    gdenomH->Write();
    gnpassH->Write();
  }//end if WHICH==9 or WHICH==10

  fthresholds->Write();  
  fthresholds->Close();

  //  double ses;                          // single-event sensitivity
  double km2sr;                        // aperture km**2-sr

  foutput << "\n";
  
  // write out summary
  foutput << "Generated " << NNU << " neutrinos with energy " << pnu << "\n";
  foutput << "Number of unweighted direct,  reflected that pass is: " << count1->npass[0] << "\t" << count1->npass[1] << "\n";
  foutput << "Number of (weighted) neutrinos that pass is: " << eventsfound << "\n";
  foutput << "Number of (weighted) neutrinos that pass,  multiplied by prob. of interacting in the ice,  is: " << eventsfound_prob << "\n";
  foutput << "Number of weighted direct,  reflected that pass is: " << allcuts_weighted[0] << "\t" << allcuts_weighted[1] << "\n";
  foutput << "Number of (weighted) neutrinos that pass (with weight>0.001) is: " << eventsfound_weightgt01 << "\n";
  foutput << "Number of (weighted) neutrinos that only traverse the crust is " << eventsfound_crust << " -> " << eventsfound_crust/eventsfound*100 << "%\n\n";
  foutput << "Number of (weighted) neutrinos that pass only VPOL trigger is: " << allcuts_weighted_polarization[0] << "\n";
  foutput << "Number of (weighted) neutrinos that pass only HPOL trigger is: " << allcuts_weighted_polarization[1] << "\n";
  foutput << "Number of (weighted) neutrinos that pass both pol triggers is: " << allcuts_weighted_polarization[2] << "\n\n";
	       
  foutput << "Volume of ice is " << volume << "\n";
  foutput << "Value of 4*pi*pi*r_earth*r_earth in km^2 " << 4*PI*PI*(EarthModel::R_EARTH*EarthModel::R_EARTH/1.E6) << "\n";
  
  
  // single event sensitivity for 150 MHz array per year
  //  everything is normalized to the active volume
  //  assume dF/dE propto E**-2
  //  assume 2pi sr or 4pi sr depending on whether below horizon is used
  //  assume 3.16e7 seconds/year
  //------------------------------------------------
  
  // ses=1/(A_eff * sr * t)
  //ses=1.0/(sigma*volume*RHOSALT*(1./M_NUCL)*3.16E7);  // 1 per year
  
  
  double nevents = 0;
  //double nweights=0;
  double nevents_db=0;
  double nevents_nfb=0;
  
  
  nevents+=eventsfound;
  nevents_db+=eventsfound_db;
  nevents_nfb+=eventsfound_nfb;
  error_plus=0;
  error_e_plus=0;
  error_mu_plus=0;
  error_tau_plus=0;
  error_minus=0;
  error_e_minus=0;
  error_mu_minus=0;
  error_tau_minus=0;
  
  error_nfb=0;
  error_km3sr_nfb=0;
  error_percent_increase_nfb=0;
  
  
  // loop over bins in weights
  for (int i=0;i<NBINS;i++) {
    if (eventsfound_binned[i]<=20) {
      error_plus+=pow(poissonerror_plus[(int)eventsfound_binned[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_e_plus+=pow(poissonerror_plus[(int)eventsfound_binned_e[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_mu_plus+=pow(poissonerror_plus[(int)eventsfound_binned_mu[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_tau_plus+=pow(poissonerror_plus[(int)eventsfound_binned_tau[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_minus+=pow(poissonerror_minus[(int)eventsfound_binned[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_e_minus+=pow(poissonerror_minus[(int)eventsfound_binned_e[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_mu_minus+=pow(poissonerror_minus[(int)eventsfound_binned_mu[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_tau_minus+=pow(poissonerror_minus[(int)eventsfound_binned_tau[i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
    } //if
    else {
      error_plus+=eventsfound_binned[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_e_plus+=eventsfound_binned_e[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_mu_plus+=eventsfound_binned_mu[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_tau_plus+=eventsfound_binned_tau[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      error_minus=error_plus;
      error_e_minus=error_e_plus;
      error_mu_minus=error_mu_plus;
      error_tau_minus=error_tau_plus;
    } //else

    error_nfb+=eventsfound_nfb_binned[i]*pow(pow(10., ((double)i+0.5)/(double)NBINS*4.+-5.), 2);
    for (int j=0;j<NBINS_DISTANCE;j++) {
      if (eventsfound_binned_distance_forerror[j][i]<=20) {
        error_distance_plus[j]+=pow(poissonerror_plus[(int)eventsfound_binned_distance_forerror[j][i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
        error_distance_minus[j]+=pow(poissonerror_minus[(int)eventsfound_binned_distance_forerror[j][i]]*pow(10., ((double)i+0.5)/(double)NBINS*(MAX_LOGWEIGHT-MIN_LOGWEIGHT)+MIN_LOGWEIGHT), 2);
      }//end if forerror
    }//end for NBINS_DISTANCE
  }//end for NBINS
  
  
  error_plus=sqrt(error_plus);
  error_e_plus=sqrt(error_e_plus);
  error_mu_plus=sqrt(error_mu_plus);
  error_tau_plus=sqrt(error_tau_plus);
  
  error_minus=sqrt(error_minus);
  error_e_minus=sqrt(error_e_minus);
  error_mu_minus=sqrt(error_mu_minus);
  error_tau_minus=sqrt(error_tau_minus);
  
  
  error_nfb=sqrt(error_nfb);
  for (int j=0;j<NBINS_DISTANCE;j++) {
    error_distance_plus[j]=sqrt(error_distance_plus[j]);
    error_distance_minus[j]=sqrt(error_distance_minus[j]);
  }
  
  
  // account for efficiency
  if (NNU != 0 && nevents!=0) {
    km3sr=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*nevents/(double)NNU/settings1->SIGMA_FACTOR;
    
    cout << nevents << " events passed out of " << NNU << "\n";
    
    error_plus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_plus/(double)NNU/settings1->SIGMA_FACTOR;
    error_minus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_minus/(double)NNU/settings1->SIGMA_FACTOR;
    
    percent_increase_db=km3sr_db/km3sr*100;
    percent_increase_nfb=km3sr_nfb/km3sr*100;
    
    error_percent_increase_nfb=sqrt(pow(error_km3sr_nfb/km3sr_nfb, 2)+pow(error_plus/km3sr, 2))*percent_increase_nfb;
    
    percent_increase_total=percent_increase_db+percent_increase_nfb;
    
    km3sr_e = (pow(1.e-3, 3))*volume*sr*(sum[0]/(double)count1->nnu_e)*sig1->RHOMEDIUM/sig1->RHOH20/settings1->SIGMA_FACTOR;
    
    cout << sum[0]/(double)nevents*100. << "% are electron neutrinos\n";
    
    error_e_plus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_e_plus/(double)count1->nnu_e/settings1->SIGMA_FACTOR;
    error_e_minus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_e_minus/(double)count1->nnu_e/settings1->SIGMA_FACTOR;
    
    km3sr_mu = (pow(1.e-3, 3))*volume*sr*(sum[1]/(double)count1->nnu_mu)*sig1->RHOMEDIUM/sig1->RHOH20/settings1->SIGMA_FACTOR;
    
    cout << sum[1]/(double)nevents*100. << "% are muon neutrinos\n";
    
    error_mu_plus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_mu_plus/(double)count1->nnu_mu/settings1->SIGMA_FACTOR;
    error_mu_minus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_mu_minus/(double)count1->nnu_mu/settings1->SIGMA_FACTOR;

    km3sr_tau = (pow(1.e-3, 3))*volume*sr*(sum[2]/(double)count1->nnu_tau)*sig1->RHOMEDIUM/sig1->RHOH20/settings1->SIGMA_FACTOR;

    cout << sum[2]/(double)nevents*100. << "% are tau neutrinos\n";

    error_tau_plus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_tau_plus/(double)count1->nnu_tau/settings1->SIGMA_FACTOR;
    error_tau_minus=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/sig1->RHOH20*sr*error_tau_minus/(double)count1->nnu_tau/settings1->SIGMA_FACTOR;

    double sum_km3sr=0;
    double sum_km3sr_error_plus=0;
    double sum_km3sr_error_minus=0;
    for (int i=0;i<NBINS_DISTANCE;i++) {
      sum_km3sr+= (pow(1.e-3, 3))*volume*sr*eventsfound_binned_distance[i]*sig1->RHOMEDIUM/sig1->RHOH20/(double)NNU/settings1->SIGMA_FACTOR;
      km3sr_distance[i]=sum_km3sr;
      sum_km3sr_error_plus+=pow(pow(1.E-3, 3)*volume*error_distance_plus[i]*sig1->RHOMEDIUM/sig1->RHOH20*sr/(double)NNU, 2)/settings1->SIGMA_FACTOR;
      error_distance_plus[i]=sqrt(sum_km3sr_error_plus);
      sum_km3sr_error_minus+=pow(pow(1.E-3, 3)*volume*error_distance_minus[i]*sig1->RHOMEDIUM/sig1->RHOH20*sr/(double)NNU, 2)/settings1->SIGMA_FACTOR;
      error_distance_minus[i]=sqrt(sum_km3sr_error_minus);
      distanceout << 700000./(double)NBINS_DISTANCE*(double)i << "\t" << km3sr_distance[i] << "\t" << error_distance_plus[i] << "\t" << error_distance_minus[i] << "\n";
    } //for

    foutput << "Total volume * solid angle is \t\t\t\t" << km3sr << " + " << error_plus << " - " << error_minus << " km^3 str\n";
    foutput << "Total volume * solid angle for electron neutrinos is \t" << km3sr_e << " + " << error_e_plus << " - " << error_e_minus << " km^3 str\n";
    foutput << "Total volume * solid angle for muon neutrinos is \t" << km3sr_mu << " + " << error_mu_plus << " - " << error_mu_minus << " km^3 str\n";
    foutput << "Total volume * solid angle for tau neutrinos is \t" << km3sr_tau << " + " << error_tau_plus << " - " << error_tau_minus << " km^3 str\n";
    cout << "Total volume * solid angle is \t\t\t\t" << km3sr << " + " << error_plus << " - " << error_minus << " km^3 str\n";
    cout << "Total volume * solid angle for electron neutrinos is \t" << km3sr_e << " + " << error_e_plus << " - " << error_e_minus << " km^3 str\n";
    cout << "Total volume * solid angle for muon neutrinos is \t" << km3sr_mu << " + " << error_mu_plus << " - " << error_mu_minus << " km^3 str\n";
    cout << "Total volume * solid angle for tau neutrinos is \t" << km3sr_tau << " + " << error_tau_plus << " - " << error_tau_minus << " km^3 str\n";

    if ( spectra1->IsMonoenergetic() )  {
      cout << "Cross section is " << sigma << "m^2\n";
      double len_int=1.0/(sigma*sig1->RHOH20*(1./M_NUCL)*1000);
      cout << "Interaction length is " << len_int << " km\n\n";
      foutput << "Interaction length is " << len_int << "\n";
      km2sr=km3sr/len_int;
      foutput << "Total area x steradians using km3sr/len_int is \t\t\t\t" << km2sr << " km^2 str\n\n";
      cout << "Total area x steradians is \t\t" << km2sr << " km^2 str\n\n";
      km2sr=ice_area/(1.E6)*PI*eventsfound_prob/(double)NNU;
      foutput << "Total area x steradians using 4*PI*R_EARTH^2*eff. is \t" << km2sr << " km^2 str\n\n";
      foutput << "These are not the same because we are not throwing all directions on all points of the surface.  Believe the first one as an approximation,  we are working on this for high cross sections.\n";
      // ses=(pnu/1.E9)/(km2sr*3.16E7);
    }//end if IsMonoenergetic

  }//end if NNU!=0 and nevents!=0
  
  
  foutput << "\t\t\t\t\t\t\tprobability of passing \t# events passing\n";
  
  foutput.precision(4);
  foutput << "No way this neutrino will see any ice \t\t\t" << (double)count1->noway[0]/(double)count_total << "\t" <<
  (double)count1->noway[1]/(double)count_total << "\t\t" <<
  count1->noway[0] << "\t" << count1->noway[1] << "\n";
  
  foutput.precision(4);
  foutput << "Wheredoesitleave in PickUnbiased gives an error \t" << (double)count1->wheredoesitleave_err[0]/(double)count1->noway[0] << "\t" <<
  (double)count1->wheredoesitleave_err[1]/(double)count1->noway[1] << "\t\t" <<
  count1->wheredoesitleave_err[0] << "\t" << count1->wheredoesitleave_err[1] << "\n";
  
  foutput.precision(4);
  foutput << "This neutrino direction never sees ice \t\t\t" << (double)count1->neverseesice[0]/(double)count1->wheredoesitleave_err[0] << "\t" <<
  (double)count1->neverseesice[1]/(double)count1->wheredoesitleave_err[1] << "\t\t" <<
  count1->neverseesice[0] << "\t" << count1->neverseesice[1] << "\n";
  
  
  foutput.precision(4);
  foutput << "WhereDoesItEnterIce in PickUnbiased gives an error \t\t\t" << (double)count1->wheredoesitenterice_err[0]/(double)count1->neverseesice[0] << "\t" <<
  (double)count1->wheredoesitenterice_err[1]/(double)count1->neverseesice[1] << "\t\t" <<
  count1->wheredoesitenterice_err[0] << "\t" << count1->wheredoesitenterice_err[1] << "\n";
  
  foutput.precision(4);
  foutput << "Interaction point too high \t\t\t\t" << (double)count1->toohigh[0]/(double)count1->wheredoesitenterice_err[0] << "\t" <<
  (double)count1->toohigh[1]/(double)count1->wheredoesitenterice_err[1] << "\t\t" <<
  count1->toohigh[0] << "\t" << count1->toohigh[1] << "\n";
  
  foutput.precision(4);
  foutput << "Interaction point too low \t\t\t\t" << (double)count1->toolow[0]/(double)count1->toohigh[0] << "\t" <<
  (double)count1->toolow[1]/(double)count1->toohigh[1] << "\t\t" <<
  count1->toolow[0] << "\t" << count1->toolow[1] << "\n";
  
  
  
  
  foutput.precision(4);
  foutput << "There is an interaction in ice \t\t\t\t" << (double)count1->iceinteraction[0]/(double)count1->toolow[0] << "\t" <<
  (double)count1->iceinteraction[1]/(double)count1->toolow[1] << "\t\t" <<
  count1->iceinteraction[0] << "\t" << count1->iceinteraction[1] << "\n";
  
  foutput.precision(4);
  foutput << "In horizon \t\t\t\t\t\t" << (double)count1->inhorizon[0]/(double)count1->iceinteraction[0] << "\t" <<
  (double)count1->inhorizon[1]/(double)count1->iceinteraction[1] << "\t\t" <<
  count1->inhorizon[0] << "\t" << count1->inhorizon[1] << "\n";
  
  
  
  foutput.precision(4);
  foutput << "From surface to balloon,  ray not intersected by earth \t" << (double)count1->nraypointsup1[0]/(double)count1->inhorizon[0] << "\t" <<
  (double)count1->nraypointsup1[1]/(double)count1->inhorizon[1] << "\t\t" <<
  count1->nraypointsup1[0] << "\t" << count1->nraypointsup1[1] << "\n";
  
  foutput.precision(4);
  foutput << "After 1/r scaling and best case attenuation, \n\tMaximum signal is detectable\t\t\t" << (double)count1->nnottoosmall[0]/(double)count1->nraypointsup1[0] << "\t" <<
  (double)count1->nnottoosmall[1]/(double)count1->nraypointsup1[1] << "\t\t"
  << count1->nnottoosmall[0] << "\t" << count1->nnottoosmall[1] << "\n";
  
  
  foutput.precision(4);
  foutput << "Viewing angle lt 90 degrees\t\t\t" << (double)count1->nviewangle_lt_90[0]/(double)count1->nnottoosmall[0] << "\t" <<
  (double)count1->nviewangle_lt_90[1]/(double)count1->nnottoosmall[1] << "\t\t"
  << count1->nviewangle_lt_90[0] << "\t" << count1->nviewangle_lt_90[1] << "\n";
  
  
  foutput.precision(4);
  foutput << "Reality check:  EM and hadronic fractions both nonzero\t" << (double)count1->ngoodfracs[0]/(double)count1->nviewangle_lt_90[0] << "\t" <<
  (double)count1->ngoodfracs[1]/(double)count1->nviewangle_lt_90[1] << "\t\t"
  << count1->ngoodfracs[0] << "\t" << count1->ngoodfracs[1] << "\n";
  foutput.precision(4);
  foutput << "\tBoth EM and hadronic fractions are zero\t\t" << (double)count1->nbadfracs[0]/(double)count1->nviewangle_lt_90[0] << "\t" <<
  (double)count1->nbadfracs[1]/(double)count1->nviewangle_lt_90[1] << "\t\t" <<
  count1->nbadfracs[0] << "\t" << count1->nbadfracs[1] << "\n";
  
  foutput.precision(4);
  foutput << "After finding neutrino direction,  \n\tchance of making through Earth\t\t\t" << (double)count_chanceofsurviving/(double)count1->ngoodfracs[0] << "\t\t\t";
  foutput.precision(10);
  foutput << count_chanceofsurviving << "\n";
  foutput.precision(4);
  foutput << "Neutrino enters ice south of 60deg S latitude\t\t" << (double)count1->nentersice[0]/(double)count_chanceofsurviving << "\t" <<
  (double)count1->nentersice[1]/(double)count_chanceofsurviving <<
  "\t\t" <<
  count1->nentersice[0] << "\t" << count1->nentersice[1] << "\n";
  
  foutput.precision(4);
  foutput << "Neutrino reasonably likely to survive trip through Earth " << (double)count1->nabsorbed[0]/(double)count1->nentersice[0] << "\t" <<
  (double)count1->nabsorbed[1]/(double)count1->nentersice[1] << "\t\t"
  << count1->nabsorbed[0] << "\t" << count1->nabsorbed[1] << "\n";
  foutput.precision(4);
  foutput << "Ray leaves the ice south of 60deg S latitude\t\t" << (double)count1->nraywithincontinent1[0]/(double)count1->nabsorbed[0] << "\t" <<
  (double)count1->nraywithincontinent1[1]/(double)count1->nabsorbed[1] << "\t" <<
  count1->nraywithincontinent1[0] << "\t" <<
  count1->nraywithincontinent1[1] << "\n";
  
  foutput.precision(4);
  foutput << "After 1/r,  best guess ice attenuation,  \n\tmaximum signal is detectable\t\t\t" << (double)count_chanceinhell0/(double)count1->nraywithincontinent1[0] << "\t\t\t";
  foutput.precision(10);
  foutput <<count_chanceinhell0 << "\n";
  
  foutput.precision(4);
  foutput << "Ray is not totally internally reflected\t\t\t" << (double)count1->nnottir[0]/(double)count_chanceinhell0 <<  "\t" <<
  (double)count1->nnottir[1]/(double)count_chanceinhell0 <<  "\t\t" <<
  count1->nnottir[0] << "\t" << count1->nnottir[1] << "\n";
  
  foutput.precision(4);
  foutput << "From surface to balloon,  ray not intersected by earth \t" << (double)count1->nraypointsup2[0]/(double)count1->nnottir[0] << "\t" <<
  (double)count1->nraypointsup2[1]/(double)count1->nnottir[1] <<
  "\t\t"
  << count1->nraypointsup2[0] << "\t" << count1->nraypointsup2[1] << "\n";
  foutput.precision(4);
  
  foutput << "Ray leaves the ice south of 60deg S latitude\t\t" << (double)count1->nraywithincontinent2[0]/(double)count1->nraypointsup2[0] <<"\t" <<
  (double)count1->nraywithincontinent2[0]/(double)count1->nraypointsup2[1] <<
  "\t\t" << count1->nraywithincontinent2[0] << "\t" <<
  count1->nraywithincontinent2[1] 	  << "\n";
  foutput.precision(4);
  
  foutput << "Ray leaves where there is ice\t\t\t\t" << (double)count1->nacceptablerf[0]/(double)count1->nraywithincontinent2[0] << "\t" <<
  (double)count1->nacceptablerf[1]/(double)count1->nraywithincontinent2[1] << "\t\t"
  << count1->nacceptablerf[0] << "\t" << count1->nacceptablerf[1] << "\n";
  foutput.precision(4);
  
  foutput << "Ray tracing converges to within 10 m\t\t\t" << (double)count1->nconverges[0]/(double)count1->nacceptablerf[0] << "\t" <<
  (double)count1->nconverges[1]/(double)count1->nacceptablerf[1] <<
  "\t\t" << count1->nconverges[0] << "\t" << count1->nconverges[1] << "\n";
  foutput.precision(4);
  
  foutput << "After fresnel coefficient,  \n\tmaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell_fresnel[0]/(double)count1->nconverges[0] << "\t" << (double)count1->nchanceinhell_fresnel[1]/(double)count1->nconverges[1] <<
  "\t\t" <<count1->nchanceinhell_fresnel[0] << "\t" << count1->nchanceinhell_fresnel[1] << "\n";
  foutput.precision(4);
  foutput << "After 1/r,  \n\tmaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell_1overr[0]/(double)count1->nchanceinhell_fresnel[0] << "\t" <<  (double)count1->nchanceinhell_1overr[1]/(double)count1->nchanceinhell_fresnel[1] << "\t\t" <<count1->nchanceinhell_1overr[0] << "\t" << count1->nchanceinhell_1overr[1] << "\n";
  foutput.precision(4);
  
  foutput << "After ice attenuation,  \n\tmaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell[0]/(double)count1->nchanceinhell_1overr[0] << "\t" <<
  (double)count1->nchanceinhell[1]/(double)count1->nchanceinhell_1overr[1] << "\t" <<count1->nchanceinhell[0] << "\t" << count1->nchanceinhell[1] << "\n";
  foutput.precision(4);
  
  foutput << "After viewing angle cut, \t\t\t\t" << (double)count1->nviewanglecut[0]/(double)count1->nchanceinhell[0] << "\t" << (double)count1->nviewanglecut[1]/(double)count1->nchanceinhell[1] << "\t\t" << count1->nviewanglecut[0] << " " << count1->nviewanglecut[1] << "\n";
  
  foutput.precision(4);
  foutput << "After factoring in off-Cerenkov cone tapering, \n\tMaximum signal is detectable\t\t\t" << (double)count1->nchanceinhell2[0]/(double)count1->nviewanglecut[0] << "\t" << (double)count1->nchanceinhell2[1]/(double)count1->nviewanglecut[1] << "\t\t" << count1->nchanceinhell2[0] << " " << count1->nchanceinhell2[1] << "\n";
  
  foutput << "Passes trigger\t\t\t\t\t\t" << (double)count1->npassestrigger[0]/(double)count1->nchanceinhell2[0] << "\t" << (double)count1->npassestrigger[1]/(double)count1->nchanceinhell2[1] << "\t\t" << count1->npassestrigger[0] << "\t" << count1->npassestrigger[1] << "\n";
  foutput << "Number of l1 triggers\t\t\t\t\t\t" << (double)count1->nl1triggers[0][0] << "\t" << (double)count1->nl1triggers[1][0] << "\n";

  foutput << "Chord is good length\t\t\t\t\t" << (double)count_chordgoodlength/(double)count1->npassestrigger[0] << "\t\t\t";
  foutput.precision(10);
  foutput <<count_chordgoodlength << "\n";
  foutput.precision(4);
  foutput << "Neutrino's path in ice more than 1m \t\t\t" << (double)count_d2goodlength/(double)count_chordgoodlength << "\t\t\t";
  foutput.precision(10);
  foutput << count_d2goodlength << "\n";
  foutput.precision(4);
  foutput << "Events that pass all cuts\t\t\t\t" << (double)count1->npass[0]/(double)count_d2goodlength << "\t" << (double)count1->npass[1]/(double)count_d2goodlength << "\t\t";
  foutput.precision(10);
  foutput <<count1->npass[0] << "\t" << count1->npass[1] << "\n";
  
  
  //  if (EXPONENT<=10||EXPONENT>100) {
  if ( spectra1->IsSpectrum() ) {
    double sum_events=0.;
    double thisenergy=0.;
    double thislen_int_kgm2=0.;
    // double *energy=spectra1->GetEnergyArray();
    // double *EdNdEdAdt=spectra1->GetEdNdEdAdt();

    // for models which don't have even spaced energy bin, 
    double even_E;
    int N_even_E = 12;
    even_E = ( spectra1->Getenergy()[spectra1->GetE_bin() - 1] - spectra1->Getenergy()[0] ) / ( (double) N_even_E );
    for (int i=0;i<N_even_E;i++) {
      thisenergy=pow(10., (spectra1->Getenergy())[0]+((double)i)*even_E);
      primary1->GetSigma(thisenergy, sigma, thislen_int_kgm2, settings1, 0, 0);
      
      // EdNdEdAdt is in #/cm^2/s
      // can also be written dN/d(lnE)dAdt
      // = dN*log(10)/d(log E)dAdt
      // the bin spacing is 0.5
      // so # events ~ dN*log(10)*0.5/d(log E)dAdt
      sum_events+=even_E*log(10.)*( spectra1->GetEdNdEdAdt(log10(thisenergy)) )/(thislen_int_kgm2/sig1->RHOH20);
      //cout << "thisenergy,  EdNdEdAdt is " << thisenergy << " " << spectra1->EdNdEdAdt[i] << "\n";
      //foutput << "interaction length is " << thislen_int_kgm2/RHOH20 << "\n";
    }//end for N_even_E
    for (int i=0;i<12;i++) {
      thisenergy=pow(10., (spectra1->Getenergy())[0]+((double)i)*0.5);
      primary1->GetSigma(thisenergy, sigma, thislen_int_kgm2, settings1, 0, 0);
      // EdNdEdAdt is in #/cm^2/s
      // can also be written dN/d(lnE)dAdt
      // = dN*log(10)/d(log E)dAdt
      // the bin spacing is 0.5
      // so # events ~ dN*log(10)*0.5/d(log E)dAdt
      sum_events+=0.5*log(10.)*(spectra1->GetEdNdEdAdt())[i]/(thislen_int_kgm2/sig1->RHOH20);
      //cout << "thisenergy,  EdNdEdAdt is " << thisenergy << " " << spectra1->EdNdEdAdt[i] << "\n";
      //foutput << "interaction length is " << thislen_int_kgm2/RHOH20 << "\n";
    } //end for i
    //km3sr=volume*pow(1.E-3, 3)*sig1->RHOMEDIUM/RHOH20*sr*nevents/(double)NNU;
    sum_events*=volume*anita1->LIVETIME*sig1->RHOMEDIUM/sig1->RHOH20*nevents/(double)NNU*sr*pow(10., 4);
    foutput << "volume,  LIVETIME,  sig1->RHOMEDIUM,  RHOH20,  nevents,  NNU,  sr are " << volume << " " << anita1->LIVETIME << " " << sig1->RHOMEDIUM << " " << sig1->RHOH20 << " " << nevents << " " << NNU << " " << sr << "\n";
    foutput << "Total events observed is " << sum_events << "\n";
  } //end if IsSpectrum

}
//end Summarize()


double GetAirDistance(double altitude_bn, double beta) { // given beta=angle wrt horizontal that the ray hits the balloon,  calculate distance that the ray traveled in air,  including curvature of earth
  return EarthModel::R_EARTH*acos((altitude_bn+EarthModel::R_EARTH)/EarthModel::R_EARTH*(1-sin(beta)*sin(beta))+1/EarthModel::R_EARTH*sin(beta)*sqrt((altitude_bn+EarthModel::R_EARTH)*(altitude_bn+EarthModel::R_EARTH)*sin(beta)*sin(beta)-2*EarthModel::R_EARTH*altitude_bn-altitude_bn*altitude_bn));
}
//end GetAirDistance()


double GetAverageVoltageFromAntennasHit(Settings *settings1, int *nchannels_perrx_triggered, double *voltagearray, double& volts_rx_sum) {
  double sum=0;
  int count_hitantennas=0;
  for (int i=0;i<settings1->NANTENNAS;i++) {
    if (nchannels_perrx_triggered[i]>=3) {
      sum+=voltagearray[i];
      count_hitantennas++;
    } //if
  } //for
  volts_rx_sum = sum;
  sum = sum/(double)count_hitantennas;
  return sum;
}
//end GetAverageVoltageFromAntennasHit()


Vector GetPolarization(const Vector &nnu, const Vector &nrf2_iceside) {
  // Want to find a unit vector in the same plane as nnu and n_refr,
  // but perpendicular to n_refr,  pointing away from nnu.
  
  // cross nnu with n_refr to get the direction of the B field.
  Vector n_bfield = nnu.Cross(nrf2_iceside);
  // cross b-field with nrf2_iceside to get the polarization vector.
  Vector n_pol = n_bfield.Cross(nrf2_iceside);
  n_pol = n_pol.Unit();
  // check and make sure E-field is pointing in the right direction.
  if (nnu*nrf2_iceside>0 && n_pol*nnu>0){
    cout << "error in GetPolarization.  Event is " << inu << "\n";
  }
  return n_pol;
}
//end GetPolarization()



void CloseTFile(TFile *hfile) {
  hfile->cd();
  hfile->Write();
  hfile->Close();
}
//end CloseTFile()


double IsItDoubleBang(double exitlength, double plepton) {
  double gamma=plepton/MTAU;
  return 1-exp(-1*exitlength/(TAUDECAY_TIME*CLIGHT*gamma));
}
//end IsItDoubleBang()


int WhereIsSecondBang(const Position &posnu, const Vector &nnu, double nuexitlength, double pnu, IceModel *antarctica1, const Position &r_bn, Position &posnu2, Position &rfexit_db, Vector &n_exit2bn_db) {
  double rnd1=0;
  double rnd2=2;
  double gamma=pnu/MTAU;
  
  if (exp(-1*nuexitlength/(TAUDECAY_TIME*CLIGHT*gamma))>0.999){
    rnd1=gRandom->Rndm()*nuexitlength;
  }
  else {
    while (rnd2>1-exp(-1*rnd1/(TAUDECAY_TIME*CLIGHT*gamma))) {
      rnd1=gRandom->Rndm()*nuexitlength;
      rnd2=gRandom->Rndm();
    } //while
  } //else
  posnu2 = posnu + rnd1*nnu;
  rfexit_db = antarctica1->Surface(posnu2)*posnu2.Unit();
  
  // unit vector pointing to antenna from exit point.
  n_exit2bn_db = (r_bn - rfexit_db) / r_bn.Distance(rfexit_db);
  
  double cosangle=(n_exit2bn_db * posnu2) / posnu2.Mag();
  if (cosangle<0){
    return 0;
  }
  return 1;
}
//end WhereIsSecondBang()


//the following is  a new function only for reflected case.
void Attenuate_down(IceModel *antarctica1, Settings *settings1, double& vmmhz_max, const Position &rfexit2, const Position &posnu, const Position &posnu_down) {
    double ATTENLENGTH=700;
    if(!settings1->VARIABLE_ATTEN){
      ATTENLENGTH=antarctica1->EffectiveAttenuationLength(settings1, posnu, 1);
    }

    int position_in_iceshelves=antarctica1->IceOnWater(posnu);
    int position_in_rossexcept=antarctica1->RossExcept(posnu);
    // int position_in_ross = antarctica->RossIceShelf(posnu);
    // int position_in_ronne = antarctica->RonneIceShelf(posnu);
    double dtemp=posnu_down.Distance(rfexit2)/ATTENLENGTH;
    
    if (dtemp<20) {
      // if(position_in_ross || position_in_ronne) {
      if(position_in_iceshelves && (!position_in_rossexcept)){
        // scalefactor_attenuation=0.310227766*exp(-dtemp);
        // vmmhz_max=vmmhz_max*exp(-dtemp)*0.310227766;//10% of power reflected
        scalefactor_attenuation=0.71*exp(-dtemp);
        vmmhz_max=vmmhz_max*0.71*exp(-dtemp);//50% of power reflected. -3dB
      } //end if
      else if(position_in_rossexcept){
        scalefactor_attenuation=0.1*exp(-dtemp);
        vmmhz_max=0.1*vmmhz_max*exp(-dtemp);//1% of power reflected. -20dB
      }//end else if
      else {
        scalefactor_attenuation=sqrt(0.001)*exp(-dtemp);
        vmmhz_max=sqrt(0.001)*vmmhz_max*exp(-dtemp);//0.1% of power reflected.-30dB
      } //else
    } //if
    else {
      scalefactor_attenuation=0;
      vmmhz_max=0;
    } //else
}
//end Attenuate_down()


void Attenuate(IceModel *antarctica1, Settings *settings1, double& vmmhz_max,  double rflength, const Position &posnu) {
    double ATTENLENGTH=700;  // constant attenuation length for now.
    if (!settings1->VARIABLE_ATTEN){
      ATTENLENGTH = antarctica1->EffectiveAttenuationLength(settings1, posnu, 0);
    }
    
    double dtemp=(rflength/ATTENLENGTH);
    if(!settings1->ROUGHNESS){
      if (dtemp<20) {
        scalefactor_attenuation=exp(-dtemp);
        vmmhz_max=vmmhz_max*exp(-dtemp);
      } //if
      else {
        scalefactor_attenuation=0;
        vmmhz_max=0;
      } //else
    }
    else{       // use a larger allowable value in case of roughness
      if (dtemp<10000) {
        scalefactor_attenuation=exp(-dtemp);
        vmmhz_max=vmmhz_max*exp(-dtemp);
      } //if
      else {
        scalefactor_attenuation=0;
        vmmhz_max=0;
      } //else
    }
}
//end Attenuate()


void IsAbsorbed(double chord_kgm2, double len_int_kgm2, double &weight1) {
    // see if neutrino is absorbed
    //  weighting works,  but not to much purpose since nu's always
    //   interact at these energies.
    double rtemp;
    
    rtemp=chord_kgm2/len_int_kgm2;
    if (rtemp<=20)
      weight1=exp(-rtemp);
    else
      weight1=0;
}
//end IsAbsorbed()


/*Position WhereDoesItEnter(IceModel *antarctica1, const Position &posnu, const Vector &nnu) {
    // now get neutrino entry point...
    double p = posnu.Mag(); // radius of interaction
    double costheta = (nnu*posnu) / p; // theta of neutrino at interaction position
    double sintheta = sqrt(1-costheta*costheta);
    
    double lon = posnu.Lon();
    double lat = posnu.Lat();
    
    double a=0; // length of chord
    
    double R = antarctica1->Surface(lon, lat);
    double delta = R - p; // depth of the interaction
    // if interaction occurs below surface,  as it should
    
    if (delta>-0.001) {
		a=p*costheta+sqrt(R*R*costheta*costheta+2*delta*R*sintheta*sintheta); // chord length
		if (a<0) {
			cout << "Negative chord length: " << a << "\n";
		} //end if
    } //end if (interaction below surface)
    else if (delta<=-0.001) {
		
		cout << "Error in interaction position.  whichray is " << whichray << "\n";
		cout << "inu,  lon,  lat from WhereDoesItEnter is " << inu << " " << lon << " " << lat << "\n";
		cout << "geoid,  surface,  p,  surface-p are " << antarctica1->Geoid(lat) << " " << antarctica1->Surface(lon, lat) << " " << p << " ,  "<<(antarctica1->Surface(lon, lat)-p)<<"\n";
		
    } //else if: error: interaction takes place above the surface
    
    Position r_in = posnu - a*nnu;
    
    lon = r_in.Lon();
    lat = r_in.Lat();
    
    r_in = antarctica1->Surface(lon, lat) * r_in.Unit();
    
    return r_in;
} //method WhereDoesItEnter
*/


void GetSmearedIncidentAngle(Vector &specular, Vector &nrf_iceside, Vector &n_exit2bn, double SMEARINCIDENTANGLE){
  //  void GetSmearedIncidentAngle(Vector &specular, Vector &nsurf_rfexit, Vector &nrf_iceside, Vector &n_exit2bn, double SMEARINCIDENTANGLE, double theta_inc_smeared) {
    // Smear the incident angle for roughness studies
    specular+=nrf_iceside; // specular is the ray that we got from Snell's law
    Vector parallel_to_surface; // find vector parallel to surface to rotate the vector around
    parallel_to_surface+=n_exit2bn; // want to cross specular with n_exit2bn
    parallel_to_surface.Cross(specular);
    nrf_iceside.Rotate(SMEARINCIDENTANGLE*(2*gRandom->Rndm()-1.), parallel_to_surface); // smear the incident ray
    //   theta_inc_smeared=acos(nrf_iceside.Dot(nsurf_rfexit));
}
//end GetSmearedIncidentAngle()


int GetRayIceSide(const Vector &n_exit2rx,  const Vector &nsurf_rfexit, double nexit,  double nenter,  Vector &nrf2_iceside) {
    // this function performs snell's law in three dimensions
    double costh=0;
    double NRATIO=nexit/nenter;
    costh=(n_exit2rx*nsurf_rfexit)/(n_exit2rx.Mag() * nsurf_rfexit.Mag()); // cos(theta) of the transmission angle
    
    if (costh<0) {
      //cout << "returning 0.  inu is " << inu << "\n";
      return 0;
    }
    double sinth=sqrt(1 - costh*costh);
    double factor=NRATIO*costh-sqrt(1-(NRATIO*sinth*NRATIO*sinth));
    nrf2_iceside = -factor*nsurf_rfexit + NRATIO*n_exit2rx;
    nrf2_iceside = nrf2_iceside.Unit(); // normalize
    return 1;
}
//end GetRayIceSide()


int GetDirection(Settings *settings1, Interaction *interaction1, const Vector &refr,  double deltheta_em,  double deltheta_had, double emfrac,  double hadfrac,  double vmmhz1m_max,  double r_fromballoon,  Ray *ray1,  Signal *sig1,  Position posnu,  Anita *anita1,  Balloon *bn1, Vector &nnu,  double& costhetanu,  double& theta_threshold) {
  int dont_count=0;
  double theta_test=0;
  double vmmhz1m_test=0;
  double costhetanu1 = 0;
  double costhetanu2 = 0;
  
  if (bn1->WHICHPATH==3) { //To make a banana plot,  force neutrino direction
    nnu = interaction1->nnu_banana;
    theta_threshold = 0; //not used for anything in banana plots
    return 1;
  } //if (make banana plot)
  
  if (settings1->SKIPCUTS || !settings1->USEDIRECTIONWEIGHTS) { // this is a setting that allows all neutrino angles,  no restriction.  Makes the code slower.
    costhetanu2=1.;
    costhetanu1=-1.;
    theta_threshold=1;
  } //end if (settings1->SKIPCUTS || !USEWEIGHTS)
  else {
    if (emfrac<=1.E-10 && deltheta_had >1.E-10) {
      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(hadfrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(sig1->changle)>1)
        //if (Tools::dMax(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(hadfrac*vmmhz1m_max*heff_max*bw/1.E6)*sin(sig1->changle)>1)
        theta_threshold=-1;
      else {
        //theta_threshold=sqrt(-1*deltheta_had*deltheta_had*log(Tools::dMax(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(hadfrac*vmmhz1m_max*heff_max*bw/1.E6)*sin(sig1->changle))/ALOG2);
        theta_threshold=sqrt(-1*deltheta_had*deltheta_had*log(anita1->VNOISE[0]/10.*anita1->maxthreshold/(hadfrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(sig1->changle))/ALOG2);
        averaging_thetas1+=theta_threshold;
      } //else
      count_inthisloop1++;
    } //if

    if (emfrac>1.E-10 && deltheta_had <=1.E-10) {
      dont_count++;
      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(emfrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(sig1->changle)>1)
      //if (Tools::dMax(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(emfrac*vmmhz1m_max*heff_max*bw/1.E6)*sin(sig1->changle)>1)
        theta_threshold=-1;
      else {
        //theta_threshold=sqrt(-1*deltheta_em*deltheta_em*log(Tools::dMax(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(emfrac*vmmhz1m_max*heff_max*bw/1.E6)*sin(sig1->changle))/0.5);
        theta_threshold=sqrt(-1*deltheta_em*deltheta_em*log(anita1->VNOISE[0]/10.*anita1->maxthreshold/(emfrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(sig1->changle))/0.5);
        averaging_thetas2+=theta_threshold;
      } //else
      count_inthisloop2++;
    } //if


    //start big code block of ifs/elses 
    if (emfrac>1.E-10 && deltheta_had>1.E-10) {
    // if the electromagnetic and hadronic components of the shower are both non-negligible
    // then theta_threshold cannot be determined analytically so we step away from the cerenkov angle in steps equal to 1/2 * deltheta_em
      if (anita1->VNOISE[0]/10.*anita1->maxthreshold/((hadfrac+emfrac)*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.) {
      //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/((hadfrac+emfrac)*vmmhz1m_max*heff_max*bw/1.E6)>1.) {
        theta_threshold=-1.; // if it's not detectable at all
      }
      else { // otherwise,  start stepping.
        theta_test=deltheta_em; // this is the angle we start stepping at
        vmmhz1m_test=vmmhz1m_max; // this will be the magnitude of the signal at theta_test away from the cerenkov cone.
        // find the magnitude of the signal at theta_test away from the cerenkov cone.
        sig1->TaperVmMHz(sig1->changle+theta_test, deltheta_em, deltheta_had, emfrac, hadfrac, vmmhz1m_test, djunk);
        //	if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.) { // is this electric field already too low to have a chance of passing the trigger threshold?
        if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.) { // is this electric field already too low to have a chance of passing the trigger threshold?
          theta_threshold=theta_test; // then that is the maximum angular deviation
        }
        else { // otherwise increment by the step size and check again.
          theta_test=1.5*deltheta_em;
          vmmhz1m_test=vmmhz1m_max;
          sig1->TaperVmMHz(sig1->changle+theta_test, deltheta_em, deltheta_had, emfrac, hadfrac, vmmhz1m_test, djunk);

          if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.) {
          //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.) {
            theta_threshold=theta_test;
          }
          else { // otherwise increment by the step size and check again.
            theta_test=2*deltheta_em;
            vmmhz1m_test=vmmhz1m_max;
            sig1->TaperVmMHz(sig1->changle+theta_test, deltheta_em, deltheta_had, emfrac, hadfrac, vmmhz1m_test, djunk);
            //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.)
            if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.)
              theta_threshold=theta_test;
            else { // otherwise increment by the step size and check again.
              theta_test=3*deltheta_em;
              vmmhz1m_test=vmmhz1m_max;
              sig1->TaperVmMHz(sig1->changle+theta_test, deltheta_em, deltheta_had, emfrac, hadfrac, vmmhz1m_test, djunk);
              //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.)

              if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.)
                theta_threshold=theta_test;
              else { // otherwise,  set is the the width of the hadronic component (much wider than the electromagnetic component)
                theta_test=deltheta_had;
                vmmhz1m_test=vmmhz1m_max;
                sig1->TaperVmMHz(sig1->changle+theta_test, deltheta_em, deltheta_had, emfrac, hadfrac, vmmhz1m_test, djunk);
                // if at the hadronic width,  you're below the threshold
                if (anita1->VNOISE[0]/10.*anita1->maxthreshold/(vmmhz1m_test/r_fromballoon*heff_max*anita1->bwmin/1.E6)>1.)
                //if (Tools::dMin(VNOISE, settings1->NLAYERS)*anita1->maxthreshold/(vmmhz1m_test*heff_max*bw/1.E6)>1.) // if at the hadronic width,  you're below the threshold
                  theta_threshold=theta_test; // set theta_threshold
                else { // otherwise,  find theta_threshold considering the hadronic component alone.  This is conservative-- an electromagnetic component would only make it narrower.
                  theta_threshold=sqrt(-1*deltheta_had*deltheta_had*log(anita1->VNOISE[0]/10.*anita1->maxthreshold/(hadfrac*vmmhz1m_max/r_fromballoon*heff_max*anita1->bwmin/1.E6)*sin(sig1->changle))/0.5);
                } // else: not below threshold at deltheta_had
              } // else: not below threshold at 3*deltheta_em

            } // else: not below threshold at 2.5*deltheta_em
          } // else: not below threshold at 2.0*deltheta_em

        } // else: not below threshold at 1.5*deltheta_em

      } // not below threshold at 1.0*deltheta_em
      count_inthisloop3++;
      averaging_thetas3+=theta_threshold;

    } // end if both the em and hadronic components are non-negligible.
    //end big code block of ifs/elses

        
    theta_threshold*=settings1->THETA_TH_FACTOR; // multiply theta_threshold by scale factor if requested,  for testing purposes.
    if (theta_threshold>0) { // we only pick the angle between 0 and pi so set the upper and lower limits accordingly.
      if (sig1->changle-theta_threshold<0 && sig1->changle+theta_threshold> PI) {
        costhetanu2=1.;
        costhetanu1=-1.;
      } //if
      else if (sig1->changle-theta_threshold>0 && sig1->changle+theta_threshold> PI) {
        costhetanu2=cos(sig1->changle-theta_threshold);
        costhetanu1=-1.;
      } //else if
      else if (sig1->changle-theta_threshold<0 && sig1->changle+theta_threshold< PI) {
        costhetanu2=1.;
        costhetanu1=cos(sig1->changle+theta_threshold);
      } //else if
      else if (sig1->changle-theta_threshold>0 && sig1->changle+theta_threshold< PI) {
        costhetanu2=cos(sig1->changle-theta_threshold);
        costhetanu1=cos(sig1->changle+theta_threshold);
      } //else if
    } // end if theta_threshold>0


  } // if SKIP_CUTS !=0
    
  if (theta_threshold>0) {
    // pick the neutrino direction,  in a coordinate system where the z axis lies along the cerenkov cone.
    costhetanu=costhetanu1+gRandom->Rndm()*(costhetanu2-costhetanu1);

    double phinu=TWOPI*gRandom->Rndm(); // pick the phi of the neutrino direction,  in the same coordinate system.
    double sinthetanu=sqrt(1-costhetanu*costhetanu);
    // 3-vector of neutrino direction,  at that same coordinate system.
    nnu = Vector(sinthetanu*cos(phinu), sinthetanu*sin(phinu), costhetanu);
    nnu = nnu.ChangeCoord(refr); // rotate so it's in our normal coordinate system.
    // now the ray is aligned along the cerenkov cone and
    // the neutrino is rotated by that same angle

    //dtryingdirection+=4*PI/(2.*theta_threshold*sin(sig1->changle)*2*PI);
    interaction1->dtryingdirection=1/((costhetanu2-costhetanu1)/2.);
    if (bn1->WHICHPATH==4) {
      //double angle=(PI/2.-sig1->changle)-ray1->rfexit[0].Angle(ray1->nrf_iceside[4])+1.*RADDEG;
      double angle=(PI/2.-sig1->changle)-ray1->rfexit[0].Angle(ray1->nrf_iceside[4]);      // this will put the viewing angle at the cerenkov angle
      double thetaposnu=posnu.Theta();
      //double phiposnu=posnu.Phi();
      costhetanu=cos(PI/2+(thetaposnu-angle));
      sinthetanu=sqrt(1-costhetanu*costhetanu);
      //phinu=0.95993;
      phinu=-1.339; // this is the phi where it's coming *from.*
      // we want the neutrino to be headed north
      nnu = Vector(-1.*sinthetanu*cos(phinu), -1.*sinthetanu*sin(phinu), -1.*costhetanu);// 3-vector of neutrino direction,  at that same coordinate system.
    }
    return 1;
  } //end if theta_threshold
  else if (theta_threshold==-1.) {
    cout << "theta_threshold is " << theta_threshold << "\n";
    return 0;
  }
  else if (emfrac<=1.E-10 && deltheta_had <= 1.E-10) {
    cout << "Error:  emfrac, hadfrac are (1st place)" << emfrac << " " << hadfrac << " " << "\n";
    return 0;
  } //else if
  
  return 0;
}
//end GetDirection()


double ScaleVmMHz(double vmmhz1m_max, const Position &posnu1, const Position &r_bn) {
  double dtemp= r_bn.Distance(posnu1);
  vmmhz1m_max= vmmhz1m_max/dtemp;
  scalefactor_distance=1/dtemp;
  //cout << "dtemp is " << dtemp << "\n";
  return vmmhz1m_max;
}
//end ScaleVmMHz()


void SetupViewangles(Signal *sig1) {
  double viewangle_max=90.*RADDEG;
  double viewangle_min=30.*RADDEG;
  for (int i=0;i<NVIEWANGLE-2;i++) {
    viewangles[i]=viewangle_max-(viewangle_max-viewangle_min)/(double)(NVIEWANGLE-2)*(double)i;
  }
  viewangles[NVIEWANGLE-2]=acos(1/sig1->N_DEPTH);
  viewangles[NVIEWANGLE-1]=90.*RADDEG;
}
//end SetupViewAngles()


double GetThisAirColumn(Settings* settings1,  Position r_in, Vector nnu, Position posnu,  double *col1,  double& cosalpha, double& mytheta, double& cosbeta0, double& mybeta) {
  double myair=0; // this is the output
  // it is the column of air in kg/m^2
  cosalpha=(r_in * nnu) / r_in.Mag(); // cosangle that the neutrino enters the earth wrt surface normal at its entrry point
  mytheta=(double)(acos(cosalpha)*DEGRAD)-90.; // turn this into an angle
  
  //------------------added on Dec 8------------------------
  if (settings1->ATMOSPHERE) {
    int index11=int(mytheta*10.); // which index this theta corresponds to
    int index12=index11+1;
    // find column of air at this theta
    myair=(col1[index11]+(col1[index12]-col1[index11])*(mytheta*10.-double(index11)))*10.;//unit is kg/m^2
  }
  else
    myair=0.;//don't include effect of atmosphere
  
  //cout<<"mytheta="<<mytheta<<"; myair="<<myair<<endl;
  //------------------added on Dec 8------------------------
  cosbeta0= (posnu * nnu) / posnu.Mag(); // cos angle of neutrino wrt person standing over the interaction point
  mybeta=(double)(acos(cosbeta0)*DEGRAD)-90.; // turn that into a theta
  return myair;
}
//end GetThisAirColumn()


void GetAir(double *col1) {
  double nothing;
  ifstream air1("data/atmosphere.dat"); // length of chord in air vs. theta (deg)
  //where theta is respect to "up"
  // binned in 0.1 degrees
  for(int iii=0;iii<900;iii++) {
    air1>>nothing>>col1[iii];
  } // read in chord lengths
}
//end GetAir()


int TIR(const Vector &n_surf, const Vector &nrf2_iceside,  double N_IN, double N_OUT) {
  double test=sin(nrf2_iceside.Angle(n_surf))*N_IN/N_OUT;
  if(test>=1)
    return 1;
  else
    return 0;
  
  return 0;
}
//end TIR()


double GetViewAngle(const Vector &nrf2_iceside, const Vector &nnu) {
  // get viewing angle of shower
  double dtemp=nrf2_iceside*nnu;
  if (dtemp>=1 && dtemp<1.02)
    dtemp=0.999999;
  if (dtemp<=-1 && dtemp>-1.02)
    dtemp=-0.9999999;
  
  return acos(dtemp);
}
//end ())etViewAngle


TStyle* RootStyle() {
  TStyle *RootStyle = new TStyle("Root-Style", "The Perfect Style for Plots ;-)");
  
#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                  // save the global style reference
  gStyle = RootStyle;
#endif
  // otherwise you need to call TROOT::SetStyle("Root-Style")
  
  // Paper size
  RootStyle->SetPaperSize(TStyle::kUSLetter);
  
  // Canvas
  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);
  
  // Pads
  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (10);
  RootStyle->SetPadBorderMode  (0);
  //  RootStyle->SetPadBottomMargin(0.13);
  RootStyle->SetPadBottomMargin(0.16);
  RootStyle->SetPadTopMargin   (0.08);
  RootStyle->SetPadLeftMargin  (0.18);
  RootStyle->SetPadRightMargin (0.05);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);
  
  // Frames
  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);
  
  
  // Histograms
  RootStyle->SetHistFillColor(0);
  RootStyle->SetHistFillStyle(1);
  RootStyle->SetHistLineColor(1);
  RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(2);
  
  // Functions
  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(1);
  
  //Legends
  RootStyle->SetStatBorderSize(2);
  RootStyle->SetStatFont      (42);
  //  RootStyle->SetOptStat       (111111);
  RootStyle->SetOptStat       (0);
  RootStyle->SetStatColor     (0);
  RootStyle->SetStatX         (0.93);
  RootStyle->SetStatY         (0.90);
  RootStyle->SetStatFontSize  (0.07);
  //  RootStyle->SetStatW         (0.2);
  //  RootStyle->SetStatH         (0.15);
  
  // Labels,  Ticks,  and Titles
  RootStyle->SetTickLength ( 0.015, "X");
  RootStyle->SetTitleSize  ( 0.10, "X");
  RootStyle->SetTitleOffset( 1.20, "X");
  RootStyle->SetTitleBorderSize(0);
  //  RootStyle->SetTitleFontSize((double)3.);
  RootStyle->SetLabelOffset( 0.015, "X");
  RootStyle->SetLabelSize  ( 0.050, "X");
  RootStyle->SetLabelFont  ( 42   , "X");
  RootStyle->SetTickLength ( 0.015, "Y");
  RootStyle->SetTitleSize  ( 0.10, "Y");
  RootStyle->SetTitleOffset( 0.600, "Y");
  RootStyle->SetLabelOffset( 0.015, "Y");
  RootStyle->SetLabelSize  ( 0.050, "Y");
  RootStyle->SetLabelFont  ( 42   , "Y");
  RootStyle->SetTitleFont  (42, "XY");
  RootStyle->SetTitleColor  (1);
  
  // Options
  RootStyle->SetOptFit     (0);
  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(0.4);
  
  //  cout << ">> Style initialized with the Root Style!" << endl;
  //  cout << ">> " << modified << endl << endl;
  return RootStyle;
}
//end RootStyle()


void GetFresnel(Roughness *rough1, int ROUGHNESS_SETTING, const Vector &surface_normal, 
  const Vector &air_rf, 
  Vector &pol, 
  const Vector &firn_rf, 
  double efield, 
  double emfrac, double hadfrac, double deltheta_em_max, double deltheta_had_max, 

  double &t_coeff_pokey, double &t_coeff_slappy, 
  double &fresnel, 
  double &mag) {

  // find angle of incidence and angle of transmission
  double incident_angle = surface_normal.Angle(firn_rf);
  double transmitted_angle = surface_normal.Angle(air_rf);
  
  //  double t_coeff_pokey, t_coeff_slappy;
  
  // this is perp the surface normal and transmitted ray,  parallel to surface
  Vector perp = air_rf.Cross(surface_normal).Unit();
  // this is in the bending plane
  Vector air_parallel = perp.Cross(air_rf).Unit();
  // this is in the bending plane
  Vector firn_parallel = perp.Cross(firn_rf).Unit();
  
  // component of polarization (in the air) perp to surface normal
  double pol_perp_firn = pol*perp; // this is the slappy component in the firn
  double pol_parallel_firn = pol*firn_parallel; // this is the pokey component in the firn
  double pol_perp_air=0, pol_parallel_air=0;
  
  if (ROUGHNESS_SETTING==0) {
    double r_coeff_pokey =  tan(incident_angle - transmitted_angle) / tan(incident_angle + transmitted_angle);
    t_coeff_pokey = sqrt((1. - r_coeff_pokey*r_coeff_pokey));
    pol_parallel_air = t_coeff_pokey * pol_parallel_firn; // find pokey component in the air

    double r_coeff_slappy = sin(incident_angle - transmitted_angle) / sin(incident_angle + transmitted_angle);
    t_coeff_slappy = sqrt((1. - r_coeff_slappy*r_coeff_slappy));
    pol_perp_air = t_coeff_slappy * pol_perp_firn; // find slappy component in the firn

    mag=sqrt( tan(incident_angle) / tan(transmitted_angle) );

    fresnel = sqrt( pow(efield * pol_perp_air, 2) + pow(efield * pol_parallel_air, 2)) / efield;

    //cout << "incident_angle,  transmitted_angle,  mag,  fresnel,  and from functions " << incident_angle << " " << transmitted_angle << " " << mag*fresnel << " " << 1. << " " << (rough1->GetPokey(incident_angle, transmitted_angle) + rough1->GetSlappy(incident_angle, transmitted_angle))/2. << "\n";

    //if (WHICHPATH==4) {
    //  cout << "t_coeff_pokey*mag is " << t_coeff_pokey*mag << "\n";
    //  cout << "t_coeff_slappy*mag is " << t_coeff_slappy*mag << "\n";
    //       }
  }//end roughness setting==0
  //else if (ROUGHNESS_SETTING==2) {
  //  mag=sqrt(tan(asin(sin(transmitted_angle)/Signal::NICE))/tan(transmitted_angle));
  //  t_coeff_pokey=rough1->GetPokey(incident_angle, transmitted_angle, emfrac, hadfrac, deltheta_em_max, deltheta_had_max);
  //  pol_parallel_air = t_coeff_pokey * pol_parallel_firn; // find pokey component in the air
  //  t_coeff_slappy=rough1->GetSlappy(incident_angle, transmitted_angle, emfrac, hadfrac, deltheta_em_max, deltheta_had_max);
  //  pol_perp_air = t_coeff_slappy * pol_perp_firn; // find slappy component in the firn
  //  fresnel = sqrt( pow(efield * pol_perp_air, 2) + pow(efield * pol_parallel_air, 2)) / efield;
  //}
  //std::cout << "In fren : " << fresnel << std::endl;
  pol = (pol_perp_air * perp + pol_parallel_air * air_parallel).Unit();
}
//end GetFresnel()



void GetBalloonLocation(Interaction *interaction1,Ray *ray1,Balloon *bn1,IceModel *antarctica) {
    // brian enter function to calculate balloon position on your map.
    // use interaction1->posnu // location of neutrino interaction
    // coordinate system:  +z="up" at the south pole
    // bn1->r_bn
    // nnu
    // ray1->nsurf_rfexit
    
    
    // brian enter function to calculate balloon position on your map.
    // use interaction1->posnu // location of neutrino interaction
    // coordinate system:  +z="up" at the south pole
    // bn1->r_bn
    // nnu
    
    
    // balloonvector = balloonvector - nuvector;//change origin to the nuetrino interaction point
    
    const Vector nuvector = interaction1->nnu;
    // double interactiondepth = nuvector[2];//NOT CORRECT! need depth BELOW the ice. this is height above center of earth.
    
    Vector zcoordvector = ray1->nsurf_rfexit;
    zcoordvector=zcoordvector.Unit();
    
    // double thetainc =acos(zcoordvector.Dot(nuvector))*180/PI;
    //nsurf_rfexit is z direction for new coordinate system. Need to make sure the n-vector is in x-z plane.
    
    Vector xcoordvector = nuvector-(zcoordvector.Dot(nuvector))*zcoordvector;//xcoordvector is such that nnu lies in the x-z plane
    xcoordvector = xcoordvector.Unit();
    
    const Vector ycoordvector = zcoordvector.Cross(xcoordvector);//Need this for ChangeCoord. 
    
    
    Vector origin_brian_tmp;
    if (interaction1->nnu.Dot(zcoordvector)>0) // up  going
  origin_brian_tmp=interaction1->nuexit; // the origin is the neutrino exit point 
    else {     
  Vector nnu_flipped=interaction1->nnu;
  nnu_flipped=nnu_flipped-2.*nnu_flipped.Dot(zcoordvector)*zcoordvector; // take it's upgoing reflection from surface
  
  Position nuexit_flipped;
  if (Ray::WhereDoesItLeave(interaction1->posnu,nnu_flipped,antarctica,nuexit_flipped))
      origin_brian_tmp=nuexit_flipped;
  
  
    }
    Vector r_bn_tmp=bn1->r_bn-origin_brian_tmp;
    r_bn_tmp=r_bn_tmp.ChangeCoord(xcoordvector,ycoordvector);//change coordinates
    
    
    double balloondist =r_bn_tmp.Mag();//this is above center of earth, if i understand correctly. Need above the surface of the earth. 
    double balloonphi = r_bn_tmp.Phi(); //phi position of the balloon
    if (balloonphi>PI)
      balloonphi=balloonphi-2*PI;

    double balloontheta = r_bn_tmp.Theta();// costheta position of the baloon
    // get this by dotting ray1->nsurf_rfexit with nnu?     
    // double thetainc = acos(interaction1->nnu[2])*180/PI; //nnu is unit vector; cos(thetainc) = z/r
    balloontheta = PI-balloontheta;//walter.cc uses a pos z as down. this corrects for that.
    
    // define a coordinate system with ray1->nsurf_rfexit defining +z
    // nnu direction defines the x-z plane
    // find balloon position in that coordinate system
    //to get the values from walter.cc we need : E_shower, correlation length, rms height and the em_frac and had_frac. the last
    // two are so we can multiply the number from sky maps by the correct frac and then add the em and hadronic portion together
    // to get the total.
    
}


void interrupt_signal_handler(int sig){
  signal(sig,  SIG_IGN);
  ABORT_EARLY = true;
  return;
}
//end interrupt_signal_handler()


#ifdef ANITA_UTIL_EXISTS
//int GetIceMCAntfromUsefulEventAnt(Anita *anita1,  AnitaGeomTool *AnitaGeom1,  int UsefulEventAnt){
int GetIceMCAntfromUsefulEventAnt(Settings *settings1,  int UsefulEventAnt){
  //int layer_temp = IceMCLayerPosition[UsefulEventIndex][0];
  //int position_temp = IceMCLayerPosition[UsefulEventIndex][1];
  //int IceMCIndex = anita1->GetRx(layer_temp,  position_temp);
  int IceMCAnt = UsefulEventAnt;
  if ((settings1->WHICH==9 || settings1->WHICH==10) && UsefulEventAnt<16) {
    IceMCAnt = (UsefulEventAnt%2==0)*UsefulEventAnt/2 + (UsefulEventAnt%2==1)*(UsefulEventAnt/2+8);
  }
  
  return IceMCAnt;
}
//end GetIceMCAntfromUsefulEventAnt()


#endif
