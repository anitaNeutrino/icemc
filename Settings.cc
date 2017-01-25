#include "TF1.h"
#include <array>
#include "position.hh"
#include "Constants.h"
#include "Settings.h"
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

Settings::Settings() {
  Initialize();

  // Bunch of variables which were global in icemc.cc but are settings:
  FILLRAYTREES=1; // fill tree for each ray in roughness simulation
  SEED=65540;      // random number seed.
  THETA_TH_FACTOR=1.0; // factor to multiply theta_th to check code is working properly
  CHANCEINHELL_FACTOR=1.0; // loosen chance in hell cuts to check code is working properly
  CONSTANTY=0; // whether or not to set y to a constant=0.2
  VARIABLE_ATTEN=0; //  0=depth dependent attenuation length, 1=fixed
  TRIGTYPE=1; //1=Trigger scheme as in the SMEX proposal  or 0= Just at least 8 channels pass with 2.3 sigma signal
  ATMOSPHERE=1;// include atmosphere
  SCALEDOWNLCPRX1=1; // scale down lcp voltage of antenna 1 by sqrt(2)
  SCALEDOWNEPOLRX1=1; // scale down power of e pol. of antenna 1 by factor of 2
  SCALEDOWNHPOLRX1=1; // scale down power of h pol. of antenna 1 by factor of 2
  SCALEDOWNEPOLRX2=1; // scale down power of e pol. of antenna 2 by factor of user's choice
  SCALEFACTOREPOLRX2=0.05; // scale power of e pol. of antenna 2 by this factor
  SCALEDOWNHPOLRX2=1; // scale down power of h pol. of antenna 2 by factor of 2
  EPOLRX2ZERO=1; // lcp channel on anita-lite is not considered for triggering.
  HPOLRX2ZERO=1; // h pol. of 2nd antenna set to zero.
  RCPRX2ZERO=1; // rcp of 2nd antenna set to zero.
  LCPRX2ZERO=1; // lcp of 2nd antenna set to zero.
  FLATSURFACE=0; // Normals of all positions on the surface are straight up.
  WRITEPOSFILE=0; //Write neutrino position information to file
  SKIPCUTS=0; //See every neutrino through to the end - don't make any of the various cuts designed to speed up the program.  (For checking distributions.)
  USEDIRECTIONWEIGHTS=1;// whether or not to restrict the neutrino angle so that the viewing angle is near the cerenkov cone
  SHOWERTYPE=0; // Type of shower for previous option
  // End of the once-global varibles.
  taumodes = 1; //Taumodes =1, taucreated in the rock.
}

void Settings::Initialize() {
  NDISCONES_PASS=3;
  DEBUG=false;                   // debugging option
  outputdir="outputs"; // directory where outputs go
  FREQ_LOW_SEAVEYS=200.E6;
  FREQ_HIGH_SEAVEYS=1200.E6;
  BW_SEAVEYS=FREQ_HIGH_SEAVEYS-FREQ_LOW_SEAVEYS;
  SIGMAPARAM=1;  // Connolly et al. 2011 default cross section parametrization
  YPARAM=1;  // Connolly et al. 2011 default y parametrization
  UNBIASED_SELECTION=1.; // (0) pick neutrino interaction in the ice and neutrino from any direction or (1) choose neutrino interaction point in the horizon on the balloon in the ice and neutrino direction on the cerenkov cone
  SIGMA_FACTOR=1;
}

void Settings::ReadInputs(ifstream &inputsfile, ofstream &foutput, Anita* anita1, Secondaries* sec1, Signal* sig1, Balloon* bn1, Ray* ray1, int& NNU, double& RANDOMISEPOL) {
  // read inputs to the code.
  // See comments in input file

  // extern int NNU;
  // extern double RANDOMISEPOL;

  string number;
  string junk;

  getline(inputsfile,junk);

  foutput << "\n\n";

  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  foutput << "Current date and time are: " << asctime (timeinfo) << "\n";


  getline(inputsfile,junk);
  foutput << junk << "\n";
  // The following inputs have to do with event input/output

  //Fenfang's livetime
  //The LIVETIME I used is 34.78days*0.75

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  NNU=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  EXPONENT=(double)atof(number.c_str());

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  UNBIASED_SELECTION=(double)atof(number.c_str());

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  HIST=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  FILLRAYTREES=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  ONLYFINAL=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  HIST_MAX_ENTRIES=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SEED=(int)atoi(number.c_str());
  cout << "SEED is " << SEED << "\n";
  gRandom->SetSeed(SEED);

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  WRITEPOSFILE=(int)atof(number.c_str());

  if (WRITEPOSFILE==1)
    cout << "Non-default setting:  WRITEPOSFILE= " << WRITEPOSFILE << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  EVENTSMAP=atoi(number.c_str());//draw the events map or not

  getline(inputsfile,junk);
  foutput << junk << "\n";
  // The following inputs have to do with the payload and balloon

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  WHICH=(int)atoi(number.c_str());


  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  NLAYERS=(int)atoi(number.c_str()); // this is number of layers, counting the upper 16 as 2 layers

  if (((WHICH==1 || WHICH==6) && NLAYERS!=4) ||
      (WHICH==0 && NLAYERS!=1) ||
      (WHICH==7 && NLAYERS!=1))
    cout << "Non-default setting:  WHICH= " << WHICH << " and NLAYERS= " << NLAYERS << "\n";

  //When you look at the Anita payload there are 4 layers, with 8,8,16 and 8 antennas each.  But in the trigger, the top two become one layer of 16 antennas.  So that means for Anita 1 and Anita 2, there are one fewer trigger layers than physical layers.
  // anything anita like
  // anita 1 simple, anita 1 kurt, anita 2 kurt, anita 3, satellite
  if (WHICH==2 || WHICH==6 || WHICH==8 || WHICH==9 || WHICH==10)
    anita1->NTRIGGERLAYERS = NLAYERS - 1;
  else
    anita1->NTRIGGERLAYERS=NLAYERS;

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  anita1->INCLINE_TOPTHREE=(double)atof(number.c_str());

  if (anita1->INCLINE_TOPTHREE!=10)
    cout << "Non-default setting:  INCLINE_TOPTHREE= " << anita1->INCLINE_TOPTHREE << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  anita1->INCLINE_NADIR=(double)atof(number.c_str());

  if (anita1->INCLINE_NADIR!=10)
    cout << "Non-default setting:  INCLINE_NADIR= " << anita1->INCLINE_NADIR << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  bn1->WHICHPATH=(int)atoi(number.c_str());

  if ((WHICH==0 && bn1->WHICHPATH!=2) || (WHICH==2 && bn1->WHICHPATH!=6))
    cout << "Non-default setting:  bn1->WHICHPATH= " << bn1->WHICHPATH << " and WHICH=" << WHICH << "\n";

  if (bn1->WHICHPATH==2)
    anita1->LIVETIME=45.*24.*3600.*0.75; // 45 days for anita-lite
  else if (bn1->WHICHPATH==0)
    anita1->LIVETIME=6.02*24.*3600.; // anita-lite
  else if (bn1->WHICHPATH==6) // kim's livetime for anita
    anita1->LIVETIME=17.*24.*3600.; // for anita, take 34.78 days * 0.75 efficiency
  else
    anita1->LIVETIME=14.*24.*3600.; // otherwise use 2 weeks by default

  if (WHICH==7) // EeVEX
    anita1->LIVETIME=100.*24.*3600.; // ultra-long duration balloon flight of 100 days

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  bn1->BN_LATITUDE=(double)atof(number.c_str());

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  bn1->BN_LONGITUDE=(double)atof(number.c_str());

  if((bn1->BN_LONGITUDE!=999 || bn1->BN_LATITUDE!=999) && bn1->WHICHPATH==0)
    cout<<"BN_LATITUDE: "<<bn1->BN_LATITUDE<<", BN_LONGITUDE: "<<bn1->BN_LONGITUDE<<endl;

  if (bn1->BN_LONGITUDE>180. && bn1->BN_LONGITUDE!=999)
    cout << "Entered balloon longitude wrong!  Should be between -180 and 180 degrees.\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  bn1->RANDOMIZE_BN_ORIENTATION=(int)atoi(number.c_str());

  if (bn1->RANDOMIZE_BN_ORIENTATION==1 && (bn1->WHICHPATH==2 || bn1->WHICHPATH==6 || bn1->WHICHPATH==7 || bn1->WHICHPATH==8))
    cout << "Warning:: Strangely you asked for a real flight path but a randomized balloon orientation.  WILL BE OVERRIDDEN.\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  bn1->BN_ALTITUDE=(double)atof(number.c_str());

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  // whether to use constant gains as entered in GetBeamWidths (0) or to use Ped's measurements as entered in ReadGains (1)
  // GAINS is actually an int, not a double...
  anita1->GAINS=(int)atof(number.c_str());
  getline(inputsfile,junk);
  foutput << junk << "\n";
  // The following inputs have to do with event input/output

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  TRIGGERSCHEME=atoi(number.c_str()); // whether it's frequency domain (0) or time domain
  TRIGTYPE=1; // ANITA, not anita-lite.  But the anita-lite code back in later

  vector<string> vnumber;
  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);

  for (int n=0;n<5;n++) {
    anita1->bwslice_thresholds[n]=(double)atof(vnumber[n].c_str());
  }
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->BANDING=atoi(number.c_str()); // whether you use anita-1 banding (0) or choose-your-own

  if (anita1->BANDING !=0 && anita1->BANDING!= 1 && anita1->BANDING!=2 && anita1->BANDING!=4) {
    cout << "Banding should be set to 0 (Anita 1), 1 (custum), 2 (Anita 2), 3 (Satellite) or 4 (Anita 3).\n";
    exit(1);
  }

  if ((TRIGGERSCHEME==0 || TRIGGERSCHEME==1) && anita1->BANDING!=1) {
    cout << "Frequency domain trigger schemes can only be used with user-set sub-bands.\n";
    exit(1);
  }

  if (TRIGGERSCHEME==2 && anita1->BANDING==1) {
    cout << "Time domain trigger scheme only works with Anita 1, Anita 2 or Anita 3 banding data, you can't set your own bands.\n";
    exit(1);
  }


  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);
  vector<string> vnumber2;
  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber2,5);

  //  double mintemp,maxtemp;
  for (int n=0;n<5;n++) {
    anita1->bwslice_min[n]=(double)atof(vnumber[n].c_str())*1.E6;
    anita1->bwslice_max[n]=(double)atof(vnumber2[n].c_str())*1.E6;
    anita1->bwslice_center[n]=(anita1->bwslice_min[n]+anita1->bwslice_max[n])/2.;
    anita1->bwslice_width[n]=(anita1->bwslice_max[n]-anita1->bwslice_min[n]);
    //cout << "center, width are " << anita1->bwslice_center[n] << " " << anita1->bwslice_width[n] << "\n";
  }

  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);

  for (int n=0;n<5;n++) {
    //anita1->bwslice is actaully an int
    anita1->bwslice_required[n]=(int)atof(vnumber[n].c_str());
  }

  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,5);

  for (int n=0;n<5;n++) {
    //anita1->bwslice_allowed[n] is still an int, not a double
    anita1->bwslice_allowed[n]=(int)atof(vnumber[n].c_str());
  }

  anita1->maxthreshold=0.;
  anita1->bwmin=1.E10;
  if (anita1->BANDING!=1)
    anita1->bwmin=200.E6;

  for (int n=0;n<5;n++) {
    if (anita1->bwslice_thresholds[n]>anita1->maxthreshold && anita1->bwslice_allowed[n]==1)
      anita1->maxthreshold=anita1->bwslice_thresholds[n];

    if (anita1->BANDING==1) {
      if ((anita1->bwslice_max[n]-anita1->bwslice_min[n])<anita1->bwmin && anita1->bwslice_allowed[n]==1)
	anita1->bwmin=anita1->bwslice_max[n]-anita1->bwslice_min[n];
    }
  }

  // now get number of bands for banding option 3 (satellite)
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->NBANDS=atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->PERCENTBW=atoi(number.c_str());
  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,2);

  anita1->NOTCH_MIN=(double)atof(vnumber[0].c_str())*1.E6; // min and max frequencies for notch filter
  anita1->NOTCH_MAX=(double)atof(vnumber[1].c_str())*1.E6;

  if (anita1->NOTCH_MIN>anita1->NOTCH_MAX) {
    cout << "Min of notch filter is greater than max.  Try again.\n";
  }
  if (anita1->NOTCH_MIN!=0 || anita1->NOTCH_MAX!=0)
    cout << "Applying a notch filter from " << anita1->NOTCH_MIN << " Hz to " << anita1->NOTCH_MAX << " Hz\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->trigRequirements[0]=atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->trigRequirements[1]=atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  antennaclump=atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->REQUIRE_CENTRE=atoi(number.c_str()); // require centre antenna in clump is one of those his
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->trigRequirements[2]=atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  LCPRCP=atoi(number.c_str());


  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,2);

  for (int n=0;n<2;n++) {
    //anita1->bwslice is actaully an int
    anita1->pol_required[n]=(int)atof(vnumber[n].c_str());
  }

  Tools::GetNumbersAsStringArray(inputsfile,foutput,vnumber,2);

  for (int n=0;n<2;n++) {
    //anita1->bwslice_allowed[n] is still an int, not a double
    anita1->pol_allowed[n]=(int)atof(vnumber[n].c_str());
  }


  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  DISCONES=(int)atof(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->INCLUDE_NADIRONLY=(double)atof(number.c_str());

  if (anita1->INCLUDE_NADIRONLY!=0)
    cout << "Non-default setting:  INCLUDE_NADIRONLY= " << anita1->INCLUDE_NADIRONLY << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  CHMASKING=atoi(number.c_str());

  if (bn1->WHICHPATH!=6 && CHMASKING==1) {
    cout << "Cannot include masking for flights other than the ANITA-1 flight. For the ANITA-3 channel masking, it is implemented together with the phi masking and it's turned on whenever the PHIMASKING is ON. CHMASKING set to 0.\n";
    CHMASKING=0;
  }

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  PHIMASKING=atoi(number.c_str());

  getline(inputsfile,junk);
  foutput << junk << "\n";
  // The following are variables that have been used for modeling anita-lite

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNLCPRX1=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNEPOLRX1=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNHPOLRX1=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNEPOLRX2=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SCALEFACTOREPOLRX2=(double)atof(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SCALEDOWNHPOLRX2=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  EPOLRX2ZERO=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  HPOLRX2ZERO=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  RCPRX2ZERO=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  LCPRX2ZERO=(int)atoi(number.c_str());

  if (WHICH==0 && !(SCALEDOWNEPOLRX1==1 && RCPRX2ZERO==1))
    cout << "Non-default setting:  WHICH= " << WHICH << " and EPOLRX2ZERO= " << EPOLRX2ZERO << "\n";

  getline(inputsfile,junk);
  foutput << junk << "\n";
  // Modify the following settings to make changes to the signal and noise
  // for testing

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SIGNAL_FLUCT=(int)atoi(number.c_str());

  if (SIGNAL_FLUCT!=1)
    cout << "Non-default setting:  SIGNAL_FLUCT= " << SIGNAL_FLUCT << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  ZEROSIGNAL=atoi(number.c_str()); // whether it's frequency domain (0) or time domain

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  RANDOMISEPOL=atoi(number.c_str());

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  sig1->SetLPM((int)atoi(number.c_str()));

  if (sig1->GetLPM()!=1)
    cout << "Non-default setting:  LPM= " << sig1->GetLPM() << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  sig1->SetJaime_Factor((double)atof(number.c_str()));


  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  THERMALNOISE_FACTOR=(double)atof(number.c_str());

  if (THERMALNOISE_FACTOR!=1)
    cout << "Non-default setting:  THERMALNOISE_FACTOR= " << THERMALNOISE_FACTOR << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);

  REMOVEPOLARIZATION=(int)atof(number.c_str());

  if (REMOVEPOLARIZATION==1)
    cout << "Non-default setting:  Polarizations turned off!\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->PULSER=atoi(number.c_str());

  if (anita1->PULSER!=0)
    cout << "Warning!  Injecting a pulser spectrum- not simulating neutrinos!  PULSER = " << anita1->PULSER << "\n";


  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  bn1->CENTER=atoi(number.c_str());

  if (bn1->CENTER!=0)
    cout << "WARNING!!  Rotating payload to center one phi sector on the incoming signal for each event.\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  ray1->MAKEVERTICAL=atoi(number.c_str());

  if (ray1->MAKEVERTICAL!=0)
    cout << "WARNING!!  Rotating polarization so it is always vertical approaching the payload.\n";

  getline(inputsfile,junk);
  foutput << junk << "\n";
  // The following are variables that have been used for modeling anita-lite

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SLOPEYSIZE=(double)atof(number.c_str());

  if (SLOPEYSIZE!=0.012)
    cout << "Non-default setting:  SLOPEYSIZE= " << SLOPEYSIZE << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SLOPEY=(int)atoi(number.c_str());

  if (SLOPEY!=1)
    cout << "Non-default setting:  SLOPEY= " << SLOPEY << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  NOFZ=(int)atoi(number.c_str());

  if (NOFZ!=1)
    cout << "Non-default setting:  NOFZ= " << NOFZ << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  VARIABLE_ATTEN=(int)atoi(number.c_str());

  if (VARIABLE_ATTEN!=0)
    cout << "Non-default setting:  VARIABLE_ATTEN= " << VARIABLE_ATTEN << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  CONSTANTICETHICKNESS=(int)atof(number.c_str());

  if (CONSTANTICETHICKNESS==1)
    cout << "Non-default setting:  CONSTANTICETHICKNESS= " << CONSTANTICETHICKNESS << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  ICE_MODEL=(int)atof(number.c_str());

  if ((CONSTANTICETHICKNESS || FIXEDELEVATION) && ICE_MODEL != 0) {
    ICE_MODEL=0;
    cout<<"Constant ice thickness and/or fixed elevation requested.  Using Crust 2.0 ice model.\n";
  } //use the Crust 2.0 data if set to constant icethickness or ground elevation

  if (ICE_MODEL==0)
    cout << "Using Crust 2.0 ice model.\n";
  else if (ICE_MODEL==1)
    cout << "Using BEDMAP ice model.\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  FLATSURFACE=(int)atof(number.c_str());

  if (FLATSURFACE==1)
    cout << "Non-default setting: all surface segments are flat.\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  FIXEDELEVATION=(int)atof(number.c_str());

  if (FIXEDELEVATION==1)
    cout << "Non-default setting:  FIXEDELEVATION= " << FIXEDELEVATION << "\n";

  //ice or salt
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  sig1->SetMedium(atoi(number.c_str()));

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  ROUGHNESS=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  ROUGHSIZE=(double)atof(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  FIRN=atoi(number.c_str());
  if (FIRN==0)
    cout << "Warning!  Non-standard parameter setting.  FIRN = " << FIRN << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  MOOREBAY=atoi(number.c_str());

  getline(inputsfile,junk);
  foutput << junk << "\n";
  // The following are variables that have been used for modeling anita-lite

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SIGMA_FACTOR=(double)atof(number.c_str());

  if (SIGMA_FACTOR!=1)
    cout << "Non-default setting:  settings->SIGMA_FACTOR= " << SIGMA_FACTOR << "\n";


  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  THETA_TH_FACTOR=(double)atof(number.c_str());

  if (THETA_TH_FACTOR!=1)
    cout << "Non-default setting:  THETA_TH_FACTOR= " << THETA_TH_FACTOR << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  CHANCEINHELL_FACTOR=(double)atof(number.c_str());

  if (CHANCEINHELL_FACTOR!=1)
    cout << "Non-default setting:  CHANCEINHELL_FACTOR= " << CHANCEINHELL_FACTOR << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SKIPCUTS=(int)atof(number.c_str());

  if (SKIPCUTS==1)
    cout << "Non-default setting:  Skipping all cuts!\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  USEDIRECTIONWEIGHTS=(int)atof(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  USEPOSITIONWEIGHTS=(int)atof(number.c_str());

  if (USEPOSITIONWEIGHTS==0)
    cout << "Non-default setting:  Not selecting events within the horizon.\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  WEIGHTABSORPTION=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  horizontal_banana_points=(int)atof(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  vertical_banana_points=(int)atof(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  FORSECKEL=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SHOWERTYPE=(int)atoi(number.c_str());
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  BORESIGHTS=atoi(number.c_str());

  if (BORESIGHTS==0)
    cout << "Warning!  Non-standard parameter setting.  BORESIGHTS = " << BORESIGHTS << "\n";

  getline(inputsfile,junk);
  foutput << junk << "\n";
  // Interactions

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  sig1->SetParameterization(atoi(number.c_str()));

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SIGMAPARAM=(int)atoi(number.c_str());

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  YPARAM=(int)atoi(number.c_str());

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  sec1->SECONDARIES=(int)atoi(number.c_str());

  if (sec1->SECONDARIES!=1)
    cout << "Non-default setting:  SECONDARIES= " << sec1->SECONDARIES << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  sec1->TAUDECAY=(int)atoi(number.c_str());

  if (sec1->TAUDECAY!=1)
    cout << "Non-default setting:  TAUDECAY= " << sec1->TAUDECAY << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  ATMOSPHERE=(int)atoi(number.c_str());

  if (ATMOSPHERE!=1)
    cout << "Non-default setting:  ATMOSPHERE= " << ATMOSPHERE << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  CONSTANTCRUST=(int)atof(number.c_str());

  if (CONSTANTCRUST==1)
    cout << "Non-default setting:  CONSTANTCRUST= " << CONSTANTCRUST << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  CONSTANTY=(int)atoi(number.c_str()); // whether to use contant of 0.2 for y (1) yes or (0) no
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  bn1->MAXHORIZON=(int)atof(number.c_str()); // max distance from interaction point to horizon

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  taumodes=(int)atoi(number.c_str()); // tau modes (1 for flat distribution for y)

  getline(inputsfile,junk);
  foutput << junk << "\n";
  // General settings

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  WHICHRAYS=(int)atoi(number.c_str());

  if (WHICHRAYS!=1)
    cout << "Non-default setting:  WHICHRAYS= " << WHICHRAYS << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  WRITE_FILE=(int)atof(number.c_str());

  if (WRITE_FILE) cout<<"Writing CreateHorizons input file.\n";

  //  cout << "WRITE_FILE is " << WRITE_FILE << "\n";

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->SIGMA_THETA=(double)atof(number.c_str());

  if (anita1->SIGMA_THETA==1)
    cout << "Non-default setting:  SIGMA_THETA = 1\n";

  anita1->SIGMA_THETA*=RADDEG; // immediately convert degrees to radians
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->FREQ_LOW=(double)atof(number.c_str());

  if (FREQ_LOW_SEAVEYS>anita1->FREQ_LOW)
    FREQ_LOW_SEAVEYS=anita1->FREQ_LOW;

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  anita1->FREQ_HIGH=(double)atof(number.c_str());
  BW=anita1->FREQ_HIGH-anita1->FREQ_LOW; // total bandwidth of simulation
  getline(inputsfile,junk);
  foutput << junk << "\n";
  // Slac

  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SLAC=atoi(number.c_str());
  // this rotates surface slope 10 degrees away from south pole

  if (SLAC==1)
    cout << "Warning!  Non-standard parameter setting.  SLAC = " << SLAC << "\n";

  if (SLAC) {
    foutput << "!!!!SLAC setting causes some settings to be superseded:\n";
    FIRN=0; // no firn
    foutput << "FIRN=0\n";
    SLOPEY=0; // slopeyness off
    foutput << "SLOPEY=0\n";
    BORESIGHTS=1; // loop over boresights
    foutput << "BORESIGHTS=1\n";
    bn1->BN_ALTITUDE=4.22/0.3; // balloon altitude in ft.!!
    foutput << "BN_ALTITUDE=4.22/0.3\n";
    bn1->RANDOMIZE_BN_ORIENTATION=0; // don't randomize the balloon orientation
    foutput << "RANDOMIZE_BN_ORIENTATION=0\n";
    SKIPCUTS=1; // don't make chance in hell cuts
    foutput << "SKIPCUTS=1\n";
    SLACSLOPE=5.8; // slope of the ice in degrees
    foutput << "SLACSLOPE=5.8\n";
    SLACICELENGTH=5.02; // length of the block of ice
    foutput << "SLACICELENGTH=5.02\n";
  }



  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SLAC_HORIZDIST=(double)atof(number.c_str());
  // horizontal distance from interaction point to payload
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SLACSLOPE=(double)atof(number.c_str());
  // slope of the ice surface
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SLACICELENGTH=(double)atof(number.c_str());
  // length of the block of ice
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  SLAC_HORIZ_DEPTH=(double)atof(number.c_str());
  // horizontal distance from interaction point to surface

  SLAC_DEPTH=tan(SLACSLOPE*RADDEG)*(SLACICELENGTH-SLAC_HORIZ_DEPTH) // height from lowest point of ice
    +21.375*CMINCH/100.; // height from beam to lowest point of ice


  getline(inputsfile,junk);
  foutput << junk << "\n";
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  COHERENT_THRESHOLD = double (atof(number.c_str()));

  // default values are 0
  APPLYIMPULSERESPONSEDIGITIZER=0;
  APPLYIMPULSERESPONSETRIGGER=0;
  USETIMEDEPENDENTTHRESHOLDS=0;

  getline(inputsfile,junk);
  foutput << junk << "\n";
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  APPLYIMPULSERESPONSEDIGITIZER=atoi(number.c_str());
  std::cout << "Apply impulse response to digitizer path: " << APPLYIMPULSERESPONSEDIGITIZER << std::endl;
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  APPLYIMPULSERESPONSETRIGGER=atoi(number.c_str());
  std::cout << "Apply impulse response to trigger path: " << APPLYIMPULSERESPONSETRIGGER << std::endl;

#ifdef ANITA_UTIL_EXISTS
  if ( (APPLYIMPULSERESPONSEDIGITIZER || APPLYIMPULSERESPONSETRIGGER) && WHICH!=8 && WHICH!=9) {
    cout << "Signal chain impulse response is only available for anita-2 and anita-3.\n";
    exit(1);
  }
#endif
#ifndef ANITA_UTIL_EXISTS
  if (APPLYIMPULSERESPONSEDIGITIZER || APPLYIMPULSERESPONSETRIGGER){
    cout << "Signal chain impulse response can only be applied when the Anita tools are sourced.\n";
    exit(1);
  }
#endif
  getline(inputsfile,junk);
  foutput << junk << "\n";
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  USETIMEDEPENDENTTHRESHOLDS=atoi(number.c_str());
  std::cout << "Use time-dependent thresholds: " << USETIMEDEPENDENTTHRESHOLDS << std::endl;

  if ( USETIMEDEPENDENTTHRESHOLDS && WHICH!=9) {
    cout << "Time-dependent thresholds are only available for anita-3.\n";
    exit(1);
  }
  getline(inputsfile,junk);
  foutput << junk << "\n";
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  NOISEFROMFLIGHTDIGITIZER=atoi(number.c_str());
  std::cout << "Use noise from flight for digitizer path: " << NOISEFROMFLIGHTDIGITIZER << std::endl;
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  NOISEFROMFLIGHTTRIGGER=atoi(number.c_str());
  std::cout << "Use noise from flight for trigger path: " << NOISEFROMFLIGHTTRIGGER << std::endl;
#ifdef ANITA_UTIL_EXISTS
  if ( (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER) && WHICH!=9) {
    cout << "Noise from flight only available for anita-3.\n";
    exit(1);
  }
  if (!APPLYIMPULSERESPONSETRIGGER && NOISEFROMFLIGHTTRIGGER ){
    cout << "Noise from flight can only be applied to trigger path if impulse reponse is also used \n";
    exit(1);
  }
#endif
#ifndef ANITA_UTIL_EXISTS
  if (NOISEFROMFLIGHTDIGITIZER || NOISEFROMFLIGHTTRIGGER){
    cout << "Noise from flight can only be applied when the Anita tools are sourced.\n";
    exit(1);
  }
#endif
  Tools::GetNextNumberAsString(inputsfile,foutput,number);
  MINBIAS=atoi(number.c_str());
  if (MINBIAS) std::cout << "Generate Minimum Bias sample: " << MINBIAS << std::endl;




} //method ReadInputs
