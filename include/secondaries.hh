////////////////////////////////////////////////////////////////////////////////////////////////
//class Secondaries:
////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef ICEMC_SECONDARIES_HH
#define ICEMC_SECONDARIES_HH

#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "Primaries.h"
#include "TRandom3.h"

#include <vector>
#include "vector.hh"

class TH1F;

namespace icemc{
  class Settings;
  class Primaries;

  //! Secondary interactions
  class Secondaries {

  private:
    TRandom3 Rand3;

    void Picky(double *y_cumulative,int NPROB_MAX,double rnd,double& y);
    // secondaries
    static const int NPROB_MAX=100; // this sets the maximum length of the arrays
    int NPROB; // this counts how many of those array elements are actually used
    double plepton;  // energy of charged lepton after interaction(s)
    double em_secondaries_max;  // em energy due to secondary interactions
    double had_secondaries_max; // had energy due to secondary interactions
    // first index=energy from 10^18 to 10^21 in half-decades
    // second index=y (100 bins)
    double dsdy_muon_brems[7][NPROB_MAX]; // probability distribution vs. y for brem
    double dsdy_muon_epair[7][NPROB_MAX]; // prob. dist. vs. y for pair production
    double dsdy_muon_pn[7][NPROB_MAX]; // photonuclear interactions

    double y_muon_brems[7][NPROB_MAX]; // inelasticity corresponding to each bin,brem
    double y_muon_epair[7][NPROB_MAX]; // same for pair production
    double y_muon_pn[7][NPROB_MAX]; // same for photonuclear

    std::vector<double> vy_muon_brems[7]; // vectors containing inelasticities distributed such
    //that they follow the dsdy curves above.
    std::vector<double> vy_muon_epair[7];
    std::vector<double> vy_muon_pn[7];

    double y_cumulative_muon_brems[7][NPROB_MAX];
    double y_cumulative_muon_epair[7][NPROB_MAX];
    double y_cumulative_muon_pn[7][NPROB_MAX];

    double int_muon_brems[7];  // integral of prob. dist.=# of expected interactions.
    // for each energy
    double int_muon_epair[7]; // same for pair prod.
    double int_muon_pn[7]; // same for photnuclear

    double max_muon_brems; // maximum value of prob. dist., brems
    double max_muon_epair;  // max value of prob. dist., pair prod.
    double max_muon_pn; // same for photonucl.

    double min_muon_brems; // minimum value of prob. dist., brems
    double min_muon_epair;  // min value of prob. dist., pair prod.
    double min_muon_pn; // same for photonucl.

    double dsdy_tauon_brems[7][NPROB_MAX];  // same as above, but for taus
    double dsdy_tauon_epair[7][NPROB_MAX];
    double dsdy_tauon_pn[7][NPROB_MAX];
    double dsdy_tauon_hadrdecay[7][NPROB_MAX]; // hadronic decay of taus
    double dsdy_tauon_edecay[7][NPROB_MAX]; // tau decay to electrons
    double dsdy_tauon_mudecay[7][NPROB_MAX]; // tau decay to muons.

    double y_tauon_brems[7][NPROB_MAX];
    double y_tauon_epair[7][NPROB_MAX];
    double y_tauon_pn[7][NPROB_MAX];
    double y_tauon_hadrdecay[7][NPROB_MAX];
    double y_tauon_edecay[7][NPROB_MAX];
    double y_tauon_mudecay[7][NPROB_MAX];

    double y_cumulative_tauon_brems[7][NPROB_MAX];
    double y_cumulative_tauon_epair[7][NPROB_MAX];
    double y_cumulative_tauon_pn[7][NPROB_MAX];
    double y_cumulative_tauon_hadrdecay[7][NPROB_MAX];
    double y_cumulative_tauon_edecay[7][NPROB_MAX];
    double y_cumulative_tauon_mudecay[7][NPROB_MAX];

    std::vector<double> vy_tauon_brems[7]; // vectors containing inelasticities distributed such
    //that they follow the dsdy curves above.
    std::vector<double> vy_tauon_epair[7];
    std::vector<double> vy_tauon_pn[7];
    std::vector<double> vy_tauon_hadrdecay[7];
    std::vector<double> vy_tauon_edecay[7];
    std::vector<double> vy_tauon_mudecay[7];

    double int_tauon_brems[7];
    double int_tauon_epair[7];
    double int_tauon_pn[7];
    double int_tauon_hadrdecay[7];
    double int_tauon_edecay[7];
    double int_tauon_mudecay[7];

    double max_tauon_brems;
    double max_tauon_epair;
    double max_tauon_pn;
    double max_tauon_hadrdecay;
    double max_tauon_edecay;
    double max_tauon_mudecay;

    double min_tauon_brems; // minimum value of prob. dist., brems
    double min_tauon_epair;  // min value of prob. dist., pair prod.
    double min_tauon_pn; // same for photonucl.
    double min_tauon_hadrdecay; // minimum value of prob. dist., brems
    double min_tauon_edecay;  // min value of prob. dist., pair prod.
    double min_tauon_mudecay; // same for photonucl.


    //   double weight_doublebang=0;
    //   double dtime_doublebang=0;


    static const int N_TAUOLA=10001;
    double tauola[N_TAUOLA][6];//tau decay energy fractions, [0]=tau neutrino, [1]=mu neutrino, [2]=Electron neutrino, [3]=hadrons, [4]=electromagnetic, [5]=muons
    std::ifstream tauolainfile;


    void readData(std::string,std::string,double (*)[NPROB_MAX], double (*)[NPROB_MAX]);
    void ReadSecondaries(); // read in prob. dist. for secondary interactions

    //For Tau Weight & Survival probability equations.
    //n.b. not in SI units.
    //from Tau neutrino propagaiton and tau energy loss 2005 Dutta, Huang, & Reno.
    //Equation 16  &  used in Equation 30.
    double B0,B1,E0;//parameterization using a logarithmic dependence on energy
    //for B, the tau elecromagnetic energy loss parameter.
    double mT;//mass of the Tau in Gev
    double cT;//Tau Decay length in cm
    // double p;//Density of Standard Rock. g/cm^3

    //Used in Connolly Calc 2011.(d_dzPsurvNu())
    //p, the Density of Standard Rock. g/cm^3
    double A; ////constant that sets the total probability to unity
    double Mn;// nucleon/ proton mass in grams,also equal to 0.938 GeV.


  public:
    Secondaries();


    int SECONDARIES;
    int TAUDECAY; // is tau decay counted as a secondary interaction
    double TAUFRAC; //fraction of tau neutrino-cc current events where the primare interaction point is the first bang
    int count_nfb;
    int secondary_e_noncons;

    void GetSecondaries(const Settings *settings1, const std::string& ,double, double&, double&, int&, TH1F*);

    void InitTauola();
    void GetTauDecay(std::string nuflavor,std::string current,std::string& taudecay, double& emfrac_db, double& hadfrac_db);

    void GetEMFracDB(double& emfrac_db, double& hadfrac_db);
    double GetDBViewAngle(const icemc::Vector &refr, const icemc::Vector &nnu);
    //void GetFirstBang(const Position &r_in, const icemc::Vector &nnu, Position &posnu, double len_int_kgm2, double d1, double &nuentrancelength);
    double NFBWeight(double ptau, double taulength);

    int GetEMFrac(const Settings *settings1,
		  std::string nuflavor,
		  std::string current,
		  std::string taudecay,
		  double y,
		  TH1F *hy,
		  double pnu,
		  int inu,

		  double& emfrac,
		  double& hadfrac,
		  int& n_interactions, int taumodes1);

    bool secondbang;
    static const bool interestedintaus=false;

    //string flavors[3]={"nue","numu","nutau"}; // the gps path of the anita-lite flight
    std::string flavors[3]; // the gps path of the anita-lite flight

  }; //class Secondaries
}

#endif
