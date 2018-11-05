////////////////////////////////////////////////////////////////////////////////////////////////
//class ShowerModel:
////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef ICEMC_SHOWER_MODEL_H
#define ICEMC_SHOWER_MODEL_H

#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "ConnollyEtAl2011.h"
#include "Interaction.h"

#include <vector>
#include "TVector3.h"
#include "RNG.h"
// #include "TGraph2D.h"
#include "TH2F.h"

class TH1F;

namespace icemc{
  class Settings;

  namespace ConnollyEtAl2011 {class CrossSectionModel;}

  /**
   * @class Shower
   * @brief Container for properties that determine the Askaryan frequency content
   */
  class Shower {
  public:
    Shower()
      : // emDeltaThetaMax(0), hadDeltaThetaMax(0),
	emFrac(0), hadFrac(0),
	nInteractions(1)
    {;}

    double sumFrac() const {return emFrac + hadFrac;}
    // double emDeltaThetaMax;
    // double hadDeltaThetaMax;
    double emFrac;
    double hadFrac;
    int nInteractions;
    Energy pnu;
    TVector3 axis;
    ClassDef(Shower, 1)
  };
  

  //! Secondary interactions
  class ShowerModel : public RNG {
  public:
    void Draw(Option_t* opt = "colz");
  private:
    const Settings* fSettings;

    void Picky(const double *y_cumulative, int NPROB_MAX, double rnd, double& y) const;
    
    static const int nP=100; // this sets the maximum length of the arrays
    static const int nE = 7;
    int NPROB; // this counts how many of those array elements are actually used
    // Energy plepton;  // energy of charged lepton after interaction(s)
    // Energy em_secondaries_max;  // em energy due to secondary interactions
    // Energy had_secondaries_max; // had energy due to secondary interactions
    // first index=energy from 10^18 to 10^21 in half-decades
    // second index=y (100 bins)
    std::array<std::array<double, nP>, nE> dsdy_muon_brems; // probability distribution vs. y for brem
    std::array<std::array<double, nP>, nE> dsdy_muon_epair; // prob. dist. vs. y for pair production
    std::array<std::array<double, nP>, nE> dsdy_muon_pn; // photonuclear interactions

    std::array<std::array<double, nP>, nE> y_muon_brems; // inelasticity corresponding to each bin,brem
    std::array<std::array<double, nP>, nE> y_muon_epair; // same for pair production
    std::array<std::array<double, nP>, nE> y_muon_pn; // same for photonuclear

    std::vector<double> vy_muon_brems[nE]; // vectors containing inelasticities distributed such that they follow the dsdy curves above.
    std::vector<double> vy_muon_epair[nE];
    std::vector<double> vy_muon_pn[nE];

    std::array<std::array<double, nP>, nE> y_cumulative_muon_brems;
    std::array<std::array<double, nP>, nE> y_cumulative_muon_epair;
    std::array<std::array<double, nP>, nE> y_cumulative_muon_pn;
    
    std::array<double, nE> int_muon_brems;  // integral of prob. dist.=# of expected interactions.
    // for each energy
    std::array<double, nE> int_muon_epair; // same for pair prod.
    std::array<double, nE> int_muon_pn; // same for photnuclear

    double max_muon_brems; // maximum value of prob. dist., brems
    double max_muon_epair;  // max value of prob. dist., pair prod.
    double max_muon_pn; // same for photonucl.

    double min_muon_brems; // minimum value of prob. dist., brems
    double min_muon_epair;  // min value of prob. dist., pair prod.
    double min_muon_pn; // same for photonucl.

    std::array<std::array<double, nP>, nE> dsdy_tauon_brems;  // same as above, but for taus
    std::array<std::array<double, nP>, nE> dsdy_tauon_epair;
    std::array<std::array<double, nP>, nE> dsdy_tauon_pn;
    std::array<std::array<double, nP>, nE> dsdy_tauon_hadrdecay; // hadronic decay of taus
    std::array<std::array<double, nP>, nE> dsdy_tauon_edecay; // tau decay to electrons
    std::array<std::array<double, nP>, nE> dsdy_tauon_mudecay; // tau decay to muons.

    // brems, epair, pn, hadrdecay, edecay, mudecay
    // brems, epair, pn

    enum class Secondary {brems,
			  epair,
			  pn,
			  hadrdecay,
			  edecay,
			  mudecay
    };

    // enum class Property {y,
    // 			 dsdy
    // };
    

    // typedef std::map<std::map<std::pair<Neutrino::Flavor,  Secondary>,  std::array<std::array<double, nP>, nE> thingy;
    typedef std::pair<Neutrino::Flavor, Secondary> DataKey;
    std::map<DataKey, TH2F> fE_YCumulative_dsdy;
    std::map<DataKey, TH2F> fE_Y_dsdy;
    std::map<DataKey, TH1F> fE_dsdy;
    std::map<DataKey,  std::array<std::array<double, nP>, nE> > fY;
    std::map<DataKey,  std::array<std::array<double, nP>, nE> > fDsdy;

    std::array<std::array<double, nP>, nE> y_tauon_brems;
    std::array<std::array<double, nP>, nE> y_tauon_epair;
    std::array<std::array<double, nP>, nE> y_tauon_pn;
    std::array<std::array<double, nP>, nE> y_tauon_hadrdecay;
    std::array<std::array<double, nP>, nE> y_tauon_edecay;
    std::array<std::array<double, nP>, nE> y_tauon_mudecay;

    std::array<std::array<double, nP>, nE> y_cumulative_tauon_brems;
    std::array<std::array<double, nP>, nE> y_cumulative_tauon_epair;
    std::array<std::array<double, nP>, nE> y_cumulative_tauon_pn;
    std::array<std::array<double, nP>, nE> y_cumulative_tauon_hadrdecay;
    std::array<std::array<double, nP>, nE> y_cumulative_tauon_edecay;
    std::array<std::array<double, nP>, nE> y_cumulative_tauon_mudecay;

    // // vectors containing inelasticities distributed such that they follow the dsdy curves above.
    // std::vector<double> vy_tauon_brems[nE];
    // std::vector<double> vy_tauon_epair[nE];
    // std::vector<double> vy_tauon_pn[nE];
    // std::vector<double> vy_tauon_hadrdecay[nE];
    // std::vector<double> vy_tauon_edecay[nE];
    // std::vector<double> vy_tauon_mudecay[nE];

    std::array<double, nE> int_tauon_brems;
    std::array<double, nE> int_tauon_epair;
    std::array<double, nE> int_tauon_pn;
    std::array<double, nE> int_tauon_hadrdecay;
    std::array<double, nE> int_tauon_edecay;
    std::array<double, nE> int_tauon_mudecay;

    // double max_tauon_brems;
    // double max_tauon_epair;
    // double max_tauon_pn;
    // double max_tauon_hadrdecay;
    // double max_tauon_edecay;
    // double max_tauon_mudecay;

    // double min_tauon_brems; // minimum value of prob. dist., brems
    // double min_tauon_epair;  // min value of prob. dist., pair prod.
    // double min_tauon_pn; // same for photonucl.
    // double min_tauon_hadrdecay; // minimum value of prob. dist., brems
    // double min_tauon_edecay;  // min value of prob. dist., pair prod.
    // double min_tauon_mudecay; // same for photonucl.


    //   double weight_doublebang=0;
    //   double dtime_doublebang=0;


    static const int N_TAUOLA=10001;
    double tauola[N_TAUOLA][6];//tau decay energy fractions, [0]=tau neutrino, [1]=mu neutrino, [2]=Electron neutrino, [3]=hadrons, [4]=electromagnetic, [5]=muons
    std::ifstream tauolainfile;

    void readData(Neutrino::Flavor f, Secondary s);
    // void readData(const std::string&,const std::string&,
    // 		  std::array<std::array<double, nP>, nE>& y,
    // 		  std::array<std::array<double, nP>, nE>& dsdy);
    
		  // double (*)[nP], double (*)[nP]);
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


  private:
    void doShower(Neutrino::Flavor ,Energy, Energy&, Energy&, int&);
    

  public:
    ShowerModel(const Settings* settings);
    virtual ~ShowerModel();

    
    Shower generate(const Neutrino& n, const Interaction& i);

    int SECONDARIES;
    int TAUDECAY; // is tau decay counted as a secondary interaction
    double TAUFRAC; //fraction of tau neutrino-cc current events where the primare interaction point is the first bang
    int count_nfb;
    int secondary_e_noncons;


    void InitTauola();
    // void GetTauDecay(const std::string& nuflavor, const std::string& current, std::string& taudecay, double& emfrac_db, double& hadfrac_db);
    void GetTauDecay(Neutrino::Flavor nuflavor, Interaction::Current current, std::string& taudecay, double& emfrac_db, double& hadfrac_db);    

    void pickEMFracDB(double& emfrac_db, double& hadfrac_db);
    double GetDBViewAngle(const TVector3 &refr, const TVector3 &nnu);
    //void GetFirstBang(const Geoid::Position &r_in, const TVector3 &nnu, Geoid::Position &posnu, double len_int_kgm2, double d1, double &nuentrancelength);
    double NFBWeight(double ptau, double taulength);

    Shower GetEMFrac(Neutrino::Flavor nuflavor,
		     Interaction::Current current,
		     double y,
		     Energy pnu);
    bool secondbang;
    static const bool interestedintaus=false;

    ClassDef(ShowerModel, 0);
  }; //class ShowerModel
}

#endif ///SHOWER_MODEL_H
