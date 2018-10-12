#ifndef ICEMC_EVENT_GENERATOR_H
#define ICEMC_EVENT_GENERATOR_H

#include "anita.hh"
#include "TVector3.h"

#include "Settings.h"

#include "Constants.h"
#include "CommandLineOptions.h"
#include "ShowerGenerator.h"
#include "RNG.h"

namespace icemc {

  class Taumodel;
  class AskaryanFactory;
  class Interaction;
  class RayTracer;
  class Roughness;
  class Screen;
  class ConnollyEtAl2011;
  class ANITA;

  namespace Source {
    class Spectra;
  }

  class Detector;

  /**
   * @class EventGenerator
   * @brief The glue that brings the simulation together
   * 
   * Contains the main neutrino generation loop in generateNeutrinos()
   */
  class EventGenerator : public RNG {
  private:
    const Settings* fSettings;
    double fStartTime;
    double fEndTime;
  public:

    // enum RayDirection {
    //   direct = 0,
    //   downgoing = 1,
    // };

    EventGenerator(const Settings* settings);
    virtual ~EventGenerator();
    static void interrupt_signal_handler(int);  // This catches ctrl-C interrupt (SIGINT)


    
    // @todo constify... if I can const this, then we're probably near the end of the refactor...
    void applyRoughness(const Settings& settings1, const int& inu, Interaction* interaction1,  RayTracer* ray1, Screen* panel1, Antarctica* antarctica1, Balloon* bn1, const AskaryanFactory* askFreqGen, Anita* anita1, const Shower& showerProps);


    
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
    int GetDirection(const Settings *settings1,  Interaction *interaction1,  const TVector3 &refr,  double deltheta_em,  double deltheta_had,  const Shower& sp,  double vmmhz1m_max,  double r_fromballoon,  RayTracer *ray1,  const AskaryanFactory* askFreqGen,  Geoid::Position posnu,  Anita *anita1,  Balloon *bn1,  TVector3 &nnu,  double& costhetanu,  double& theta_threshold) ;
    void GetFresnel(Roughness *rough1,  int ROUGHNESS_SETTING,  const TVector3 &nsurf_rfexit,  const TVector3 &n_exit2rx,  TVector3 &n_pol,  const TVector3 &nrf2_iceside,  double efield, double deltheta_em, double deltheta_had,  double &t_coeff_pokey,  double &t_coeff_slappy,  double &fresnel,  double &mag) const;
    // void GetFresnel(Roughness *rough1,  int ROUGHNESS_SETTING,  const TVector3 &nsurf_rfexit,  const TVector3 &n_exit2rx,  TVector3 &n_pol,  const TVector3 &nrf2_iceside,  double efield,  const Shower& ,  double deltheta_em, double deltheta_had,  double &t_coeff_pokey,  double &t_coeff_slappy,  double &fresnel,  double &mag) const;    
    double GetViewAngle(const TVector3 &nrf2_iceside,  const TVector3 &nnu) const;
    int TIR(const TVector3 &n_surf,  const TVector3 &nrf2_iceside,  double N_IN,  double N_OUT) const;
    // void IntegrateBands(Anita *anita1,  int k,  Screen *panel1,  double *freq,  double scalefactor,  double *sumsignal) const;

    // // @todo constify... needs some love to constify
    // void Summarize(const Settings *settings1,  Anita* anita1,  Source::Spectra *spectra1, const AskaryanFactory* askFreqGen,  ConnollyEtAl2011 *primary1,  double,  double eventsfound,  double,  double,  double,  double*,  double,  double,  double&,  double&,  double&,  double&, TString);
    // void WriteNeutrinoInfo(const int& inu, const Geoid::Position&,  const TVector3&,  const Geoid::Position&,  double,  Neutrino::Flavor,  Neutrino::Interaction::Current,  double,  std::ofstream &nu_out) const;

    /** 
     * @brief Run the neutrino generation
     * 
     * This function does the stuff that used to be the main in the icemc executable
     * @todo needs some love to constify... probably not doable...
     */    
    void generateNeutrinos(Detector& detector);


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
    Detector* fDetector = nullptr;
    Taumodel* fTauPtr = nullptr;
  private:

    

  };
}



#endif //ICEMC_EVENT_GENERATOR_H
