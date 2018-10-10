#ifndef ASKARYAN_FACTORY_H
#define ASKARYAN_FACTORY_H

#include <cmath>
#include <iostream>
#include "AskaryanFreqs.h"
#include "ShowerGenerator.h"
#include "Report.h"
#include "FTPair.h"

namespace icemc {

  class Shower;
  
  /**
   * @class AskaryanFactory
   * @brief Generate Askaryan radiation from a neutrino of a given energy in a given medium.
   * 
   * In theory one should set the medium and give it an energy and generate a set of frequencies.
   */
  class AskaryanFactory {
  public:

    // AskaryanFactory(); ///< Default constructor
    AskaryanFactory(int n, double dt); ///< Default constructor

    void TaperVmMHz(double viewangle, double deltheta_em, double deltheta_had,
		    double emfrac, double hadfrac, double& vmmhz1m, double& vmmhz_em) const;

    void TaperVmMHz(double viewangle, double deltheta_em, double deltheta_had,
		    const Shower& sp, double& vmmhz1m, double& vmmhz_em) const {
      TaperVmMHz(viewangle, deltheta_em, deltheta_had, sp.emFrac, sp.hadFrac, vmmhz1m, vmmhz_em);
    }

    ///@todo make this more elegent once you understand it better, (move to AskaryanFreqs class and maybe put the loop over k inside the function)
    void TaperVmMHz(double viewangle, double deltheta_em, double deltheta_had, const Shower& sp, AskaryanFreqs& radioSignal, int k) const {
      double dummyVariable;
      TaperVmMHz(viewangle,  deltheta_em, deltheta_had,  sp, radioSignal.vmmhz[k], dummyVariable);
    }

    void GetSpread(Energy pnu, double emfrac, double hadfrac, double freq,
		   double& deltheta_em_max, double& deltheta_had_max) const;

    void GetSpread(Energy pnu, const Shower& sp, double freq,
		   double& deltheta_em_max, double& deltheta_had_max) const {
      GetSpread(pnu, sp.emFrac,  sp.hadFrac, freq, deltheta_em_max, deltheta_had_max);
    }
    
    
    double GetVmMHz1m(double pnu, double freq) const;
    AskaryanFreqs generateAskaryanFreqs(double vmmhz_max, double vmmhz1m_max, double pnu, int numFreqs, const double *freq_Hz, double notch_min, double notch_max, const Shower* sp) const;
    // FTPair generate(double pnu, int numFreqs, const double *freq_Hz, const Shower* sp) const;

    void GetVmMHz(double vmmhz_max, double vmmhz1m_max, double pnu, const double *freq,
		  double notch_min,double notch_max, double *vmmhz, int nfreq) const;
    
    int GetLPM() const;
    Energy GetELPM() const;

    void Initialize();

    void SetParameterization(int whichparameterization);


    void SetMedium(int medium) {
      MEDIUM = medium;
      if (MEDIUM!=0) {
	icemc::report() << severity::info << "Medium is " << MEDIUM << ", which is a non-default setting:  Not ice!\n";
      }
      InitializeMedium();
    }
    
    void InitializeMedium();

    void SetNMediumReceiver(double nmedium_receiver) {
      NMEDIUM_RECEIVER=nmedium_receiver;
    }
    void SetLPM(double lpm) {
      LPM=(int)lpm;
    }
    void SetKelvins(double kelvins) {
      KELVINS=kelvins;
    }
    void SetBetaMedium(double betamedium) {
      BETAMEDIUM=betamedium;
    }

    void SetRhoMedium(double rhomedium) {
      RHOMEDIUM=rhomedium;
    }
    void SetKrMedium(double kr_medium) {
      KR_MEDIUM=kr_medium;
    }
    void SetKlMedium(double kl_medium) {
      KL_MEDIUM=kl_medium;
    }

    void SetRmMedium(double rm_medium) {
      RM_MEDIUM=rm_medium;
    }
    void SetNDepth(double n_depth) {
      N_DEPTH=n_depth;
      SetChangle(acos(1/N_DEPTH));
      SetrhoDepth((N_DEPTH-1.)/0.86*1000.);
      SetX0Depth(X0MEDIUM); // no dependence on rho
    }
    
    void SetX0Depth(double x0_depth) {
      X0_DEPTH=x0_depth;
    }
    void SetrhoDepth(double rho_depth) {
      RHO_DEPTH=rho_depth;
    }
    void SetKeMedium(double ke_medium) {
      KE_MEDIUM=ke_medium;
    }
    void SetEcMedium(double ecmedium) {
      ECMEDIUM=ecmedium;
    }
    void SetX0Medium(double x0medium) {
      X0MEDIUM=x0medium;
    }

    double GetChangle() const {return changle;}
    
    void SetChangle(double thischangle) {
      changle=thischangle;
    }
    void SetAlphaMedium(double alphamedium) {
      ALPHAMEDIUM=alphamedium;
    }
    void SetAexMedium(double aexmedium) {
      AEXMEDIUM=aexmedium;
    }
    void SetKdelta_Medium(double kdelta_medium) {
      KDELTA_MEDIUM=kdelta_medium;
    }
    void SetJaime_Factor(double jaime_factor) {
      JAIME_FACTOR=jaime_factor;
      if (JAIME_FACTOR!=1){
	icemc::report() << severity::info << "Non-default setting: JAIME_FACTOR = " << JAIME_FACTOR << "\n";
      }
    }







    // double vmmhz1m_max; // V/m/MHz at 1m
    double X0MEDIUM;                  // radiation length of medium
    double KELVINS;                   // temperature of medium + system
    double N_DEPTH;                   // index of refraction at the interaction depth
    double RHO_DEPTH;                 // density at the interaction depth
    double X0_DEPTH;                  // density at the interaction depth
    double NMEDIUM_RECEIVER;          // index of refraction at receiver
    double RHOMEDIUM;                 // density of medium

    static const double RHOICE;       // density of ice (kg/m**3)
    static const double RHOAIR;       // density of air (kg/m**3)
    static const double RHOH20;       // density of water (kg/m**3) 
    static const double N_AIR;        // index of refr for air
    static const double NICE;         // index of refraction of ice
    static const double NSALT;        // index of refraction of salt
    static const double CHANGLE_ICE;  // cherenkov angle in ice
    static const double VIEWANGLE_CUT;



  protected:
    const int fNumFreqs;
    const double fDeltaF_Hz;    
    const std::vector<double> fFreqs_Hz;
    
    // properties of ice
    double x0ice; 
    double ecice;                     // critical energy in ice (MeV)
    double nice;                      // index of refraction of ice
    double nfirn;                     // index of refraction at the very surface - Peter
    double invnfirn; 
    double invnice;
    double rhoice;                    // density of ice (kg/m**3)
    double kelvins_ice;               // temperature in Kelvin (ice+system)
    double changle_ice;
    double aex_ice;                   // efficiency for producing charge asymmetry relative to ice.  1 by definition
                                      // double n_depth;  // index of refraction at the interaction depth
    double alphaice;                  // exponent that goes into cutting off the spectrum at high frequencies
    double rm_ice;                    // moliere radius, in g/cm^2
    double ke_ice;                    // const staticant in jaime's parameterization, in V/cm/MHz
    double kl_ice;                    // const staticant in jaime's parameterization
    double kdelta_ice;                // const staticant in jaime's parameterization
    double kr_ice;                    // const staticant in jaime's parameterization
    double betaice;                   // exponent, in jaime's parameterization
    double nu0_modified;              // nu_0 modified for a specific medium
    double nu_r;                      // parameter for signal parameterization
    int WHICHPARAMETERIZATION;
    double vmmhz1m_reference;         // reference value for V/m/MHz at f=1 MHz and pnu=10^18 eV
    double freq_reference;            // reference frequency in MHz
    double pnu_reference;             // reference energy in eV
    double KR_MEDIUM;                 // constant in jaime's parameterization
    double RM_MEDIUM;                 // moliere radius, in g/cm^2
    double KL_MEDIUM;                 // constant in jaime's parameterization  
    double KE_MEDIUM;                 // constant in jaime's parameterization, in V/cm/MHz
    double ECMEDIUM;                  // radiation length of medium
    double ALPHAMEDIUM;               // exponent that goes into cutting off the spectrum at high frequencies
    double AEXMEDIUM;                 // efficiency for making charge asymmetry
    double KDELTA_MEDIUM;             // constant in jaime's parameterization
    double BETAMEDIUM;                // exponent, in jaime's parameterization
    double JAIME_FACTOR;              // factor to multiply Jaime's parameterization for error analysis
    int MEDIUM;
    int LPM;

    double changle;                   // Cherenkov angle (radians) property of choice of medium.

    static const double RHOSALT;      // density of salt (kg/m**3)
    static const double RM_ICE;       // moliere radius, in g/cm^2
    static const double RM_SALT;      // moliere radius, in g/cm^2
    static const double KR_SALT;      // constant in jaime's parameterization
    static const double KR_ICE;       // constant in jaime's parameterization
    static const double X0SALT;       // radiation length of salt (meters)
    static const double ECSALT;       // critical energy in salt (MeV)
    static const double X0ICE; 
    static const double ECICE;        // critical energy in ice (MeV)
    static const double AEX_ICE;      // efficiency for producing charge asymmetry relative to ice.  1 by definition
    static const double ALPHAICE;     // exponent that goes into cutting off the spectrum at high frequencies
    static const double AEX_SALT;     // efficiency for producing charge asymmetry relative to ice
    static const double ALPHASALT;    // exponent that goes into cutting off the spectrum at high frequencies
    static const double KE_SALT;      // constant in jaime's parameterization, in V/cm/MHz
    static const double KL_SALT;      // constant in jaime's parameterization
    static const double KDELTA_SALT;  // constant in jaime's parameterization
    static const double KE_ICE;       // constant in jaime's parameterization, in V/cm/MHz
    static const double KL_ICE;       // constant in jaime's parameterization
    static const double KDELTA_ICE;   // constant in jaime's parameterization
    static const double KELVINS_ICE;  // temperature in Kelvin (ice+system)
    static const double KELVINS_SALT; // temperature in salt (350) + receiver temp (150)
    static const double BETAICE;      // exponent, in jaime's parameterization
    static const double BETASALT;     // exponent, in jaime's parameterization   
    
  };
}
#endif // ASKARYAN_FREQS_GENERATOR_H
