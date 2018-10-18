#ifndef ASKARYAN_RADIATION_MODEL_H
#define ASKARYAN_RADIATION_MODEL_H

#include <cmath>
#include <iostream>
#include "ShowerModel.h"
#include "Report.h"
#include "FTPair.h"
#include "Detector.h"

namespace icemc {

  class Shower;
  class Settings;
  
  /**
   * @class AskaryanRadiationModel
   * @brief Generate Askaryan radiation from a neutrino of a given energy in a given medium.
   * 
   * In theory one should set the medium and give it an energy and generate a set of frequencies.
   */
  class AskaryanRadiationModel {
  public:

    // AskaryanRadiationModel(); ///< Default constructor
    AskaryanRadiationModel(const Settings* settings, int n, double dt); ///< Default constructor

    
    
    /** 
     * Generates Askaryan signal 1m along shower axis from interaction, finds component pointing along outgoing RF direction
     * 
     * @param nu is the neutrino which triggered the interaction
     * @param shower contains details of the simulated shower
     * @param opticalPath is the contains direction the RF signal travels way from the interaction
     * 
     * @return the component of the Askaryan signal pointing along outgoingRfDirection
     */
    PropagatingSignal generate(const Neutrino& nu, const Shower& shower, const OpticalPath& opticalPath) const;

  private:
    const Settings* fSettings;

    double GetVmMHz1m(Energy pnu, double freq) const;
    FTPair generateOnAxisAt1m(Energy energy) const;
    void taperWaveform(FTPair& waveform/*modified*/, double viewAngleRadians, Energy energy, const Shower& shower) const;
    void GetSpread(Energy pnu, const Shower& sp, double freq, double& coneWidthEm, double& coneWidthHad) const;
    TVector3 getPolarizationVector(const TVector3& rfDir, const TVector3& showerAxis) const;
    
  public:
    
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


    template<class T>
    T taperSingleFreq(T amplitude,
		      double viewangle,
		      double coneWidthEm,
		      double coneWidthHad,
		      double emfrac,
		      double hadfrac) const {

      auto taper_component = [&](double threshold, double coneWidth, double frac){
			       double rtemp = (viewangle-changle)*(viewangle-changle)/(coneWidth*coneWidth);
			       T amp = 0;
			       // the power goes like exp(-(theta_v-theta_c)^2/Delta^2)
			       // so the e-field is the same with a 1/2 in the exponential
			       if (frac>threshold) { // if there is an em component
				 if (rtemp<=20) {
				   // if the viewing angle is less than 20 sigma away from the cerankov angle
				   // this is the effect of the em width on the signal
				   amp = amplitude*exp(-rtemp);
				 }
				 else{
				   // if it's more than 20 sigma just set it to zero 
				   amp=0.;
				 }
			       }
			       else{ // if the em component is essentially zero than set this to zero
				 amp = 0;
			       }
			       return amp;
			     };

      // std::cout << emfrac << "\t" << hadfrac << std::endl;
      T amplitude_em  = taper_component(1e-10, coneWidthEm,  emfrac);
      T amplitude_had = taper_component(1e-10, coneWidthHad, hadfrac);

      // std::cout << __PRETTY_FUNCTION__ << "\t" << amplitude_em << "\t" << amplitude_had << "\t" << std::endl;
      amplitude = sin(viewangle)*(emfrac*amplitude_em + hadfrac*amplitude_had);

      return amplitude;
      
      // //--EM 
      // T amplitude_em=0;  // V/m/MHz at 1m due to EM component of shower

      // // this is the number that get exponentiated
      // //  double rtemp=0.5*(viewangle-changle)*(viewangle-changle)/(deltheta_em*deltheta_em);
      // double rtemp=(viewangle-changle)*(viewangle-changle)/(deltheta_em*deltheta_em);

      // // the power goes like exp(-(theta_v-theta_c)^2/Delta^2)
      // // so the e-field is the same with a 1/2 in the exponential
      // if (emfrac>pow(10.,-10.)) { // if there is an em component
      // 	if (rtemp<=20) {
      // 	  // if the viewing angle is less than 20 sigma away from the cerankov angle
      // 	  // this is the effect of the em width on the signal
      // 	  amplitude_em=amplitude*exp(-rtemp);
      // 	}
      // 	else{
      // 	  // if it's more than 20 sigma just set it to zero 
      // 	  amplitude_em=0.;
      // 	}
      // }
      // else{ // if the em component is essentially zero than set this to zero
      // 	amplitude_em=0;
      // }


      // T amplitude_had=0; // V/m/MHz at 1m due to HAD component of shower  
      // //--HAD
      // // this is the quantity that gets exponentiated
      // rtemp=(viewangle-changle)*(viewangle-changle)/(deltheta_had*deltheta_had);

      // //std::cout << "rtemp (had) is " << rtemp << "\n";

      // if (hadfrac!=0) { // if there is a hadronic fraction
      // 	if (rtemp<20) { // if we are less than 20 sigma from cerenkov angle
      // 	  amplitude_had=amplitude*exp(-rtemp); // this is the effect of the hadronic width of the signal
      // 	}
      // 	else{ // if we're more than 20 sigma from cerenkov angle
      // 	  amplitude_had=0.; // just set it to zero
      // 	}
      // }
      // else {
      // 	amplitude_had=0.;
      // }
      // amplitude=sin(viewangle)*(emfrac*amplitude_em+hadfrac*amplitude_had);

      // return amplitude;

    } //TaperVmMHz
    
    
  };
}
#endif // ASKARYAN_RADIATION_MODEL_H
