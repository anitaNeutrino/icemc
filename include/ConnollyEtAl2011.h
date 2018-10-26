#ifndef CONNOLLY_ET_AL_2011_H
#define CONNOLLY_ET_AL_2011_H

#include <iostream>
#include "CrossSectionModel.h"

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

#include "Geoid.h"
#include "TVector3.h"
#include "Neutrino.h"
#include "Inelasticity.h"

namespace icemc {
  class Interaction;
  class Settings;  

  /**
   * @namespace ConnollyEtAl2011
   * @brief Contains things defined in https://arxiv.org/pdf/1102.0691.pdf
   */
  namespace ConnollyEtAl2011 {

    const double A0_low = 0.0;    ///< Table V of Connolly et al. for use in Eq. 16.
    const double A1_low = 0.0941; ///< Table V of Connolly et al. for use in Eq. 16.
    const double A2_low = 4.72;   ///< Table V of Connolly et al. for use in Eq. 16.
    const double A3_low = 0.456;  ///< Table V of Connolly et al. for use in Eq. 16.

    typedef std::pair<Neutrino::L, Neutrino::Interaction::Current> LC; // pair representing nu/nubar & neutral/charged current interactions
    typedef std::map<LC, double> Map_LC_Double; // stores a double for each combo of nu/nubar + neutral/charged current interactions

    const Map_LC_Double A0_high = {{{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged}, -0.0026},
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged}, -0.008 },
				   {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral}, -0.005 },   
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral}, -0.005 }};

    const Map_LC_Double A1_high = {{{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged},  0.085},
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged},  0.26 },
				   {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral},  0.23 },
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral},  0.23 }};

    const Map_LC_Double A2_high = {{{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged},  4.1 },
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged},  3.0 },
				   {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral},  3.0 },
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral},  3.0 }};

    const Map_LC_Double A3_high = {{{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged},  1.7 },
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged},  1.7 },
				   {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral},  1.7 },
				   {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral},  1.7 }};

    
    const double B0 = 2.55;            ///<  Eq. 17 of Connolly et al.
    const double B1 = -0.0949;         ///<  Eq. 17 of Connolly et al.

    ///< Table III in Connolly et al.
    const Map_LC_Double C0 = {{{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral}, -1.826 },
			      {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged}, -1.826 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral}, -1.033 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged}, -1.033 }};

    const Map_LC_Double C1 = {{{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral}, -17.31 },
			      {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged}, -17.31 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral}, -15.95 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged}, -15.95 }};

    const Map_LC_Double C2 = {{{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral}, -6.448 },
			      {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged}, -6.406 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral}, -7.296 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged}, -7.247 }};

    const Map_LC_Double C3 = {{{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral}, 1.431 },
			      {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged}, 1.431 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral}, 1.569 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged}, 1.569 }};

    const Map_LC_Double C4 = {{{Neutrino::L::Matter,     Neutrino::Interaction::Current::Neutral}, -18.61 },
			      {{Neutrino::L::Matter,     Neutrino::Interaction::Current::Charged}, -17.91 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Neutral}, -18.30 },
			      {{Neutrino::L::AntiMatter, Neutrino::Interaction::Current::Charged}, -17.72 }};


    // Parameters for equation18, given in the text just below
    const double F0 = 0.128;
    const double F1 = -0.197;
    const double F2 = 21.8;

    const double ymin_low  = 0.00002; ///< Minimum y in low-y region, Connolly et al.
    const double ymax_low  = 0.001;   ///< Maximum y in low-y region, Connolly et al.
    const double ymin_high = 0.001;   ///< Minimum y in high-y region, Connolly et al.
    const double ymax_high = 1.;      ///< Maximum y in high-y region, Connolly et al.
    // const double dy_low    = 0.00002; ///< y step in low region, Connolly et al.
    // const double dy_high   = 0.001;   ///< y step in high region, Connolly et al.    
    
    const std::string equation14 = "x+pow( z*pow([1]-x, (-1./y + 1)) + (1-z)*pow([0]-x,-1./y + 1) , y/(y-1))"; // x = C1, y = C2, z = R, [0] = y_min, [1] = y_max
    const std::string equation15 = "(pow([1]-x,y) / pow([0]-x, y-1)) + x"; // x = C1', y = R,  [0] = y_min, [1] = y_max    
    const std::string equation16 = "[0]-[1]*(-exp(-(x-[2])/[3]))"; // x = epsilon, [i] is A0, A1, A2, A3
    const std::string equation17 = "[0]+[1]*x"; // x = epsilon, [i] are the B0, B1
    const std::string equation18 = "[0]*sin([1]*(x-[2]))"; //x = epsilon, [i] are the F0, F1, F2




    
    /**
     * @class CrossSectionModel
     * @brief Implements the cross section model given in the paper
     */
    class CrossSectionModel : public icemc::CrossSectionModel {
    public:
      CrossSectionModel(const Settings* settings); ///< Constructor 
      virtual ~CrossSectionModel(){;}
      virtual double getSigma(Energy pnu, Neutrino::L leptonNumber, Neutrino::Interaction::Current currentint) const override;      

    private:
      std::map<std::pair<Neutrino::L, Neutrino::Interaction::Current>, TF1> fSigma;
    };



    class YGenerator :  public icemc::YGenerator, public RNG {

    private:
      const Settings* fSettings;
      const double epsLow  = 7;
      const double epsHigh = 12;
      TF1 fFracLowHigh;          ///< This is the fraction of the distribution in the low y region given by Equation 18.
      std::map<LC, TF1> fC1_high;       ///< parameterization of parameter C1 in the high y region according to Equation 16
      TF1 fC1_low;        ///< parameterization of parameter C1 in the low y region according to Equation 16.
      TF1 fC2;            ///< parameterization of parameter C2 in the low y region according to Equation 17.
      TF3 fy0_low;        ///< For picking inelasticity in low y region according to Equation 14.
      TF2 fy0_high;       ///< For picking inelasticity in high y region according to Equation 15.

      enum class YRegion {High = 0,Low = 1};
      
      
      

    public:
      YGenerator(const Settings* settings);
      virtual double pickY(Energy E,Neutrino::L l, Neutrino::Interaction::Current c) override;
      double Getyweight(Energy E, double y, Neutrino::L l, Neutrino::Interaction::Current c); ///@todo ?
      TH1D* plot(Energy energy = 1e18*Energy::Unit::eV,
		 Neutrino::L l = Neutrino::L::Matter,
		 Neutrino::Interaction::Current c = Neutrino::Interaction::Current::Neutral,
		 int n = 10000);
    };
    
  } //ConnollyEtAl2011
    
} /// icemc 



#endif
