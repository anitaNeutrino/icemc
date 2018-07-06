////////////////////////////////////////////////////////////////////////////////////////////////
//class Primaries:
////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PRIMARIES_H_
#define PRIMARIES_H_

#include "TRandom3.h" 
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TF3.h"
#include "TH2D.h"

#include "position.hh"
#include "vector.hh"

namespace icemc {

  class Interaction;
  class Primaries;
  class IceModel;
  class Counting;
  class Settings;

  /**
   * @class NuFlavor
   * @brief enum for neutrino flavour
   */
  enum class NuFlavor : int {
			   e = 1,
			   mu = 2,
			   tau = 3
  };


  /**
   * @class CurrentType
   * @brief enum for type of interaction
   */
  enum class CurrentType : int {
				Charged = 0,
				Neutral = 1
  };

  /**
   * @class Y
   * @brief Inelasticity distributions: stores parametrizations and picks inelasticities
   */
  class Y {  

  private:
    TF1* ffrac;							///< This is the fraction of the distribution in the low y region given by Equation 18. 
    TF1* fC1_high[2][2];					///< parameterization of parameter C1 in the high y region according to Equation 16
    TF1 *fC1_low;						///< parameterization of parameter C1 in the low y region according to Equation 16.
    TF1 *fC2;							///< parameterization of parameter C2 in the low y region according to Equation 17.
    TF3* fy0_low;						///< For picking inelasticity in low y region according to Equation 14.
    TF2* fy0_high;						///< For picking inelasticity in high y region according to Equation 15.
    TRandom3 Rand3;

    double pickYConnollyetal2011(int NU,CurrentType CURRENT,double e);	///< pick an inelasticity using recipe in Connolly et al. (2011)
    // NU=0: nubar, NU=1: nu
    // CURRENT=0: CC, CURRENT-1: NC

    static constexpr double ymin_low  = 0.00002;		///< Minimum y in low-y region, Connolly et al.
    static constexpr double ymax_low  = 0.001;			///< Maximum y in low-y region, Connolly et al.
    static constexpr double ymin_high = 0.001;			///< Minimum y in high-y region, Connolly et al.
    static constexpr double ymax_high = 1.;			///< Maximum y in high-y region, Connolly et al.
    static constexpr double dy_low    = 0.00002;		///< y step in low region, Connolly et al.
    static constexpr double dy_high   = 0.001;			///< y step in high region, Connolly et al.

    /** 
     * @brief THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM Ghandhi,Reno,Quigg,Sarcevic hep-ph/9512364
     * 
     * The curves are not in their later article.
     * There is also a slow energy dependence.
     * 
     * @return a parameterized y value
     */
    double pickYGandhietal();

    static constexpr double R1 = 0.36787944;			///< 1/e, used for Gandhi et al.
    static constexpr double R2 = 0.63212056;			///< 1-R1, used for Gandhi et al.

    
  public:
    Y();

    /** 
     * Pick inelasticity y according to chosen model
     * 
     * @param settings1 are the icemc settings 
     * @param pnu 
     * @param nu_nubar 
     * @param currentint 
     * 
     * @return 
     */
    double pickY(const Settings *settings1, double pnu, int nu_nubar, CurrentType currentint);

    /** 
     * @brief If you want to choose y from a flat distribution this is the weight it should have according to Connolly et al. (2011)
     * 
     * @param pnu 
     * @param y 
     * @param nu_nubar 
     * @param currentint 
     * 
     * @return 
     */
    double Getyweight(double pnu, double y, int nu_nubar, CurrentType currentint);
    
  };//Y


  /**
   * @class Primaries
   * @brief Functions you need to generate a primary interaction including cross sections and picking charged current/neutral current and flavor
   */
  class Primaries {
    
  private:
    TRandom3 Rand3;
    
    TH2D* m_hsigma;	///< plot of cross section vs. log(E/GeV)
    TCanvas* m_csigma;	///< canvas
    Y* m_myY; 
    int run_old_code;
    
  public:
    double pickY(const Settings *settings1,double pnu,int nu_nubar, CurrentType currentint);///<pick inelasticity y according to chosen model
    double Getyweight(double pnu,double y,int nu_nubar,CurrentType currentint);///< in case you choose y from a flat distribution, this is the weight you should give it according to Connolly et al. (2011)


    double A_low[4];      ///< Table V of Connolly et al. for use in Eq. 16.  Same for any nu_nubar and current type.
    double A0_high[2][2]; ///< Table V of Connolly et al. for use in Eq. 16.  
    double A1_high[2][2]; ///< Table V of Connolly et al. for use in Eq. 16.  
    double A2_high[2][2]; ///< Table V of Connolly et al. for use in Eq. 16.  
    double A3_high[2][2]; ///< Table V of Connolly et al. for use in Eq. 16.  
    double b0;            ///<  Eq. 17 of Connolly et al.
    double b1;            ///<  Eq. 17 of Connolly et al.
    
    TF1* m_fy[2][2];
    TF1* m_fsigma[2][2];
    
    double c0[2][2];      ///< Table V of Connolly et al. for Eq. 7
    double c1[2][2];      ///< Table V of Connolly et al. for Eq. 7
    double c2[2][2];      ///< Table V of Connolly et al. for Eq. 7
    double c3[2][2];      ///< Table V of Connolly et al. for Eq. 7
    double c4[2][2];      ///< Table V of Connolly et al. for Eq. 7
    
    static constexpr int NSIGMAS=2;///< number of possible cross section models
    ///< 0=Gandhi et al.
    ///< 1=Connolly et al. 2011
    double mine[NSIGMAS]; ///< minimum energy for cross section parametrizations, in eV
    double maxe[NSIGMAS]; ///< maximum energy for cross section parametrizations, in eV
    
    Primaries(); ///< Constructor 
    ~Primaries();///< Destructor 
    
    /// Neutrino-nucleon cross-sections using model chosen
    int GetSigma(double pnu,double& sigma,double &len_int_kgm2,const Settings *settings1,int nu_nubar,CurrentType currentint);



    // string GetCurrent();
    NuFlavor GetNuFlavor() const;    
  protected:
  };///<Primaries



  class RayTracer;
  class Antarctica;
  
  /**
   * @class Interaction
   * @brief Stores everything about a particular neutrino interaction.
   */
  class Interaction  {
    
  private:
        

  public:
    
    /** 
     * Constructor
     */
    Interaction(Primaries *primary1, const Settings *settings1); //, int whichray); //, Counting *count1);

    int PickDownwardInteractionPoint(int ibnposition, const Position& r_bn, const Settings *settings1, const Antarctica *antarctica1, RayTracer *ray1);
    
    void PickAnyDirection();
    
    int noway;
    int wheredoesitleave_err;
    int neverseesice;
    int wheredoesitenterice_err;
    int toohigh;
    int toolow;
    
    double pathlength_inice;
    
    Vector nnu;                 ///< direction of neutrino (+z in south pole direction)
    double costheta_nutraject;  ///< theta of nnu with earth center to balloon as z axis 
    double phi_nutraject;       ///< phi of nnu with earth center to balloon as z axis

    double weight_nu;           ///< Weight for neutrino that survives to posnu
    double weight_nu_prob;      ///< Weight for neutrino that survives to posnu and interacts in the ice
    
    Position r_in;              ///< position where neutrino enters the earth
    Position r_enterice;        ///< position where neutrino enters the ice
    Position nuexit;            ///< place where neutrino would have left the earth
    Position nuexitice;         ///< place where neutrino would have left the ice
    double chord;               ///< chord in m from earth entrance to rock-ice boundary
    double logchord;            ///< log_10 of chord length earth entrance to where it enters ice
    double weight_bestcase;     ///< what weight1 would be if whole earth had density of crust - for quick and dirty calculation of best case scenario
    double chord_kgm2_bestcase; ///< the chord the neutrino would traverse if it all was crust density
    double chord_kgm2_ice;      ///< from ice entrance to interaction point
    double d1;                  ///< same as chord in m (earth entrance to rock-ice boundary)
    double d2;                  ///< ice-rock boundary to interaction point in m
    
    
    void setNuFlavor(const Primaries *primary1, const Settings *settings1);//, int whichray, Counting *count1);
    CurrentType GetCurrent();
    void setCurrent();
    int getPdgCode() const;

    Position posnu;
    Position posnu_down;
    NuFlavor nuflavor;	  ///< neutrino flavor    
    CurrentType current;	  ///< CC or NC?
    

    double dtryingdirection;	  ///< weighting factor: how many equivalent tries each neutrino counts for after having reduced angular phase space for possibly detectable events
    double dnutries;		  ///< product of dtryingdirection and dtryingposition
    
    double altitude_int;	  ///< depth of interaction
    double altitude_int_mirror;	  ///< depth of the mirror point of interaction.
    
    double r_fromballoon[2];	  ///< distance from interaction to balloon for each ray
    
    double r_fromballoon_db;	  ///< same, for double bangs
    double r_exit2bn;		  ///< exit to balloon
    double r_exit2bn_measured;	  ///< exit to balloon deduced from measured theta
    int iceinteraction;		  ///< whether or not there is an interaction in the ice
    
    
  };

}


/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param f is the Flavor class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::NuFlavor& f);

/** 
 * For a nice cout/cerr/logging experience
 * 
 * @param os is a output string stream
 * @param c is the Current class enum
 * 
 * @return the updated output string stream
 */
std::ostream& operator<<(std::ostream& os, const icemc::CurrentType& c);


#endif
