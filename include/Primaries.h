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

#include "Geoid.h"
#include "TVector3.h"
#include "Neutrino.h"



namespace icemc {

  class Interaction;
  class Primaries;
  class IceModel;
  class Counting;
  class Settings;
  class Y;

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
    double pickY(const Settings *settings1,double pnu,int nu_nubar, Neutrino::CurrentType currentint);///<pick inelasticity y according to chosen model
    double Getyweight(double pnu,double y,int nu_nubar,Neutrino::CurrentType currentint);///< in case you choose y from a flat distribution, this is the weight you should give it according to Connolly et al. (2011)


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
    int GetSigma(double pnu,double& sigma,double &len_int_kgm2,const Settings *settings1,int nu_nubar,Neutrino::CurrentType currentint);



    // string GetCurrent();
    Neutrino::Flavor GetNuFlavor() const;    
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

    int PickDownwardInteractionPoint(const Geoid::Position& r_bn, const Settings *settings1, const Antarctica *antarctica1);
    
    void PickAnyDirection();
    
    int noway;
    int wheredoesitleave_err;
    int neverseesice;
    int wheredoesitenterice_err;
    int toohigh;
    int toolow;
    
    double pathlength_inice;
    
    TVector3 nnu;                 ///< direction of neutrino (+z in south pole direction)
    double costheta_nutraject;  ///< theta of nnu with earth center to balloon as z axis 
    double phi_nutraject;       ///< phi of nnu with earth center to balloon as z axis

    double weight_nu;           ///< Weight for neutrino that survives to posnu
    double weight_nu_prob;      ///< Weight for neutrino that survives to posnu and interacts in the ice
    
    Geoid::Position r_in;              ///< position where neutrino enters the earth
    Geoid::Position r_enterice;        ///< position where neutrino enters the ice
    Geoid::Position nuexit;            ///< place where neutrino would have left the earth
    Geoid::Position nuexitice;         ///< place where neutrino would have left the ice
    double chord;               ///< chord in m from earth entrance to rock-ice boundary
    double logchord;            ///< log_10 of chord length earth entrance to where it enters ice
    double weight_bestcase;     ///< what weight1 would be if whole earth had density of crust - for quick and dirty calculation of best case scenario
    double chord_kgm2_bestcase; ///< the chord the neutrino would traverse if it all was crust density
    double chord_kgm2_ice;      ///< from ice entrance to interaction point
    double d1;                  ///< same as chord in m (earth entrance to rock-ice boundary)
    double d2;                  ///< ice-rock boundary to interaction point in m
    
    
    void setNuFlavor(const Primaries *primary1, const Settings *settings1);//, int whichray, Counting *count1);
    Neutrino::CurrentType GetCurrent();
    void setCurrent();
    int getPdgCode() const;

    Geoid::Position posnu;
    Geoid::Position posnu_down;
    Neutrino::Flavor nuflavor;	  ///< neutrino flavor    
    Neutrino::CurrentType current;	  ///< CC or NC?
    

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



#endif
