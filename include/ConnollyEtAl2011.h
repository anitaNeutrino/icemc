#ifndef CONNOLLY_ET_AL_2011_H
#define CONNOLLY_ET_AL_2011_H

#include <iostream>
#include "CrossSectionModel.h"
#include "RNG.h"

#include "TF1.h"
#include "TF2.h"
#include "TF3.h"

#include "Geoid.h"
#include "TVector3.h"
#include "Neutrino.h"
#include "Inelasticity.h"

// std::vector;
namespace icemc {  
  
  class Interaction;
  class Settings;

  /**
   * @class ConnollyEtAl2011
   * @brief To compute Cross sections given in https://arxiv.org/pdf/1102.0691.pdf
   */
  class ConnollyEtAl2011 : public CrossSectionModel, public RNG {
    
    static constexpr int NSIGMAS=2;///< number of possible cross section models
    ///< 0=Gandhi et al.
    ///< 1=Connolly et al. 2011
    // double mine[NSIGMAS]; ///< minimum energy for cross section parametrizations, in eV
    // double maxe[NSIGMAS]; ///< maximum energy for cross section parametrizations, in eV
    
  public:
    // double pickY(const Settings *settings1,double pnu,int nu_nubar, Neutrino::Current currentint);///<pick inelasticity y according to chosen model
    double pickY(double pnu, Neutrino::L leptonNumber, Neutrino::Current currentint);///<pick inelasticity y according to chosen model    
    double Getyweight(double pnu,double y,Neutrino::L leptonNumber,Neutrino::Current currentint);///< in case you choose y from a flat distribution, this is the weight you should give it according to Connolly et al. (2011)

    ConnollyEtAl2011(const Settings* settings); ///< Constructor 
    
    /// Neutrino-nucleon cross-sections using model chosen
    virtual double getSigma(double pnu, Neutrino::L leptonNumber, Neutrino::Current currentint) const override;

    Neutrino::Flavor pickFlavor();    
  protected:
    const Settings* fSettings;

    
  private:
    
    Y m_myY; ///< defined in Inelasticity.h
    
    typedef std::map<std::pair<Neutrino::L, Neutrino::Current>, double> DoubleLCC; // represent a double for nu/nubar, neutral/charged currents
    std::array<double, 4> A_low; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A0_high; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A1_high; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A2_high; ///< Table V of Connolly et al. for use in Eq. 16.
    DoubleLCC A3_high; ///< Table V of Connolly et al. for use in Eq. 16.
    double b0;         ///<  Eq. 17 of Connolly et al.
    double b1;         ///<  Eq. 17 of Connolly et al.
    
    // TF1* m_fy[2][2];
    

    std::map<std::pair<Neutrino::L, Neutrino::Current>, TF1> fSigma;
    
    DoubleLCC c0;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c1;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c2;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c3;      ///< Table V of Connolly et al. for Eq. 7
    DoubleLCC c4;      ///< Table V of Connolly et al. for Eq. 7
    
  };///<ConnollyEtAl2011
  



  class Antarctica;
  
  /**
   * @class Interaction
   * @brief Stores everything about a particular neutrino interaction.
   */
  class Interaction : public RNG {
    
  private:
        

  public:
    
    /** 
     * Constructor
     */
    Interaction(ConnollyEtAl2011 *primary1, const Settings *settings1); //, int whichray); //, Counting *count1);

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
    
    // void setNuFlavor(const ConnollyEtAl2011 *primary1, const Settings *settings1);//, int whichray, Counting *count1);
    
    Neutrino::Current pickCurrent();
    int getPdgCode() const;

    Geoid::Position posnu;
    Geoid::Position posnu_down;
    Neutrino::Flavor nuflavor;	  ///< neutrino flavor
    Neutrino::Current current;	  ///< CC or NC?

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
