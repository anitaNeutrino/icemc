#ifndef INTERACTION_H
#define INTERACTION_H

#include "Geoid.h"
#include "RNG.h"
#include "Neutrino.h"

namespace icemc {

  class Settings;
  class Antarctica;
  
  /**
   * @class Interaction
   * @brief Stores everything about a particular neutrino interaction.
   */
  class Interaction : public RNG {

  public:    
    /** 
     * Constructor
     */
    Interaction(const Settings *settings1); //, int whichray); //, Counting *count1);

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

#endif //INTERACTION_H
