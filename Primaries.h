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
#include "vector.hh"
#include "position.hh"

class Vector;
class Position;
class Interaction;
class Primaries;
class IceModel;
class Counting;
class Settings;

//using namespace std;
using std::string;
//
using std::cout;




//! Inelasticity distributions: stores parametrizations and picks inelasticities
class Y {

  

private:
    TF1* ffrac; //!< This is the fraction of the distribution in the low y region given by Equation 18. 
    

    TF1* fC1_high[2][2]; //!< parameterization of parameter C1 in the high y region according to Equation 16
    
    TF1 *fC1_low; //!< parameterization of parameter C1 in the low y region according to Equation 16.
    
    TF1 *fC2; //!< parameterization of parameter C2 in the low y region according to Equation 17.

    
    TF3* fy0_low;   //!< For picking inelasticity in low y region according to Equation 14.
    TF2* fy0_high;     //!< For picking inelasticity in high y region according to Equation 15.
    TRandom3 Rand3;

    double pickYConnollyetal2011(int NU,int CURRENT,double e); //!< pick an inelasticity using recipe in Connolly et al. (2011)
    // NU=0: nubar, NU=1: nu
    // CURRENT=0: CC, CURRENT-1: NC

    static constexpr double ymin_low = 0.00002;//!< Minimum y in low-y region, Connolly et al.
    constexpr static double ymax_low=0.001;//!< Maximum y in low-y region, Connolly et al.
    constexpr static double ymin_high=0.001;//!< Minimum y in high-y region, Connolly et al.
    constexpr static double ymax_high=1.;//!< Maximum y in high-y region, Connolly et al.
    constexpr static double dy_low=0.00002;//!< y step in low region, Connolly et al.
    constexpr static double dy_high=0.001;//!< y step in high region, Connolly et al.

    double pickYGandhietal(); //!< THIS IS A ROUGH PARAMETRIZATION OF PLOT 6 FROM Ghandhi,Reno,Quigg,Sarcevic  hep-ph/9512364 (the curves are not in their later article).  There is also a slow energy dependence.

    
    constexpr static double R1=0.36787944;  //!< 1/e, used for Gandhi et al.
    constexpr static double R2=0.63212056;  //!< 1-R1, used for Gandhi et al.

    
public:
    Y();
    double pickY(Settings *settings1,double pnu,int nu_nubar,int currentint);//!< pick inelasticity y according to chosen model
    double Getyweight(double pnu,double y,int nu_nubar,int currentint);
    //!< If you want to choose y from a flat distribution this is the weight it should have according to Connolly et al. (2011)
    
};//Y
//! Functions you need to generate a primary interaction including cross sections and picking charged current/neutral current and flavor
class Primaries {
    
private:
    TH2D *m_hsigma; //!< plot of cross section vs. log(E/GeV)
    TCanvas *m_csigma;//!< canvas
    Y *m_myY; 
    int run_old_code;
    
public:
    double pickY(Settings *settings1,double pnu,int nu_nubar,int currentint);//!<pick inelasticity y according to chosen model
    double Getyweight(double pnu,double y,int nu_nubar,int currentint);//!< in case you choose y from a flat distribution, this is the weight you should give it according to Connolly et al. (2011)


    double A_low[4];//!< Table V of Connolly et al. for use in Eq. 16.  Same for any nu_nubar and current type.
    double A0_high[2][2];//!< Table V of Connolly et al. for use in Eq. 16.  
    double A1_high[2][2];//!< Table V of Connolly et al. for use in Eq. 16.  
    double A2_high[2][2];//!< Table V of Connolly et al. for use in Eq. 16.  
    double A3_high[2][2];//!< Table V of Connolly et al. for use in Eq. 16.  
    double b0; //!<  Eq. 17 of Connolly et al.
    double b1; //!<  Eq. 17 of Connolly et al.
    
    TF1* m_fy[2][2];
    TF1* m_fsigma[2][2];
    
    double c0[2][2];//!< Table V of Connolly et al. for Eq. 7
    double c1[2][2];//!< Table V of Connolly et al. for Eq. 7
    double c2[2][2];//!< Table V of Connolly et al. for Eq. 7
    double c3[2][2];//!< Table V of Connolly et al. for Eq. 7
    double c4[2][2];//!< Table V of Connolly et al. for Eq. 7
    
    static constexpr int NSIGMAS=2;//!< number of possible cross section models
    //!< 0=Gandhi et al.
    //!< 1=Connolly et al. 2011
    double mine[NSIGMAS];//!< minimum energy for cross section parametrizations, in eV
    double maxe[NSIGMAS]; //!<minimum energy for cross section parametrizations, in eV
    
    Primaries();//!<constructor 
    ~Primaries();//!<destructor 
    //!<*primary1 must be manually deleted in icemc for deconstructor to actually be called.
    
//! Neutrino-nucleon cross-sections using model chosen
    int GetSigma(double pnu,double& sigma,double &len_int_kgm2,Settings *settings1,int nu_nubar,int currentint);



    // string GetCurrent();
    string GetNuFlavor();
protected:
};//!<Primaries

//! Stores everything about a particular neutrino interaction.  
class Interaction  {
    
private:
    
    
    
    Vector tmp_banana; //!<Intermediate vector
    
    //!<For banana plot
    
    //!< static const double RADDEG_TMP=3.14159/180.;
    static constexpr double nu_banana_theta_angle=-0.413 * 3.14159/180.;//!< don't let me use RADDEG which is annoying 
    
    
    static constexpr double altitude_nu_banana=-400.;//!<Depth of interaction of banana neutrino
    
    
    static constexpr double lat_nu_banana=0.; 
    static constexpr double lon_nu_banana=0.;
    
    
    static constexpr double banana_slopey=0.;//!<Turn slopyness off for banana plots (SLOPEY)
    static constexpr double nu_banana_phi_angle=0. * 3.14159/180.; 
    
    
    
public:
    
    static constexpr double phi_nu_banana=3.14159/4; //!<Location in phi
    
    static constexpr double banana_observation_distance=600000.;//!<How far from the surface above the interaction are we when we measure the voltages? (meters) Note: Should be at least 100000 for best results.
    static constexpr double theta_nu_banana=170.*3.14159/180.;//!<Location of banana neutrino in theta

    static constexpr int kcc=0;
    static constexpr int knc=1;

    double banana_phi_obs;
    Vector banana_obs; //!<Vector from the neutrino interaction to the observation point
    Interaction(string inttype,Primaries *primary1,Settings *settings1,int whichray,Counting *count1);//! Constructor

    void PickAnyDirection();
    
    int noway;
    int wheredoesitleave_err;
    int neverseesice;
    int wheredoesitenterice_err;
    int toohigh;
    int toolow;
    
    double pathlength_inice;
    
    Vector nnu;  //!< direction of neutrino (+z in south pole direction)
    double costheta_nutraject; //!<theta of nnu with earth center to balloon as z axis 
    double phi_nutraject; //!<phi of nnu with earth center to balloon as z axis

    double weight_nu; //!< Weight for neutrino that survives to posnu
    double weight_nu_prob; //!< Weight for neutrino that survives to posnu and interacts in the ice
    
    Position r_in; //!< position where neutrino enters the earth
    Position r_enterice; //!< position where neutrino enters the ice
    Position nuexit; //!< place where neutrino would have left the earth
    Position nuexitice; //!< place where neutrino would have left the ice
    double chord;  //!< chord in m from earth entrance to rock-ice boundary
    double logchord; //!< log_10 of chord length earth entrance to where it enters ice
    double weight_bestcase; //!< what weight1 would be if whole earth had density of crust - for quick and dirty calculation of best case scenario
    double chord_kgm2_bestcase; //!< the chord the neutrino would traverse if it all was crust density
    double chord_kgm2_ice; //!< from ice entrance to interaction point
    double d1;  //!<same as chord in m (earth entrance to rock-ice boundary)
    double d2;  //!< ice-rock boundary to interaction point in m
    double nearthlayers; //! number of earth layers traversed
    double total_kgm2; // the total kgm2 traversed
    int crust_entered; 
    int mantle_entered; 
    int core_entered;
    
    
    static constexpr double pnu_banana=2.00E19;
    static constexpr double banana_y=0.2;//!<Elasticity.  0.2 is an average number.
    double banana_weight;//!<Weight measurement locations to account for phase space
    double banana_theta_obs;
    double banana_volts;//!<Total voltage measured at a spot on the sky
    static constexpr double banana_signal_fluct=0.;//!<Turn off noise for banana plots (settings1->SIGNAL_FLUCT) (shouldn't matter)
    static constexpr double banana_sigma=0.;//!<NSIGMA in the case of a banana plot
    
    
    void  setNuFlavor(Primaries *primary1,Settings *settings1,int whichray,Counting *count1);
    string GetCurrent();
    void setCurrent();
    Position posnu;
    Position posnu_down;
    string  nuflavor;                   //!< neutrino flavor
    string  current;                    //!<  CC or NC?
    int nuflavorint;                //!< Added by Stephen for output purposes
    int currentint;                 //!< Ditto - Stephen
    
    
    double surface_over_banana_nu;
    
    string banana_flavor; //!<Force interaction to be a muon neutrino
    string banana_current;  //!<Force interaction to be a neutral current
    
    Vector nnu_banana; //!<Forced neutrino direction 
    
    Position nu_banana;  //!<The forced interaction point of the neutrino for the banana plots
    Position nu_banana_surface; //!<The location of the surface above the forced neutrino interaction point  
    
    //!< phase space weighting
    double dtryingdirection; //!<weighting factor: how many equivalent tries each neutrino counts for after having reduced angular phase space for possibly detectable events
    double dnutries; //!<product of dtryingdirection and dtryingposition

    double d_effective_area; //!< In unbaised mode, projection of surface onto nnu
    
    double altitude_int;//!< depth of interaction
    double altitude_int_mirror;//!<depth of the mirror point of interaction.
    
    double r_fromballoon[2]; //!< distance from interaction to balloon for each ray
    
    double r_fromballoon_db; //!< same, for double bangs
    double r_exit2bn; //!< exit to balloon
    double r_exit2bn_measured; //!< exit to balloon deduced from measured theta
    int iceinteraction;//!< whether or not there is an interaction in the ice
    
    
protected:
};//!<Interaction
#endif
