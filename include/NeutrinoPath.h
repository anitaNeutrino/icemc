#ifndef NEUTRINO_PATH_H
#define NEUTRINO_PATH_H

#include "Geoid.h"

namespace icemc {


  class WorldModel;
  
  /**
   * @class NeutrinoPath
   * @brief Holds all path related information, for use in icemc::EventGenerator
   * 
   * Also holds an individual event's weight information, since once you have 
   * an event position and direction you have its weight.
   */
  class NeutrinoPath {
  private:
    const WorldModel* fW = nullptr;
    const Geoid::Position fInteractionPos;
    TVector3 fNeutrinoDir;    
  public:

    /** 
     * Default constructor, zeros all member variables
     */
    NeutrinoPath();

    /** 
     * From these two
     * 
     * @param interaction 
     * @param rfDir 
     */
    NeutrinoPath(const Geoid::Position& interaction, const TVector3& rfDir, const WorldModel* m);
    

    /** 
     * Default destructor, currently does nothing
     */
    virtual ~NeutrinoPath() {;}


    void project();


    /** 
     * @brief Sets all member variables to zero
     * 
     * Must be updated by hand as new variables are added.
     */
    void reset();

    double theta_in;	///< theta where neutrino enters earth (radians, south pole=0)
    double lat_in;	///< latitude where neutrino enters earth (degrees, south pole=-90)
    double nearthlayers;///< how many layers (core, mantle, crust) does nnu traverse
    double weight_prob;	///< event weight,  including probability it interacts somewhere in ice along its path
    double weight1;	///< event weight,  just based on absorption in earth,  see note
    double weight;	///< total event weight (product of the previous 2)
    double logweight;	///< log of the previous number
    double len_int;	///< interaction length in m
    double pieceofkm2sr;///< Use this for making plots comparing different cross sections.  The integral of a plot from a run will be the total Area*sr of the detector.  That way it is proportional to the cross section and the integral is something meaningful to people.

    // ClassDef(NeutrinoPath, 1);
  };
}

#endif //NEUTRINO_PATH_H
