#ifndef ICEMC_SEAVEY_H
#define ICEMC_SEAVEY_H

#include "TVector3.h"
#include "Detector.h" // for PropagatingSignal
#include "Constants.h"
#include "AskaryanFactory.h"// for N_AIR, although I'm going to move it

class TCanvas;

namespace icemc {

  class Settings;
  class ANITA;

  /**
   * @class Seavey
   * @brief A lightweight class to store the ANITA antenna information
   */

  class Seavey {
  public:

    /**
     * The Pol enum class is the opposite (!) convention from eventReaderRoot and all downstream ANITA tools!
     * You have been warned!
     */
    enum class Pol : int      {V,       H};
    enum class XPol : int      {VtoH,    HtoV};
    enum class AngleDir : int  {Azimuth, Elevation};

    Seavey(const TVector3& positionV, const TVector3& positionH,
	   const TVector3& ePlane, const TVector3& hPlane, const TVector3& normal,
	   const Settings* settings = NULL, double refractiveIndexOfMedium = icemc::AskaryanFactory::N_AIR);

    
    void addSignal(const PropagatingSignal& incomingSignal);
    
    double getHeight(Pol pol, double freqHz) const;
    double getHeight(XPol xPol, double freqHz) const;
    double getOffAxisResponse(Pol pol,  AngleDir dir, double freqHz, double angleRad) const;


    /** 
     * Is this frequency allowed by the passbands?
     * @see #fPassBandsHz
     * 
     * @param freqHz is the frequency in question
     * 
     * @return true if allowed by any (or no passbands are specified), false otherwise
     */
    bool freqAllowedByPassBands(double freqHz) const;

    const FTPair& get(Pol pol) const;
    

    /** 
     * Where is the Seavey (after accounting for payload movement and position)
     * 
     * @param pol is the desired polarization (defaults to VPol)
     * 
     * @return the position of the phase center
     */    
    const TVector3& getPosition(Pol pol = Pol::V) const {
      switch(pol){
      case Pol::H: return fPositionH.global;
      case Pol::V: // is also the default...
      default:
	return fPositionV.global;
      }
    }
    
    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompEvector(const TVector3& n_eplane, const TVector3& n_hplane, const TVector3& n_pol,
				     double& e_component, double& h_component, double& n_component);
    

    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompkvector(const TVector3& n_eplane, const TVector3& n_hplane, const TVector3& n_normal, const TVector3 n_exit2bn,
				     double& e_component_kvector, double& h_component_kvector, double& n_component_kvector);

    
    ///@todo make this private when the refactor is complete?
    static void GetHitAngles(double e_component_kvector,double h_component_kvector,double n_component_kvector, 
			     double& hitangle_e, double& hitangle_h);

    

    /** 
     * Draws the gains data. The canvas will delete the graphs when you delete it (which is up to you)
     * 
     * @return The canvas on which the gains data are plotted (your responsibility to delete)
     */
    static TCanvas* plotGains();

    /** 
     * Draw the antenna gains for a refractive index, n
     * The canvas will delete the graphs when you delete it (which is up to you)
     * 
     * @param n is the refractiveIndex
     * 
     * @return The canvas on which the height data are plotted (your responsibility to delete)
     */    
    static TCanvas* plotHeights(double n);
    
    /** 
     * Draws the gains data. The canvas will delete the graphs when you delete it (which is up to you)
     * 
     * @return The canvas on which the gains data has been plotted
     */
    static TCanvas* plotInterpolatedGains(double freq, int nBins = 50);


    /** 
     * Turn on some verbose debugging output
     * 
     * @param b true to switch on, false to switch off
     */
    void setDebug(bool b){
      fDebug = b;
    }

    /**
     * Is the debugging turned on?
     * 
     * @return the current state of #fDebug
     */
    bool getDebug() const {return fDebug;}


    /** 
     * @todo finish this
     * Move everything representing the Seavey into a new position
     * 
     * @param anitaPos 
     * @param heading 
     * @param pitch 
     * @param roll 
     */
    void updatePosition(const Geoid::Position& anitaPos, double heading, double pitch, double roll);

  private:
    
    /**
     * @struct VectorPair, a pair of vectors
     * 
     * This is just to group the two representations of the Seavey positions together
     * (payload coordinates and Earth centered (global) coordinates)
     */
    struct VectorPair {
      const TVector3 payload;
      TVector3 global;
      VectorPair(const TVector3& v) : payload(v){
	global = payload;
      }
    };

    VectorPair fPositionV; ///< VPol phase center position 
    VectorPair fPositionH; ///< HPol phase center position
    VectorPair fEPlane; ///< Seavey E-plane 
    VectorPair fHPlane; ///< Seavey H-plane
    VectorPair fNormal; ///< Normal to the antenna
    
    FTPair fVPol;
    FTPair fHPol;

    double fRefractiveIndex = icemc::AskaryanFactory::N_AIR; ///< This is the refractive index at the antenna (formerly known as nmedium_receiver)
    bool fDebug = false;
    std::vector<std::pair<double, double> > fPassBandsHz; ///< Passbands frequencies (Hz) pairs go(low, high), filled in constructor if icemc::Settings are passed
  };
}

std::ostream& operator<<(std::ostream& os, icemc::Seavey::Pol pol);

#endif // ICEMC_SEAVEY_H
