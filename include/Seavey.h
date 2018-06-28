
#ifndef ICEMC_SEAVEY_H
#define ICEMC_SEAVEY_H

#include "vector.hh"
#include "Detector.h" // for PropagatingSignal
#include "Constants.h"
#include "AskaryanFreqsGenerator.h"// for N_AIR, although I'm going to move it

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

    Seavey(const Settings* settings = NULL, double refractiveIndexOfMedium = icemc::AskaryanFreqsGenerator::N_AIR);

    
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
     * Since ANITA moves and freely rotates, the antenna position needs to be updated.
     * 
     * @param position Put the antenna in this position.
     */
    void setPosition(const icemc::Vector& position){ fPosition = position;}

    /** 
     * What position is the antenna in?
     * @return const reference to #fPosition
     */
    const icemc::Vector& getPosition() const {return fPosition;}
    
    const icemc::Vector& getEPlane() const {return fEPlane;}
    const icemc::Vector& getHPlane() const {return fHPlane;}
    const icemc::Vector& getNormal() const {return fNormal;}
    

    
    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompEvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_pol,
				     double& e_component, double& h_component, double& n_component);
    

    ///@todo  make this private when the refactor is complete?
    static void GetEcompHcompkvector(const Vector& n_eplane, const Vector& n_hplane, const Vector& n_normal, const Vector n_exit2bn,
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

    ///@todo make these all private
    icemc::Vector fPosition; ///< Position in payload centered coordinates
    icemc::Vector fEPlane; ///< Seavey E-plane
    icemc::Vector fHPlane; ///< Seavey H-plane
    icemc::Vector fNormal; ///< Normal to the antenna
    
  private:

    FTPair fVPol;
    FTPair fHPol;

    double fRefractiveIndex = icemc::AskaryanFreqsGenerator::N_AIR; ///< This is the refractive index at the antenna (formerly known as nmedium_receiver)
    bool fDebug = false;
    std::vector<std::pair<double, double> > fPassBandsHz; ///< Passbands frequencies (Hz) pairs go(low, high), filled in constructor if icemc::Settings are passed
  };
}

#endif // ICEMC_SEAVEY_H
