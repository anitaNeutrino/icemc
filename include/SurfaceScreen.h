#ifndef ICEMC_SURFACE_SCREEN_H
#define ICEMC_SURFACE_SCREEN_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdexcept>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <vector>

#include "TVector3.h"
#include "Geoid.h"


namespace icemc{
  class Detector;
  class Settings;

  
  class SurfaceScreen {
  private:
    double fedgeLength;               ///< the full length of one side
    Geoid::Position fcentralPoint;           ///< coordinates of screen center
    double fcosineProjectionFactor;   ///< cosine projection factor of the screen onto the ground, corrects for the long extension so sampling is faster; cos(angle between local normal at RF exit and vector to balloon)
  
    TVector3 fnormal;                   ///< screen orientation, '+' = pointing back to balloon
    TVector3 funit_x;                   ///< X unit vector in screen (parallel to ground surface, perp. to screen normal)
    TVector3 funit_y;                   ///< Y unit vector in screen (~ perp. to ground surface, perp. to screen normal)

    int fNsamples;                    ///< number of samples in X-direction (and Y-, assuming symmetry)
    int fNvalidpoints;                ///< total number of points on the screen

    std::vector<double> fVmmhz_freq;  ///< container for the valid screen points giving the frequency dependence magnitude for each point; every anita::NFREQ will be each screen point; final size will be (anita::NFREQ * fNsamples)
    std::vector<double> fVmmhz0;      ///< container for vmmhz[0]
    std::vector<double> fViewangle;
    std::vector<double> fDelays;      ///< container for the relative propagation phase delays for each frequency and screen point; final size will be (anita::NFREQ *fNsamples) after the push_backs
    std::vector<TVector3> fVec2blns;    ///< container of 'vector to balloon'
    std::vector<TVector3> fPols;        ///< container of transmitted polarizations
    std::vector<Geoid::Position> fImpactPt;  ///< container of ground impact points
    std::vector<double> fWeight;      ///< container for weight of a screen point ( == area of screen element), normalized when used
    std::vector<double> fIncAngles;   ///< container for incidence angles
    std::vector<double> fTransAngles; ///< container for transmission angle
    double fWeightNorm;               ///< normalization of the weights == simple weight sum
    std::vector<double> fFacetLength; ///< edge length [m] of individual contributing facet
    std::vector<double> fTcoeff_parl_polparl; ///< transmission coefficients of parallel components for parallel polarization
    std::vector<double> fTcoeff_perp_polparl; ///< transmission coefficients of perpendicular components for parallel polarization
    std::vector<double> fTcoeff_parl_polperp; ///< transmission coefficients of parallel components for perpendicular polarization
    std::vector<double> fTcoeff_perp_polperp; ///< transmission coefficients of perpendicular components for perpendicular polarization

  public:
    //! Creates an instance of a screen
    /**
     * @param a - unused (needs removed)
     */
    SurfaceScreen(int a);

    //! Sets number of grid divisions for the base screen
    /**
     * @param i - number of grid point
     */
    void SetNsamples(int i);

    //! Gets number of grid divisions for the base screen
    /**
     * @return integer
     */
    int GetNsamples() const;

    //! Sets the physical length of a side of the screen
    /**
     * @param a - size (meters)
     */
    void SetEdgeLength(double a);

    //! Sets the position of the central point of the screen
    /**
     * @param a - position vector
     */
    void SetCentralPoint(Geoid::Position a);

    //! Sets the projection factor of the screen relative to the specular RF exit point
    /**
     * @param a - cosine of the viewing angle w.r.t. to screen normal
     */
    void SetCosineProjectionFactor(double a);

    //! Gets the projection factor
    /**
     * @return double
     */
    double GetCosineProjectionFactor() const;

    //! Sets the screen normal
    /**
     * @param a - normal vector
     */
    void SetNormal(TVector3 a);

    //! Sets an orientation vector of the screen
    /**
     * @param a - vector
     */
    void SetUnitX(TVector3 a);

    //! Sets another orientation vector of the screen
    /**
     * @param a - vector
     */
    void SetUnitY(TVector3 a);

    //! Gets the screen length
    /**
     * @return double
     */
    double GetEdgeLength() const ;

    //! Gets the position of the screen's central point
    /**
     * @return Position
     */
    Geoid::Position GetCentralPoint() const;

    //! Gets the screen normal
    /**
     * @return Vector
     */
    TVector3 GetNormal() const;

    //! Gets an orientation vector
    /**
     * @return Vector
     */
    TVector3 GetUnitX() const;

    //! Gets another orientation vector
    /**
     * @return Vector
     */
    TVector3 GetUnitY() const;

    //! Calculates the X index of the screen corresponding to the specified counter value
    /**
     * @param i - index
     * @return double
     */
    double CalcXindex(int i) const ;

    //! Calculates the Y index of the screen corresponding to the specified counter value
    /**
     * @param i - index
     * @return double
     */
    double CalcYindex(int i) const ;

    //! Calculates the physical position of the screen corresponding to the specified counter value
    /**
     * @param i - index
     * @return double
     */
    Geoid::Position GetPosition(int i, int j) const;

    //! Appends a Vmmhz value to the fVmmhz_freq array
    /**
     * @param A - Vmmhz value
     */
    void AddVmmhz_freq(double A);

    //! Get the Vmmhz value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetVmmhz_freq(int i) const;

    //! Appends a Vmmhz value (for the lowest Anita frequency) to the fVmmhz0 array
    /**
     * @param A - Vmmhz value
     */
    void AddVmmhz0(double A);

    //! Gets the Vmmhz0 value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetVmmhz0(int i) const;

    //! Appends a viewangle value to the fViewangle array
    /**
     * @param A - viewangle
     */
    void AddViewangle(double A);

    //! Get the viewangle value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetViewangle(int i) const;

    //! Appends a delay value to the fDelays array
    /**
     * @param A - delay
     */
    void AddDelay(double A);

    //! Get the delay value stores at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetDelay(int i) const;

    //! Sets the total number of points on the screen
    /**
     * @param i - number of points
     */
    void SetNvalidPoints(int i);

    //! Gets the total number of points
    /**
     * @return int
     */
    double GetNvalidPoints() const;

    //! Appends a vector to the fVec2blns array
    /**
     * @param v - Vector
     */
    void AddVec2bln(TVector3 v);

    //! Gets the to-balloon vector at the specified index
    /**
     * @param i - index
     * @return Vector
     */
    TVector3 GetVec2bln(int i) const;

    //! Appends a vector to the fPols array
    /**
     * @param v - Vector
     */
    void AddPol(TVector3 v);

    //! Gets the polarization vector at the specified index
    /**
     * @param i - index
     * @return Vector
     */
    TVector3 GetPol(int i) const;

    //! Appends a vector to the fImpactPt array
    /**
     * @param p - Position
     */
    void AddImpactPt(Geoid::Position p);

    //! Gets the position at the specified index
    /**
     * @param i - index
     * @return Position
     */
    Geoid::Position GetImpactPt(int i) const;

    //! Appends a weight value to the fWeight array
    /**
     * @param a - double
     */
    void AddWeight(double a);

    //! Gets the weight value at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetWeight(int i) const;

    //! Sets the normalization factor for the weights (so they sum to 1)
    /**
     * @param a - double
     */
    void SetWeightNorm(double a);

    //! Gets the weight normalization factor
    /**
     * @return double
     */
    double GetWeightNorm() const;

    //! Appends an incidence angle value to the fIncAngles array
    /**
     * @param A - viewangle
     */
    void AddIncidenceAngle(double A);

    //! Get the incidence angle value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetIncidenceAngle(int i) const;

    //! Appends a transmission angle value to the fTransAngles array
    /**
     * @param A - viewangle
     */
    void AddTransmissionAngle(double A);

    //! Get the transmission angle value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetTransmissionAngle(int i) const;

    //! Appends a facet edge length value to the fFacetLength array
    /**
     * @param A - viewangle
     */
    void AddFacetLength(double A);

    //! Get the facet edge length value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetFacetLength(int i) const;

    //! Appends a parallel transmission coefficient value to the fTcoeff_parl array
    /**
     * @param A - coefficient
     */
    void AddTparallel_polParallel(double A);

    //! Get the parallel transmission coefficient value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetTparallel_polParallel(int i) const;

    //! Appends a perpendicular transmission coefficient value to the fTcoeff_perp array
    /**
     * @param A - coefficient
     */
    void AddTperpendicular_polParallel(double A);

    //! Get the perpendicular transmission coefficient value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetTperpendicular_polParallel(int i) const;

    //! Appends a parallel transmission coefficient value to the fTcoeff_parl array
    /**
     * @param A - coefficient
     */
    void AddTparallel_polPerpendicular(double A);

    //! Get the parallel transmission coefficient value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetTparallel_polPerpendicular(int i) const;

    //! Appends a perpendicular transmission coefficient value to the fTcoeff_perp array
    /**
     * @param A - coefficient
     */
    void AddTperpendicular_polPerpendicular(double A);

    //! Get the perpendicular transmission coefficient value stored at the specified index
    /**
     * @param i - index
     * @return double
     */
    double GetTperpendicular_polPerpendicular(int i) const;

    //! Resets the following screen parameters (fNvalidpoints,fVmmhz_freq,fVmmhz0,fViewangle,fDelays,fVec2blns,fPols,fImpactPt,fWeight,fWeightNorm)
    void ResetParameters();



  };
}

#endif
