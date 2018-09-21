#ifndef SPECTRA_H_
#define SPECTRA_H_

#include "RNG.h"
#include "TSpline.h"
#include <string>

namespace icemc{

  //! Neutrino spectra
  class Spectra : public RNG {

  private:
    double maxflux;   // max flux value
    //  static const int NSPECTRA_MAX=300;  // why need this??
    static const int E_bin_max = 50;
    int E_bin;   // initialize # of energy bins (max = 50)
  
    //  double energy[E_bin_max]; // energies that correspond to the fluxes in the previous array  
    //  double EdNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
    //  double E2dNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
  
    void GetFlux(std::string filename);    // read neutrino flux EdNdEdAdt (in GeV) from filename file
  
    TGraph *gEdNdEdAdt;   //graph for EdNdEdAdt flux
    TGraph *gE2dNdEdAdt;  //graph for E2dNdEdAdt flux

    TGraph *CDF;
    TGraph *inverse_CDF;

    TSpline3 *sEdNdEdAdt; //spline of EdNdEdAdt
    TSpline3 *sE2dNdEdAdt;    //spline of E2dNdEdAdt
    int EXPONENT; // set flux model

  public:  

    double energy[E_bin_max]; // energies that correspond to the fluxes in the previous array  
    double EdNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
    double E2dNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
  
    Spectra(int EXPONENT); // constructor  
  
    double GetNuEnergy(); // get the neutrino energy which follows neutrino flux. 
    double GetCDFEnergy();//get Energy from 'CDF'
    void GetCDF();//set up CDF and inverse CDF;
    TGraph *GetGEdNdEdAdt();
    TGraph *GetGE2dNdEdAdt();

    TSpline3 *GetSEdNdEdAdt();
    TSpline3 *GetSE2dNdEdAdt();

    double *Getenergy();
    double *GetEdNdEdAdt();
    double *GetE2dNdEdAdt();
    double GetEdNdEdAdt(double E_val);    // return flux value from TSpline
    double GetE2dNdEdAdt(double E_val);   // return flux value from TSpline

    double Getmaxflux();
  
    int GetE_bin();   // return energy bin number


    int IsSpectrum(); // return 1 or 0 depend on EXPONENT value
    int IsMonoenergetic();    // return 1 or 0 depend of EXPONENT value

    /** 
     * Save ROOT TGraphs and png images of the energy spectra used to generate neutrinos
     * 
     * @param fileNameNoSuffix is the base name WITHOUT any suffix, e.g. "temp" will produce temp.root, temp1.png, temp2.png
     */
    void savePlots(const TString& fileNameNoSuffix);

    /** 
     * Saves
     * 
     * @param fileName 
     */
    void savePlots2(const TString& fileName);

    // destructor

  }; //class Spectra
}


#endif
