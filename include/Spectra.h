#ifndef SPECTRA_H_
#define SPECTRA_H_

#include "RNG.h"
#include "AbstractSources.h"
#include "TSpline.h"
#include <string>
#include "Neutrino.h"

namespace icemc{

  class Settings;

  namespace Source {
  
    //! Neutrino spectra
    class Spectra : public EnergyModel, public RNG {

    private:
      const Settings* fSettings;
      double maxflux;   // max flux value
      //  static const int NSPECTRA_MAX=300;  // why need this??
      // static const int E_bin_max = 50;
      // int E_bin;   // initialize # of energy bins (max = 50)
  
      //  double energy[E_bin_max]; // energies that correspond to the fluxes in the previous array  
      //  double EdNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
      //  double E2dNdEdAdt[E_bin_max]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt
  
      void GetFlux(std::string filename);    // read neutrino flux EdNdEdAdt (in GeV) from filename file
  
      TGraph *gEdNdEdAdt = nullptr;   //graph for EdNdEdAdt flux
      TGraph *gE2dNdEdAdt = nullptr;  //graph for E2dNdEdAdt flux

      TGraph *CDF = nullptr;
      TGraph *inverse_CDF = nullptr;

      TSpline3 *sEdNdEdAdt = nullptr; //spline of EdNdEdAdt
      TSpline3 *sE2dNdEdAdt = nullptr;    //spline of E2dNdEdAdt
      int EXPONENT; // set flux model


      std::vector<double> energy;//[E_bin_max]; // energies that correspond to the fluxes in the previous array  
      std::vector<double> EdNdEdAdt;//[E_bin_max]; //flux of incident neutrinos vs. energy E*dN/dE/dA/dt
      std::vector<double> E2dNdEdAdt;//[E_bin_max]; //flux of incident neutrinos vs. energy E^2*dN/dE/dA/dt

      inline void resizeFluxVectors(int n){
	energy.resize(n, 0);
	EdNdEdAdt.resize(n ,0);
	E2dNdEdAdt.resize(n, 0);
      }
      
    public:
      Spectra(const Settings* settings); // constructor
    
      virtual Energy pickNeutrinoEnergy();
    
      Energy GetNuEnergy(); // get the neutrino energy which follows neutrino flux. 
      Energy GetCDFEnergy();//get Energy from 'CDF'
      void GetCDF();//set up CDF and inverse CDF;
      TGraph* GetGEdNdEdAdt();
      TGraph* GetGE2dNdEdAdt();

      TSpline3* GetSEdNdEdAdt();
      TSpline3* GetSE2dNdEdAdt();

      // double* Getenergy();
      // double* GetEdNdEdAdt();
      // double* GetE2dNdEdAdt();
      double GetEdNdEdAdt(double E_val);    // return flux value from TSpline
      double GetE2dNdEdAdt(double E_val);   // return flux value from TSpline

      double Getmaxflux();
  
      // int GetE_bin();   // return energy bin number


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

  } // namespace Source
} // namespace icemc


#endif
