#ifndef _icemc_source_hh 
#define _icemc_source_hh 

/** For point-source analyses, we want to draw an energy and direction
 * from the sources we're considering. 
 *
 *  Flux units are now in GeV/cm^2/s (no steradian anymore) 
 *
 *   ALL ENERGY UNITS ARE IN GeV  NOT log eV. That was easiest for me to figure out :)  
 *
 *   
 *
*/

#include <vector> 
#include "vector.hh"  // only two kinds of vectors? we need more! 
#include "TRandom3.h" 
#include "TMath.h" 
#include "TF1.h" 



/** A source model holds a collection of sources */ 

class Source; 

class SourceModel 
{

  public: 

    SourceModel(const char * model_name, unsigned seed  = 0); 

/** We will label our source models by a string so they're easy to pick. 
 *  
 *  We will document the available models here
 *
 *  Since it may make sense to stack some source models, you can select multiple models with a + 
 *
 *
 */ 

    static SourceModel * getSourceModel(const char * key, unsigned seed = 0); 

    /** Add a source to our model. 
     * This class will then own the source (it will release its memory). 
     **/ 
    void addSource(Source * source) { sources.push_back(source) ; }
    
    const char * getName() const { return name; } 
    int getDirectionAndEnergy( Vector * nudir, double t, double & nuE, double minE = 1e9, double maxE = 1e13) ; 
    int getDirection( Vector &nudir, double t, double nuE = 1e10) { return getDirectionAndEnergy( &nudir, t, nuE, nuE, nuE); }
    TH1 * estimateFlux (double tmin, double tmax, double Emin, double Emax, int nbins = 100, int Ntrials = 1e6); 
    unsigned getNSources() const { return sources.size(); } 
    virtual ~SourceModel(); 
  private:
    std::vector<Source*> sources; 
    const char * name; 
    TRandom3 rng; 
}; 


/** Virtual source flux class */
class SourceFlux
{
  public: 
    virtual double getFlux(double E, double t) const = 0; 
    virtual double getFluxBetween(double Emin, double Emax, double t) const = 0; 
    virtual double pickEnergy(double Emin, double Emax, double t, TRandom * rng = gRandom) const = 0; 
    virtual ~SourceFlux() { ; }
}; 




/*
 * Source class... for now assume point sources. Can add extended sources later. 
 * */ 
class Source
{

  public: 
    /* The source will own the flux */ 
    Source (const char * name, double RA, double dec, SourceFlux * flux); 
    Vector getDirection( double t) const; 
    const SourceFlux * getFlux() const { return flux; } 
    virtual ~Source() { delete flux; } 

  protected:
    const char * name; 
    SourceFlux * flux; 
    double RA, dec; 
};

/** A time invariant flux with an exponential distribution */ 
class ConstantExponentialSourceFlux : public SourceFlux
{

  public: 
    //gamma is the spectral index (so it's positive). norm is the normalization (in units of GeV / cm^2 / s) at normE, where normE is in GeV
    ConstantExponentialSourceFlux(double gamma, double norm, double normE=0.1); 
    virtual double getFlux(double E, double t) const 
    { (void) t; return A * TMath::Power(E,-gamma) ; } 
    virtual double getFluxBetween(double Emin, double Emax, double t) const 
    { (void) t; return A * (  TMath::Power(Emin,-gamma+1) / (gamma-1)  - TMath::Power(Emax,-gamma+1) / (gamma-1)); }
    virtual double pickEnergy(double Emin, double Emax, double t, TRandom * rng = gRandom) const; 
    virtual ~ConstantExponentialSourceFlux() { ; } 
 
  private: 
    double gamma; 
    double A; 
    mutable TF1 f; 
}; 




/** uses fermi light curves from known TeV blazars 
 *  but converting these to a flux is... non-trivial. 
 *
 **/
class FermiFlux : public SourceFlux
{

  public: 

//    FermiFlux(const TGraph * flux, const TGraph * bg,  const TGraph * sigma, double minsig, double gamma ); 


  private: 


}; 

#endif 
