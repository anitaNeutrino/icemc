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
#include "TMath.h" 
#include "TF1.h" 
#include <stdint.h> 
#include "TRandom.h" 


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

  //a way to further restrict source types. 
  struct Restriction
  {
    const char * whichSubtype; 
    const char * whichSources; 
    double whichStartTime; 
    double whichEndTime; 
    double maxDec; 
    static double fromString(const char * time); 
    Restriction(double max_dec = 30, const char * sources = "All", const char * subtype = "All",
        double startTime = 0, double endTime= 2147483647) 
      : whichSubtype(subtype), whichSources(sources), whichStartTime(startTime), whichEndTime(endTime), maxDec(max_dec) {; }
  };

  static SourceModel * getSourceModel(const char * key, unsigned seed = 0, Restriction r = Restriction()); 

  /** this must be called before asking for a time weight */ 
  void setUpWeights(double t0, double t1, double minE = 1e9, double maxE=1e12, int N = 1e6); 

  /** Add a source to our model. 
   * This class will then own the source (it will release its memory). 
   **/ 
  void addSource(Source * source) { sources.push_back(source) ; }
  //this is the flux at this time divided by the average flux (computed by setUpWeights). 
  double getTimeWeight(double t) const; 
    
  const char * getName() const { return name; } 
  /** Returns the index of the source used ! */ 
  int getDirectionAndEnergy( Vector * nudir, double t, double & nuE, double minE = 1e9, double maxE = 1e12) ; 
  /** Returns the index of the source used ! */ 
  int getDirection(Vector &nudir, double t, double nuE = 1e10) { return getDirectionAndEnergy( &nudir, t, nuE, nuE, nuE); }
  TH1 * estimateFlux (double tmin, double tmax, double Emin, double Emax, int nbins = 100, int Ntrials = 1e6); 
  const Source * getSource(int i ) const { return  sources[i]; } 
  unsigned getNSources() const { return sources.size(); } 
  virtual ~SourceModel(); 

  /** fills a vector with the times that sources turn on and off */ 
  void computeFluxTimeChanges(std::vector<double> * changes) const; 
private:
  std::vector<Source*> sources; 
  double weight_Emin;
  double weight_Emax;
  double average_flux; 
  const char * name;
}; 


/** Virtual source flux class */
class SourceFlux
{
public: 
  virtual double getFlux(double E, double t) const = 0; 
  virtual double getFluxBetween(double Emin, double Emax, double t) const = 0; 
  virtual double pickEnergy(double Emin, double Emax, double t, TRandom * rng = gRandom) const = 0; 
  virtual void getFluxTimeChanges(std::vector<double> * changes) const { (void) changes ; }
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
  const char * getName() const { return name.c_str(); }
  double getRA() const { return RA; }
  double getDec() const { return dec; } 
  Vector getDirection( double t) const; 
  const SourceFlux * getFlux() const { return flux; } 
  virtual ~Source() { delete flux; } 

protected:
  std::string name; 
  SourceFlux * flux; 
  double RA, dec; 
};

/** A time invariant flux with an exponential distribution */ 
class ConstantExponentialSourceFlux : public SourceFlux
{

public: 
  //gamma is the spectral index (so it's positive). norm is the normalization (in units of GeV / cm^2 / s) at normE, where normE is in GeV
  ConstantExponentialSourceFlux(double gamma, double norm, double normE=1e5); 
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





/** Same as above, but only between time t0 and t1 */
class TimeWindowedExponentialSourceFlux : public SourceFlux
{

public: 

  TimeWindowedExponentialSourceFlux(double t0, double t1, double gamma, double norm, double normE = 0.1, double cutoff = 0 ); 

  virtual double getFlux(double E, double t) const 
  { 
    if (t < t0 || t > t1){ return 0; }
    if (cutoff && E > cutoff) return 0; 
    else return A * TMath::Power(E,-gamma) ;
  } 
  virtual double getFluxBetween(double Emin, double Emax, double t) const 
  { 
    if (t < t0 || t > t1){ return 0; }
    if (cutoff && Emin > cutoff) return 0; 
    if (cutoff && Emax > cutoff) Emax = cutoff; 
    return A * (  TMath::Power(Emin,-gamma+1) / (gamma-1)  - TMath::Power(Emax,-gamma+1) / (gamma-1));
  }
  virtual double pickEnergy(double Emin, double Emax, double t, TRandom * rng = gRandom) const; 
  virtual ~TimeWindowedExponentialSourceFlux() { ; } 
 
  virtual void getFluxTimeChanges(std::vector<double> * changes) const { changes->push_back(t0); changes->push_back(t1); }
private: 

  double gamma; 
  double A; 
  mutable TF1 f; 
  double t0, t1; 
  double cutoff; 
}; 

#endif 
