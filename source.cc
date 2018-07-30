#include "source.hh" 
#include "TMath.h" 
#include <map>
#include "TObjString.h" 
#include "TTimeStamp.h" 



// SourceModel 


SourceModel::SourceModel(const char * n, unsigned seed) 
  : name(n), rng(seed) { }

SourceModel::~SourceModel() 
{
  for (unsigned i = 0; i < sources.size(); i++) delete sources[i]; 
}

static std::map<std::string, SourceModel * > models; 

SourceModel * SourceModel::getSourceModel(const char * key, unsigned seed) 
{

  if (!key || !strcasecmp(key,"NONE")) return NULL; 


  if (models.count(key)) 
    return models[key]; 

  SourceModel * m = new SourceModel(key,seed); 

  TString str(key); 
  TObjArray * tokens = str.Tokenize("+"); 

  for (int i =0; i < tokens->GetEntries(); i++)
  {
    TObjString * tok = (TObjString*) tokens->At(i); 

    TString stripped( tok->String().Strip( TString::kBoth)); 

    // the flux from TXS 
    if (stripped == "TXS15")
    {

      m->addSource(new Source("TXS 0506+056", 5 + 9/60. +  25.9637 / 3600, 
                                              5 + 41/60. + 35.3279/3600,
          new ConstantExponentialSourceFlux(2,8e4) //from IceCube paper , assume it's constant over flight for now
          )) ; 
    }

    else if (stripped == "GRB15") 
    {
      //add ANITA3 GRB's


    }
    else if (stripped == "TEVBLAZARS") 
    {
      //

    }


  }

  delete tokens; 

  
  if (!m->getNSources())
  {
    fprintf(stderr,"WARNING: no sources added for key %s\n", key); 
    delete m; 
    m = NULL; 
  }

  models[key] = m; 

  return m; 
}

int SourceModel::getDirectionAndEnergy( Vector & nudir, double t, double  & nuE, double minE, double maxE)
{
  if (!sources.size()) return 1; 

  bool fixedEnergy = minE==maxE;
  if (fixedEnergy) nuE = minE; 

  //pick a random source, weighted by the flux 
  double total_flux = 0 ;
  std::vector<double> fluxes(sources.size()); 

  for (unsigned i = 0; i < sources.size(); i++) 
  {
     double f = fixedEnergy ? sources[i]->getFlux()->getFlux(nuE,t) : sources[i]->getFlux()->getFluxBetween(minE,maxE,t); 
     total_flux += f; 
     fluxes[i] = total_flux; 
  }

  double random = rng.Uniform(0, total_flux); 
  const Source * which = sources[std::upper_bound(fluxes.begin(), fluxes.end(), random) - fluxes.begin()]; 

  nudir = which->getDirection(t); 

  if (!fixedEnergy) nuE =  which->getFlux()->pickEnergy(minE,maxE,t,&rng); 

  return 0; 
}


Source::Source(const char * nm, double ra, double dec, SourceFlux * f) 
: name(nm), flux(f) , RA(ra), z(sin(dec*TMath::DegToRad())) 
{
}

Vector Source::getDirection(double t) const 
{
  time_t secs = t; 
  int nsecs = 1e9 *(t-secs); 
  TTimeStamp ts(secs, nsecs); 

  double lst = ts.AsGMST();
  double h = lst - RA; 
  double phi = h  *15 * TMath::DegToRad();  

#ifdef _NICEMC_
  return Vector( cos(phi), sin(phi), z); 
#else
  return Vector( sin(phi), cos(phi), -z); 
#endif
}


ConstantExponentialSourceFlux::ConstantExponentialSourceFlux(double e, double norm, double normE)
:  gamma(e) , f("f", "[0] * x^-[1]",1e18,1e22) 
{
  //convert the normalization to the normalization at 0. 
  A = norm * TMath::Power(normE,gamma); 
  f.SetParameters(A,gamma); 
}

double ConstantExponentialSourceFlux::pickEnergy(double Emin, double Emax, double t,  TRandom * rng)  const
{
  (void) t; 
  TRandom * old = gRandom; 
  gRandom = rng; 

  //I'm too bad at math to figure out the proper quantile function here. It's probably trivial, but this works and isn't THAT slow. 
  if (f.GetXmin() != Emin && f.GetXmax() !=Emax) f.SetRange(Emin,Emax); 

  double E = f.GetRandom(); 
  gRandom = old; 
  return E; 
}



