#include "source.hh" 
#include "TMath.h" 
#include <map>
#include "TObjString.h" 
#include "TTimeStamp.h" 
#include "TH1.h" 




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
          new ConstantExponentialSourceFlux(2,2.2e-8,0.1) //from IceCube paper , assume it's constant over flight for now
          )) ; 
    }

    else if (stripped == "GRB15") 
    {
      //add ANITA3 GRB's


    }
    else if (stripped == "A3_FAVA") 
    {
      // Let's load all the blazars from FAVA that occurred during the A3 flight


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

TH1 * SourceModel::estimateFlux(double tmin, double tmax, double Emin, double Emax, int nbins, int N) 
{
  double from = log10(Emin); 
  double to = log10(Emax); 

  TH1D * spectrum = new TH1D("spec", name, nbins, from, to); 
  spectrum->SetDirectory(0); 
  spectrum->GetXaxis()->SetTitle("log10 (E GeV)"); 
  for (int i = 0; i < N; i++) 
  {
    double t = rng.Uniform(tmin,tmax); 
    for (unsigned j = 0; j < sources.size(); j++) 
    {
      //figure out the flux between each energy bin 
      for (int E = 1; E<= nbins; E++) 
      {
        double l = TMath::Power(10,spectrum->GetBinLowEdge(E)); 
        double h = TMath::Power(10,spectrum->GetBinLowEdge(E) + spectrum->GetBinWidth(E)); 
        spectrum->Fill(spectrum->GetBinCenter(E), sources[j]->getFlux()->getFluxBetween(l,h,t));
      }
    }
  }

  spectrum->Scale(1./N); 
  return spectrum; 
}


int SourceModel::getDirectionAndEnergy( Vector * nudir, double t, double  & nuE, double minE, double maxE)
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

  if (nudir) *nudir = which->getDirection(t); 

  if (!fixedEnergy) nuE =  which->getFlux()->pickEnergy(minE,maxE,t,&rng); 

  return 0; 
}


Source::Source(const char * nm, double ra, double dc, SourceFlux * f) 
: name(nm), flux(f) , RA(ra), dec(dc*TMath::DegToRad()) 
{
}

Vector Source::getDirection( double t) const 
{

#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
  std::cerr << "ROOT 5 doesn't have AsLMST. Returning nonsense." << std::endl;
  return Vector(1,0,0); 
#else


  time_t secs = t; 
  int nsecs = 1e9 *(t-secs); 
  TTimeStamp ts(secs, nsecs); 
  double lst = ts.AsLMST(-90);
  double h = lst - RA; 
  h *= (TMath::Pi()/12); 
  return Vector(cos(h) * cos(dec), sin(h) * cos(dec), sin(dec)); 
#endif
}


ConstantExponentialSourceFlux::ConstantExponentialSourceFlux(double e, double norm, double normE)
:  gamma(e) , f("f", "[0] * x^-[1]",1e9,1e13) 
{
  //figure out what A is 
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






