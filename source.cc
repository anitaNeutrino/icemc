#include <map>
#include <string>
#include <iostream>
#include "TMath.h" 
#include "TObjString.h" 
#include "TTimeStamp.h" 
#include "TH1.h" 
#include "TFile.h" 
#include "TTree.h"
#include "EnvironmentVariable.h" 
#include "blazars/fava.h"
#include "source.hh"

// SourceModel 


SourceModel::SourceModel(const char * n, unsigned seed) 
  : name(n), rng(seed) { }

SourceModel::~SourceModel() 
{
  for (unsigned i = 0; i < sources.size(); i++) delete sources[i]; 
}

static std::map<std::string, SourceModel * > models; 

// TXS blazar (from icecube)
const double txs_index = 2;
const double txs_norm = 2.2e-8; 
const double txs_E0 = 0.1; 
const double txs_flux = 3.8; 
const double txs_z = 0.3365;

// SNe
const double gamma_index = 2;

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
          new ConstantExponentialSourceFlux(txs_index,txs_flux,txs_E0) //from IceCube paper, assume it's constant over flight for now
          )) ; 
    }
    
    else if(stripped=="SNA4")
      {
	
	TFile *file = new TFile("./supernovae/data/supernovaeTNS.root"); 
	TTree *tree = (TTree*) file->Get("tnsTree");

	// Tree vars
	float ra, dec;
	std::string *name = new std::string;
	std::string *objType = new std::string;
	int discoveryUnixTime = 0;

	// New vars
	int a4_tmin = 1480707643;
	int a4_tmax = 1483004563;
	
	std::string desiredType;
	size_t found;
	
	tree->SetBranchAddress("ra",&ra); 
	tree->SetBranchAddress("dec",&dec); 
	tree->SetBranchAddress("name",&name);
	tree->SetBranchAddress("objType",&objType);
	tree->SetBranchAddress("discoveryUnixTime",&discoveryUnixTime);

	for (unsigned int i = 0; i < tree->GetEntries(); i++) 
	  {
	    tree->GetEntry(i);

	    // Only look for SNs detected within A4
	    if( (discoveryUnixTime < a4_tmin) || (discoveryUnixTime > a4_tmax) ){continue;}
	    // Declination cut (ANITA won't see neutrinos from these sources)
	    if( abs(dec)>30 ){continue;}
	    // Only look at SNe of type II for now
	    // Core collapse SNe, associated with type II, can accelerate CRs to high energies
	    // Thus, only look for those beginning with SN II
	    desiredType = "SN II"; // <- DEFAULT
	    found = objType->find(desiredType);
	    if(found != 0){continue;}
	    // Only look at a single source for now
	    if(*name != "SN 2016jby"){continue;} 
	    
	    //std::cout << "ra = " << ra << std::endl;
	    //std::cout << "name = " << *name << std::endl;
	    //std::cout << name->c_str() << std::endl;
	    
	    m->addSource(new Source(name->c_str(), ra, dec,
				    new ConstantExponentialSourceFlux(gamma_index,txs_flux,txs_E0) // Just use these params as the first step
				    )) ; 
	  }
	
      }

    else if (stripped == "GRB15") 
    {
      //add ANITA3 GRB's


    }
    else if (stripped.BeginsWith("A3_FAVA"))
    {
      // Let's load all the blazars from FAVA that occurred during the A3 flight

      bool flux_weighted = strcasestr(stripped.Data(),"FLUX"); 
      bool z_weighted = strcasestr(stripped.Data(),"REDSHIFT"); 

      int a3_tmin = 1418938406; 
      int a3_tmax = 1420777814; 

      const char * dir = EnvironmentVariable::ICEMC_SRC_DIR() ?: "."; 

      TFile ffava(Form("%s/blazars/fava.root", dir)); 
      TTree * fava_tree = (TTree*) ffava.Get("fava"); 
      FAVAEntry *fava = 0; 

      fava_tree->SetBranchAddress("fava",&fava); 

      for (int i = 0; i < fava_tree->GetEntries(); i++) 
      {
        fava_tree->GetEntry(i); 

//        fava_tree->Show(i); 
       
        printf("%s %s %d %d %g %g %g\n",fava->association.GetString().Data(), fava->source_class.GetString().Data(), 
            fava->unix_tmin, fava->unix_tmax, fava->dec, fava->he_sigma, fava->he_flux); 
        if (fava->unix_tmax < a3_tmin || fava->unix_tmin > a3_tmax) continue;
        printf("passed time cut\n");


        //say no to |dec| >30 
        if (fabs(fava->dec) > 30) continue; 
        printf("passed dec cut\n");

        if (fava->he_sigma < 4 && fava->he_flux < 0) continue; //say no to low HE flux 
        printf("passed flux cut \n");

        //only blazars
        if (fava->source_class.GetString() ==  "bcu" || fava->source_class.GetString() == "fsrq" || fava->source_class.GetString() == "bll")
        {
          m->addSource(new Source(fava->association.GetString().Data(), fava->ra, fava->dec, 
                       new TimeWindowedExponentialSourceFlux( fava->unix_tmin, fava->unix_tmax, txs_index, 
                        txs_norm * (flux_weighted ? fava->he_flux / txs_flux : 
                                    z_weighted ?   
                                    txs_flux * txs_z/fava->z: 1), txs_E0))
                       ); 

        }

      }
    }
  }

  delete tokens; 

  std::cout << "----------------------" << std::endl;
  
  if (!m->getNSources())
  {
    fprintf(stderr,"WARNING: no sources added for key %s\n", key); 
    delete m; 
    m = NULL; 
  }
  else
    {
      std::cout << "Sources added: " << m->getNSources() << std::endl;
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
  unsigned index = std::upper_bound(fluxes.begin(), fluxes.end(), random) - fluxes.begin(); 
//  printf("random: %g total_flux%g, index:%u \n",random, total_flux, index); 
  if (total_flux == 0) 
  {
    nuE = minE; // do something... 
    return 1; 
  }
  const Source * which = sources[index]; 

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
:  gamma(e) , f("f", "[0] * x^-[1]",1e9,1e12) 
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
  if (f.GetXmin() != Emin || f.GetXmax() !=Emax) f.SetRange(Emin,Emax); 

  double E = f.GetRandom(); 
  gRandom = old; 
  return E; 
}



TimeWindowedExponentialSourceFlux::TimeWindowedExponentialSourceFlux(double t0, double t1, double e, double norm, double normE)
:  gamma(e) , f("f", "[0] * x^-[1]",1e9,1e12) , t0(t0), t1(t1) 
{
  //figure out what A is 
  A = norm * TMath::Power(normE,gamma); 
  f.SetParameters(A,gamma); 
}

double TimeWindowedExponentialSourceFlux::pickEnergy(double Emin, double Emax, double t,  TRandom * rng)  const
{
  if (t < t0 || t > t1) return 0; 
  TRandom * old = gRandom; 
  gRandom = rng; 

  //I'm too bad at math to figure out the proper quantile function here. It's probably trivial, but this works and isn't THAT slow. 
  if (f.GetXmin() != Emin || f.GetXmax() !=Emax) f.SetRange(Emin,Emax); 

  double E = f.GetRandom(); 
  gRandom = old; 
  return E; 
}
