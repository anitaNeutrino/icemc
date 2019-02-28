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


SourceModel::SourceModel(const char * n, unsigned seed, std::string whichSources, std::string whichSubtype, std::string whichStartTime, std::string whichEndTime, int decCutLimit, std::string customName, double customRA, double customDec, double customGamma)
  : name(n), rng(seed), whichSources(whichSources), whichSubtype(whichSubtype), whichStartTime(whichStartTime), whichEndTime(whichEndTime), decCutLimit(decCutLimit), customName(customName), customRA(customRA), customDec(customDec), customGamma(customGamma) { }

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

SourceModel * SourceModel::getSourceModel(const char * key, unsigned seed, std::string whichSources, std::string whichSubtype, std::string whichStartTime, std::string whichEndTime, int decCutLimit, std::string customName, double customRA, double customDec, double customGamma) 
{

  if (!key || !strcasecmp(key,"NONE")) return NULL; 


  if (models.count(key)) 
    return models[key]; 

  SourceModel * m = new SourceModel(key,seed,whichSources, whichSubtype, whichStartTime, whichEndTime, decCutLimit, customName, customRA, customDec, customGamma); 

  TString str(key); 
  TObjArray * tokens = str.Tokenize("+"); 

  ///// Assign times in which we look for sources
  int TMIN = 0;
  int TMAX = 0;
  // Start times
  if(whichStartTime == "A4_MIN")
    {
      TMIN = 1480707643;
    }
  else if(whichStartTime == "A3_MIN")
    {
      TMIN = 1418938406;
    }
  else
    {
      TMIN = stoi(whichStartTime);
    }
  // End times
  if(whichEndTime == "A4_MAX")
    {
      TMAX = 1483004563;
    }
  else if(whichEndTime == "A3_MAX")
    {
      TMAX = 1420777814;
    }
  else if(whichEndTime == "MAX")
    {
      TMAX = 2147483646;
    }
  else
    {
      TMAX = stoi(whichEndTime);
    }
  /////
  
  for (int i =0; i < tokens->GetEntries(); i++)
    {
      TObjString * tok = (TObjString*) tokens->At(i); 

      TString stripped( tok->String().Strip( TString::kBoth)); 

      // The flux from TXS
      if (stripped == "TXS")
	{

	  std::string *tns_name = new std::string;
	  *tns_name = "TXS 0506+056";
      
	  m->addSource(new Source(*tns_name, 5 + 9/60. +  25.9637 / 3600, 
				  5 + 41/60. + 35.3279/3600,
				  new ConstantExponentialSourceFlux(txs_index,txs_flux,txs_E0) //from IceCube paper, assume it's constant over flight for now
				  )) ; 
	}
    
      // A custom source 
      else if (stripped=="CUSTOM")
	{

	  m->addSource(new Source(customName, customRA, 
				  customDec,
				  new ConstantExponentialSourceFlux(customGamma,txs_flux,txs_E0) // for now just use these flux constants
				  )) ; 
	
	}

      // Supernovae
      else if(stripped=="SN")
	{	
	  TFile *file = new TFile("./supernovae/data/supernovaeTNS.root"); 
	  TTree *tree = (TTree*) file->Get("tnsTree");

	  // Tree vars
	  float ra, dec;
	  std::string *objName = new std::string;
	  std::string *fullObjType = new std::string;
	  std::string *objSubtype = new std::string;
	  int discoveryUnixTime = 0;
	
	  tree->SetBranchAddress("ra",&ra); 
	  tree->SetBranchAddress("dec",&dec); 
	  tree->SetBranchAddress("name",&objName);
	  tree->SetBranchAddress("fullObjType",&fullObjType);
	  tree->SetBranchAddress("objSubtype",&objSubtype);
	  tree->SetBranchAddress("discoveryUnixTime",&discoveryUnixTime);
	
	  for (unsigned int i = 0; i < tree->GetEntries(); i++) 
	    {
	      tree->GetEntry(i);
	     
	      //////////// Cuts
	      // Time cut for finding sources within a certain time period
	      if( (discoveryUnixTime < TMIN) || (discoveryUnixTime > TMAX) ){continue;}
	      // Declination cut (ANITA won't see neutrinos from these sources)
	      if( abs(dec)>decCutLimit){continue;}
	      // Search for specific subtype
	      if(whichSources!="All")
		{
		  // Account for the chosen one
		  if(*objName==whichSources)
		    {
		      std::cout << "Using the specified source: " << whichSources << std::endl;
		    }
		  else{continue;}
		}	    
	      // Only look at a single subtype (REMEMBER this searches for a complete string
	      // but SNII can take the formats of: SN IIa ,SN IIb, SN IIa-91bg-like, etc)
	      if(whichSubtype!="All")
		{
		  if(*objSubtype==whichSubtype)
		    {
		      //std::cout << "Using the specified subtype: " << whichSubtype << std::endl;
		    }
		  else{continue;}
		}
	      ////////////
	    
	      // Legacy stuff
	      //if( abs(dec)>10){continue;}	    
	      // Only look at SNe of type II for now
	      // Core collapse SNe, associated with type II, can accelerate CRs to high energies
	      // Thus, only look for those beginning with SN II
	      //desiredType = "SN II"; //
	      //found = objType->find(desiredType);
	      //if(found != 0){continue;}
	      // Only look at 2 sources for now
	      //if(*objName != "SN 2017jy" && *objName != "SN 2017dhd"){continue;}
	    
	      //std::cout << "ra = " << ra << std::endl;
	      //std::cout << "name = " << *name << std::endl;
	      //std::cout << objName->c_str() << std::endl;
	      //std::cout << *objName << std::endl;
	      ////////////
	    
	      m->addSource(new Source(*objName, ra, dec,
				      new ConstantExponentialSourceFlux(gamma_index,txs_flux,txs_E0) // Just use these params as the first step
				      )) ; 
	    }
	
	}

      else if (stripped == "GRB") 
	{
	  TFile *file = new TFile("./GRBs/data/GRB.root"); 
	  TTree *tree = (TTree*) file->Get("grbTree");

	  // Tree vars
	  float RA, dec;
	  std::string *objName = new std::string;
	  int unixTime = 0;
	
	  tree->SetBranchAddress("RA",&RA); 
	  tree->SetBranchAddress("dec",&dec); 
	  tree->SetBranchAddress("name",&objName);
	  tree->SetBranchAddress("unixTime",&unixTime);
	
	  for (unsigned int i = 0; i < tree->GetEntries(); i++) 
	    {
	      tree->GetEntry(i);
	     
	      //////////// Cuts
	      // Time cut for finding sources within a certain time period
	      if( (unixTime < TMIN) || (unixTime > TMAX) ){continue;}
	      // Declination cut (ANITA won't see neutrinos from these sources)
	      if( abs(dec)>decCutLimit){continue;}
	      // Search for specific subtype
	      if(whichSources!="All")
		{
		  // Account for the chosen one
		  if(*objName==whichSources)
		    {
		      std::cout << "Using the specified source: " << whichSources << std::endl;
		    }
		  else{continue;}
		}
	    
	      m->addSource(new Source(*objName, RA, dec,
				      new ConstantExponentialSourceFlux(gamma_index,txs_flux,txs_E0) // Just use these params as the first step
				      )) ; 
	    }
	}
    
      else if (stripped=="FAVA")
	{
	  // Let's load all the blazars from FAVA that occurred during a flight

	  bool flux_weighted = strcasestr(stripped.Data(),"FLUX"); 
	  bool z_weighted = strcasestr(stripped.Data(),"REDSHIFT"); 

	  const char * dir = EnvironmentVariable::ICEMC_SRC_DIR() ?: "."; 

	  TFile ffava(Form("%s/blazars/fava.root", dir)); 
	  TTree * fava_tree = (TTree*) ffava.Get("fava"); 
	  FAVAEntry *fava = 0;
      
	  fava_tree->SetBranchAddress("fava",&fava);      
      
	  for (int i = 0; i < fava_tree->GetEntries(); i++) 
	    {
	      fava_tree->GetEntry(i);

	      //std::cout << "dec = " << fava->dec << std::endl;
	
	      //printf("%s %s %d %d %g %g %g\n",fava->association.GetString().Data(), fava->source_class.GetString().Data(), fava->unix_tmin, fava->unix_tmax, fava->dec, fava->he_sigma, fava->he_flux);
	
	      // time cut
	      if( (fava->unix_tmax < TMIN) || (fava->unix_tmin > TMAX) ){continue;}
	      //printf("passed time cut\n");

	      std::string *fava_name = new std::string;	
	      *fava_name = fava->association.GetString().Data();
	      //std::cout << "name = " << *fava_name << std::endl;

	      // If we didn't specific all sources
	
	      if(whichSources!="All")
		{
		  // Account for the chosen one
		  if(*fava_name==whichSources)
		    {
		      std::cout << "Using the specified source: " << whichSources << std::endl;
		    }
		  else{continue;}
		}
	
	
	      //say no to |dec| >decCutLimit
	      if (fabs(fava->dec) > decCutLimit) continue; 
	      //printf("passed dec cut\n");

	      //say no to low HE flux 
	      if (fava->he_sigma < 4 && fava->he_flux < 0) continue;
	      //printf("passed flux cut \n");

	      //only blazars
	      if (fava->source_class.GetString() ==  "bcu" || fava->source_class.GetString() == "fsrq" || fava->source_class.GetString() == "bll")
		{
      
		  m->addSource(new Source(*fava_name, fava->ra, fava->dec, 
					  new TimeWindowedExponentialSourceFlux( fava->unix_tmin, fava->unix_tmax, txs_index, 
										 txs_norm * (flux_weighted ? fava->he_flux / txs_flux : 
											     z_weighted ?   
											     txs_flux * txs_z/fava->z: 1), txs_E0))
			       ); 

		}

	    }
	}
      else
	{
	  std::cerr << "Make sure you have specified a source!" << std::endl;
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


std::string SourceModel::getDirectionAndEnergyAndOriginInfo(std::string objName, double & RA, double & dec, Vector * nudir, double t, double  & nuE, double minE, double maxE)
{
  if (!sources.size()) return "noObject";

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
      return "noObject"; 
    }
  // This is the chosen source, so we should retain some info about it
  // Then we can easily access info about *each* simulated neutrino's origin
  const Source * which = sources[index];
  objName = which->getObjName();
  RA = which->getObjRA();
  dec = which->getObjDEC();
  //std::cout << objName << std::endl;
  //std::cout << RA << std::endl;
  //std::cout << dec << std::endl;
  
  if (nudir) *nudir = which->getDirection(t); 

  if (!fixedEnergy) nuE =  which->getFlux()->pickEnergy(minE,maxE,t,&rng);
      
  return objName; 
}


Source::Source(std::string nm, double ra, double dc, SourceFlux * f) 
  : objName(nm), flux(f) , RA(ra), dec(dc*TMath::DegToRad()) 
{
}

Vector Source::getDirection(double t) const 
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
