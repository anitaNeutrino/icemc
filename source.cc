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
#include "icemc_random.h" 

// SourceModel 


SourceModel::SourceModel(const char * n, unsigned seed)
  : name(n)
{
  weight_Emin = 0;
  weight_Emax = 0;
  average_flux = -1; 
}

SourceModel::~SourceModel() 
{
  for (unsigned i = 0; i < sources.size(); i++) delete sources[i]; 
}


// convert redshift to MPC
static double luminosity_distance(double z) 
{
  static TF1 * f_comoving = 0; 
  if (!f_comoving) 
  {
    f_comoving = new TF1("comoving_distance", "[0]  / sqrt([1]*(1+x)^4 + [2] * (1+x)^3 + [3] * (1+x)^2 + [4])", 0,10); 
    f_comoving->SetParameters(67.74, 5e-5, 0.3089, 0, 0.6911); 
  }

  double d_M = f_comoving->Integral(0,z); 
  return (1+z) * d_M; 
} 


// TXS blazar (from icecube)
const double txs_index = 2;
const double txs_norm = 2.2e-8; 
const double txs_E0 = 100e3; 
const double txs_flux = 3.8;
const double txs_z = 0.3365;
const double txs_D = luminosity_distance(txs_z); 

// SNe
const double gamma_index = 2;

double SourceModel::Restriction::fromString(const char * the_time) 
{
  printf("PARSING \"%s\"\n",the_time); 
  // Start times
  if(strcasestr(the_time, "A4_MIN"))
    {
      return 1480707643;
    }
    else if(strcasestr(the_time, "A3_MIN"))
    {
      return 1418938406;
    }
    else if(strcasestr(the_time,"A4_MAX"))
    {
      return 1483004563;
    }
    else if(strcasestr(the_time, "A3_MAX"))
    {
      return 1420777814;
    }
   else if(strcasestr(the_time,"MAX"))
    {
      return 2147483646;
    }

  return atof(the_time); 
}
 


SourceModel * SourceModel::getSourceModel(const char * key, unsigned seed, SourceModel::Restriction r) 
{

  if (!key || !strcasecmp(key,"NONE")) return NULL; 

  SourceModel * m = new SourceModel(key,seed);

  TString str(key); 
  TObjArray * tokens = str.Tokenize("+"); 

 
  for (int i =0; i < tokens->GetEntries(); i++)
  {
     TObjString * tok = (TObjString*) tokens->At(i); 

     TString stripped( tok->String().Strip( TString::kBoth)); 

     // The flux from TXS
     if (stripped == "TXS")
     {

       double txs_ra = 77.3582; // decimal degrees
       double txs_dec = 5.69314; // decimal degrees

       m->addSource(new Source("TXS 0506+056", txs_ra/15., txs_dec,
          new ConstantExponentialSourceFlux(txs_index,txs_flux,txs_E0) //from IceCube paper, assume it's constant over flight for now
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
      
      int two_weeks = 24*3600*14; 
      int two_days = 24*3600*2; 
      for (unsigned int i = 0; i < tree->GetEntries(); i++) 
      {
         tree->GetEntry(i);
              
         //////////// Cuts
         // Time cut for finding sources within a certain time period
         if( (discoveryUnixTime < r.whichStartTime - two_weeks) || (discoveryUnixTime > r.whichEndTime) ){continue;}
             // Declination cut (ANITA won't see neutrinos from these sources)
             if( abs(dec)>r.maxDec){continue;}
             // Search for specific subtype
             if(strcasecmp(r.whichSources,"All"))
             {
                 // Account for the chosen one
                 if(strcasestr(r.whichSources, objName->c_str()))
                 {
                   std::cout << "Using the specified source: " << objName << std::endl;
                 }
                 else{continue;}
             }            
             // Only look at a single subtype (REMEMBER this searches for a complete string
             // but SNII can take the formats of: SN IIa ,SN IIb, SN IIa-91bg-like, etc)
             if(strcasecmp(r.whichSubtype,"All"))
             {
               if(strcasestr(r.whichSubtype,objSubtype->c_str()))
               {
                     //std::cout << "Using the specified subtype: " << whichSubtype << std::endl;
               }
               else{continue;}
             }
            
             m->addSource(new Source(objName->c_str(), ra/=15., dec,
                          new TimeWindowedExponentialSourceFlux(discoveryUnixTime-two_days, discoveryUnixTime + two_weeks - two_days, 2,txs_flux,txs_E0) // Just use these params as the first step
                         )) ; 
       }
     
    }

    else if (stripped.BeginsWith("GRB"))
    {
      TFile *file = new TFile("./GRBs/data/GRBIceCube.root"); 
      TTree *tree = (TTree*) file->Get("iceCubeTree");

      // Tree vars
      float RA, dec;
      std::string *objName = new std::string;
      int unixTriggerTime = 0;
      float t90,  fluence; 
      
      tree->SetBranchAddress("RA",&RA); 
      tree->SetBranchAddress("dec",&dec); 
      tree->SetBranchAddress("name",&objName);
      tree->SetBranchAddress("unixTriggerTime",&unixTriggerTime);
      tree->SetBranchAddress("T90",&t90);
      tree->SetBranchAddress("fluence",&fluence);

      for (unsigned int i = 0; i < tree->GetEntries(); i++) 
      {
        tree->GetEntry(i);
             
        if (fluence < 0) fluence = 1e-6; 

        //////////// Cuts
        // Time cut for finding sources within a certain time period
        if( (unixTriggerTime < r.whichStartTime) || (unixTriggerTime > r.whichEndTime) ){continue;}
        // Declination cut (ANITA won't see neutrinos from these sources)
        if( abs(dec)>r.maxDec){continue;}

        printf("%d\n", unixTriggerTime); 

        // Search for specific subtype
        if(strcasecmp(r.whichSources,"All"))
        {
            // Account for the chosen one
            if(strcasestr(r.whichSources,objName->c_str()))
            {
                      std::cout << "Using the specified source: " << objName << std::endl;
            }
            else{continue;}
        }



        if (strcasestr(stripped.Data(),"AFTERGLOW") )
        {
          m->addSource(new Source(objName->c_str(), RA/=15., dec,
                  new TimeWindowedExponentialSourceFlux(unixTriggerTime,
                                                                      unixTriggerTime + 3600,1.5,
                                                                      fluence,1e9, 1e10) 
                ));
         
        }
        else if (strcasestr(stripped.Data(),"PROMPT") )
        {
          m->addSource(new Source(objName->c_str(), RA/=15., dec,
              //fireball model for E > 1e18 eV, I think
                                new TimeWindowedExponentialSourceFlux(unixTriggerTime,
                                                                      unixTriggerTime + t90,4,
                                                                      fluence,1e9) 
                                )) ; 

        }

        // conservative case... 5 minutes before, one hour after, use afterglow flux. 
        else
        {
          m->addSource(new Source(objName->c_str(), RA/=15., dec,
              //fireball model for E > 1e18 eV, I think
                                new TimeWindowedExponentialSourceFlux(unixTriggerTime-300,
                                                                      unixTriggerTime + 3600,1.5,
                                                                      fluence,1e9, 1e10) 
                                )) ; 

        }
 
      }
    }
     
    else if (stripped.BeginsWith("FAVA"))
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
         
        // time cut
        if( (fava->unix_tmax < r.whichStartTime) || (fava->unix_tmin > r.whichEndTime) ){continue;}

        const char * fava_name = fava->association.GetString().Data();

	// Skip flares that are not associated with a source
	if(strcasestr(fava_name,"None"))
	  {
	    continue;
	  }

	// If we didn't specific all sources
            
        if(strcasecmp(r.whichSources,"All"))
        {
           // Account for the chosen one
           if(strcasestr(r.whichSources,fava_name))
           {
             std::cout << "Using the specified source: " << fava_name << std::endl;
           }
           else{continue;}
        }
            
            
        //say no to |dec| >30 
        if (fabs(fava->dec) > r.maxDec) continue; 
              
        //say no to low HE flux 
             
        if (fava->he_sigma < 4 && fava->he_flux < 0) continue;
        //printf("passed flux cut \n");

        //only blazars
        if (fava->source_class.GetString() ==  "bcu" || fava->source_class.GetString() == "fsrq" || fava->source_class.GetString() == "bll")
        {
       
	  m->addSource(new Source(fava_name, (fava->ra)/15., fava->dec, 
                         new TimeWindowedExponentialSourceFlux( fava->unix_tmin, fava->unix_tmax, txs_index, 
                          txs_norm * (flux_weighted ? fava->he_flux / txs_flux : z_weighted ?  txs_flux * pow(txs_D/luminosity_distance(fava->z),2): 1), txs_E0))
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
  
  return m; 
}

TH1 * SourceModel::estimateFlux(double tmin, double tmax, double Emin, double Emax, int nbins, int N) 
{
  double from = log10(Emin); 
  double to = log10(Emax); 

  TH1D * spectrum = new TH1D("spec", name, nbins, from, to); 
  spectrum->SetDirectory(0); 
  spectrum->GetXaxis()->SetTitle("log10 (E GeV)"); 
  TRandom * rng = getRNG(RNG_SOURCE);
  for (int i = 0; i < N; i++) 
    {
      double t = rng->Uniform(tmin,tmax); 
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


void SourceModel::computeFluxTimeChanges(std::vector<double> * changes) const 
{

  changes->clear(); 

  for (unsigned i = 0; i < sources.size(); i++) 
  {
    sources[i]->getFlux()->getFluxTimeChanges(changes); 
  }

  //sort and remove dupes 
  std::sort(changes->begin(), changes->end()); 
  std::vector<double>::iterator it = std::unique(changes->begin(),changes->end()); 
  changes->resize(std::distance(changes->begin(),it));
}



int SourceModel::getDirectionAndEnergy( Vector * nudir, double t, double  & nuE, double minE, double maxE)
{
  if (!sources.size()) return -1; 

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

  TRandom * rng = getRNG(RNG_SOURCE);
  double random = rng->Uniform(0, total_flux); 
  unsigned index = std::upper_bound(fluxes.begin(), fluxes.end(), random) - fluxes.begin();

  
  //  printf("random: %g total_flux%g, index:%u \n",random, total_flux, index); 
  if (total_flux == 0) 
  {
    nuE = minE; // do something... 
    return -1; 
  }
  // This is the chosen source, so we should retain some info about it
  // Then we can easily access info about *each* simulated neutrino's origin
  const Source * which = sources[index];
  //std::cout << objName << std::endl;
  //std::cout << RA << std::endl;
  //std::cout << dec << std::endl;
  
  if (nudir) *nudir = which->getDirection(t); 

  if (!fixedEnergy) nuE =  which->getFlux()->pickEnergy(minE,maxE,t,rng);
      
  return index; 
}


Source::Source(const char * nm, double ra, double dc, SourceFlux * f) 
: name(nm), flux(f) , RA(ra), dec(dc*TMath::DegToRad()) 
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



TimeWindowedExponentialSourceFlux::TimeWindowedExponentialSourceFlux(double t0, double t1, double e, double norm, double normE, double cutoff)
  :  gamma(e) , f("f", "[0] * x^-[1]",1e9,cutoff ? cutoff : 1e12) , t0(t0), t1(t1) , cutoff(cutoff) 
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
  if (cutoff && Emax > cutoff) Emax = cutoff; 

  //I'm too bad at math to figure out the proper quantile function here. It's probably trivial, but this works and isn't THAT slow. 
  if (f.GetXmin() != Emin || f.GetXmax() !=Emax) f.SetRange(Emin,Emax); 

  double E = f.GetRandom(); 
  gRandom = old; 
  return E; 
}


void SourceModel::setUpWeights(double t0, double t1, double minE, double maxE, int N) 
{

  average_flux = 0; 
  weight_Emin = minE; 
  weight_Emax = maxE; 

  TRandom * rng = getRNG(RNG_SOURCE);
  for (int i =0; i < N; i++) 
  {
    double t = rng->Uniform(t0,t1); 
    for (unsigned j = 0; j < sources.size(); j++) 
    {
      average_flux += sources[j]->getFlux()->getFluxBetween(weight_Emin,weight_Emax,t)/N; 
    }
  }
  printf("Set up weights: average flux is %g\n", average_flux); 
}


double SourceModel::getTimeWeight(double t) const
{
  static int nnag = 0; 
  if (average_flux < 0) 
  {
    if (nnag++ < 100) fprintf(stderr,"WARNING: setUpWeights() hasn't been called yet. Will return crazy value!\n"); 
  }
  double weight = 0;
  for (unsigned i = 0; i < sources.size(); i++) 
  {
    weight += sources[i]->getFlux()->getFluxBetween(weight_Emin,weight_Emax,t); 
  }
  return weight/average_flux; 
}

