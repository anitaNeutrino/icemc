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
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

#include "blazars/fava.h"
#include "source.hh"
#include "icemc_random.h" 

// SourceModel 


SourceModel::SourceModel(const char * n)
  : name(n)
{
  weight_Emin = 0;
  weight_Emax = 0;
  average_flux = -1; 
  average_nonzero_flux= -1; 
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


// M77 
const double m77_index = 3.16; 



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
 


SourceModel * SourceModel::getSourceModel(const char * key, SourceModel::Restriction r) 
{

  if (!key || !strcasecmp(key,"NONE")) return NULL; 

  SourceModel * m = new SourceModel(key);

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

     else if (stripped.BeginsWith("M77")) 
     {

       double index = m77_index; 
       const char * index_string = strcasestr(stripped.Data(),"INDEX="); 
       if (index_string)
       {
         sscanf(index_string,"HOURS=%lf%*s", &index); 


       }
       m->addSource(new Source("M77", 2 + 42./60 + 40.771/3600, -47.84/3600, new ConstantExponentialSourceFlux(index, txs_flux, txs_E0))); //flux and E0 don't matter 
     }

    //A3 SN
    else if(stripped=="A3SN")
    {
      TFile f("./supernovae/data/supernovaA3.root"); 
      TTree *tree = (TTree*) f.Get("sn");
      std::string *objName = new std::string;
      std::string *objType = new std::string;
      int tdiscovery;
      int texplode; 
      double ra, dec; 

      tree->SetBranchAddress("ra",&ra); 
      tree->SetBranchAddress("dec",&dec); 
      tree->SetBranchAddress("tdiscovery",&tdiscovery); 
      tree->SetBranchAddress("texplode",&texplode); 
      tree->SetBranchAddress("name",&objName); 
      tree->SetBranchAddress("type",&objType); 
      int two_weeks = 24*3600*14; 

      for (int i = 0; i < tree->GetEntries(); i++) 
      {
        tree->GetEntry(i); 
        if( abs(dec)>r.maxDec){continue;}
        m->addSource(new Source(objName->c_str(), ra, dec,
         new TimeWindowedExponentialSourceFlux(texplode, texplode+two_weeks, 2,txs_flux,txs_E0) // Just use these params as the first step
       )) ; 
      }

      delete objName; 
      delete objType; 
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
      //int two_days = 24*3600*2; 
      for (unsigned int i = 0; i < tree->GetEntries(); i++) 
      {
         tree->GetEntry(i);
              
         //////////// Cuts
         // Time cut for finding sources within a certain time period
	 // NOTE we are working with _neutrinos_ instead of photons. Neutrinos reach Earth before photons because they escape the star before the shockwaves of the supernova explosion.
	 // assume photons arrive, at max, 2 weeks after the neutrinos
	 // discoveryUnixTime is the discovery time of the supernova using photons, not neutrinos, thus:
         if( (discoveryUnixTime < r.whichStartTime) || (discoveryUnixTime > r.whichEndTime + two_weeks) ){continue;}
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
                          new TimeWindowedExponentialSourceFlux(discoveryUnixTime - two_weeks, discoveryUnixTime, 2,txs_flux,txs_E0) // Just use these params as the first step
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
      float t90, fluence, redshift; 
      
      tree->SetBranchAddress("RA",&RA); 
      tree->SetBranchAddress("dec",&dec); 
      tree->SetBranchAddress("name",&objName);
      tree->SetBranchAddress("unixTriggerTime",&unixTriggerTime);
      tree->SetBranchAddress("T90",&t90);
      tree->SetBranchAddress("fluence",&fluence);
      tree->SetBranchAddress("redshift",&redshift);

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

	// Calculate the explicit maximum neutrino energy via fireball model
	bool modelCalculation = true;
	double energyCutoff = 1e10;
	if(modelCalculation==true)
	  {

	    // Cosmology
	    double omegaRad = 9.236e-5;
	    double omegaM = 0.315;
	    double omegaLambda = 0.679;
	    double H0 = 67.4 * pow(10,3); // We want this in ms^-1/Mpc, not the standard kms^-1/Mpc
	    double c = 299792458; // ms^-1, as we we calculate c/H0
	    double EGammaIso = 0.;
	    double dL = 0.; // we will calculate in Mpc
	    double MpcToCm = 3.08e+24;

	    // GRB model vars
	    double r = 0.;
  
	    // Model assumptions
	    double lorentzFactor = 300.;
	    double eB = 0.1;
	    double eE = 0.1;
	    /// Afterglow
	    double nEx = 1.; // cm^-3
	    double EKinIso = 0.;
	    double shockLorentzFactor = 0.;
	    double energyCGamma = 0.;
	    double fPiAG = 0.;
	    double energyBreakNeutrinoAG = 0.;
	    double shockRadius = 0.;
	    double massProton = 0.0015; // erg/c^2
	    double magFieldShock = 0.;
	    double maxShockNuE = 0.;

	    // GRB spectra calculation starts
	    // Some T90s are not recorded, set to 20s.
	    if(t90==-999){t90 = 20;}
	    // Some redshifts are not calculated, set to 2.
	    if(redshift==-999){redshift = 2;}
	    // Just set the fluence to this for now
	    if(fluence==-999){fluence=1e-6;}
            
	    // Cosmology stuff
	    // Assume flat Universe
	    // Curvature is k = 0, omegaK = 0      
	    TF1 Ez("E(z)", "1/(pow([0] * pow((1+x),3) + [1] + [2] * pow((1+x),4),0.5))", 0, redshift);
	    //TF1 Ez("E(z)", "1/(pow([0] * pow((1+x),3) + [1],0.5))", 0, redshift);
	    Ez.SetParameter(0,omegaM); Ez.SetParName(0,"omegaM");
	    Ez.SetParameter(1,omegaLambda); Ez.SetParName(1,"omegaLambda");
	    Ez.SetParameter(2,omegaRad); Ez.SetParName(2,"omegaRad");
	    ROOT::Math::WrappedTF1 wEz(Ez);
	    ROOT::Math::GaussIntegrator ig;
	    ig.SetFunction(wEz);
	    ig.SetRelTolerance(0.001);
	    double integral =  ig.Integral(0, redshift);
      
	    // Now calc lum dist
	    dL = (1+redshift) * c/H0 * integral; // Mpc
      
	    // Calc gamma-ray bolometric energy
	    EGammaIso = 4 * TMath::Pi() * dL*dL * fluence * 1/(1+redshift); // Mpc^2 * erg cm^-2
	    EGammaIso*=pow(MpcToCm,2); // erg
      
	    //////////////////////////////
	    // Afterglow specific calcs
	    /////////////////////////////

	    // Calculate total jet kinetic energy and bulk Lorentz factor of GRB afterglow
	    EKinIso = EGammaIso/eE; //erg
      
	    // Radii, depths
	    shockRadius = pow(((3*EKinIso)/(4*TMath::Pi() * nEx * massProton * c*c * lorentzFactor*lorentzFactor)),(1./3.)); // in cm
	    shockRadius*=100; // in meters, to compare with internal rad
	    if(shockRadius < r)
	      {
		std::cout << "Afterglow shock radius < burst radius. Something might be wrong." << std::endl;
	      }

	    shockRadius/=100;
	    magFieldShock = pow(8 * TMath::Pi() * eB * nEx * massProton,0.5);
      
	    shockLorentzFactor = 195 * pow(((EKinIso)/(pow(10,54))),0.125) * pow(t90/10,-0.375) * pow(nEx/1.,-0.125);

	    // Peak sync gamma energy radiated by electrons in B field:
	    energyCGamma = 0.4 * pow(((EKinIso)/(pow(10,54))),(-25./24.)) * pow(((lorentzFactor)/(300)),(4./3.)) * pow(((t90)/(10)),(9./8.)) * pow(((nEx)/(1)),(-11./24.));
      
	    // Calculate proton-to-pion conversion factor
	    fPiAG = 0.2;
	    //fPiAG = 0.2 * pow(((EKinIso)/(pow(10,54))),(33./24.)) * pow(t90/10,(-9./8.)) * pow(nEx/1,(9./8.));

	    // convert to GeV
	    energyCGamma/=1e+9;
	    // break energies
	    energyBreakNeutrinoAG = 0.015 * shockLorentzFactor * 1/(1+redshift) * pow(energyCGamma,-1); // GeV
	    // max shock neutrino energy
	    maxShockNuE = fPiAG/(4*(1+redshift)) * eE * shockLorentzFactor * magFieldShock * shockRadius;
      
	    // normal case
	    if(maxShockNuE > energyBreakNeutrinoAG)
	      {
		// All good
	      }
	    else if(maxShockNuE < 1e9) // account for minimum possible energy
	      {
		maxShockNuE = 1e9;
	      }
	    else
	      {
		maxShockNuE = energyBreakNeutrinoAG+energyBreakNeutrinoAG*0.01;
	      }
	    
	    energyCutoff = maxShockNuE;

	  }

        double afterglow_hrs = 1; 

        const char * hours_string = strcasestr(stripped.Data(),"HOURS="); 
        if (hours_string)
        {
          sscanf(hours_string, "HOURS=%lf%*s", &afterglow_hrs); 
          printf("AFTERGLOW HOURS is %f\n", afterglow_hrs); 
        }

      
        if (strcasestr(stripped.Data(),"AFTERGLOW") )
        {
          m->addSource(new Source(objName->c_str(), RA/=15., dec,
                  new TimeWindowedExponentialSourceFlux(unixTriggerTime,
                                                                      unixTriggerTime + 3600*afterglow_hrs,1.5,
                                                                      fluence,1e9, energyCutoff) 
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

        // conservative case... 5 minutes before, one hour after, use afterglow flux up to 1e10 GeV. 
        else
        {
          m->addSource(new Source(objName->c_str(), RA/=15., dec,
              //fireball model for E > 1e18 eV, I think
                                new TimeWindowedExponentialSourceFlux(unixTriggerTime-300,
                                                                      unixTriggerTime + 3600*afterglow_hrs,1.5,
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
        //if (fabs(fava->dec) > r.maxDec) continue; 
              
        //say no to low HE flux 
             
        //if (fava->he_sigma < 4 && fava->he_flux < 0) continue;
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
  if (Emin == Emax) return Emin; 

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
  average_nonzero_flux = 0;

  per_source_average_flux.clear();
  per_source_average_nonzero_flux.clear();
  per_source_average_flux.resize(sources.size());
  per_source_average_nonzero_flux.resize(sources.size());
  std::vector<int> Nnonzero_per_source(sources.size());
  int Nnonzero = 0; 

  weight_Emin = minE; 
  weight_Emax = maxE; 

  TRandom * rng = getRNG(RNG_SOURCE);
  for (int i =0; i < N; i++) 
  {
    double t = rng->Uniform(t0,t1); 
    double sumFlux = 0;
    for (unsigned j = 0; j < sources.size(); j++) 
    {
      double dFlux =  minE == maxE ? sources[j]->getFlux()->getFlux(minE, t) : sources[j]->getFlux()->getFluxBetween(weight_Emin,weight_Emax,t); 
      sumFlux += dFlux; 

      average_flux += dFlux/N; 
      per_source_average_flux[j]+=dFlux/N; 
      average_nonzero_flux += dFlux; 

      if (dFlux!=0)
      {
        Nnonzero_per_source[j]++; 
        per_source_average_nonzero_flux[j]+= dFlux; 
      }
    }
    if (sumFlux != 0) 
    {
      Nnonzero++;
    }
  }

  for (unsigned j = 0; j < sources.size(); j++) 
  {
   if (per_source_average_nonzero_flux[j])  per_source_average_nonzero_flux[j] /= Nnonzero_per_source[j]; 
  }

  if (average_nonzero_flux) average_nonzero_flux /= Nnonzero; 
  printf("Set up weights: average flux is %g (average non-zero flux: %g)\n", average_flux, average_nonzero_flux); 
}

double SourceModel::getPerSourceTimeWeight(double t, int i, bool use_average_nonzero_flux) const
{
  double weight = weight_Emin == weight_Emax ? sources[i]->getFlux()->getFlux(weight_Emin,t)
                                             : sources[i]->getFlux()->getFluxBetween(weight_Emin,weight_Emax,t); 

  return weight / (use_average_nonzero_flux ? per_source_average_nonzero_flux[i] : per_source_average_flux[i]); 
}

double SourceModel::getTimeWeight(double t, bool use_average_nonzero_flux) const
{
  static int nnag = 0; 
  if (average_flux < 0) 
  {
    if (nnag++ < 100) fprintf(stderr,"WARNING: setUpWeights() hasn't been called yet. Will return crazy value!\n"); 
  }
  double weight = 0;
  for (unsigned i = 0; i < sources.size(); i++) 
  {
    weight += weight_Emin == weight_Emax ? sources[i]->getFlux()->getFlux(weight_Emin,t) : sources[i]->getFlux()->getFluxBetween(weight_Emin,weight_Emax,t); 
  }
  return weight / (use_average_nonzero_flux ? average_nonzero_flux : average_flux); 
}

