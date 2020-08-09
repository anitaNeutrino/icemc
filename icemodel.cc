#include "vector.hh"
#include "TH1.h" 
#include "TFile.h" 
#include "TChain.h"
#include "Constants.h"
#include "Settings.h"
#include "earthmodel.hh"
#include "icemodel.hh"
#include <assert.h> 


#include "signal.hh"
#include "position.hh"
#include "Primaries.h"
#include "anita.hh"
#include "ray.hh"
#include "balloon.hh"
#include "Settings.h"
#include "EnvironmentVariable.h"

#include "icemc_random.h" 

#include <fstream>
#include <iostream>


using std::endl;
using std::cout;
using std::string;
using std::ifstream;
using std::cerr;

const string ICEMC_SRC_DIR = EnvironmentVariable::ICEMC_SRC_DIR();
const string ICEMC_DATA_DIR = ICEMC_SRC_DIR+"/data/";

  
#ifdef ICEMODEL_DEBUG_TREE
ClassImp(icemodel_debug); 
#include "TTree.h" 
static TFile * debugfile = 0; 
static TTree * debugtree = 0; 

static icemodel_debug * dbg = new icemodel_debug; 


static void write_the_tree() 
{
  if (debugtree)
  {
    debugfile->cd(); 
    debugtree->Write(); 
    delete debugfile ;
  }
}

static void setupTree() 
{
  if (!debugtree) 
  {

    debugfile = new TFile(Form("icemodel_debug_%d.root",getpid()),"RECREATE"); 
    debugtree = new TTree("debug","debug"); 
    debugtree->Branch("dbg",dbg); 
    atexit(write_the_tree); 
  }
}





class FillTreeOnReturn
{

  public: 
    FillTreeOnReturn(TFile * cdto, TTree * tree)  : f(cdto), t(tree) { ; } 
    ~FillTreeOnReturn() {
      TDirectory * tmp = gDirectory; 
      f->cd(); 
      t->Fill();
      tmp->cd(); 
    } 
  private: 
    TFile * f; 
    TTree * t; 
};

#endif


static bool inHistBounds(double x, double y, const TH2 & h)
{
  if (x < h.GetXaxis()->GetXmin()) return false;
  if (y < h.GetYaxis()->GetXmin()) return false;
  if (x > h.GetXaxis()->GetXmax()) return false;
  if (y > h.GetYaxis()->GetXmax()) return false;

  return true; 
}




IceModel::IceModel(int model,int earth_model,int WEIGHTABSORPTION_SETTING) : EarthModel(earth_model,WEIGHTABSORPTION_SETTING) 
{
    
    bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(71*RADDEG)) / (1 - eccentricity*sin(71*RADDEG)) ),eccentricity/2) * tan((PI/4) - (71*RADDEG)/2); //varies with latitude, defined here for 71 deg S latitude
    bedmap_nu = bedmap_R / cos(71*RADDEG);
    
    //Parameters of the BEDMAP ice model. (See http://www.antarctica.ac.uk/aedc/bedmap/download/)
    nCols_ice=1200; //number of columns in data, set by header file (should be 1200)
    nRows_ice=1000; //number of rows in data, set by header file (should be 1000)
    cellSize=5000; //in meters, set by header file (should be 5000) - same for both files
    xLowerLeft_ice=-3000000; 
    yLowerLeft_ice=-2500000;
    nCols_ground=1068;
    nRows_ground=869;
    xLowerLeft_ground=-2661600;
    yLowerLeft_ground=-2149967;
    nCols_water=1200;
    nRows_water=1000;
    xLowerLeft_water=-3000000;
    yLowerLeft_water=-2500000;
    NODATA=-9999;
    
    h_ice_thickness.SetTitle("BEDMAP Ice Thickness; Easting (m), Northing (m)"); 

    h_ground_elevation.SetTitle("BEDMAP Ground Elevation; Easting (m), Northing (m)"); 

    h_water_depth.SetTitle("BEDMAP Water Depth; Easting (m), Northing (m)"); 
    
    sample_x = 0; 
    sample_y = 0; 
    
    ice_model=model;
    DEPTH_DEPENDENT_N = (int) (model / 10);
    ice_model -= DEPTH_DEPENDENT_N * 10;
    
    
    if (ice_model != 0 && ice_model != 1) {
	cout<<"Error!  Unknown ice model requested!  Defaulting to Crust 2.0.\n";
	ice_model = 0;
    } //if
    else if (ice_model==1) {
	ReadIceThickness();
	ReadWaterDepth();
	ReadGroundBed();
    } //else if (BEDMAP)
    //read in attenuation length data for direct signals
    int i=0;
    
    ifstream sheetup((ICEMC_DATA_DIR+"/icesheet_attenlength_up.txt").c_str());
    if(sheetup.fail())
    {
	cerr << "Failed to open icesheet_attenlength_up.txt" << endl;
	exit(1);
    }
    
    i=0;
    while(sheetup>>d_sheetup[i]>>l_sheetup[i])
    {
	i++;
    }
    sheetup.close();
    
    ifstream shelfup((ICEMC_DATA_DIR+"/iceshelf_attenlength_up.txt").c_str());
    if(shelfup.fail())
    {
	cerr << "Failed to open iceshelf_attenlength_up.txt" << endl;
	exit(1);
    }
    
    i=0;
    while(shelfup>>d_shelfup[i]>>l_shelfup[i])
    {
	i++;
    }
    shelfup.close();
    
    ifstream westlandup((ICEMC_DATA_DIR+"/westland_attenlength_up.txt").c_str());
    
    if(westlandup.fail())
    {cerr << "Failed to open westland_attenlength_up.txt";
	exit(1);
    }
    i=0;
    while(westlandup>>d_westlandup[i]>>l_westlandup[i])
    {
	i++;
    }
    westlandup.close();

    
    //read in attenuation length for downgoing signals
    ifstream sheetdown((ICEMC_DATA_DIR+"/icesheet_attenlength_down.txt").c_str());
    if(sheetdown.fail())
    {
	cerr << "Failed to open icesheet_attenlength_down.txt" << endl;
	exit(1);
    }
    
    i=0;
    while(sheetdown>>d_sheetdown[i]>>l_sheetdown[i])
    {
	i++;
    }
    sheetdown.close();
    
    ifstream shelfdown((ICEMC_DATA_DIR+"/iceshelf_attenlength_down.txt").c_str());
    if(shelfdown.fail())
    {
	cerr << "Failed to open iceshelf_attenlength_down.txt" << endl;
	exit(1);
    }
    
    i=0;
    while(shelfdown>>d_shelfdown[i]>>l_shelfdown[i])
    {
	i++;
    }
    shelfdown.close();
    
    ifstream westlanddown((ICEMC_DATA_DIR+"/westland_attenlength_down.txt").c_str());
    if(westlanddown.fail())
    {cerr << "Failed to open westland_attenlength_down.txt";
	exit(1);
    }
    i=0;
    while(westlanddown>>d_westlanddown[i]>>l_westlanddown[i])
    {
	i++;
    }
    westlanddown.close();


    cart_ice_top = 0; 
    cart_ice_bot = 0; 
    cart_ice_file = 0; 
    cart_resolution = 0; 
    cart_min_z = -1; 
    cart_max_z = -1; 
}
//constructor IceModel(int model)
//
IceModel::~IceModel() 
{
  if (cart_ice_top) delete cart_ice_top; 
  if (cart_ice_bot) delete cart_ice_bot; 
  if (cart_ice_file) delete cart_ice_file; 

}


Position IceModel::PickBalloonPosition() {
    Vector temp;
    return temp;
    
}

Position IceModel::PickInteractionLocation(int ibnposition, Settings *settings1, const Position &rbn, Interaction *interaction1) {
    
    // random numbers used
    double rnd1=0;
    double rnd3=1.;
    
    double vol_bybin=0.; // volume of ice in each bin
    int whichbin_forcrust20=0; // choice of a bin of ice for the interaction to occur in, for Crust 2.0
    int which_coord=0;// choice of coordinates for the interaction, for Bedmap
    double phi=0,theta=0; // phi and theta for interaction location
    double lon=0,lat=0; // longitude and latitude corresponding to those phi and theta
    int ilat=0,ilon=0; // ilat and ilon bins where interaction occurs
    int e_coord=0; //east coordinate of interaction, from Bedmap
    int n_coord=0; // north coordinate of interaction, from Bedmap
    
    TRandom * rng = getRNG(RNG_INTERACTION_LOCATION); 
    
    (void) settings1; 
    (void) rbn; 
    (void) interaction1; 
    
    //in case of roughness, create an array of allowable indices for the lookup (so this stay local to this function and we don't modify *horizon[ibnposition] and fark up other things downstream)
    /*if(settings1->ROUGHNESS){

      interaction1->noway=0;
      interaction1->wheredoesitleave_err=0;
      interaction1->neverseesice=0;
      interaction1->wheredoesitenterice_err=0;
      interaction1->toohigh=0;
      interaction1->toolow=0;

      bool bl_foundit = false;
      Position temppos;

      while (!bl_foundit){
        // go through selecting a shift off the balloon position to get lon,lat
        Vector blnormal = GetSurfaceNormal(rbn).Unit();
        Vector unitx = Vector();
        Vector unity;

        unitx = unitx - unitx.Dot(blnormal)*blnormal;
        unity = blnormal.Cross(unitx);

        temppos = rbn + rng->Rndm() * settings1->ROUGH_INTPOS_SHIFT * unitx + rng->Rndm() * settings1->ROUGH_INTPOS_SHIFT * unity;
        lon = temppos.Lon();
        lat = temppos.Lat();
        //if( !IceThickness(lon,lat)){   //ignore if not thick enough
        //  continue;                       // PickPosnuForaLonLat below complained too much
        //}
        bl_foundit = true;
      }

      theta = lat*RADDEG;
      phi=LongtoPhi_0is180thMeridian(lon); // convert longitude to phi
      //cout << rbn.Lon() << "  "<< rbn.Lat() << "  " <<temppos.Lon()<< "  "<<temppos.Lat()<< "  :: ";
    }
    else{ // NO roughness case*/
      //cout << "in pickinteractionlocation, size of ilat_inhorizon is " << ilat_inhorizon[ibnposition].size() << "\n";
      if (ice_model == 0) { // this is Crust 2.0
        //cout << "Inside Crust 2.0 if statement.\n";
        // vol_bybin is initialized to 0
        while(vol_bybin/maxvol_inhorizon[ibnposition]<rnd3) {
          rnd3=rng->Rndm(); // pick random numbers between 0 and 1
          rnd1=rng->Rndm();
	  
          whichbin_forcrust20=(int)(rnd1*(double)ilat_inhorizon[ibnposition].size());
          ilat=ilat_inhorizon[ibnposition][whichbin_forcrust20];

          ilon=ilon_inhorizon[ibnposition][whichbin_forcrust20];

          vol_bybin=icethkarray[ilon][ilat]*1000.*area[ilat];
        } //while
        phi=SmearPhi(ilon, rng->Rndm());
        theta=SmearTheta(ilat, rng->Rndm());
        lon = GetLon(phi);
        lat = GetLat(theta);
	//        cout << "ibnposition, phi, theta, lon, lat are " << ibnposition << " " << phi << " " << theta << " " << lon << " " << lat << "\n";
      } //end if(Crust 2.0)
      else if (ice_model==1) { // this is Bedmap
        //cout << "Inside Bedmap if statement.\n";
        while(vol_bybin/maxvol_inhorizon[ibnposition]<rnd3) {
          rnd3=rng->Rndm();
          rnd1=rng->Rndm();

          which_coord=(int)(rnd1*(double)easting_inhorizon[ibnposition].size());

          e_coord=easting_inhorizon[ibnposition][which_coord];
          n_coord=northing_inhorizon[ibnposition][which_coord];

          //      cout << "e_coord, n_coord are " << e_coord << " " << n_coord << "\n";

          GroundENtoLonLat(e_coord,n_coord,lon,lat); //Recall that the e / n coordinates in horizon were picked from the ground bed array.
          //lon+=180.;
          //cout << "lon, lat are " << lon << " " << lat << "\n";

          if (e_coord > 1068 || e_coord < 0 || n_coord < 0 || n_coord > 869)
          cout<<"Problem in PickDownward: e_coord, n_coord : "<<e_coord<<" , "<<n_coord<<endl;
          vol_bybin=IceThickness(lon,lat)*Area(lat);
        } //while
        theta = lat*RADDEG;
        phi=LongtoPhi_0is180thMeridian(lon); // convert longitude to phi
      } //end if(BEDMAP)
    //}

    //all routines (roughness or no) \it{should} deliver a lon, lat, theta, phi
    //roughness may not sometimes, should add checks or something
    Vector posnu=PickPosnuForaLonLat(lon,lat,theta,phi);
    //Vector blnormal = GetSurfaceNormal(rbn).Unit();
    //cout << blnormal.Angle(Vector(rbn[0]-posnu[0],rbn[1]-posnu[1],rbn[2]-posnu[2]))*180./PI<<"\n";
    return posnu;
} //PickInteractionLocation




static int inu = 0; 

const int PICK_RANDOM_SEGMENT = 1; 
static void sampleLocation( double len_int_kgm2, 
                           std::vector<std::pair<double,double> > & ice_intersections, 
                           Interaction *  interaction1, const Vector & x,  TRandom * rng, IceModel * model ) 

{

  //bah, this is dumb. But we don't want to have to call Getchord twice for the segments. 
    std::vector<double> interaction_weights(ice_intersections.size()); 
    std::vector<double> absorption_weights(ice_intersections.size()); 
    double cumulative_prob = 0; 
    std::vector<double> lengths(ice_intersections.size()); 
    std::vector<double> cumulative_probs(ice_intersections.size()); 
    std::vector<double> chords(ice_intersections.size()); 
    std::vector<double> nearthlayers(ice_intersections.size()); 
    std::vector<double> total_kgm2(ice_intersections.size()); 
    std::vector<int> crust_entered(ice_intersections.size()); 
    std::vector<int> mantle_entered(ice_intersections.size()); 
    std::vector<int> core_entered(ice_intersections.size()); 
    std::vector<Vector> enterices; 
    std::vector<Vector> exitices; 


    //we need a point on the trajectory for this (doesn't have to be the final location!) Use halfway within the first ice intersection as trial point. 
   
    interaction1->r_in = model->WhereDoesItEnter(x + 0.5 * (ice_intersections[0].second + ice_intersections[0].first) * interaction1->nnu,  interaction1->nnu);

	  double L_ice=len_int_kgm2/Signal::RHOICE;
    interaction1->pathlength_inice = 0;

    // we need to find the relative probability of having an interaction in each ice segment. 
    // use the probability of interacting in the ice segment (interaction weight) multiplied by the probability of not having absorbed prior to do this 
    // since these are calculated by Getchord, call Getchord once for each ice segment, using the ice entrance point. 
    // Save all the values so that they can be set when the ice interaction point is selected. 
    for (unsigned i = 0; i < ice_intersections.size(); i++) 
    {
      double this_length = ice_intersections[i].second - ice_intersections[i].first; 
      lengths[i] = this_length; 

      //get absorption weight to enter point
      enterices.push_back(x +  ice_intersections[i].first * interaction1->nnu);
      exitices.push_back( x +  ice_intersections[i].second * interaction1->nnu);
 

  
      model->Getchord(true, 
                      len_int_kgm2, 
                      interaction1->r_in, 
                      this_length, true, //include absorption from OTHER ice. This ice shouldn't be included since using entrance point
                      enterices[i], 
                      inu, //may not be perfecly synchronized... for debug only 
                      chords[i], 
                      interaction_weights[i], 
                      absorption_weights[i], 
                      nearthlayers[i], 
                      0, //TODO myair, not supported for point sources right now
                      total_kgm2[i], crust_entered[i], mantle_entered[i], core_entered[i]); 


      interaction1->pathlength_inice += this_length;

      cumulative_prob += interaction_weights[i]; 
      cumulative_probs[i] = cumulative_prob; 
    }

    int which = -1; 
    if (PICK_RANDOM_SEGMENT)
    {
      which = rng->Integer(ice_intersections.size()); 
    }
    else 
    {
      //now select the ice we interact in 
      double p = rng->Uniform(0, cumulative_prob); 
      which = std::upper_bound(cumulative_probs.begin(), cumulative_probs.end(), p) - cumulative_probs.begin(); 
    }

    interaction1->r_enterice = enterices[which]; 
    interaction1->nuexitice = exitices[which]; 
    interaction1->total_kgm2 = total_kgm2[which];
    interaction1->nearthlayers = nearthlayers[which];
    interaction1->crust_entered = crust_entered[which];
    interaction1->mantle_entered = mantle_entered[which];
    interaction1->core_entered = core_entered[which];
    interaction1->weight_nu = absorption_weights[which];
    if (PICK_RANDOM_SEGMENT ) 
    {
      interaction1->weight_nu_prob = interaction_weights[which] * ice_intersections.size(); 

    }
    else
    {
      interaction1->weight_nu_prob = cumulative_prob;
    }
    

    //let's exponentially attenuate within the ice (ignoring any gaps in ice... those will be covered in principle by the attenuation weight)
    double this_length = lengths[which]; 
    double distance = 0; 

    //If we have a very short path length in ice, it's ok to just assume it's uniform
    if (this_length/ L_ice < 1e-3) 
    {
       distance = rng->Rndm() * this_length; 
    }
    else
    {
       distance = -log(rng->Uniform(exp(-this_length/L_ice),1))*L_ice; 
    }
    



    //finally assign the position
   interaction1->posnu = enterices[which] + distance * interaction1->nnu; 
   inu++; 


#ifdef ICEMODEL_DEBUG_TREE
   dbg->sum_weight_probs = cumulative_prob; 
   for (unsigned i = 0; i < enterices.size(); i++) 
   {
     dbg->abs_weights.push_back(absorption_weights[i]); 
     dbg->weight_probs.push_back(interaction_weights[i]); 
     dbg->lengths.push_back(lengths[i]); 
   }


#endif
}

int IceModel::PickUnbiasedPointSourceNearBalloon(Interaction * interaction1, const Position * balloon_position, 
                                      double maxdist, double chord_step,  double len_int_kgm2, const Vector * force_dir) 
{

#ifdef ICEMODEL_DEBUG_TREE
    setupTree(); 
    dbg->reset() ; 
    FillTreeOnReturn fill (debugfile, debugtree); 
#endif

    // first pick the neutrino direction
    if (!force_dir) 
    {
      interaction1->PickAnyDirection();
    }
    else
    {
      interaction1->nnu = *force_dir; 
    }
    
    TRandom * rng = getRNG(RNG_INTERACTION_LOCATION); 


    //now we pick a point in the plane normal to the neutrino direction
    // We will use the point directly below ANITA with an altitude at 0 km as a reference then consider offsets
    // relative to that in that plane. 
    
    //piece of shit! 
    Position ref(balloon_position->Lon()+180, balloon_position->Lat(), R_EARTH); 


    // For now, randomly pick an offset up to maxdist km away 
    // This can be optimized by precalculating the polygon 
    // visible from this angle or finding the bounds in some other way... 
    // Especially for nearly orthogonal events this will be far too big... 
    double x= rng->Uniform(-maxdist*1e3,maxdist*1e3); 
    double y= rng->Uniform(-maxdist*1e3,maxdist*1e3); 
    sample_x =x; sample_y = y; 
     
    // get nnu as a TVector3 since I like those better
    Vector nudir = interaction1->nnu.Unit(); //is probably already a unit vector? 
    Vector orth = nudir.Orthogonal(); 
    Vector orth2 = nudir.Cross(orth);

    Vector p0 = ref + x*orth + y*orth2; 

    //Where do we intersect the Earth? 
    Position int1; 
    Position int2; 

    //check if we intersect with the (super) geoid. 
    int nintersect = GeoidIntersection(p0, nudir, &int1, &int2); 


#ifdef ICEMODEL_DEBUG_TREE
    dbg->balloon.SetXYZ(balloon_position->X(), balloon_position->Y(), balloon_position->Z()); 
    dbg->ref.SetXYZ(ref.X(), ref.Y(), ref.Z()); 
    dbg->nudir.SetXYZ(nudir.X(), nudir.Y(), nudir.Z()); 
    dbg->nintersect = nintersect; 
    dbg->x = x; 
    dbg->y = y; 
    dbg->int1.SetXYZ(int1.X(), int1.Y(), int1.Z()); 
    dbg->int2.SetXYZ(int2.X(), int2.Y(), int2.Z()); 
#endif

    if (nintersect == 0) 
    {

#ifdef ICEMODEL_DEBUG_TREE
      dbg->noway = true;
#endif
      interaction1->neverseesice=1;
      interaction1->noway = 1; 
      return 0; 
    }





    std::vector<std::pair<double,double> > ice_intersections; 
    int n_intersections = GetIceIntersectionsCartesian(int1, nudir, ice_intersections,chord_step); 

    if (n_intersections == 0) 
    {
 #ifdef ICEMODEL_DEBUG_TREE
      dbg->noway = true;
#endif
      interaction1->pathlength_inice = 0;
      interaction1->neverseesice = 1;
      return 0; 
    }
//    printf("Hits the ice %d times!\n", n_intersections); 

    sampleLocation(len_int_kgm2, ice_intersections, interaction1,int1, rng, this);



#ifdef ICEMODEL_DEBUG_TREE
    dbg->exitice.SetXYZ(interaction1->nuexitice.X(),interaction1->nuexitice.Y(),interaction1->nuexitice.Z()); 
    dbg->enterice.SetXYZ(interaction1->r_enterice.X(),interaction1->r_enterice.Y(),interaction1->r_enterice.Z()); 
    dbg->nupos.SetXYZ(interaction1->posnu.X(),interaction1->posnu.Y(),interaction1->posnu.Z()); 
#endif



    if (interaction1->posnu.Mag()-Surface(interaction1->posnu)>0) {
      interaction1->toohigh=1;
      return 0;
    }
    if (interaction1->posnu.Mag()-Surface(interaction1->posnu)+IceThickness(interaction1->posnu)<0) {
      interaction1->toolow=1;
      return 0; 
    }    


    return 1;
}

     

int IceModel::PickUnbiased(Interaction *interaction1, double len_int_kgm2, double & position_weight, double chord_step, Vector * force_dir) {
    
    TRandom * rng = getRNG(RNG_INTERACTION_LOCATION); 
    if (!force_dir) 
    {
      interaction1->PickAnyDirection(); // first pick the neutrino direction
    }
    else
    {
      interaction1->nnu = *force_dir; 
    }
    

    //now we pick a point somewhere on the surface of the superellipsoid. 
    //up to the coastline 
    //TODO, can pick closer to payload but hard to do while keeping ellipsoidal geometry... 
    //instead, will just reject points far away... 

    double sinth = rng->Uniform(0, sin(COASTLINE*RADDEG)); 
    double phi = rng->Uniform(0,2*TMath::Pi()); 

    double costh = sqrt(1-sinth*sinth); 

    const double RMAX = GEOID_MAX + 5500; 
    const double RMIN = GEOID_MIN + 5500; 
    Vector p0(RMAX * cos(phi) * sinth, RMAX * sin(phi) * sinth, RMIN*costh); 


    // The position weight is (maybe) the dot product of the position with the neutrino direction
    // this is actually the inverse weight
    position_weight = 1./fabs(p0.Unit() * interaction1->nnu); 

    //but actually, we can be double counting if the other superellipsoid intesrection is above coastline because
    //we can sample the same trajectory twice. So let's divide position_weight by two if that is indeed the case. 

    Position ints[2]; 
    GeoidIntersection(p0, interaction1->nnu, &ints[0],&ints[2]); 

    int nabovecoastline = 0; 
    for (int i =0; i < 2; i++) 
    {
      if (ints[i].Lat() < COASTLINE) nabovecoastline++; 
    }
    assert(nabovecoastline>0); 
    position_weight /= nabovecoastline; 
    interaction1->r_enterice = p0; 

    //now let's set the ice intersactions
    std::vector<std::pair<double,double> > ice_intersections; 
    int n_intersections = GetIceIntersectionsCartesian(p0, interaction1->nnu, ice_intersections,chord_step); 


    if (!n_intersections)
    {
      interaction1->noway = 1; 
      interaction1->neverseesice=1;
      interaction1->pathlength_inice = 0;
      return 0; 
    }

    sampleLocation(len_int_kgm2, ice_intersections, interaction1,p0, rng,this);
    
    
    if (interaction1->posnu.Mag()-Surface(interaction1->posnu)>0) {
        interaction1->toohigh=1;
        //cout << "inu, toohigh is " << inu << " " << interaction1->toohigh << "\n";
        return 0;
    }
    if (interaction1->posnu.Mag()-Surface(interaction1->posnu)+IceThickness(interaction1->posnu)<0) {
      interaction1->toolow=1;
      //cout << "inu, toolow is " << inu << " " << interaction1->toolow << "\n";
      return 0;
    }    

    return 1;
    
}

Vector IceModel::GetSurfaceNormal(const Position &r_out) {
    Vector n_surf = r_out.Unit();
    if (FLATSURFACE) 
	return n_surf;
    
    if (ice_model==0) {
	double theta=r_out.Theta();
	
	int ilon,ilat;
	GetILonILat(r_out,ilon,ilat);
	
	int ilon_previous=ilon-1;
	if (ilon_previous<0)
	    ilon_previous=NLON-1;
	
	int ilon_next=ilon+1;
	if (ilon_next==NLON)
	    ilon_next=0;
	
	double r=(geoid[ilat]+surfacer[ilon][ilat])*sin(theta);
	
	double slope_phi=(surfacer[ilon_next][ilat]-surfacer[ilon_previous][ilat])/(r*2*phistep);
	
	int ilat_previous=ilat-1;
	if (ilat_previous<0)
	    ilat_previous=0;
	
	int ilat_next=ilat+1;
	if (ilat_next==NLAT)
	    ilat_next=NLAT-1;
	
	double slope_costheta=(surfacer[ilon][ilat_next]-surfacer[ilon][ilat_previous])/((geoid[ilat]+surfacer[ilon][ilat])*2*thetastep);
	
	// first rotate n_surf according to tilt in costheta and position on continent - rotate around the y axis.
	double angle=atan(slope_costheta);
	
	n_surf = n_surf.RotateY(angle);
	
	// now rotate n_surf according to tilt in phi - rotate around the z axis.
	angle=atan(slope_phi);
	
	n_surf = n_surf.RotateZ(angle);
    } //end if(Crust 2.0)
    else if (ice_model==1) {
	double dist_to_check = 7500; //check elevations at this distance north, south, east and west of event
	double lon,lat;
	double lon_prev,lon_next;
	double lat_prev,lat_next;
	lon = r_out.Lon();
	lat = r_out.Lat(); //longitude and latitude of interaction
	double local_surface_elevation = Surface(lon,lat);
	
	lat_next = lat + dist_to_check * (180 / (local_surface_elevation * PI)); //the latitude 7.5 km south of the interaction
	lat_prev = lat - dist_to_check * (180 / (local_surface_elevation * PI)); //the latitude 7.5 km north of the interaction
	
	lon_next = lon + dist_to_check * (180 / (sin(lat*RADDEG) * local_surface_elevation * PI)); 
	lon_prev = lon - dist_to_check * (180 / (sin(lat*RADDEG) * local_surface_elevation * PI)); 
	
	if (lat_next > 90) {
	    //cout<<"lat_next is > 90"<<endl;
	    lat_next = 90 - (lat_next - 90);  //if we went past the pole, set coordinates for the other side
	    lon_next += 180;
	    lon_prev += 180;
	} //end if
	//cout<<"lon, lat: "<<lon<<" , "<<lat<<endl;
	//correct any out of range longitudes
	if (lon_next > 360) {
	    //cout<<"lon_next > 360\n";
	    lon_next -= 360;
	}
	else if (lon_next < 0) {
	    //cout<<"lon_next < 0\n";
	    lon_next += 360;
	}
	if (lon_prev > 360) {
	    //cout<<"lon_prev > 360\n";
	    lon_prev -= 360;
	}
	else if (lon_prev < 0) {
	    //cout << "lon_prev < 0";
	    lon_prev += 360;
	}
	
	double slope_phi=(SurfaceAboveGeoid(lon_next,lat)-SurfaceAboveGeoid(lon_prev,lat))/(2*dist_to_check);
	
	double slope_costheta=(SurfaceAboveGeoid(lon,lat_next)-SurfaceAboveGeoid(lon,lat_prev))/(2*dist_to_check);
	
	// first rotate n_surf according to tilt in costheta - rotate around the y axis.
	double angle=atan(slope_costheta);
	
	n_surf = n_surf.RotateY(angle);
	
	// now rotate n_surf according to tilt in phi - rotate around the z axis.
	angle=atan(slope_phi);
	
	n_surf = n_surf.RotateZ(angle);
    } //end if(BEDMAP)
    
    return n_surf;
    
} //method GetSurfaceNormal

int IceModel::WhereDoesItEnterIce(const Position &posnu,
				  const Vector &nnu,
				  double stepsize,
				  Position &r_enterice) {
    // now get exit point...
    //   see my geometry notes.
    // parameterize the neutrino trajectory and just see where it
    // crosses the earth radius.
    
    //  Position r_enterice;
    double distance=0;
    int left_edge=0;
    Position x = posnu;
    double x2;
    
    Position x_previous = posnu;
    
    double x_previous2= x_previous * x_previous;
    x2=x_previous2;
    
    double lon = x.Lon(),lat = x.Lat();
    double lon_old = lon,lat_old = lat;
    double local_surface = Surface(lon,lat);
    double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
    double surface_previous2=pow(local_surface,2);
    
    double rock2=rock_previous2;
    double surface2=surface_previous2;
    int foundit=0;  // keeps track of whether you found an ice entrance point
    
    //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
    //cout << "x2 at start is " << x2 << "\n";
    while (distance<2*local_surface+1000) {
	
	distance+=stepsize;
	
	x -= stepsize*nnu;
	x2=x*x;
	//cout << "x2 is " << x2 << "\n";
	lon = x.Lon();
	lat = x.Lat();
	
	double ice_thickness=IceThickness(lon,lat);
	if (lon!=lon_old || lat!=lat_old) {
	    local_surface = Surface(lon,lat);
	    
	    //if (lat>COASTLINE) 
	    //left_edge=1;
	    
	    rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
	    surface2=pow(local_surface,2);    
	    
	    if (ice_model==0) {
		if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
		    left_edge=1;
	    } //if (Crust 2.0)
	} //if (neutrino has stepped into new lon/lat bin)
	
	if ((((x_previous2>rock_previous2 && x2<rock2) // crosses rock boundary from above
	      || (x_previous2<surface_previous2 && x2>surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from below
	    || left_edge) {
	    //  cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
	    //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
	    r_enterice = x;
	    // this gets you out of the loop.
	    //continue;
	    distance=3*Geoid(lat);
	    foundit=1;
	    //cout << "foundit is " << foundit << "\n";
	    //cout << "r_enterice is ";r_enterice.Print();
	    //continue;
	} //if
	
	x_previous = x;
	x_previous2 = x2;
	//cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
	
	if (lon!=lon_old || lat!=lat_old) {
	    rock_previous2 = rock2;
	    surface_previous2 = surface2;
	    lat_old = lat;
	    lon_old = lon;
	} //if
	
    } //while
    
    return foundit;
}//WhereDoesItEnterIce

int IceModel::WhereDoesItExitIce(const Position &posnu,
				 const Vector &nnu,
				 double stepsize,
				 Position &r_enterice) {
    // now get exit point...
    //   see my geometry notes.
    // parameterize the neutrino trajectory and just see where it
    // crosses the earth radius.
    
    //  Position r_enterice;
    double distance=0;
    int left_edge=0;
    Position x = posnu;
    double x2;   
    
    Position x_previous = posnu;
    
    double x_previous2= x_previous * x_previous;
    x2=x_previous2;
    
    double lon = x.Lon(),lat = x.Lat();
    double lon_old = lon,lat_old = lat;
    double local_surface = Surface(lon,lat);
    double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
    double surface_previous2=pow(local_surface,2);
    
    double rock2=rock_previous2;
    double surface2=surface_previous2;
    int foundit=0;  // keeps track of whether you found an ice entrance point
    
    
    
    //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
    //cout << "x2 at start is " << x2 << "\n";
    int nsteps=0;
    while (distance<2*local_surface+1000) {
	//cout << "another step.\n";
	distance+=stepsize;
	nsteps++;
	//    cout << "inu, nsteps is " << inu << " " << nsteps << "\n";
	x -= stepsize*nnu;
	x2=x*x;
	//cout << "x2 is " << x2 << "\n";
	lon = x.Lon();
	lat = x.Lat();
	
	double ice_thickness=IceThickness(lon,lat);
	if (lon!=lon_old || lat!=lat_old) {
	    local_surface = Surface(lon,lat);
	    
	    //if (lat>COASTLINE) 
	    //left_edge=1;
	    
	    rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
	    surface2=pow(local_surface,2);    
	    
	    if (ice_model==0) {
		if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
		    left_edge=1;
	    } //if (Crust 2.0)
	} //if (neutrino has stepped into new lon/lat bin)
	
	if ((((x_previous2<rock_previous2 && x2>rock2) // crosses rock boundary from above
	      || (x_previous2>surface_previous2 && x2<surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from above
	    || left_edge) {
	    //  cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
	    //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
	    r_enterice = x;
	    // this gets you out of the loop.
	    //continue;
	    distance=3*Geoid(lat);
	    foundit=1;
	    //cout << "foundit is " << foundit << "\n";
	    //continue;
	} //if
	
	x_previous = x;
	x_previous2 = x2;
	//cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
	
	if (lon!=lon_old || lat!=lat_old) {
	    rock_previous2 = rock2;
	    surface_previous2 = surface2;
	    lat_old = lat;
	    lon_old = lon;
	} //if
	
    } //while
    
    return foundit;
}//WhereDoesItExitIce
int IceModel::WhereDoesItExitIceForward(const Position &posnu,
					const Vector &nnu,
					double stepsize,
					Position &r_enterice) {
    // now get exit point...
    //   see my geometry notes.
    // parameterize the neutrino trajectory and just see where it
    // crosses the earth radius.
    
    //  Position r_enterice;
    double distance=0;
    int left_edge=0;
    Position x = posnu;
    double x2;
    
    Position x_previous = posnu;
    
    double x_previous2= x_previous * x_previous;
    x2=x_previous2;
    
    double lon = x.Lon(),lat = x.Lat();
    double lon_old = lon,lat_old = lat;
    double local_surface = Surface(lon,lat);
    double rock_previous2= pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
    double surface_previous2=pow(local_surface,2);
    
    double rock2=rock_previous2;
    double surface2=surface_previous2;
    int foundit=0;  // keeps track of whether you found an ice entrance point
    
    //  cout << "lon, lat are " << posnu.Lon() << " " << posnu.Lat() << "\n";
    //cout << "x2 at start is " << x2 << "\n";
    while (distance<2*local_surface+1000) {
	
	distance+=stepsize;
	
	x += stepsize*nnu;
	x2=x*x;
	//cout << "x2 is " << x2 << "\n";
	lon = x.Lon();
	lat = x.Lat();
	
	double ice_thickness=IceThickness(lon,lat);
	if (lon!=lon_old || lat!=lat_old) {
	    local_surface = Surface(lon,lat);
	    
	    //if (lat>COASTLINE) 
	    //left_edge=1;
	    
	    rock2=pow((local_surface - IceThickness(lon,lat) - WaterDepth(lon,lat)),2);
	    surface2=pow(local_surface,2);    
	    
	    if (ice_model==0) {
		if ((int)(lat)==COASTLINE && rock_previous2 < x2 && surface2 > x2)
		    left_edge=1;
	    } //if (Crust 2.0)
	} //if (neutrino has stepped into new lon/lat bin)
	
	if ((((x_previous2<rock_previous2 && x2>rock2) // enters rock boundary from above
	      || (x_previous2>surface_previous2 && x2<surface2)) && ice_thickness>0 && lat<COASTLINE) // crosses surface boundary from below
	    || left_edge) {
	    //  cout << "lat, COASTLINE, left_edge is " << lat << " " << COASTLINE<< " " << left_edge << "\n";
	    //cout << "x_previous2, surface_previous, x2, surface2 are " << x_previous2 << " " << surface_previous2 << " " << x2 << " " << surface2 << "\n";
	    r_enterice = x;
	    // this gets you out of the loop.
	    //continue;
	    distance=3*Geoid(lat);
	    foundit=1;
	    //cout << "foundit is " << foundit << "\n";
	    //continue;
	} //if
	
	x_previous = x;
	x_previous2 = x2;
	//cout << "x_previous, x_previous2 " << x << " " << x2 << "\n";
	
	if (lon!=lon_old || lat!=lat_old) {
	    rock_previous2 = rock2;
	    surface_previous2 = surface2;
	    lat_old = lat;
	    lon_old = lon;
	} //if
	
    } //while
    
    return foundit;
}//WhereDoesItExitIceForward


double IceModel::IceThickness(double lon, double lat) {
    //This method returns the thickness of the ice in meters at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
    double ice_thickness=0;
    //cout << "ice_model is " << ice_model << "\n";
    //cout << "icethkarray is " << icethkarray[(int)(lon/2)][(int)(lat/2)]*1000. << "\n";
    if (ice_model==1) {
	double E=0;
	double N=0;
	LonLattoEN(lon,lat,E,N);
        if (inHistBounds(E,N,h_ice_thickness))
        {
	    ice_thickness = h_ice_thickness.Interpolate(E,N); //if this region has BEDMAP data, use it.
        }
	else
        {
	    ice_thickness = icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.; //if the location given is not covered by BEDMAP, use Crust 2.0 data
        }
    } //BEDMAP ice thickness
    else if (ice_model==0) {
	ice_thickness = icethkarray[(int)(lon/2)][(int)(lat/2)]*1000.;
	//cout << "ilon, ilat are " << (int)(lon/2) << " " << (int)(lat/2) << "\n";
    } //Crust 2.0 ice thickness
    
    return ice_thickness;
} //method IceThickness
double IceModel::IceThickness(const Position &pos) {
    //This method returns the thickness of the ice in meters at a location under a given position vector.  Code by Stephen Hoover.
    
    return IceThickness(pos.Lon(),pos.Lat());
} //method IceThickness(position)

double IceModel::SurfaceAboveGeoid(double lon, double lat)  {
    //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a latitude and longitude (in degrees).  In areas covered by water where no ice present, the method returns 0.  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
    // lon must be 0 to 360
    double surface=0;
    
    if (ice_model==1) {
        double E, N ; 
	LonLattoEN(lon,lat,E,N);

        if (inHistBounds(E,N, h_ground_elevation))
        {
	    surface = h_ground_elevation.Interpolate(E,N) + h_ice_thickness.Interpolate(E,N) + h_water_depth.Interpolate(E,N);
        }
	else
        {
	    surface = surfacer[(int)(lon/2)][(int)(lat/2)]; //If the position requested is outside the bounds of the BEDMAP data, use the Crust 2.0 data, regardless of the ice_model flag.
        }
    } //Elevation of surface above geoid according to BEDMAP
    else if (ice_model==0) {
	surface = surfacer[(int)(lon/2)][(int)(lat/2)];
    } //Elevation of surface above geoid according to Crust 2.0
    
    return surface;
} //method SurfaceAboveGeoid

double IceModel::SurfaceAboveGeoid(const Position &pos)  {
    //This method returns the elevation above the geoid of the surface of the ice (or bare ground, if no ice is present) in meters, at a location specified by a position vector.  Code by Stephen Hoover.
    
    return SurfaceAboveGeoid(pos.Lon(),pos.Lat());
} //method SurfaceAboveGeoid(position)

double IceModel::Surface(double lon,double lat)  {
    return (SurfaceAboveGeoid(lon,lat) + Geoid(lat)); // distance from center of the earth to surface
} //Surface

double IceModel::Surface(const Position& pos)  {
    return Surface(pos.Lon(),pos.Lat());
} //Surface

double IceModel::WaterDepth(double lon, double lat)  {
    //This method returns the depth of water beneath ice shelves in meters, at a location specified by a latitude and longitude (in degrees).  A switch in the input file can be set to determine whether the Crust 2.0 model or the BEDMAP model is used to find the ice depth.  Code by Stephen Hoover.
    double water_depth_value=0;
    
    if (ice_model==0) {
	water_depth_value = waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
    } //if(Crust 2.0)
    else if (ice_model==1) {
        double E, N; 
	LonLattoEN(lon,lat,E,N);
        if (inHistBounds(E,N,h_water_depth))
        {
	    water_depth_value = h_water_depth.Interpolate(E,N); 
        }
	else
	    water_depth_value = waterthkarray[(int)(lon/2)][(int)(lat/2)]*1000;
    } //else if(BEDMAP)
    
    return water_depth_value;
} //method WaterDepth(longitude, latitude)
double IceModel::WaterDepth(const Position &pos) {
    //This method returns the depth of water beneath ice shelves in meters, at a location specified by a position vector.  Code by Stephen Hoover.
    
    return WaterDepth(pos.Lon(),pos.Lat());
} //method WaterDepth(position)

int IceModel::IceOnWater(const Position &pos) {
    if(IceThickness(pos)>0.&&WaterDepth(pos)>0.)
	return 1;
    else return 0;
    
}
int IceModel::RossIceShelf(const Position &pos) {
    int ilon,ilat;
    
    GetILonILat(pos,ilon,ilat);
    
    if ((ilat==2 && ilon>=5 && ilon<=14) ||
	(ilat==3 && (ilon>=168 || ilon<=14)) ||
	(ilat==4 && (ilon>=168 || ilon<=13)) ||
	(ilat==5 && (ilon>=168 || ilon<=14)))
	return 1;
    else
	return 0;
}//RossIceShelf

int IceModel::RossExcept(const Position &pos) {
    int ilon,ilat;
    GetILonILat(pos,ilon,ilat);
    if(ilon<=178&&ilon>=174&&ilat>=4&&ilat<=5)
	return 1;
    else 
	return 0;
}


int IceModel::RonneIceShelf(const Position &pos) {
    int ilon,ilat;
    
    GetILonILat(pos,ilon,ilat);
    
    if ((ilat==4 && ilon>=52 && ilon<=74) ||
	(ilat==5 && ilon>=50 && ilon<=71) ||
	(ilat==6 && ilon>=55 && ilon<=64))
	return 1;
    else
	return 0;
    
}//RonneIceShelf

int IceModel::WestLand(const Position &pos) {
    double lon = pos.Lon() , lat = pos.Lat();
    
    if((lat>=4&&lat<=26)&&((lon>=0&&lon<=180)||lon>=336))
	return 1;
    else return 0;
    
}//WestLand

double IceModel::GetBalloonPositionWeight(int ibnpos) {
  //  cout << "ibnpos, volume_inhorizon, volume are " << ibnpos << " " << volume_inhorizon[ibnpos] << " " << volume << "\n";
    if (volume_inhorizon[ibnpos]==0) {
      return 0; //we will be skipped anyway 
    }
    
    return volume/volume_inhorizon[ibnpos];
} //GetBalloonPositionWeight

int IceModel::OutsideAntarctica(const Position &pos) {
    return (pos.Lat() >= COASTLINE);
} //OutsideAntarctica(Position)

int IceModel::OutsideAntarctica(double lat) {
    return (lat >= COASTLINE);
} //OutsideAntarctica(double lat)

int IceModel::AcceptableRfexit(const Vector &nsurf_rfexit,const Position &rfexit,const Vector &n_exit2rx) {
    
    //Make sure there's actually ice where the ray leaves
    if (rfexit.Lat()>COASTLINE || IceThickness(rfexit)<0.0001) {
	cout << "latitude is " << rfexit.Lat() << " compared to COASTLINE at " << COASTLINE << "\n";
	cout << "ice thickness is " << IceThickness(rfexit) << "\n";
	return 0;
	
    } //if
    
    if (nsurf_rfexit*n_exit2rx<0) {
	//cout << "dot product is " << nsurf_rfexit*n_exit2rx << "\n";
	return 0;
    } //if
    
    return 1;
} //AcceptableRfexit

double IceModel::GetN(double altitude) {
    // these are Peter's fit parameters
    double a1=0.463251;
    double b1=0.0140157;
    double n=0;
    
    if (altitude < FIRNDEPTH) 
	n=Signal::NICE;
    else if (altitude >= FIRNDEPTH && altitude <=0 && DEPTH_DEPENDENT_N) 
	//    N_DEPTH=NFIRN-(4.6198+13.62*(altitude_int/1000.))*
	//(altitude_int/1000.);   // Besson's equation for n(z)
	n=NFIRN+a1*(1.0-exp(b1*altitude));   // Peter's equation for n(z)
    else if (altitude > 0)
	cout<<"Error!  N requested for position in air!\n";
    else if (!DEPTH_DEPENDENT_N)
	n = NFIRN;
    
    return n;
} //GetN(altitude)

double IceModel::GetN(const Position &pos) {
    return GetN(pos.Mag() - Surface(pos.Lon(),pos.Lat()));
} //GetN(Position)

double IceModel::EffectiveAttenuationLength(Settings *settings1,const Position &pos,const int &whichray) {
    double localmaxdepth = IceThickness(pos);
    double depth = Surface(pos) - pos.Mag();
    //cout << "depth is " << depth << "\n";
    int depth_index=0;
    double attenuation_length=0.0;
    
    if(WestLand(pos) && !CONSTANTICETHICKNESS) 
    {
	depth_index=int(depth*419.9/localmaxdepth);//use 420 m ice shelf attenuation length data as the standard, squeeze or stretch if localmaxdepth is longer or shorter than 420m.
	if(RossIceShelf(pos) || RonneIceShelf(pos)) 
	{	  
	    if(whichray==0)
		attenuation_length=l_shelfup[depth_index];
	    else if(whichray==1)
		attenuation_length=l_shelfdown[depth_index];
	    else
		cerr << " wrong attenuation length " <<endl;
	    
	    //for sanity check
	    if((depth_index+0.5)!=d_shelfup[depth_index])
	    {
		cerr << "the index of the array l_iceshelfup is wrong!" << endl;
		exit(1);
	    }
	}
	else //in ice sheet of westland
	{
	    if(whichray==0)
		attenuation_length=l_westlandup[depth_index]; 
	    else if(whichray==1)
		attenuation_length=l_westlanddown[depth_index];
	    else
		cerr << " wrong attenuation length " <<endl;
      	}
	
	if(settings1->MOOREBAY)//if use Moore's Bay measured data for the west land
	    attenuation_length*=1.717557; //about 450 m (field attenuation length) for one whole way when assuming -3dB for the power loss at the bottom
    }
    else //in east antarctica or constant ice thickness
    { 
	
	depth_index =int(depth*(2809.9/localmaxdepth));
	
	
	if(whichray==0)
	    attenuation_length =l_sheetup[depth_index];
	else if(whichray==1)
	    attenuation_length =l_sheetdown[depth_index];
	else
	    cerr << " wrong attenuation length " <<endl;
    } //else
    
    return attenuation_length;
} //EffectiveAttenuationLengthUp

double IceModel::Area(double latitude) {
    //Returns the area of one square of the BEDMAP data at a given latitude. 
    double lat_rad = (90 - latitude) * RADDEG;
    
    return (pow(cellSize* ((1 + sin(71*RADDEG)) / (1 + sin(lat_rad))),2));
} //method Area

void IceModel::LonLattoEN(double lon, double lat,  double& E, double& N) {
    //takes as input a latitude and longitude (in degrees) and converts to indicies for BEDMAP matricies. Needs a location for the corner of the matrix, as not all the BEDMAP files cover the same area.  Code by Stephen Hoover.
    
    
    double lon_rad = (lon - 180) * RADDEG; //convert to radians, and shift origin to conventional spot
    double lat_rad = (90 - lat) * RADDEG;
    
    bedmap_R = scale_factor*bedmap_c_0 * pow(( (1 + eccentricity*sin(lat_rad)) / (1 - eccentricity*sin(lat_rad)) ),eccentricity/2) * tan((PI/4) - lat_rad/2);
    
    E = bedmap_R * sin(lon_rad);
    N = bedmap_R * cos(lon_rad);
    return;
} //method LonLattoEN


void IceModel::ENtoLonLat(int e_coord, int n_coord, double xLowerLeft, double yLowerLeft, double& lon, double& lat)  {
    //Takes as input the indicies from a BEDMAP data set, and turns them into latitude and longitude coordinates.  Information on which data set (surface data, ice depth, water depth) is necessary, in the form of coordinates of a corner of the map.  Code by Stephen Hoover.
    
    double isometric_lat=0;
    double easting = xLowerLeft+(cellSize*(e_coord+0.5)); //Add offset of 0.5 to get coordinates of middle of cell instead of edges.
    double northing = -(yLowerLeft+(cellSize*(n_coord+0.5)));
    
    //  cout << "easting, northing are " << easting << " " << northing << "\n";
    
    //first set longitude
    
    if (northing!=0)
	lon = atan(easting/northing);
    else
	lon = 90*RADDEG;
    
    // this puts lon between -pi and pi
    if (easting > 0 && lon < 0) //adjust sign of longitude
	lon += PI;
    else if (easting < 0 && lon > 0)
	lon -= PI;
    else if (easting == 0 && northing < 0)
	lon += PI;
    
    //  now find latitude
    
    if (easting != 0)
	bedmap_R = fabs(easting/sin(lon));
    else if (easting == 0 && northing != 0)
	bedmap_R = fabs(northing);
    else {
	lat = 0; //at the pole, set lat=0 degrees
	lon = lon*DEGRAD; // now put lon between 180 and 180 (only at pol)
	return;
    } //else
    
    isometric_lat = (PI/2) - 2*atan(bedmap_R/(scale_factor*bedmap_c_0));
    
    lat = isometric_lat + bedmap_a_bar*sin(2*isometric_lat) + bedmap_b_bar*sin(4*isometric_lat) + bedmap_c_bar*sin(6*isometric_lat) + bedmap_d_bar*sin(8*isometric_lat);
    
    lon = lon * DEGRAD + 180;  //convert to degrees, shift 0 to line up with bin 0 of Crust 2.0
    lat = 90 - lat*DEGRAD; //convert to degrees, with 0 degrees at the south pole
    
    //  if (lon>160 && lon<165)
    //cout << "e_coord, n_coord, easting, northing, lon are " << e_coord << " " << n_coord << " " << easting << " " << northing << " " << lon << "\n";
    return;
    
} //method ENtoLonLat

void IceModel::IceENtoLonLat(int e, int n, double& lon, double& lat) {
    //Converts indicies of the BEDMAP ice thickness matrix into longitude and latitude.  Code by Stephen Hoover.
    // cout << "I'm inside IceENtoLonLat.\n";
    ENtoLonLat(e,n,xLowerLeft_ice,yLowerLeft_ice,lon,lat);
}//IceENtoLonLat
void IceModel::GroundENtoLonLat(int e, int n, double& lon, double& lat) {
    //Converts indicies of the BEDMAP ground elevation matrix into longitude and latitude.  Code by Stephen Hoover.
    ENtoLonLat(e,n,xLowerLeft_ground,yLowerLeft_ground,lon,lat);
}//GroundENtoLonLat
void IceModel::WaterENtoLonLat(int e, int n, double& lon, double& lat) {
    //Converts indicies of the BEDMAP water depth matrix into longitude and latitude.  Code by Stephen Hoover.
    ENtoLonLat(e,n,xLowerLeft_water,yLowerLeft_water,lon,lat);
}//WaterENtoLonLat



void IceModel::GetMAXHORIZON(Balloon *bn1) {
    
    double altitude_inmeters=bn1->BN_ALTITUDE;
    if (bn1->BN_ALTITUDE==0.)
	bn1->MAXHORIZON=8.E5; // if it is a standard flight then use a horizon of 800 km
    else
	bn1->MAXHORIZON=(sqrt((R_EARTH+altitude_inmeters)*(R_EARTH+altitude_inmeters)-R_EARTH*R_EARTH))*1.1; // find distance from hrizon to balloon, increase it by 10% to be conservative.
    cout << "MAXHORIZON is " << bn1->MAXHORIZON << "\n";
}


void IceModel::CreateHorizons(Settings *settings1,Balloon *bn1,double theta_bn,double phi_bn,double altitude_bn,ofstream &foutput) {
    
    // add up volume of ice within horizon of payload
    // goes a little beyond horizon.
    
    // when we select a path in a circle at 80deg S latitude,
    // vectors are binned in phi (longitude).
    
    // when we select the Anita-lite path, 
    // vectors are binned in time.
    
    //for (int i=0; i<60;i++)
    //cout<<"area at lat "<<(double)i/2<<" is "<<Area((double)i/2)<<endl;
  
    volume = 0.; // initialize volume to zero
    
    double total_area=0; // initialize total area to zero
    int NBALLOONPOSITIONS; // number of balloon positions considered here
    if (bn1->WHICHPATH==2) // if anita-lite
	NBALLOONPOSITIONS=(int)((double)bn1->NPOINTS/(double)bn1->REDUCEBALLOONPOSITIONS); //only take 1/100 of the total balloon positions that we have because otherwise it's overkill.
    else if (bn1->WHICHPATH==6 || bn1->WHICHPATH==7 || bn1->WHICHPATH==8 || bn1->WHICHPATH==9) {
	NBALLOONPOSITIONS=(int)((double)bn1->flightdatachain->GetEntries()/(double)bn1->REDUCEBALLOONPOSITIONS)+1;
    }
    else if (bn1->WHICHPATH==1) // for picking random point along 80 deg south
	NBALLOONPOSITIONS=NPHI; // NPHI is the number of bins in phi for the visible ice in the horizon.  This is not the same as NLON, the number of bins in longitude for crust 2.0
    else // includes fixed position (bn1->WHICHPATH=0)
	NBALLOONPOSITIONS=1;
    
    volume_inhorizon.resize(NBALLOONPOSITIONS); 
    maxvol_inhorizon.resize(NBALLOONPOSITIONS); 
    ilon_inhorizon.resize(NBALLOONPOSITIONS); 
    ilat_inhorizon.resize(NBALLOONPOSITIONS); 
    easting_inhorizon.resize(NBALLOONPOSITIONS); 
    northing_inhorizon.resize(NBALLOONPOSITIONS); 
   


    double phi_bn_temp=0; //phi of balloon, temporary variable
    Position r_bn_temp; //position of balloon
    Position r_bin; // position of each bin
    
    double surface_elevation=0;
    double local_ice_thickness=0;
    double lat=0;
    double lon=0;
    //double r=0; // r,theta and phi are temporary variables
    double theta=0;
    double phi=0;
    int volume_found=0;
    char horizon_file[80];
    FILE *bedmap_horizons = 0;
    char line[200];
    int e_coord = 0;
    int n_coord = 0;
    
    if (bn1->WHICHPATH==0) 
    {
      sprintf(horizon_file,"bedmap_horizons_fixed_%g_%g_%g.dat",bn1->BN_LATITUDE, bn1->BN_LATITUDE, bn1->BN_ALTITUDE); 
    }
    else
    {
      sprintf(horizon_file,"bedmap_horizons_whichpath%i_weights%i.dat",bn1->WHICHPATH,settings1->USEPOSITIONWEIGHTS);
    }
    
    if (ice_model==1 && !settings1->WRITE_FILE) { // for bedmap model, need to be able to read file
	if(!(bedmap_horizons = fopen(horizon_file, "r"))) {
	    printf("Error: unable to open %s.  Creating new file.\n", horizon_file);
	    settings1->WRITE_FILE=1;
	}//if
    } //if
    if (ice_model==1 && settings1->WRITE_FILE) { // for bedmap model, need to be able to write to file
	if(!(bedmap_horizons = fopen(horizon_file, "w"))) {
	    printf("Error: unable to open %s\n", horizon_file);
	    exit(1);
	}//if
    } //else if
    
    if (bn1->WHICHPATH!=2 && bn1->WHICHPATH!=6) // not anita and not anita-lite
	lat=GetLat(theta_bn); //get index of latitude, same for all balloon positions
    
    
    
    for (int i=0;i<NBALLOONPOSITIONS;i++) { // loop over balloon positions
	
	maxvol_inhorizon[i]=-1.; // volume of bin with the most ice in the horizon
//        double the_time = 0; 
	
	if (bn1->WHICHPATH==2) { // anita or anita-lite path
	    theta_bn=(90+bn1->latitude_bn_anitalite[i*bn1->REDUCEBALLOONPOSITIONS])*RADDEG; //theta of the balloon wrt north pole
	    lat=GetLat(theta_bn); // latitude  
	    phi_bn_temp=(-1*bn1->longitude_bn_anitalite[i*bn1->REDUCEBALLOONPOSITIONS]+90.); //phi of the balloon, with 0 at +x and going counter clockwise looking down from the south pole
	    if (phi_bn_temp<0) //correct phi_bn if it's negative
		phi_bn_temp+=360.;
	    phi_bn_temp*=RADDEG;// turn it into radians
	    
	    
	    altitude_bn=bn1->altitude_bn_anitalite[i*bn1->REDUCEBALLOONPOSITIONS]*12.*CMINCH/100.; // get the altitude for this balloon posistion (ANITA-lite had values in feet)
	    
	} // end if anita-lite
	
       
	
	else if (bn1->WHICHPATH==6 || bn1->WHICHPATH==7 || bn1->WHICHPATH==8 || bn1->WHICHPATH==9) {
	    
	    bn1->flightdatachain->GetEvent(i*bn1->REDUCEBALLOONPOSITIONS);
//      the_time = bn1->realTime_flightdata_temp; //??!? 
	    
	    theta_bn=(90+(double)bn1->flatitude)*RADDEG; //theta of the balloon wrt north pole
	    lat=GetLat(theta_bn); // latitude  
	    phi_bn_temp=(-1*(double)bn1->flongitude+90.); //phi of the balloon, with 0 at +x and going counter clockwise looking down from the south pole
	    if (phi_bn_temp<0) //correct phi_bn if it's negative
		phi_bn_temp+=360.;
	    phi_bn_temp*=RADDEG;// turn it into radians
	    
	    altitude_bn=(double)bn1->faltitude; // for anita, altitude in is meters
	    
	} // end if anita or anita-lite
	else if (bn1->WHICHPATH==1) // for picking random position along 80 deg south
	    phi_bn_temp=dGetPhi(i); //get phi of the balloon just based on the balloon positions.  Remember nballoonpositions runs from 0 to 179 here.
	// the output of dGetPhi is in radians and runs from -pi/2 to 3*pi/2 relative to 
	else // includes bn1->WHICHPATH=0 
	    phi_bn_temp=phi_bn;
	
	// altitude_bn has already been set in SetDefaultBalloonPosition 
	// same for theta_bn
	// lat set above
	
	lon=GetLon(phi_bn_temp); //get longitude for this phi position
	// lon goes from 0 (at -180 longitude) to 360 (at 180 longitude)
	
	// position of balloon
	surface_elevation = this->Surface(lon,lat); // distance from center of the earth to surface at this lon and lat
	r_bn_temp = Vector(sin(theta_bn)*cos(phi_bn_temp)*(surface_elevation+altitude_bn),
			   sin(theta_bn)*sin(phi_bn_temp)*(surface_elevation+altitude_bn),
			   cos(theta_bn)*(surface_elevation+altitude_bn)); // vector from center of the earth to balloon
	
	// cout << "MAXHORIZON is " << MAXHORIZON << "\n";
	
	
	if (ice_model==0) { // Crust 2.0
	    for (int j=0;j<NLON;j++) { // loop over bins in longitude
		for (int k=0;k<ILAT_COASTLINE;k++) { // loop over bins in latitude
		    
		    // get position of each bin
		    r_bin = Vector(sin(dGetTheta(k))*cos(dGetPhi(j))*(geoid[k]+surfacer[j][k]),
				   sin(dGetTheta(k))*sin(dGetPhi(j))*(geoid[k]+surfacer[j][k]),
				   cos(dGetTheta(k))*(geoid[k]+surfacer[j][k])); // vector from center of the earth to the surface of this bin
		    
		    
		    if (!volume_found) 
			volume += icethkarray[j][k]*1000*area[k]; // add this to the total volume of ice in antarctica
		    if (!volume_found && icethkarray[j][k] > 0)
			total_area += area[k]; // add this to the total area of ice in antarctica
		    
		    // if the bin is within the maximum horizon of the balloon or if we don't care
		    //r=(geoid[k]+surfacer[j][k]);
		    //phi=dGetPhi(j);
		    //theta=dGetTheta(k);
		    
		    //	  cout << "USEWEIGHTS is " << USEWEIGHTS << "\n";
		    if ((settings1->USEPOSITIONWEIGHTS && r_bin.Distance(r_bn_temp)<bn1->MAXHORIZON)
			|| !settings1->USEPOSITIONWEIGHTS) {
			// then put this latitude and longitude in vector
			
			
			ilat_inhorizon[i].push_back(k); 
			ilon_inhorizon[i].push_back(j);
			// add up volume in horizon
			
			
			volume_inhorizon[i]+=icethkarray[j][k]*1000*area[k];
			
			
			// finding which bin in horizon has maximum volume
			if (icethkarray[j][k]*1000*area[k]>maxvol_inhorizon[i]) {
			    maxvol_inhorizon[i]=icethkarray[j][k]*1000.*area[k];
			}
		    } //end if (distance < 800 km)
		    
		    
		} //end for (k loop)
	    } //end for (j loop)   
	    
	    //      cout << "i, volume_inhorizon are " << i << " " << volume_inhorizon[i] << "\n";
	    
	    // ifi the balloon is too close to the ice, it will think there aren't any
	    // bins in the horizon, so force it the include the ones just below the balloon
	    int ilat_bn,ilon_bn;
	    GetILonILat(r_bn_temp,ilon_bn,ilat_bn); // find which longitude and latitude the balloon is above
	    
	    if (ilat_inhorizon[i].size()==0 || ilon_inhorizon[i].size()==0) {
		maxvol_inhorizon[i]=icethkarray[ilon_bn][ilat_bn]*1000.*area[ilat_bn]; // need to give it a maximum volume for a bin in horizon 
		volume_inhorizon[i]=icethkarray[ilon_bn][ilat_bn]*1000.*area[ilat_bn]; // and a total volume for the horizon
	    }
	    
	    if (ilat_inhorizon[i].size()==0) // for the ith balloon position, if it didn't find a latitude bin in horizon
		ilat_inhorizon[i].push_back(ilat_bn); // force it to be the one below the balloon       
	    
	    if (ilon_inhorizon[i].size()==0) // for the ith balloon position, if it didn't find a longitude bin in horizon
		ilon_inhorizon[i].push_back(ilon_bn); // force it to be the one below the balloon
	    
	    
	    
	} //end if (ice_model==0) Crust 2.0
	
	else if (ice_model==1 && !settings1->WRITE_FILE) { // for bedmap model
	    fgets(line,200,bedmap_horizons);
	    while(line[0] != 'X') {
		e_coord = atoi(strtok(line,","));
		n_coord = atoi(strtok(NULL,","));
		easting_inhorizon[i].push_back(e_coord);
		northing_inhorizon[i].push_back(n_coord);
		fgets(line,200,bedmap_horizons);
	    } //while there are more eastings and northings
	    strtok(line,",");
	    volume_inhorizon[i] = atof(strtok(NULL,",")); // volume in the horizon
	    maxvol_inhorizon[i] = atof(strtok(NULL,","));
	    
	    if (!volume_found) {
		total_area = atof(fgets(line,200,bedmap_horizons)); // total area on the continent
		volume = atof(fgets(line,200,bedmap_horizons)); // total volume on the continent
	    } //if
	} //end if (ice_model==1) && !settings1->WRITE_FILE
	
	else if (ice_model==1 && settings1->WRITE_FILE) { //for BEDMAP model, look through all easting and northing coordinates in the groundbed map (our smallest).  Output what we find to a file for later use.
	    
	    for (int n_coord=0;n_coord<nRows_ground;n_coord++) {
		for (int e_coord=0;e_coord<nCols_ground;e_coord++) {
		    
		    GroundENtoLonLat(e_coord,n_coord,lon,lat);
		    
		    theta = lat * RADDEG;
		    phi=LongtoPhi_0is180thMeridian(lon);
		    
		    surface_elevation = this->Surface(lon,lat);
		    local_ice_thickness = this->IceThickness(lon,lat);
		    
		    // get position of each bin
		    r_bin = Vector(sin(theta)*cos(phi)*surface_elevation,
				   sin(theta)*sin(phi)*surface_elevation,
				   cos(theta)*surface_elevation);
		    
		    if (!volume_found)
			volume += local_ice_thickness*Area(lat);
		    if (!volume_found && local_ice_thickness > 5)
			total_area += Area(lat);
		    if ((settings1->USEPOSITIONWEIGHTS && r_bn_temp.Distance(r_bin)<bn1->MAXHORIZON) || !settings1->USEPOSITIONWEIGHTS) {
			fprintf(bedmap_horizons,"%i,%i,\n",e_coord,n_coord);
			// then put this latitude and longitude in vector
			easting_inhorizon[i].push_back(e_coord);
			northing_inhorizon[i].push_back(n_coord);
			// add up volume in horizon
			
			GroundENtoLonLat(e_coord,n_coord,lon,lat);
			
			
			volume_inhorizon[i]+=local_ice_thickness*Area(lat);
			
			// finding which bin in horizon has maximum volumey
			if (local_ice_thickness*Area(lat)>maxvol_inhorizon[i]) {
			  maxvol_inhorizon[i]=local_ice_thickness*Area(lat);
			}
		    } //end if (distance < 800 km & ice present)
		} //end for (e_coord loop)
	    } //end for (n_coord loop)
	    
	    fprintf(bedmap_horizons,"X,%f,%f,\n",volume_inhorizon[i],maxvol_inhorizon[i]);
	    if (!volume_found) {
		fprintf(bedmap_horizons,"%f\n",total_area);
		fprintf(bedmap_horizons,"%f\n",volume);
	    } //if
	    
            fclose(bedmap_horizons); 

	} //end if (ice_model==1) && settings1->WRITE_FILE
	
	if (!volume_found) {
	    cout<<"Total surface area covered with ice (in m^2) is : "<<total_area<<endl;
	    volume_found=1;
	} //if

    } //end loop over balloon positions
    
    // finding average volume in horizon over all balloon positions.
    volume_inhorizon_average=0;
    
    for (int i=0;i<NBALLOONPOSITIONS;i++) {
	volume_inhorizon_average+=volume_inhorizon[i];
	
    } //for
    volume_inhorizon_average/=(double)NBALLOONPOSITIONS;
    
    //cout << "Total volume of ice in Antarctica with this ice model (m^3): " << volume << "\n";
    //cout << "Average volume of ice within a horizon is " << volume_inhorizon_average << "\n";
    
    foutput << "Average volume of ice within a horizon is " << volume_inhorizon_average << "\n";
    
    foutput << "Average thickness of ice within horizon, averages over balloon positions " << volume_inhorizon_average/PI/pow(bn1->MAXHORIZON,2) << "\n";
} //method CreateHorizons 


void IceModel::ReadIceThickness() {
    //Reads the BEDMAP ice thickness data.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
    
  ifstream IceThicknessFile((ICEMC_DATA_DIR+"/icethic.asc").c_str());
    if(!IceThicknessFile) {
	cerr << "Couldn't open: $ICEMC_DATA_DIR/icethic.asc" << endl;
	exit(1);
    }
    
    cout<<"Reading in BEDMAP data on ice thickness.\n";
    
    string tempBuf1;
    string tempBuf2;
    string tempBuf3;
    string tempBuf4;
    string tempBuf5;
    string tempBuf6;
    int temp1,temp2,temp3,temp4,temp5,temp6;
    
    IceThicknessFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
    >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
    >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
    
    if(tempBuf1 == string("ncols")) {
	nCols_ice=temp1;
    }
    if(tempBuf2 == string("nrows")) {
	nRows_ice=temp2;
    }
    if(tempBuf3 == string("xllcorner")) {
	xLowerLeft_ice=temp3;
    }
    if(tempBuf4 == string("yllcorner")) {
	yLowerLeft_ice=temp4;
    }
    if(tempBuf5 == string("cellsize")) {
	cellSize=temp5;
    }
    if(tempBuf6 == string("NODATA_value")) {
	NODATA=temp6;
    }
    //cout<<"nCols_ice, nRows_ice "<<nCols_ice<<" , "<<nRows_ice<<endl;
    //cout<<"xLL_ice, yLL_ice, cellsize "<<xLowerLeft_ice<<" , "<<yLowerLeft_ice<<" , "<<cellSize<<endl<<endl;
    
    assert(nRows_ice == 1000 && nCols_ice==1200); 

    h_ice_thickness.SetBins(nCols_ice, xLowerLeft_ice, xLowerLeft_ice + nCols_ice*cellSize, 
                           nRows_ice, yLowerLeft_ice, yLowerLeft_ice+nRows_ice*cellSize);

    double theValue;
    volume=0.;
    ice_area=0.;
    max_icevol_perbin=0.;
    max_icethk_perbin=0.;
    double lon_tmp,lat_tmp;
    for(int rowNum=0;rowNum<nRows_ice;rowNum++) {
	for(int colNum=0;colNum<nCols_ice;colNum++) {
	    IceENtoLonLat(colNum,rowNum,lon_tmp,lat_tmp);
	    IceThicknessFile >> theValue;
	    if(theValue==NODATA)
		theValue=0; //Set ice depth to 0 where we have no data.

            //note that the histogram has a backwards indexing on rows since it's defined in a different order
	    h_ice_thickness.SetBinContent(1+colNum, nRows_ice-rowNum, theValue); 
	    volume+=theValue*Area(lat_tmp);
	    if (theValue>0)
		ice_area+=Area(lat_tmp);
	    if (theValue*Area(lat_tmp)>max_icevol_perbin)
		max_icevol_perbin=theValue*Area(lat_tmp);
	    if (theValue>max_icethk_perbin)
		max_icethk_perbin=theValue;
	}//for
    }//for
    
    IceThicknessFile.close();
    return;
} //method ReadIceThickness

void IceModel::ReadGroundBed() {
    //Reads the BEDMAP data on the elevation of the ground beneath the ice.  If there is water beneath the ice, the ground elevation is given the value 0.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  ifstream GroundBedFile((ICEMC_DATA_DIR+"/groundbed.asc").c_str());
    if(!GroundBedFile) {
	cerr << "Couldn't open: ICEMC_DATA_DIR/groundbed.asc" << endl;
	exit(1);
    }
    
    cout<<"Reading in BEDMAP data on elevation of ground.\n";
    
    string tempBuf1;
    string tempBuf2;
    string tempBuf3;
    string tempBuf4;
    string tempBuf5;
    string tempBuf6;
    int temp1,temp2,temp3,temp4,temp5,temp6;
    
    GroundBedFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
    >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
    >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
    
    if(tempBuf1 == string("ncols")) {
	nCols_ground=temp1;
    }
    if(tempBuf2 == string("nrows")) {
	nRows_ground=temp2;
    }
    if(tempBuf3 == string("xllcorner")) {
	xLowerLeft_ground=temp3;
    }
    if(tempBuf4 == string("yllcorner")) {
	yLowerLeft_ground=temp4;
    }
    if(tempBuf5 == string("cellsize")) {
	cellSize=temp5;
    }
    if(tempBuf6 == string("NODATA_value")) {
	NODATA=temp6;
    }
    assert(nRows_ground == 869 && nCols_ground==1068); 
    
    h_ground_elevation.SetBins(nCols_ground, xLowerLeft_ground, xLowerLeft_ground + nCols_ground*cellSize, 
                               nRows_ground, yLowerLeft_ground, yLowerLeft_ground + nRows_ground*cellSize);

    //cout<<"nCols_ground, nRows_ground "<<nCols_ground<<" , "<<nRows_ground<<endl;
    //cout<<"xLL_ground, yLL_ground, cellsize "<<xLowerLeft_ground<<" , "<<yLowerLeft_ground<<" , "<<cellSize<<endl<<endl;
    
    double theValue;
    for(int rowNum=0;rowNum<nRows_ground;rowNum++) {
	for(int colNum=0;colNum<nCols_ground;colNum++) {
	    GroundBedFile >> theValue;
	    
	    if(theValue==NODATA)
		theValue=0; //Set elevation to 0 where we have no data.
	    h_ground_elevation.SetBinContent(colNum+1,nRows_ground-rowNum, double(theValue));
	    //if (theValue != -96 && theValue != 0)
	    //cout<<"ground_elevation: "<<theValue<<endl;
	}//for
    }//for
    
    GroundBedFile.close();
    return;
} //method ReadGroundBed

void IceModel::ReadWaterDepth() {
    //Reads BEDMAP data on the depth of water beneath the ice.  Where no water is present, the value 0 is entered.  Assumes the file is in directory "data".  Code by Ryan Nichol, added to Monte Carlo by Stephen Hoover
  ifstream WaterDepthFile((ICEMC_DATA_DIR+"/water.asc").c_str());
    if(!WaterDepthFile) {
	cerr << "Couldn't open: ICEMC_DATA_DIR/water.asc" << endl;
	exit(1);
    }
    
    cout<<"Reading in BEDMAP data on water depth.\n";
    
    string tempBuf1;
    string tempBuf2;
    string tempBuf3;
    string tempBuf4;
    string tempBuf5;
    string tempBuf6;
    int temp1,temp2,temp3,temp4,temp5,temp6;
    
    WaterDepthFile >> tempBuf1 >> temp1 >> tempBuf2 >> temp2 
    >> tempBuf3 >> temp3 >> tempBuf4 >> temp4 
    >> tempBuf5 >> temp5 >> tempBuf6 >> temp6;
    
    if(tempBuf1 == string("ncols")) {
	nCols_water=temp1;
    }
    if(tempBuf2 == string("nrows")) {
	nRows_water=temp2;
    }
    if(tempBuf3 == string("xllcorner")) {
	xLowerLeft_water=temp3;
    }
    if(tempBuf4 == string("yllcorner")) {
	yLowerLeft_water=temp4;
    }
    if(tempBuf5 == string("cellsize")) {
	cellSize=temp5;
    }
    if(tempBuf6 == string("NODATA_value")) {
	NODATA=temp6;
    }
    
    //cout<<"nCols_water, nRows_water "<<nCols_water<<" , "<<nRows_water<<endl;
    //cout<<"xLL_water, yLL_water, cellsize "<<xLowerLeft_water<<" , "<<yLowerLeft_water<<" , "<<cellSize<<endl<<endl;
    
    assert(nRows_water == 1000 && nCols_water==1200); 

    h_water_depth.SetBins(nCols_water, xLowerLeft_water, xLowerLeft_water + nCols_water*cellSize, 
                          nRows_water, yLowerLeft_water, yLowerLeft_water+nRows_water*cellSize);
    double theValue;
    for(int rowNum=0;rowNum<nRows_water;rowNum++) {
	for(int colNum=0;colNum<nCols_water;colNum++) {
	    WaterDepthFile >> theValue;
	    
	    if(theValue==NODATA)
		theValue=0; //Set depth to 0 where we have no data.
	    h_water_depth.SetBinContent(1+colNum,nRows_water-rowNum, double(theValue));
	}//for
    }//for
    
    WaterDepthFile.close();
    return;
} //method ReadWaterDepth




//Methods related to the 2D cartesian surfaces
//


void IceModel::CreateCartesianTopAndBottom(int resolution, bool force)  
{
  TString fname;

  if (resolution == cart_resolution) 
  {
    return ;
  } 

  cart_resolution = resolution; 

  if (cart_ice_file) delete cart_ice_file; 

  fname.Form("%s/cartesian_icemodel_%d_%d.root", getenv("ICEMC_SRC_DIR")?: ".", ice_model ,  resolution); 

  TDirectory * olddir = gDirectory; 

  //check if it exists
  if (!force) 
  {
    cart_ice_file = new TFile(fname); 

    if (!cart_ice_file->IsOpen())
    {
      delete cart_ice_file; 
    }
    else 
    {
      cart_ice_top = (TH2*) cart_ice_file->Get("top");
      cart_ice_bot = (TH2*) cart_ice_file->Get("bot");

      if (!cart_ice_top || !cart_ice_bot) 
      {
        delete cart_ice_file; 
      }

      else
      {
        olddir->cd();
        return; 
      }
    }
  }

  //if we got this far, we don't already have it 
  cart_ice_file = new TFile(fname,"RECREATE"); 
  //the bounds will roughly be +/- REARTH /2,so we'll convert that to our resolution

  printf("Creating new cartesian ice file at %s...\n", cart_ice_file->GetName()); 

  double bound = ceil(6371000/2/resolution)*resolution; 
  int nbins = (2*bound)/resolution; 
  
  TString hname; 
  hname.Form("icetop_%d_%d", ice_model, resolution); 
  TString htitle; 
  htitle.Form("Ice Top (model=%d, res=%d m)", ice_model, resolution); 
  cart_ice_top = new TH2D(hname.Data(), htitle.Data(), nbins,-bound,bound,nbins,-bound,bound);

  hname.Form("icebot_%d_%d", ice_model, resolution); 
  htitle.Form("Ice Bottom (model=%d, res=%d m)", ice_model, resolution); 
  cart_ice_bot= new TH2D(hname.Data(), htitle.Data(), nbins,-bound,bound,nbins,-bound,bound);

  // auxiliary hists
  hname.Form("ice_bot_r_%d_%d",ice_model, resolution);
  htitle.Form("Ice Bottom (r) (model=%d, res=%d m)", ice_model, resolution); 

  TH2 * r_bot= new TH2D(hname.Data(), htitle.Data(), nbins,-bound,bound,nbins,-bound,bound);

  hname.Form("ice_top_r_%d_%d",ice_model, resolution);
  htitle.Form("Ice top (r) (model=%d, res=%d m)", ice_model, resolution); 

  TH2 * r_top= new TH2D(hname.Data(), htitle.Data(), nbins,-bound,bound,nbins,-bound,bound);

  hname.Form("ice_bot_alt_%d_%d",ice_model, resolution);
  htitle.Form("Ice Bottom (alt) (model=%d, res=%d m)", ice_model, resolution); 

  TH2 * alt_bot= new TH2D(hname.Data(), htitle.Data(), nbins,-bound,bound,nbins,-bound,bound);

  hname.Form("ice_top_alt_%d_%d",ice_model, resolution);
  htitle.Form("Ice top (alt) (model=%d, res=%d m)", ice_model, resolution); 

  TH2 * alt_top= new TH2D(hname.Data(), htitle.Data(), nbins,-bound,bound,nbins,-bound,bound);


  for (int j = 1; j <= nbins; j++) 
  {
    double y = cart_ice_top->GetYaxis()->GetBinCenter(j);
    for (int i = 1; i <= nbins; i++) 
    {
      double x = cart_ice_top->GetXaxis()->GetBinCenter(i);

      double r2 = x*x+y*y; 
      const double a = 6378137;
      const double a2 = a*a; 
      if (r2 > a2) continue; 
      const double b = 6356752.314245; ;
      const double b2 = b*b; 
      double z_geoid = sqrt((1-r2/a2)*b2); 

      Position p(Vector(x,y,z_geoid)); 
      double lon = p.Lon(); 
      double lat = p.Lat(); 
      double ice_depth = IceThickness(lon,lat); 
      if (ice_depth) 
      {
        double surface_above_geoid = SurfaceAboveGeoid(lon,lat); 
        double geoid = Geoid(lat); 
        double surface = surface_above_geoid+geoid; 

        Position top(lon,lat,surface);
        Position bot(lon,lat,surface-ice_depth);

        cart_ice_top->SetBinContent(i,j,top.Z()); 
        cart_ice_bot->SetBinContent(i,j,bot.Z()); 

        r_top->SetBinContent(i,j,surface); 
        r_bot->SetBinContent(i,j,surface-ice_depth); 

        alt_top->SetBinContent(i,j,surface_above_geoid); 
        alt_bot->SetBinContent(i,j,surface_above_geoid-ice_depth); 
      }
    }
  }

  cart_ice_top->Write("top");
  cart_ice_bot->Write("bot");
  r_top->Write("r_top");
  r_bot->Write("r_bot");
  alt_top->Write("alt_top");
  alt_bot->Write("alt_bot");
  cart_ice_file->Flush(); 
  printf("...Done!\n"); 
  olddir->cd();

}


bool IceModel::CartesianIsInIce(double x,double y, double z) 
{
  const TH2 * htop = GetCartesianTop();
  if (z > cart_max_z) return false; 
  if (z < cart_min_z) return false; 
  if (x < htop->GetXaxis()->GetXmin()) return false;
  if (x > htop->GetXaxis()->GetXmax()) return false;
  if (y < htop->GetYaxis()->GetXmin()) return false;
  if (y > htop->GetYaxis()->GetXmax()) return false;

  double top = ((TH2*)htop)->Interpolate(x,y); 
  if (!top) return false; 

  return z < top && z > ((TH2*)GetCartesianBottom())->Interpolate(x,y);
}

int IceModel::GetIceIntersectionsCartesian(const Position & posnu, const Vector & nnu_in, 
                                           std::vector<std::pair<double,double> > & intersection_points, 
                                           double start_step, int map_resolution) 
{
  if (map_resolution!= cart_resolution) CreateCartesianTopAndBottom(map_resolution); 

  if (cart_max_z == -1)
  {
    cart_max_z =  GetCartesianTop()->GetMaximum();

    cart_min_z = cart_max_z; 

    const TH2 * hbot = GetCartesianBottom(); 

    for (int j = 1; j <= hbot->GetNbinsY(); j++) 
    {
      for (int i = 1; i <= hbot->GetNbinsX(); i++) 
      {
        double val = hbot->GetBinContent(i,j); 
        if (val > 0 && val < cart_min_z) cart_min_z = val; 
      }
    }
  }


  std::vector<double> boundaries; 

  //is our start point in the ice? 

  bool start_in_ice = CartesianIsInIce(posnu.X(), posnu.Y(), posnu.Z()); 

  Vector nnu = nnu_in.Unit(); 
  //let's start at the position and go both ways
  for (int dir = -1; dir <=1 ; dir+=2) 
  {

    Vector x = posnu; 
    bool was_in_ice = start_in_ice; 

    // go until we are way to far away from the Earth to matter
    
    double displacement = 0;
    bool jumped = false;
    bool just_jumped = false;
    while (true) //this is mean Earth radius + 5 km. 
    {
//      printf("%g  ->  %g,%g,%g\n", displacement, x.X(),x.Y(), x.Z()); 
      displacement += start_step*dir; 
      x = posnu + displacement*nnu; 
      double mag2 = x.X()*x.X() + x.Y()*x.Y() + x.Z()*x.Z(); 
 //     printf(" %g\n", sqrt(mag2)); 
      if (mag2 > (7000e3*7000e3))
      {
        printf("wtf:  [%g,%g,%g], %g\n", x.X(), x.Y(), x.Z(), mag2); 

      }
      if (mag2 > (6365e3*6365e3))
      {
        break; 
      }

      //jump to the other side of the world if we dip down too much! 
      if (mag2 < (6355.5e3*6355.5e3) && !jumped && !was_in_ice) 
      {
       // printf("Trying to go to other side of the earth! start displacement: %g, mag: %g\n", displacement,sqrt(mag2));
        double ds[2]; 
        GeoidIntersection(x,dir*nnu, 0,0,-1000,ds); 
        //now figure out which is on the opposite side of the Earth from where we are...  since we're using our current position and direction it should be the positive one! 
        displacement = (ds[0] > 0 ? ds[0] : ds[1])*dir; 
        //printf("End displacement: %g\n", displacement); 
        jumped = true; 
        just_jumped = true; 
        continue; 
      }
      
      bool is_in_ice = CartesianIsInIce(x.X(),x.Y(),x.Z()); 

      //this is no bueno...we won't find the right boundary in this case!  
      // let's back up 500 m from where we were at the start of this iteration
      if (just_jumped && is_in_ice) 
      {
        displacement -= dir*(500+start_step); 
        continue; 
      }

      just_jumped = false; 

      if (was_in_ice!=is_in_ice) 
      {
        //binary search for the boundary to 1 m resolution

        Vector xsearch = x; 
        double step = -start_step*dir/2; 
        bool search_was_in_ice = is_in_ice; 
        double boundary_displacement = displacement; 

        while (fabs(step) > 1) 
        {
          boundary_displacement+=step;
          xsearch = posnu + boundary_displacement * nnu; 
          bool search_is_in_ice = CartesianIsInIce(xsearch.X(),xsearch.Y(),xsearch.Z()); 

          step*=0.5; 
          if (search_is_in_ice!=search_was_in_ice) step*=-1; 
          search_was_in_ice = search_is_in_ice; 
        }

        //we found a boundary! (if there is more than one boundary within this segment... we probably picked a random one. Oh well. 
        boundaries.push_back(boundary_displacement); 
      }
      was_in_ice = is_in_ice; 
    }
  }

  std::sort(boundaries.begin(), boundaries.end()); 

  if (boundaries.size() % 2) 
  {
    fprintf(stderr,"Uh oh, have an odd number of boundaries (%lu) . That means there's a bug...\n", boundaries.size()); 
    for (unsigned i = 0; i < boundaries.size(); i++) printf(" %g ", boundaries[i]); 
    printf("\n"); 
    assert(0);
    return -1; 
  }

  for (unsigned i = 0; i < boundaries.size()/2; i++)
  {
    intersection_points.push_back(std::pair<double,double>(boundaries[2*i], boundaries[2*i+1])); 
  }

  return boundaries.size()/2; 
}

