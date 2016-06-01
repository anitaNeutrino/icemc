#include "vector.hh"
#include "position.hh"
#include "roughness.hh"
#include <cmath>
#include "Constants.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "signal.hh"
#include "Settings.h"
#include "Primaries.h"
#include "anita.hh"
#include "balloon.hh"
#include "earthmodel.hh"
#include "icemodel.hh"


#include "ray.hh"

Roughness::Roughness() {
    
}





void Roughness::GetBalloonLocation(Interaction *interaction1,Ray *ray1,Balloon *bn1,IceModel *antarctica) {
    // brian enter function to calculate balloon position on your map.
    // use interaction1->posnu // location of neutrino interaction
    // coordinate system:  +z="up" at the south pole
    // bn1->r_bn
    // nnu
    // ray1->nsurf_rfexit
    
    
    // brian enter function to calculate balloon position on your map.
    // use interaction1->posnu // location of neutrino interaction
    // coordinate system:  +z="up" at the south pole
    // bn1->r_bn
    // nnu
    
    
    // balloonvector = balloonvector - nuvector;//change origin to the nuetrino interaction point
    
    const Vector nuvector = interaction1->nnu;
    // double interactiondepth = nuvector[2];//NOT CORRECT! need depth BELOW the ice. this is height above center of earth.
    
    Vector zcoordvector = ray1->nsurf_rfexit;
    zcoordvector=zcoordvector.Unit();
    
    // double thetainc =acos(zcoordvector.Dot(nuvector))*180/PI;
    //nsurf_rfexit is z direction for new coordinate system. Need to make sure the n-vector is in x-z plane.
    
    Vector xcoordvector = nuvector-(zcoordvector.Dot(nuvector))*zcoordvector;//xcoordvector is such that nnu lies in the x-z plane
    xcoordvector = xcoordvector.Unit();
    
    const Vector ycoordvector = zcoordvector.Cross(xcoordvector);//Need this for ChangeCoord. 
    
    
    Vector origin_brian_tmp;
    if (interaction1->nnu.Dot(zcoordvector)>0) // up	going
	origin_brian_tmp=interaction1->nuexit; // the origin is the neutrino exit point	
    else {	   
	Vector nnu_flipped=interaction1->nnu;
	nnu_flipped=nnu_flipped-2.*nnu_flipped.Dot(zcoordvector)*zcoordvector; // take it's upgoing reflection from surface
	
	Position nuexit_flipped;
	if (Ray::WhereDoesItLeave(interaction1->posnu,nnu_flipped,antarctica,nuexit_flipped))
	    origin_brian_tmp=nuexit_flipped;
	
	
    }
    Vector r_bn_tmp=bn1->r_bn-origin_brian_tmp;
    r_bn_tmp=r_bn_tmp.ChangeCoord(xcoordvector,ycoordvector);//change coordinates
    
    
    balloondist =r_bn_tmp.Mag();//this is above center of earth, if i understand correctly. Need above the surface of the earth. 
    balloonphi = r_bn_tmp.Phi(); //phi position of the balloon
    if (balloonphi>PI)
	balloonphi=balloonphi-2*PI;
    balloontheta = r_bn_tmp.Theta();// costheta position of the baloon
    // get this by dotting ray1->nsurf_rfexit with nnu?     
    // double thetainc = acos(interaction1->nnu[2])*180/PI; //nnu is unit vector; cos(thetainc) = z/r
    balloontheta = PI-balloontheta;//walter.cc uses a pos z as down. this corrects for that.
    
    // define a coordinate system with ray1->nsurf_rfexit defining +z
    // nnu direction defines the x-z plane
    // find balloon position in that coordinate system
    //to get the values from walter.cc we need : E_shower, correlation length, rms height and the em_frac and had_frac. the last
    // two are so we can multiply the number from sky maps by the correct frac and then add the em and hadronic portion together
    // to get the total.
    
}









double Roughness::GetPokey(double incident_angle,double transmitted_angle,
			   double emfrac,double hadfrac,double deltheta_em_max,double deltheta_had_max) {
    
    // first get normalization
    double combined_deltheta=GetCombinedDeltheta(emfrac,hadfrac,deltheta_em_max,deltheta_had_max);
    double deltheta_prime=sqrt(rough_sigma*rough_sigma+combined_deltheta*combined_deltheta/2.);
    double A,a,b,c,d;
    
    //A=(a-b)/(c-d)
    //a=cos(CHANGLE_ICE-combined_deltheta)
    if (Signal::CHANGLE_ICE-combined_deltheta<0)
	a=1.;
    else
	a=cos(Signal::CHANGLE_ICE-combined_deltheta);
    
    if (Signal::CHANGLE_ICE+combined_deltheta>PI)
	b=-1.;
    else 
	b=cos(Signal::CHANGLE_ICE+combined_deltheta);
    
    if (Signal::CHANGLE_ICE-deltheta_prime<0)
	c=1.;
    else
	c=cos(Signal::CHANGLE_ICE-deltheta_prime);
    
    if (Signal::CHANGLE_ICE+deltheta_prime>PI)
	d=-1.;
    else
	d=cos(Signal::CHANGLE_ICE+deltheta_prime);
    
    
    A=(a-b)/(c-d);
    
    
    
    double t=A*fpokey2->Eval(transmitted_angle,incident_angle);
    
    
    return t;
}

double Roughness::GetSlappy(double incident_angle,double transmitted_angle,
			    double emfrac,double hadfrac,double deltheta_em_max,double deltheta_had_max) {
    
    
    // first get normalization
    double combined_deltheta=(emfrac*deltheta_em_max+hadfrac*deltheta_had_max)/(emfrac+hadfrac);
    double deltheta_prime=sqrt(rough_sigma*rough_sigma+combined_deltheta*combined_deltheta);
    
    double A,a,b,c,d;
    
    //A=(a-b)/(c-d)
    //a=cos(CHANGLE_ICE-combined_deltheta)
    if (Signal::CHANGLE_ICE-combined_deltheta<0)
	a=1.;
    else
	a=cos(Signal::CHANGLE_ICE-combined_deltheta);
    
    if (Signal::CHANGLE_ICE+combined_deltheta>PI)
	b=-1.;
    else 
	b=cos(Signal::CHANGLE_ICE+combined_deltheta);
    
    if (Signal::CHANGLE_ICE-deltheta_prime<0)
	c=1.;
    else
	c=cos(Signal::CHANGLE_ICE-deltheta_prime);
    
    if (Signal::CHANGLE_ICE+deltheta_prime>PI)
	d=-1.;
    else
	d=cos(Signal::CHANGLE_ICE+deltheta_prime);
    
    
    A=(a-b)/(c-d);
    
    
    double t=A*fslappy2->Eval(transmitted_angle,incident_angle);
    //  double t= h2_slappy_smeared->GetBinContent(h2_slappy_smeared->FindBin(incident_angle,transmitted_angle));
    
    return t;
    
}


void Roughness::GetExitforRoughness(Settings *settings1,IceModel *antarctica,double emfrac,double hadfrac,double deltheta_em_max,double deltheta_had_max,
				    Ray *ray1,
				    Vector &nnu,Vector &r_bn,
				    Vector &posnu) {
    
    // want to construct the ray that is not bent at the ice-air interface, and hits the balloon.  That is the "reference ray"
    // this function just draws a straight line from the interaction to the balloon.
    
    Vector referenceray_iceside_specular; // direction of reference ray on ice side
    Vector n_referenceray_exit2bn; // in line with incident ray
    Position referenceray_rfexit; // position where the reference ray exits the ice
    
    // in the absence of the firn, the ray goes straight from the interaction position to the balloon.
    // here is the normalized vector from exit point to balloon
    n_referenceray_exit2bn=(r_bn-posnu).Unit();
    
    // here is the normalized vector from interaction to exit point (here, same as above)
    
    Vector totherightvect,tempvect,rotationvector;
    // double angle,viewingangle_rough,viewingangle_norough,rotationangle;
    double viewingangle_rough,viewingangle_norough,rotationangle;
    if (settings1->ROUGHNESS==1)
	referenceray_iceside_specular=(r_bn-posnu).Unit();
    else if (settings1->ROUGHNESS==2) {
	// find the vector pointing to "right" from the neutrino's perspective
	// where "up" would be (nnu x rfexit) x nnu
	totherightvect=(nnu.Cross(ray1->rfexit[0])).Unit();
	
	// this tempvect will be aligned with totheright when nrf_iceside is a ray emitted at the top of the cone
	tempvect=(nnu.Cross(ray1->nrf_iceside[3])).Unit();
	
	// angle=acos(tempvect.Dot(totherightvect)); // this tells you what "phi" one the cone the ray is at.
	
	// need to rotate the ray about this vector to change its viewing angle but not its "phi" around cone
	rotationvector=nnu.Cross(ray1->nrf_iceside[3]);
	
	viewingangle_norough=acos(nnu.Dot(ray1->nrf_iceside[3]));
	
	//    viewingangle_rough=Signal::CHANGLE_ICE-(Signal::CHANGLE_ICE-viewingangle_norough)*(Signal::CHANGLE_ICE-rough_sigma)/Signal::CHANGLE_ICE;
	
	
	double combined_deltheta=GetCombinedDeltheta(emfrac,hadfrac,deltheta_em_max,deltheta_had_max);
	
	// need a factor of 2 in front of combined_deltheta because
	viewingangle_rough=(viewingangle_norough*(combined_deltheta*combined_deltheta/2.)+Signal::CHANGLE_ICE*(rough_sigma*rough_sigma))/(rough_sigma*rough_sigma+combined_deltheta*combined_deltheta/2.);
	
	rotationangle=viewingangle_rough-viewingangle_norough;
	referenceray_iceside_specular=ray1->nrf_iceside[3].Rotate(rotationangle,rotationvector);
	tempvect=(nnu.Cross(referenceray_iceside_specular)).Unit();
	// angle=acos(tempvect.Dot(totherightvect)); // this tells you what "phi" one the cone the ray is at.
	ray1->nrf_iceside[4]=referenceray_iceside_specular;
	ray1->nrf_iceside[3]=referenceray_iceside_specular;
    }
    
    // here is the position of the exit point
    Ray::WhereDoesItLeave(posnu,referenceray_iceside_specular,antarctica,
			  referenceray_rfexit);
    
    // the exit point, specular ray and the ray from the exit to the balloon are here defined by the "reference ray"
    ray1->rfexit[1]=referenceray_rfexit;
    nrf_iceside_specular=referenceray_iceside_specular; // this is a member of the roughness class
    ray1->n_exit2bn[2]=r_bn-referenceray_rfexit;
    
}

double Roughness::GetCombinedDeltheta(double emfrac,double hadfrac,double deltheta_em_max,double deltheta_had_max) {
    
    return (emfrac*deltheta_em_max+hadfrac*deltheta_had_max)/(emfrac+hadfrac);
    
}
