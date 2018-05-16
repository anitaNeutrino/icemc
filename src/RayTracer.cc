#include "TF1.h"
#include "TRandom3.h"

#include "Constants.h"
#include "TRandom3.h"
#include "Settings.h"
#include "earthmodel.hh"
#include "icemodel.hh"
#include "AskaryanFreqsGenerator.h"
#include "vector.hh"
#include "position.hh"
#include "anita.hh"
#include "RayTracer.h"
#include <cmath>


icemc::RayTracer::RayTracer() {
  sum_slopeyness=0.;
  xaxis=Vector(1.,0.,0.);
  yaxis=Vector(0.,1.,0.);
  zaxis=Vector(0.,0.,1.);
  slopeyx=0.;
  slopeyy=0.;
  slopeyz=0.; // these are a measure of how much the surface is sloped in the x,y and z directions
}

void icemc::RayTracer::PrintAnglesofIncidence() const {

  std::cout << "angle of incidence (firn-air) is " << nrf_iceside[3].Angle(nsurf_rfexit)*constants::DEGRAD << "\n";
  std::cout << "angle of incidence (ice-firn) is " << nrf_iceside[4].Angle(nsurf_rfexit)*constants::DEGRAD << "\n";
}



void icemc::RayTracer::Initialize() {
  
  for (int i=0;i<3;i++) {
    n_exit2bn[i]=Vector(0.,0.,0.); // normal vector in direction of exit point to balloon - 5 iterations, 3 directions for eac
    nrf_iceside[i]=Vector(0.,0.,0.);  // direction of rf [tries][3d]
    rfexit_db[i]=Vector(0.,0.,0.);
    rfexit[i]=Vector(0.,0.,0.); // position where the rf exits the ice- 5 iterations, 3 dimensions eac
    
    for (int j=0;j<Anita::NPHI_MAX;j++) {
      for (int k=0;k<Anita::NLAYERS_MAX;k++) {
	n_exit2bn_eachboresight[i][k][j]=Vector(0.,0.,0.); // normal vector in direction of exit point to each antenna boresight - 5 iterations
	nrf_iceside_eachboresight[i][k][j]=Vector(0.,0.,0.);  // direction of rf [tries][3d]
  	
	rfexit_eachboresight[i][k][j]=Vector(0.,0.,0.);
      }
    }
    
  }
  nsurf_rfexit=Vector(0.,0.,0.); // normal of the surface at the place where the rf leaves
  nsurf_rfexit_db=Vector(0.,0.,0.);

}


void icemc::RayTracer::GetRFExit(const Settings *settings1, Anita *anita1, int whichray, Position posnu, Position posnu_down, Position r_bn,
				 Position r_boresights[Anita::NLAYERS_MAX][Anita::NPHI_MAX], int whichtry, const IceModel *antarctica){



  if (whichray==0)
    WhereDoesItLeave(posnu,nrf_iceside[2*whichtry],antarctica, // inputs
		     rfexit[whichtry]); // output
  
  
  //******wufan******
  if (whichray==1) // reflected rays
    WhereDoesItLeave(posnu_down,nrf_iceside[2*whichtry],antarctica,rfexit[whichtry]);  //from mirror point position and the direction of signals to the ice 
  //to find the exit point at the surface of the Earth.wufan 
  
  n_exit2bn[whichtry] = (r_bn - rfexit[whichtry]).Unit();
  
  if (settings1->BORESIGHTS) { // now find rfexit and n_exit2bn for each boresight, still first iteration
    //std::cout << "first iteration.\n";
    for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
      for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
        WhereDoesItLeave(posnu,nrf_iceside_eachboresight[2*whichtry][ilayer][ifold],antarctica,
			 rfexit_eachboresight[whichtry][ilayer][ifold]);
        n_exit2bn_eachboresight[whichtry][ilayer][ifold] = (r_boresights[ilayer][ifold] - rfexit_eachboresight[whichtry][ilayer][ifold]).Unit(); 
      }
    }
  } // end if we are doing this for each boresight
  
  if (settings1->SLAC) {
    // ray comes out a little earlier because of the slope of the surface.
    // use law of cosines the get how much "distance" should be cut short
    
    double x=sin(settings1->SLACSLOPE*constants::RADDEG)*(settings1->SLACICELENGTH/2.+rfexit[0].Distance(rfexit[whichtry]))/sin(constants::PI/2.+acos(nrf_iceside[2*whichtry].Dot(nsurf_rfexit)));
    
    rfexit[whichtry]-=x*nrf_iceside[2*whichtry];
    
    if (settings1->BORESIGHTS) {
      for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
        for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
          
          x=sin(settings1->SLACSLOPE*constants::RADDEG)*(settings1->SLACICELENGTH/2.+rfexit_eachboresight[0][ilayer][ifold].Distance(rfexit_eachboresight[whichtry][ilayer][ifold]))/sin(constants::PI/2.+acos(nrf_iceside_eachboresight[2*whichtry][ilayer][ifold].Dot(nsurf_rfexit)));
          
          rfexit_eachboresight[whichtry][ilayer][ifold]-=x*nrf_iceside_eachboresight[2*whichtry][ilayer][ifold];
          
          n_exit2bn_eachboresight[whichtry][ilayer][ifold] = (r_boresights[ilayer][ifold] - rfexit_eachboresight[whichtry][ilayer][ifold]).Unit();
        } // end loop over antennas on a layer
      } // end loop over layers of the payload
    } // end if we're keeping track of all antenna boresights
  } // end if we're modeling the slac run
}

//###########
// icemc::RayTracer::WhereDoesItLeave() is defined in ray.hh since it is a statis member function // MS 2/1/2017


int icemc::RayTracer::RandomizeSurface(const Settings *settings1, Position rfexit_temp, Vector posnu, const IceModel *antarctica, double &slopeyangle, int whichtry){

  double howmuch=settings1->SLOPEYSIZE;
  Position nsurf_rfexit_temp;
  
  // rotate the surface normal according to the local slope.
  if (!settings1->SLAC){
    nsurf_rfexit_temp = antarctica->GetSurfaceNormal(rfexit_temp); // find the normal to the surface taking into account the tilt from the differential heights between neighboring bins
  }
  else if (settings1->SLAC) { // if we are simulating the slac test, rotate the surface by 10 degrees away from the south pole
    Vector zaxis(0.,0.,1.);
    nsurf_rfexit_temp=(rfexit_temp.Unit()).Rotate(-settings1->SLACSLOPE*constants::RADDEG,posnu.Cross(zaxis));
  }
  
  Position nsurf_rfexit_temp_copy=nsurf_rfexit_temp;
  // tilt local surface based on slopeyness.
  if (settings1->SLOPEY) {
    
    // randomizing surface direction
    
    // 0.10=5.4, 0.2=7.4 deg mean
    
    // 02/06/07: 0.012=0.84, 0.10=6.78,3.76 deg mean
    
    double slopeyness=0;
    if (whichtry==0) { // only reset the surface slopeyness for the first try an then repeat the same for each subsequent try
      slopeyx=howmuch*gRandom->Gaus();
      slopeyy=howmuch*gRandom->Gaus();
      slopeyz=howmuch*gRandom->Gaus();
    }
    Vector ntemp2 = nsurf_rfexit_temp + slopeyx * xaxis
      + slopeyy * yaxis
      + slopeyz * zaxis;
    
      
    ntemp2 = ntemp2 / ntemp2.Mag();
    
    double rtemp= ntemp2.Dot(nsurf_rfexit_temp);
    
    if (rtemp<=1) {
      slopeyness=acos(rtemp);
      sum_slopeyness+=slopeyness;
      nsurf_rfexit_temp = ntemp2;
    }//if
    else
      slopeyness=-999;
    
    if (settings1->DEBUG) {
      rtemp=nsurf_rfexit_temp.Mag(); 
      if (fabs(rtemp-1.0)>=0.01) {
	std::cout << "error: nx,ny,nz_surf="<<nsurf_rfexit_temp<<std::endl;
	return 0;
      }//if
    }//if
    
    
    // to avoid FP errors, fudge the rare case of going out the poles
    if (fabs(nsurf_rfexit_temp[2])>=0.99999) {  // a hack
      nsurf_rfexit_temp = Vector(sqrt(1.-0.99999*0.99999-nsurf_rfexit_temp[1]*nsurf_rfexit_temp[1]), nsurf_rfexit_temp[1], 0.99999);
    }//if
  }  
  slopeyangle=constants::DEGRAD*acos(nsurf_rfexit_temp.Dot(nsurf_rfexit_temp_copy));
  nsurf_rfexit=nsurf_rfexit_temp;  
  return 1;
  
}//RandomizeSurface


// int icemc::RayTracer::GetSurfaceNormal(IceModel *antarctica,Vector posnu,Position *rfexit) {
int icemc::RayTracer::GetSurfaceNormal(const Settings *settings1,const IceModel *antarctica,Vector posnu,double &slopeyangle,int whichtry){
      
  Position rfexit_temp;

  rfexit_temp=rfexit[whichtry];

 
  if (!RandomizeSurface(settings1,rfexit_temp,posnu,antarctica,slopeyangle,whichtry)) {

    return 0;  //can happen only in debug mode
  }






  return 1;
}



int icemc::RayTracer::TraceRay(const Settings *settings1,Anita *anita1,int iter,double n_depth) { // iter is which iteration (1 or 2)

   
  // use snell's law to get the first guess at the 
  // direction of the rf as it leaves ice surface.
  // 0th guess was simply radially outward from interaction position
  // this now takes into account balloon position and surface normal.
    
  if (settings1->FIRN) {

    
    if (!GetRayIceSide(n_exit2bn[iter-1],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,constants::NFIRN,
		       nrf_iceside[2*iter-1])) { // nrf_iceside[1] is the rf direction in the firn

      return 0; // reject if TIR.
    } 
    // This could be throwing away events that the final guess would have kept
    //  Need to look into this and update in the future
      
    // use snell's law again to get direction of rf ray emitted from
    // interaction point.  Physically, it is a continuous bending
    // through the firn but it turns out it obeys Snell's law (Seckel)
            
    if (!GetRayIceSide(nrf_iceside[2*iter-1],nsurf_rfexit,constants::NFIRN,n_depth,
		       nrf_iceside[2*iter])) { // nrf_iceside[2] is the rf direction in the ice

      return 0;   // note to self:  need to check whether rays are totally internally reflected within ice   
    }




    // find where the refracted ray leaves.
    // rejected if the rf leaves beyond the boundary of the 
    // continent. Should never happen.
      
    if (settings1->BORESIGHTS) {
	
      for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
	for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
	    
	  if (!GetRayIceSide(n_exit2bn_eachboresight[iter-1][ilayer][ifold],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,constants::NFIRN,
			     nrf_iceside_eachboresight[2*iter-1][ilayer][ifold])) // nrf_iceside[1] is the rf direction in the firn
	    return 0; // reject if TIR.  


	  //	    std::cout << "ITER " << iter-1 << " " << (n_exit2bn_eachboresight[iter-1][ilayer][ifold]) << std::endl;

	  // This could be throwing away events that the final guess would have kept
	  //  Need to look into this and update in the future
	    
	  // use snell's law again to get direction of rf ray emitted from
	  // interaction point.  Physically, it is a continuous bending
	  // through the firn but it turns out it obeys Snell's law (Seckel)
	    
	    
	  if (!GetRayIceSide(nrf_iceside_eachboresight[2*iter-1][ilayer][ifold],nsurf_rfexit,constants::NFIRN,n_depth,
			     nrf_iceside_eachboresight[2*iter][ilayer][ifold])) // nrf_iceside[2] is the rf direction in the ice
	    return 0;   // note to self:  need to check whether rays are totally internally reflected within ice   
	    
	  // find where the refracted ray leaves.
	  // rejected if the rf leaves beyond the boundary of the 
	  // continent. Should never happen.
	    
	    
	} // end phi foldings
      } // end layers
    } // end if doing this for all boresights
  } // end if we are modeling the firn
    
  else { // no firn

    if (!GetRayIceSide(n_exit2bn[iter-1],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,AskaryanFreqsGenerator::NICE,
		       nrf_iceside[2*iter-1])) // nrf_iceside[1] is the rf direction in the ice
      return 0; // reject if TIR.  
    nrf_iceside[2*iter]=nrf_iceside[2*iter-1]; // no firn so the next element is the same
      
    if (settings1->BORESIGHTS) {

      for(int ilayer=0;ilayer<settings1->NLAYERS;ilayer++) {
	for(int ifold=0;ifold<anita1->NRX_PHI[ilayer];ifold++) {
	    
	  if (!GetRayIceSide(n_exit2bn_eachboresight[iter-1][ilayer][ifold],nsurf_rfexit,AskaryanFreqsGenerator::N_AIR,AskaryanFreqsGenerator::NICE,
			     nrf_iceside_eachboresight[2*iter-1][ilayer][ifold])) // nrf_iceside[1] is the rf direction in the ice
	    return 0; // reject if TIR.  
	  nrf_iceside_eachboresight[2*iter][ilayer][ifold]=nrf_iceside_eachboresight[2*iter-1][ilayer][ifold]; // no firn so the next element is the same
	    
	    
	} // number of foldings in phi
      } // number of layers
    } // end if doing this for all boresights
      
      
      
  } // end if we are not modeling the firn
    
    

  return 1;
}


int icemc::RayTracer::GetRayIceSide(const Vector &n_exit2bn,
				    const Vector &nsurf_rfexit,
				    double nexit,
				    double nenter, 
				    Vector &nrf2_iceside) const {

  // this function performs snell's law in three dimensions

  double costh=0;

  double NRATIO=nexit/nenter;

  //  std::cout << " INPUT " << n_exit2bn << " " << nsurf_rfexit << " " << nexit << " " << nenter << std::endl;
  
  costh=(n_exit2bn.Dot(nsurf_rfexit))/(n_exit2bn.Mag() * nsurf_rfexit.Mag()); // cos(theta) of the transmission angle

  if (costh<0) {
    //    std::cout << "returning 0.  inu is " << inu << "\n";
    return 0;
  
  }
  double sinth=sqrt(1 - costh*costh);
  
  //  double sinth_i=nexit*sinth/nenter;
  
  //  double angle=asin(sinth)-asin(sinth_i);
  
  
  // std::cout << "sinth, sinth_i are " << sinth << " " << sinth_i << "\n";
  // std::cout << "th, th_i are " << asin(sinth)*DEGRAD << " " << asin(sinth_i)*DEGRAD << "\n";
  // std::cout << "angle is " << angle*DEGRAD << "\n";
  
  //  nrf2_iceside=n_exit2rx.Rotate(angle,n_exit2rx.Cross(nsurf_rfexit));


  double factor=NRATIO*costh-sqrt(1-(NRATIO*sinth*NRATIO*sinth));


  nrf2_iceside = -factor*nsurf_rfexit + NRATIO*n_exit2bn;
  nrf2_iceside = nrf2_iceside.Unit(); // normalize
  
  
    
  

  return 1;
}//GetRayIceSide



// static
int icemc::RayTracer::WhereDoesItLeave(const Position &posnu, const Vector &ntemp, const IceModel *antarctica, Position &r_out){
  
  double distance=0;
  double posnu_length=posnu.Mag(); // distance from center of earth to interaction
  double lon,lat;//,lon_old,lat_old; //latitude, longitude indices for 1st and 2nd iteration
  lon = posnu.Lon(); // what latitude, longitude does interaction occur at
  lat = posnu.Lat();
  // lon_old=lon; // save this longitude and latitude so we can refer to it later
  // lat_old=lat;

  // use law of cosines to get distance from interaction to exit point for the ray
  // need to solve for that distance using the quadratic formula

  // angle between posnu and ntemp vector for law of cosines.
  double costheta=-1*(posnu.Dot(ntemp))/posnu_length;

  // a,b,c for quadratic formula, needed to solve for 
  double a=1;
  double b=-1*2*posnu_length*costheta;
  double c=posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2);

  if (b*b-4*a*c<0.) {
    return 0;
  }

  // use the "+" solution because the other one is where the ray is headed downward toward the rock
  distance=(-1*b+sqrt(b*b-4*a*c))/2;
				
  // now here is the exit point for the ray
  r_out = posnu + distance*ntemp;
				
  lon = r_out.Lon(); // latitude and longitude of exit point
  lat = r_out.Lat();
				
  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
  // sometimes though the new surface is lower than the one over posnu which causes a problem.
  if (b*b-4*a*c<0.) {
    //std::cout << "inu is " << inu << "\n";  
    // try halving the distance
    distance=distance/2.;
    //std::cout << "bad.  distance 1/2 is " << distance << "\n";
    r_out = posnu + distance*ntemp;
    lon = r_out.Lon(); // latitude and longitude of exit point
    lat = r_out.Lat();
    c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
    if (b*b-4*a*c<0.) { // if we still have the problem back up more
      distance=distance/2.; // now we are at 1/4 the distance
      //std::cout << "bad.  distance 1/4 is " << distance << "\n";
      r_out = posnu + distance*ntemp;
      lon = r_out.Lon(); // latitude and longitude of exit point
      lat = r_out.Lat();
      c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
      if (b*b-4*a*c<0.) { // the problem is less then 1/4 of the way in
	distance=distance/2.; // now we are at 1/8 the distance
	//std::cout << "bad.  distance 1/8 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4*a*c<0.) {
	  // still have the problem so just make the distance 0
	  distance=0.; // now we are at 1/8 the distance
	  //std::cout << "bad.  distance is " << distance << "\n";
	  lon = posnu.Lon(); // latitude and longitude of exit point
	  lat = posnu.Lat();
	  r_out=antarctica->Surface(lon,lat)/posnu.Mag()*posnu;
	}
      } // now we are at 1/8 the distance
      else {// if this surface is ok problem is between 1/4 and 1/2
	distance=distance*1.5; // now we are at 3/8 the distance
	//	std::cout << "good.  distance 3/8 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4.*a*c<0.) {
	  distance=distance*2./3.; // go back to 1/4
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	  //std::cout << "good at distance 1/4 is " << distance << "\n";
	}
      } // now we are at 3/8 the distance
		
		
    } // now we are at 1/4 the distance
    else { // if this surface at 1/2 distance is ok see if we can go a little further
      distance=distance*1.5; // now we are at 3/4 the distance
      //std::cout << "good.  distance 3/4 is " << distance << "\n";
      r_out = posnu + distance*ntemp;
      lon = r_out.Lon(); // latitude and longitude of exit point
      lat = r_out.Lat();
      c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
      if (b*b-4*a*c<0.) { // the problem is between 1/2 and 3/4 of the way in
	distance=distance*5./6.; // now we are at 5/8 the distance
	//std::cout << "bad.  distance 5/8 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4*a*c<0.) {
	  distance=distance*4./5.;
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	  //std::cout << "good at distance 3/4 is " << distance << "\n";
	}
      } // now we are at 1/8 the distance
      else {// if this surface is ok problem is between 1/4 and 1/2
	distance=distance*7./6.; // now we are at 7/8 the distance
	//std::cout << "good.  distance 7/8 is " << distance << "\n";
	r_out = posnu + distance*ntemp;
	lon = r_out.Lon(); // latitude and longitude of exit point
	lat = r_out.Lat();
	c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	if (b*b-4*a*c<0) {
	  // now found the problem so go back to 3/4 distance
	  distance=distance*6./7.;
	  //std::cout << "good at  distance 3/4 is " << distance << "\n";
	  r_out = posnu + distance*ntemp;
	  lon = r_out.Lon(); // latitude and longitude of exit point
	  lat = r_out.Lat();
	  c = posnu_length*posnu_length - pow(antarctica->Surface(lon,lat),2); // redo the law of cosines
	}
      } // now we are at 3/8 the distance
    } // now we are at 3/4 distance
  } // if exit point we initially found was not ok
  else {
    distance=(-1*b+sqrt(b*b-4*a*c))/2; // and quadratic formula
    r_out = posnu + distance*ntemp;
  }
  return 1;
}




