#include <cmath>
#include <algorithm>

#include "vector.hh"
#include "position.hh"
#include "roughness.hh"
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
#include "spline.h"

Roughness::Roughness(int a){
  file_roughness = "data/roughness_full.txt";
  froughsetting = a;
  if (froughsetting==0){
    gritvalue = -1;
    Ntheta0 = 8;

    amplitude = 11533.5079801;
    x_mean = -70.7895043835*PI/180.;
    y_mean = -117.265622073*PI/180.;
    x_stddev = -62.0590662533*PI/180.;
    y_stddev = 0.278860665027*PI/180.;
    gaustheta = 1.02423010135;
    maxmeaspower = 240.000000;
    fitfuncmax = 1026.090470;
  }
  else if (froughsetting==1){
    gritvalue = 400;
    Ntheta0 = 8;

    amplitude = 1.93882105435;
    x_mean = 2.06659758103*PI/180.;
    y_mean = 2.45567914294*PI/180.;
    x_stddev = 43.1937102786*PI/180.;
    y_stddev = 4.21294249006*PI/180.;
    gaustheta = 0.967793908981;
    maxmeaspower = 2.115675;
    fitfuncmax = 1.936978;
  }
  else if (froughsetting==2){
    gritvalue = 1000;
    Ntheta0 = 9;

    amplitude = 4.98141776637;
    x_mean = 5.82182777326*PI/180.;
    y_mean = 7.98147417219*PI/180.;
    x_stddev = 33.4842444088*PI/180.;
    y_stddev = 2.3019635059*PI/180.;
    gaustheta = 1.00169997698;
    maxmeaspower = 6.014093;
    fitfuncmax = 4.972172;
  }
  else if (froughsetting==3){
    gritvalue = 1500;
    Ntheta0 = 9;

    amplitude = 10.6775743805;
    x_mean = 2.72644088639*PI/180.;
    y_mean = 4.30910349029*PI/180.;
    x_stddev = 30.5405941232*PI/180.;
    y_stddev = 1.19768826165*PI/180.;
    gaustheta = 1.01394332614;
    maxmeaspower = 11.930000;
    fitfuncmax = 10.643788;
  }
  else{
    froughsetting = 1; //default to 400 grit
    gritvalue = 400;
    Ntheta0 = 8;

    amplitude = 1.93882105435;
    x_mean = 2.06659758103*PI/180.;
    y_mean = 2.45567914294*PI/180.;
    x_stddev = 43.1937102786*PI/180.;
    y_stddev = 4.21294249006*PI/180.;
    gaustheta = 0.967793908981;
    maxmeaspower = 2.115675;
    fitfuncmax = 1.936978;
  }

  for (int ii=0; ii<8; ii++){
    if(ii==0){
      corrfactor_thetas.push_back(0.); corrfactor_fresnel.push_back(0.92); corrfactor_loss.push_back(0.91); theta_g2a.push_back(0.);
    }
    if(ii==1){
      corrfactor_thetas.push_back(10.); corrfactor_fresnel.push_back(0.92); corrfactor_loss.push_back(0.91); theta_g2a.push_back(6.6);
    }
    if(ii==2){
      corrfactor_thetas.push_back(20.); corrfactor_fresnel.push_back(0.92); corrfactor_loss.push_back(0.91); theta_g2a.push_back(13.2);
    }
    if(ii==3){
      corrfactor_thetas.push_back(30.); corrfactor_fresnel.push_back(0.91); corrfactor_loss.push_back(0.91); theta_g2a.push_back(19.5);
    }
    if(ii==4){
      corrfactor_thetas.push_back(40.); corrfactor_fresnel.push_back(0.90); corrfactor_loss.push_back(0.90); theta_g2a.push_back(25.4);
    }
    if(ii==5){
      corrfactor_thetas.push_back(50.); corrfactor_fresnel.push_back(0.87); corrfactor_loss.push_back(0.90); theta_g2a.push_back(30.7);
    }
    if(ii==6){
      corrfactor_thetas.push_back(60.); corrfactor_fresnel.push_back(0.80); corrfactor_loss.push_back(0.89); theta_g2a.push_back(35.3);
    }
    if(ii==7){
      corrfactor_thetas.push_back(70.); corrfactor_fresnel.push_back(0.65); corrfactor_loss.push_back(0.89); theta_g2a.push_back(38.8);
    }
    if(ii==8){
      corrfactor_thetas.push_back(90.); corrfactor_fresnel.push_back(0.); corrfactor_loss.push_back(0.); theta_g2a.push_back(41.81031489);
    }
  }


  std::cerr<<"Reading roughness data file:  "<< file_roughness<<std::endl;
  ReadDataFile();
  std::cerr<<"Constructing roughness splines"<<std::endl;
  ConstructSplines();

};


void Roughness::ReadDataFile(void){
  std::ifstream in;
  std::string line = "";

  in.open(file_roughness.c_str());
  // read the header
  getline(in, line);

  double pG, pT0, pT, pP;

  for (int i = 0; i < 648; i++){ // we know there are X lines in the file
    in >> pG >> pT0 >> pT >> pP;
    if (pG != gritvalue)
      continue;
    theta_0.push_back(pT0);  //stores only those lines for specified grit value
    theta.push_back(pT);
    power.push_back(pP);
    //std::cerr<<pT0<<"  "<<pT<<"  "<<pP<<std::endl;
  }

  for(int i=0; i<theta_0.size(); i++){
    if(theta_0[i]!=theta_0[i+1])
      theta_0_unique.push_back(theta_0[i]);
  }

};


void Roughness::ConstructSplines(void){
// make a spline for each theta0 and push it back onto 'spline_theta0'
// read the theta and power vector by Ntheta, make a spline and push, then repeat
  tk::spline *spl_ptr;

  for (int j=0; j<Ntheta0; j++){
    std::vector<double> X;
    std::vector<double> Y;

    
    spl_ptr = new tk::spline;

    if(froughsetting>0){
      for (int i=0; i<19; i++){
        X.push_back( theta[ j*19 + i ] );
        Y.push_back( power[ j*19 + i ] );
        //std::cerr<<theta[ j*19 + i ]<<"  "<<power[ j*19 + i ]<<std::endl;
      }
    }
    else{  //need to treat flat glass separately since some theta_0 have more than 19 measurements
      if(j<2){
        for (int i=0; i<19; i++){
          X.push_back( theta[ j*19 + i ] );
          Y.push_back( power[ j*19 + i ] );
        }
      }
      if(j==2){
        for (int i=0; i<20; i++){
          X.push_back( theta[ j*19 + i ] );
          Y.push_back( power[ j*19 + i ] );
        }
      }
      if(j==3){
        for (int i=0; i<19; i++){
          X.push_back( theta[ 2*19 + 20 + i ] );
          Y.push_back( power[ 2*19 + 20 + i ] );
        }
      }
      if(j==4){
        for (int i=0; i<20; i++){
          X.push_back( theta[ 3*19 + 20 + i ] );
          Y.push_back( power[ 3*19 + 20 + i ] );
        }
      }
      if(j>4){
        for (int i=0; i<19; i++){
          X.push_back( theta[ (j-2)*19 + 2*20 + i ] );
          Y.push_back( power[ (j-2)*19 + 2*20 + i ] );
        }
      }
    }

    // simple check that X is in increasing order (which it's not in the file),
    //  which is required for tk::spline
    if (X[1]<X[0]){
      std::reverse(X.begin(), X.end());
      std::reverse(Y.begin(), Y.end());
    }

    spl_ptr->set_points(X,Y);
    //now push onto stack
    spline_theta0.push_back( spl_ptr );
  }



  // generate the splines for the frensel and loss correction factors
  spl_cf_fresnel_ptr = new tk::spline;
  spl_cf_fresnel_ptr->set_points(corrfactor_thetas, corrfactor_fresnel);
  //
  spl_cf_loss_ptr = new tk::spline;
  spl_cf_loss_ptr->set_points(corrfactor_thetas, corrfactor_loss);


  //generate splines for converting incidence angle for air->glass interface of measurements, and back
  // these values taken from Table 5 ELOG #077
  spl_ag2ga_ptr = new tk::spline;
  spl_ag2ga_ptr->set_points(corrfactor_thetas, theta_g2a);

  spl_ga2ag_ptr = new tk::spline;
  spl_ga2ag_ptr->set_points(theta_g2a, corrfactor_thetas);

};


double Roughness::GetLaserPower(){
  return 580.; // \muW, taken from ELOG 077
};


double Roughness::InterpolatePowerValue(double T0, double T){
  // T0 [degrees] is the incident angle of the ray-in-ice with respect the the surface normal pointed into the air
  // theta [degrees] is the exiting angle from the surface towards the balloon, measured with respect to the surface normal pointing into the air
  double p;

  //++++++++++++++++++++++++++++++++++++++
  // use this for 1+1 interpolation
  // procedure: for theta, go through spline vector and get that value. write to a vector, then spline that and get value for the T0
  // use Roughness::T0 for X proxy
  //std::vector<double> spl_values; // Y proxy
  //for (int j=0; j<Ntheta0; j++){
  //  spl_values.push_back(  (*(spline_theta0[j]))(T)  );
  //}
  //tk::spline s;
  //s.set_points(theta_0_unique, spl_values);
  //p = s(T0);

  //++++++++++++++++++++++++++++++++++++++
  // use this for analytic evaluation
  // get the fit value and renormalize to the maximum measured power
  p = (evaluate2dGaussian(T0, T) * maxmeaspower / fitfuncmax);
  return p;
};


double Roughness::evaluate2dGaussian(double T0, double T){
  double x = T0*PI/180.;
  double y = T*PI/180.;

  double a = cos(gaustheta)*cos(gaustheta)/2./x_stddev/x_stddev + sin(gaustheta)*sin(gaustheta)/2./y_stddev/y_stddev;
  
  double b = sin(2.*gaustheta)/2./x_stddev/x_stddev - sin(2.*gaustheta)/2./y_stddev/y_stddev;
  
  double c = sin(gaustheta)*sin(gaustheta)/2./x_stddev/x_stddev + cos(gaustheta)*cos(gaustheta)/2./y_stddev/y_stddev;

  double expvalue = -a*(x-x_mean)*(x-x_mean)-b*(x-x_mean)*(y-y_mean)-c*(y-y_mean)*(y-y_mean);

  return (amplitude * exp(expvalue));
};


double Roughness::GetFresnelCorrectionFactor(double T0){
   return (*spl_cf_fresnel_ptr)(T0);
};


double Roughness::GetLossCorrectionFactor(double T0){
   return (*spl_cf_loss_ptr)(T0);
};


double Roughness::ConvertTheta0AirGlass_to_GlassAir(double T0){
  return (*spl_ag2ga_ptr)(T0);
};


double Roughness::ConvertTheta0GlassAir_to_AirGlass(double T1){
  return (*spl_ga2ag_ptr)(T1);
};
