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
  file_roughness = "data/roughness_full_dPdtheta.txt";
  froughsetting = a;
  if (froughsetting==0){
    gritvalue = -1;
    Ntheta0 = 8;

    amplitude = 1112.67009843;
    x_mean = -20.0115394201*PI/180.;
    y_mean = -33.7371829382*PI/180.;
    x_stddev = -34.497529739*PI/180.;
    y_stddev = 0.231052092516*PI/180.;
    gaustheta = 1.02432941913;
    maxmeaspower = 99.006880;
    fitfuncmax = 560.591438;
  }
  else if (froughsetting==1){
    gritvalue = 400;
    Ntheta0 = 8;

    amplitude = 0.798064953724;
    x_mean = -4.84261442966*PI/180.;
    y_mean = -7.05473728178*PI/180.;
    x_stddev = 40.3467981404*PI/180.;
    y_stddev = 3.49683983104*PI/180.;
    gaustheta = 0.928258773739;
    maxmeaspower = 0.874336;
    fitfuncmax = 0.781501;
  }
  else if (froughsetting==2){
    gritvalue = 1000;
    Ntheta0 = 9;

    amplitude = 1.99315873416;
    x_mean = -1.43758569242*PI/180.;
    y_mean = -2.9617313585*PI/180.;
    x_stddev = 32.2255020021*PI/180.;
    y_stddev = 2.05472924304*PI/180.;
    gaustheta = 0.97544239723;
    maxmeaspower = 2.087395;
    fitfuncmax = 1.985588;
  }
  else if (froughsetting==3){
    gritvalue = 1500;
    Ntheta0 = 9;

    amplitude = 6.39968975801;
    x_mean = -18.3706965926*PI/180.;
    y_mean = -29.1330198008*PI/180.;
    x_stddev = 41.5606379916*PI/180.;
    y_stddev = 1.09720021655*PI/180.;
    gaustheta = 1.00846567317;
    maxmeaspower = 4.921467;
    fitfuncmax = 4.539596;
  }
  else{
    froughsetting = 1; //default to 400 grit
    gritvalue = 400;
    Ntheta0 = 8;

    amplitude = 0.798064953724;
    x_mean = -4.84261442966*PI/180.;
    y_mean = -7.05473728178*PI/180.;
    x_stddev = 40.3467981404*PI/180.;
    y_stddev = 3.49683983104*PI/180.;
    gaustheta = 0.928258773739;
    maxmeaspower = 0.874336;
    fitfuncmax = 0.781501;
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

  //++++++++++++++++++++++++++++++++++++++
  // use this for 1+1 interpolation
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
