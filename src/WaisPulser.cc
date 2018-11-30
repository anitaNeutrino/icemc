#include "WaisPulser.h"
#include "EnvironmentVariable.h"
#include "Interaction.h"
#include "FTPair.h"
#include "WorldModel.h"
#include "Report.h"
#include "Settings.h"

#ifdef ANITA_UTIL_EXISTS
#include "AnitaConventions.h"
#endif

icemc::WaisPulser::WaisPulser(const Settings* settings, std::shared_ptr<const WorldModel> worldModel)
  : fSettings(settings)
{

  
  // ANITA-3 WAIS location by default
  // If ANITA_UTIL_EXISTS then the appropriate location will be loaded
  Double_t latWAIS  = - ( 79 + (27.93728/60) ) ;
  Double_t lonWAIS  = - (112 + ( 6.74974/60) ) ;
  Double_t altWAIS  = 1775.68;
  
#ifdef ANITA_UTIL_EXISTS
  latWAIS = AnitaLocations::getWaisLatitude()  ;
  lonWAIS = AnitaLocations::getWaisLongitude() ;
  altWAIS = AnitaLocations::getWaisAltitude()  ;
#endif
  
  fWaisPulserPosition.SetLonLatAlt(lonWAIS,  latWAIS, altWAIS);

  double newAlt = worldModel->SurfaceAboveGeoid(fWaisPulserPosition);
  std::cout << "WorldModel gives surface at WAIS at " << newAlt << ", old was " << fWaisPulserPosition.Altitude() << std::endl;
  ///@todo check, do we want to put this slightly below/above the surface?
  fWaisPulserPosition.SetAltitude(newAlt + 0.01); // put it slightly above the surface?

  fSurfaceNormal = worldModel->GetSurfaceNormal(fWaisPulserPosition);  
  
  readModel();
  
}

icemc::WaisPulser::~WaisPulser(){
}



icemc::Interaction icemc::WaisPulser::generateInteraction(const icemc::Neutrino& n, const Geoid::Position& detector){  
  Interaction interaction;
  interaction.position = fWaisPulserPosition;
  return interaction;
  
}


void icemc::WaisPulser::readModel(){

  if(fVoltsPerMeter.size()==0){

    std::cout << "Reading pulser model: ";
    //SAW sample rate on this model is not quite right.
    //string fname = ICEMC_SRC_DIR +"/data/var31_minialfa_icemc2icemc_input_electric_field_1m.dat";
    //
    // SAW The file below has the correct spacing for the phase and electric field in volts/m/mhz
    //
    // We need the frequency spectrum to have 256 frequency bins with a bin spacing of 5.07812 MHz
    // so that it can be properly convolved with the ANITA signal electronics.
    //
    // Since RFSignal calculates a the frequency binning from the length of the waveform and the deltaT,
    // that sets the requirements on the waveform to be
    // length = 2 x NFREQ - 1 = 511
    // deltaT = 1/(length x df ) = 3.85368436041e-10 s
  
  
    std::string fileName(EnvironmentVariable::ICEMC_SRC_DIR());
    fileName += "/data/var31_minialfa_icemc2icemc_input_electric_field_1m_511points_window.dat";

    std::ifstream input_stream(fileName);
    if(!input_stream.good()){
      icemc::report() << severity::warning << "Couldn't open " << fileName << std::endl;
      return;
    }
    else {
      icemc::report() << severity::info << "Reading " << fileName << std::endl;
    }

    std::string buffer;
    int c = 0;
    
    while(!input_stream.eof()){      
      getline(input_stream,buffer,'\n');
      //skip the header
      if( c != 0){
	//simple stringstream read which assumes each line is "index,volts,time" with no spaces
	//or any other funny business. The ignore commands skip the commans
	std::stringstream ss;
	int ind=-999;
	float v=-999., t=-999.;
	ss << buffer;
            
	ss >> ind;
	ss.ignore(1);
	ss >> v;
	ss.ignore(1);
	ss >> t;
	ss.ignore(1);
            
	if( ind != -999 ){
	  fVoltsPerMeter.push_back(v);

	  fTimes.push_back(t);

	  //cout << buffer << endl;
	  // std::cout << c-1 << "\t" << time[c-1] << "\t" << voltsperm[c-1] << std::endl;
	}
      }

      if(!input_stream.eof()){
	c = c + 1;
      }
    }
    std::cout << "Counted " << c << " rows from electric field file" << fileName << std::endl;
    
    input_stream.close();  
  }  
}



icemc::PropagatingSignal icemc::WaisPulser::generateImpulse(const OpticalPath& opticalPath, const Neutrino& nu, const Shower& shower) const {

  auto eField = fVoltsPerMeter;
  for(auto& e : eField){
    e *= fSettings->jamieFactor;
  }
  
  FTPair signal(eField, fTimes[1] - fTimes[0], fTimes[0]);

  const TVector3& rfDir = opticalPath.steps.at(0).direction().Unit();
  TVector3 polarization = rfDir.Cross(fSurfaceNormal); 
  PropagatingSignal ps(signal, opticalPath.steps.at(0).direction(), polarization.Unit());
  return ps;

};
