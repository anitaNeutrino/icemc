#include "PropagatingSignal.h"
#include "Constants.h"
#include "Report.h"

double icemc::PropagatingSignal::energy() const {

  const TGraph& gr = waveform.getTimeDomain();
  double E_squared = 0;
  double dt = gr.GetX()[1] - gr.GetX()[0];
  for(int i=0; i < gr.GetN(); i++){
    E_squared += gr.GetY()[i]*gr.GetY()[i]*dt;
  }
  return E_squared*0.5*constants::VACUUM_PERMITTIVITY;
}


double icemc::PropagatingSignal::maxEField() const {
  const TGraph& gr = waveform.getTimeDomain();
  return *std::max_element(gr.GetY(), gr.GetY()+gr.GetN());  
}

double icemc::PropagatingSignal::maxVmMHz() const {
  const std::vector<std::complex<double> >& freqs = waveform.getFreqDomain();
  
  double max = 0;

  for( auto f : freqs){
    if(std::abs(f) > max)
      max = std::abs(f); 
  }
  
  return max;
}


void icemc::PropagatingSignal::propagate(const OpticalPath& opticalPath){
  
  // 1/r loss from power intensity on spherical wavefront
  const double distanceFactor = 1./opticalPath.distance();
  
  // propagation through dielectric
  const double attenuationFactor = opticalPath.attenuation();
   
  polarization = opticalPath.polarization(polarization);
  double fresnelFactor = polarization.Mag();

  // Magnification at ice-air barrier
  double magnificationFactor = opticalPath.magnification();

  const double tolerance = 1e-10;
  if(fresnelFactor > 1 + tolerance || fresnelFactor<=0){
    icemc::report() << severity::error << "Un-physical Fresnel factor! " <<  fresnelFactor << std::endl;
  }
  
  polarization *= (1./fresnelFactor); // make unit length

  const double totalFieldReduction = distanceFactor*attenuationFactor*fresnelFactor*magnificationFactor;

  // apply the losses
  auto& amps = waveform.changeFreqDomain();
  for(auto& amp : amps){
    amp *= totalFieldReduction;
  }  
}

void icemc::SignalSummary::set(const PropagatingSignal* s){
  if(s){
    maxEField = s->maxEField();
    energy = s->energy();
    maxVmMHz = s->maxVmMHz();
  }
}



