#include "PropagatingSignal.h"
#include "Constants.h"

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


icemc::SignalSummary icemc::PropagatingSignal::propagate(const OpticalPath& opticalPath){
  
  // 1/r loss from power intensity on spherical wavefront
  const double distanceFactor = 1./opticalPath.distance();

  // propagation through dielectric
  const double attenuationFactor = opticalPath.attenuation();

  polarization = opticalPath.polarization(polarization);

  double fresnelFactor = polarization.Mag();
  assert(fresnelFactor < 1);
  
  polarization *= (1./fresnelFactor); // make unit length

  const double totalFieldReduction = distanceFactor*attenuationFactor*fresnelFactor;

  // apply the losses
  auto& amps = waveform.changeFreqDomain();
  for(auto& amp : amps){
    amp *= totalFieldReduction;
  }


  SignalSummary ss(this);
  return ss;  
}






