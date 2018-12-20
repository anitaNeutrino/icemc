#include "RicianRNG.h"
#include "TMath.h"


icemc::RicianRNG::~RicianRNG(){

}

void icemc::RicianRNG::pickRician(double signal_amplitude, double signal_phase, double sigma, double &amplitude, double &phase){
  double rand_gauss_a, rand_gauss_b;
  get_circular_bivariate_normal_random_variable(rand_gauss_a, rand_gauss_b);
    
  // Gives the gaussian-distributed random variables a standard deviation of sigma
  rand_gauss_a *= sigma;
  rand_gauss_b *= sigma;
    
  // Gives the gaussian-distributed random variables a mean of (v*cos(theta), v*sin(theta)) when v is the mean of the desired rician distribution
  rand_gauss_a += signal_amplitude * cos(signal_phase);
  rand_gauss_b += signal_amplitude * sin(signal_phase);
    
  // The Rician Distribution produces the probability of the the absolute value (radius) of a circular bivariate normal random variable:
  amplitude = sqrt(rand_gauss_a * rand_gauss_a + rand_gauss_b * rand_gauss_b);
  // Thus, the descriptor other than amplitude for the circular bivariate is given by a phase:
  phase = atan2(rand_gauss_b, rand_gauss_a);
  return;
}


void icemc::RicianRNG::get_circular_bivariate_normal_random_variable(double& rand_gauss_a, double& rand_gauss_b){
  double rand_uni_a = pickUniform(); //gRandom->Rndm() produces uniformly-distributed floating points in ]0,1]
  double rand_uni_b = pickUniform();
  double first_term = sqrt(-2. * TMath::Log(rand_uni_a));
  // Box-Muller transform from a bivariate uniform distribution from 0 to 1 to a gaussian with mean = 0 and sigma = 1
  rand_gauss_a = first_term * cos(2. * M_PI * rand_uni_b);
  rand_gauss_b = first_term * sin(2. * M_PI * rand_uni_b);
  return;
}

