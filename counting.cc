#include "vector.hh"
#include "position.hh"
#include "counting.hh"
#include "Tools.h"
#include "TMath.h"

Counting::Counting() {
  Tools::Zero(npass,2);
  Tools::Zero(npassestrigger,2);
  Tools::Zero(nchanceinhell2,2);
  Tools::Zero(nviewanglecut,2);
  Tools::Zero(nchanceinhell,2);
  Tools::Zero(nchanceinhell_1overr,2);
  Tools::Zero(nchanceinhell_fresnel,2);
  Tools::Zero(nconverges,2);
  Tools::Zero(nacceptablerf,2);
  Tools::Zero(nraywithincontinent1,2); // reality check:  exiting ray is within 30 degrees of south pole
  Tools::Zero(nraywithincontinent2,2); // same, after next iteration.
  Tools::Zero(nraypointsup1,2); // same, after next iteration.
  Tools::Zero(nraypointsup2,2); // same, after next iteration.
  Tools::Zero(nnottoosmall,2); // same, after next iteration.
  Tools::Zero(nviewangle_lt_90,2); // same, after next iteration.
  Tools::Zero(ngoodfracs,2); // same, after next iteration.
  Tools::Zero(nbadfracs,2); // same, after next iteration.
  Tools::Zero(nnottir,2); // same, after next iteration.
  Tools::Zero(nentersice,2); // same, after next iteration.
  Tools::Zero(nabsorbed,2); // same, after next iteration.
  Tools::Zero(noway,2); // same, after next iteration.
  Tools::Zero(wheredoesitleave_err,2); // same, after next iteration.
  Tools::Zero(neverseesice,2); // same, after next iteration.
Tools::Zero(iceinteraction,2); // same, after next iteration.
Tools::Zero(inhorizon,2); // same, after next iteration.
 Tools::Zero(wheredoesitenterice_err,2); // same, after next iteration.
Tools::Zero(toohigh,2); // same, after next iteration.
Tools::Zero(toolow,2); // same, after next iteration.
  for (int i=0;i<NCOSTHETA;i++) {
    Tools::Zero(weights_rin[i],NPHI); // same, after next iteration.
  }
// variables for counting neutrinos and reporting results.
  nnu_e=0;  //counting the number of e,mu,tau neutrinos
  nnu_mu=0;  
  nnu_tau=0;



}
void Counting::findCosthetaPhiBins(Position r_in,int &icostheta,int &iphi) {

  icostheta=(int)((cos(r_in.Theta())-COSTHETAMIN)/(COSTHETAMAX-COSTHETAMIN)*(double)NCOSTHETA);
  iphi=(int)((r_in.Phi()-PHIMIN)/(PHIMAX-PHIMIN)*(double)NPHI);


}
void Counting::IncrementWeights_r_in(Position r_in,double weight) {
  int iphi,icostheta;
  findCosthetaPhiBins(r_in,icostheta,iphi);
  weights_rin[icostheta][iphi]+=weight;

}
