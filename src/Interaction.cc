#include "Interaction.h"
#include "Constants.h"
#include "Settings.h"
#include "Antarctica.h"
#include "Neutrino.h"

icemc::Interaction::Interaction(const Settings *settings1) {

  noway=0;
  wheredoesitleave_err=0;
  neverseesice=0;

  wheredoesitenterice_err=0;
  toohigh=0;
  toolow=0;

  iceinteraction=0;
  dtryingdirection=0.;
  dnutries=0.;

  weight_nu=0;
  weight_nu_prob=0;

  current = pickCurrent(); //setCurrent();
}


void icemc::Interaction::PickAnyDirection() {
  double rndlist[2];
  rndlist[0] = pickUniform();
  rndlist[1] = pickUniform();
  
  costheta_nutraject=2*rndlist[0]-1;

  // pick a neutrino azimuthal angle
  phi_nutraject=2*constants::PI*rndlist[1];
  
  // check that these give the right result
  double thetanu=acos(costheta_nutraject);
  
  double sinthetanu=sin(thetanu);
  
  // find direction vector of neutrino
  // **** are cosine and sine flipped?
  nnu.SetX(sinthetanu*cos(phi_nutraject));
  nnu.SetY(sinthetanu*sin(phi_nutraject));
  nnu.SetZ(costheta_nutraject);
}




int icemc::Interaction::PickDownwardInteractionPoint(const Geoid::Position&r_bn, const Settings *settings1, const Antarctica *antarctica1) {

  if(settings1->UNBIASED_SELECTION==1){
    if(antarctica1->PickUnbiased(this)){ // pick neutrino direction and interaction point
      dtryingdirection=1.;
      iceinteraction=1;
    }
    else{
      iceinteraction=0;
    }
  }
  else{
    iceinteraction=1;

    // If we require neutrinos from a particular position
    // we generate that cartesian position here

    static Geoid::Position specific_position;

    if (settings1->SPECIFIC_NU_POSITION) 
      {
	specific_position.SetLonLatAlt(settings1->SPECIFIC_NU_POSITION_LONGITUDE,
				       settings1->SPECIFIC_NU_POSITION_LATITUDE,
				       settings1->SPECIFIC_NU_POSITION_ALTITUDE);
      }

    do{
      ///@todo ibnposition
      posnu = antarctica1->pickInteractionPosition(r_bn);
    } while(settings1->SPECIFIC_NU_POSITION &&  (posnu - specific_position).Mag() > settings1->SPECIFIC_NU_POSITION_DISTANCE);
  }
  
  return 0;
}//PickDownwardInteractionPoint





/**
 * Choose CC or NC: get from ratios in Ghandi etal paper,
 * updated for the CTEQ6-DIS parton distribution functions (M.H. Reno, personal communication).
 * Need to add capability of using ratios from Connolly et al.
 */

icemc::Neutrino::Interaction::Current icemc::Interaction::pickCurrent() {
  double rnd = pickUniform();
  if (rnd<=0.6865254){ // 10^18 eV - 10^21 eV (use this one for ANITA)
//if (rnd<=0.6893498) // 10^17 eV - 10^20 eV (use this one for SalSA)
    return Neutrino::Interaction::Current::Charged;//"cc";
  }
  else{
    return Neutrino::Interaction::Current::Neutral;//"nc";  
  }
} //GetCurrent



int icemc::Interaction::getPdgCode() const {
  int pdgcode = -1;
  if (nuflavor==Neutrino::Flavor::e){
    pdgcode = 12;
  }
  else if (nuflavor==Neutrino::Flavor::mu){
    pdgcode = 14;
  }
  else if (nuflavor==Neutrino::Flavor::tau){
    pdgcode = 16;
  }
  return pdgcode;
}


