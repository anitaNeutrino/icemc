#include "Neutrino.h"

std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::Interaction::Current& c){
  switch(c){
  case icemc::Neutrino::Interaction::Current::Charged:
    return os << "Neutrino:Interaction:Current:Charged";
  case icemc::Neutrino::Interaction::Current::Neutral:
    return os << "Neutrino:Interaction:Current:Neutral";
  default:
    return os << "Neutrino:Interaction:Current:Unknown";
  }
}

std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::Flavor& f){
  switch(f){
  case icemc::Neutrino::Flavor::e:
    return os << "Neutrino:Flavor:e";
  case icemc::Neutrino::Flavor::mu:
    return os << "Neutrino:Flavor:mu";
  case icemc::Neutrino::Flavor::tau:
    return os << "Neutrino:Flavor:tau";
  default:
    return os << "Neutrino:Flavor:Unknown";
  }
}


std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::L& l){
  switch(l){
  case icemc::Neutrino::L::Matter:
    return os << "Neutrino:L:Matter";
  case icemc::Neutrino::L::AntiMatter:
    return os << "Neutrino:L:AntiMatter";
  default:
    return os << "Neutrino:L:Unknown";
  }
}
