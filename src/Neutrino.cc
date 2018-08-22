#include "Neutrino.h"

std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::CurrentType& c){
  switch(c){
  case icemc::Neutrino::CurrentType::Charged:
    return os << "CurrentType::Charged";
  case icemc::Neutrino::CurrentType::Neutral:
    return os << "CurrentType::Neutral";
  default:
    return os << "Unknown CurrentType!";
  }
}

std::ostream& operator<<(std::ostream& os, const icemc::Neutrino::Flavor& f){
  switch(f){
  case icemc::Neutrino::Flavor::e:
    return os << "Neutrino::Flavor::e";
  case icemc::Neutrino::Flavor::mu:
    return os << "Neutrino::Flavor::mu";
  case icemc::Neutrino::Flavor::tau:
    return os << "Neutrino::Flavor::tau";
  default:
    return os << "Unknown Flavor!";
  }
}