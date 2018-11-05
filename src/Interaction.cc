#include "Interaction.h"



std::ostream& operator<<(std::ostream& os, const icemc::Interaction::Current& c){
  switch(c){
  case icemc::Interaction::Current::Charged:
    return os << "Interaction:Current:Charged";
  case icemc::Interaction::Current::Neutral:
    return os << "Interaction:Current:Neutral";
  default:
    return os << "Interaction:Current:Unknown";
  }
}
