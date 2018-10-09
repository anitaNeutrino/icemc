#include "Energy.h"

std::ostream& operator<<(std::ostream& os, const icemc::Energy& e){
  return os << e.in(icemc::Energy::Unit::eV) << " eV";
}

icemc::Energy operator* (double lhs, icemc::Energy rhs){
  return rhs *= lhs;
}
