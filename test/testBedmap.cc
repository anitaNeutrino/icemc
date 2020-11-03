#include "Crust.h"
#include "Antarctica.h"

int main(int argc, char **argv){

 
  std::cout << "Creating Antarctica..." << std::endl;

  auto antarctica = std::make_shared<icemc::Antarctica>(0);
  
  std::cout << antarctica->GetTotalIceVolume() << " m^3 of Antarctic ice" << std::endl;
  std::cout << antarctica->GetTotalIceArea() << " m^2 of Antarctic is covered by ice" << std::endl;

  std::cout << "Done!" << std::endl;
}
