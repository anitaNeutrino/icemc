#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "hot-loop.h"
#include "TApplication.h"

std::string some_imortant_str = "some_important_str value";

int main(void) {

  TApplication *theApp = new TApplication("tapp", NULL, NULL);
  theApp = theApp;

  std::cout << "Hello world!" << std::endl;
  // hot_loop("./canvas.so");
  hot_loop("/nfs/data_disks/scratch1/bugaev/PROGS/INTER_C/canvas.so");
  return 0; // Never returns.
}
