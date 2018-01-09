// #include <iostream>
// #include <fstream>
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
#include "hot-loop.h"
// #include "TApplication.h"
#include <dlfcn.h>
// #include "TROOT.h"

std::string some_imortant_str = "some_important_str value";

int main(void) {

  // TApplication *theApp = new TApplication("tapp", NULL, NULL);
  // theApp = theApp;
  // gROOT->SetStyle("Plain");

  hot_loop("/nfs/data_disks/scratch1/bugaev/PROGS/INTER_C/canvas.so");

  return 0;
}
