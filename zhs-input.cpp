#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include "bvv-macro.h"
#include "bvv-util.h"

extern "C" void atmosinit_(int *modlabel, char* atmosname, int atmosname_len);
extern "C" double zfromdepth_(double *depth, int *layer);


class TSimPar {
  char buf[256];
public: 
  double Ze;
  double Az;   // astro-ph/9911331v1: The azimuth angle phi is the angle between the horizontal projection of the shower axis and  the x-axis (0 <= phi < 360 deg).
  double HMax; // vertical, g/cm^2;
  double AntX; // Antenna X [m].
  double AntY; // Antenna Y [m].
  double AntZ; // Antenna Z [m].
  double Hg  ; // Ground level [m].

  TSimPar(double Ze = -9999, double Az = -9999, double HMax = -9999, double AntX = -9999, double AntY = -9999, double AntZ = -9999,
          double Hg = -9999): Ze(Ze), Az(Az), HMax(HMax), AntX(AntX), AntY(AntY), AntZ(AntZ), Hg(Hg) {}
  std::string str() {
    sprintf(buf, "Ze: %9.4f, Az: %6.2f, HMax: %9.2f\nAntX: %6.1f, AntY: %6.1f, AntZ: %6.1f\nHg: %6.1f", Ze, Az, HMax, AntX, AntY, AntZ, Hg);
    return std::string(buf);
  }
};

void atmosinit() {
    int modlabel = 2;
    char atmosname[256] = "never-mind-atmo-name";
    atmosinit_(&modlabel, atmosname, strlen(atmosname));
}

double AltFromDepth(double depth) {
  // Returning altitude in meters.
  double Alt = -1;
  int layer = -1;
  Alt = zfromdepth_(&depth, &layer);
  // std::cout << "layer: " << layer << std::endl;
  return Alt;
}

double ZenSry(std::vector<std::string> tokens){
  // Find shower Ze in an array of tokens.
  double result = -9999;
  if (tokens.size() == 5 && tokens.at(0) == "Primary" && tokens.at(1) == "zenith") {
    result = atof(tokens.at(3).c_str());
  }
  return result;
} 

double ZenInp(std::vector<std::string> tokens){
  // Find shower Ze in an array of tokens.
  double result = -9999;
  if (tokens.size() == 3 && tokens.at(0) == "PrimaryZenAngle") {
    result = atof(tokens.at(1).c_str());
  }
  return result;
} 

double AzSry(std::vector<std::string> tokens){
  // Find shower Az in an array of tokens.
  double result = -9999;
  if (tokens.size() == 5 && tokens.at(0) == "Primary" && tokens.at(1) == "azimuth") {
    result = atof(tokens.at(3).c_str());
  }
  return result;
} 

double HMaxSry(std::vector<std::string> tokens) {
  // Find vert. shower max in an array of tokens.
  // If match not found, return -9999.
  // If math is found, the units are g/cm^2.
  double result = -9999;
  if (tokens.size() >= 6) {
    if (tokens[0] == "Vt." && tokens[4] == "(g/cm2):" ) {
      result = atof(tokens.at(5).c_str());
    }
  }
  return result;
}

int AntXYZInp( std::vector<std::string> tokens, double& AntXCandidate, double& AntYCandidate, double& AntZCandidate){
  // Find Antenna X, Y, Z coordinates in an array of tokens.
  // If match not found, return -9999.
  // If math is found, the unit is [m].
  double result = -9999;
  if (tokens.size() == 4) {
    if (tokens[0] == "AddAntenna") {
      AntXCandidate = atof(tokens.at(1).c_str());
      AntYCandidate = atof(tokens.at(2).c_str());
      AntZCandidate = atof(tokens.at(3).c_str());
      result = 0;
    }
  }
  return result;
}


double HgInp(std::vector<std::string> tokens){
  // Find ground altitude (Hg) in an array of tokens.
  double result = -9999;
  if (tokens.size() == 3 && tokens.at(0) == "GroundAltitude") {
    result = atof(tokens.at(1).c_str());
  }
  return result;
} 

/* Find shower max ZHS (Sx, Sy, Sz) coordinates */
double ShowerXYZ(double& Sx, double& Sy, double& Sz, double Ze = 70.291, double Az = 30, double RE = 1.4828, double Hg = 1.01718, double Hmax = 4.091) {
  // Inputs:
  // Default values (except for Az) are for the test case shown in https://www.geogebra.org/m/hunra7du;
  // The above plot is 2D, therefore Az does not apply.
  // Ze [deg]: shower Ze from ZHS input file.
  // Az [deg]: shower Az from ZHS input file.
  // RE: Earth radius [m]. I am not sure what to do with the fact that IceMC RE is different from AireS RE.
  // Hg: Ground level [m].
  // Hmax: Altitude of shower max [m].
  // Outputs:
  // Sx, Sy, Sz: shower max coordinates.
  // return value: Hsl - "slant altitude" of shower max - the distance between ZHS (0, 0, 0) and shower max.
  const double pi = atan(1.0) * 4;
  const double deg = pi / 180;
  double v1 = cos(Ze * deg);
  double v2 = Hg * Hg;
  double Hsl = sqrt((RE*RE + 2.0*Hg*RE + v2)*v1*v1 +(2.0*Hmax - 2.0*Hg)*RE + Hmax*Hmax - 1.0*v2) + (-1.0*RE - Hg)*v1;
  Sx = Hsl * cos(Az * deg);
  Sy = Hsl * sin(Az * deg);
  Sz = Hsl * v1;
  return Hsl;
}

/* UnitRA: find the unit vector describing EM wave direction of propagation,
which is the vector connecting the reflection point with ANITA, renormalized to 1.
*/
void UnitRA(double& ux, double& uy, double& uz, double Sx=-5, double Sy=-3, double Sz=5, double Ax=5, double Ay=-5, double Az=10){
  // Input parameters:
  // S_[xyz]: shower core coordinates in the ZHS reference frame.
  // A_[xyz]: ANITA       coordinates in the ZHS reference frame. 
  // Default values correspond to the geometry described in https://www.geogebra.org/m/kfedka7s
  // Output parameters:
  // u_[xyz]: The unit vector coordinates.
  double v1 = Az*Az;
  double v2 = Sz*Sz;
  double v3 = sqrt(v2 + 2*Az*Sz + v1);
  double v4 = 1./sqrt(v1*v2 + 2*v1*Az*Sz + v1*Sy*Sy - 2*Ay*v1*Sy + v1*Sx*Sx - 2*Ax*v1*Sx + v1*v1 + (Ay*Ay + Ax*Ax)*v1);
  double v5 = -Sx;
  double v6 = v5 + Ax;
  double v7 = -Sy;
  double v8 = v7 + Ay;
  double v9 = sqrt(v8*v8 +v6*v6);
  double v10 = 1/(Sz + Az);
  double v11 = v10/(v9*v10);
  double v12 = v10*v11;
  double v13 = v3*v4;
  ux = v13*(-v6*v9*Sz*v12 + v5 + Ax);
  uy = v13*(-v9*v8*Sz*v12 + v7 + Ay);
  uz = Az*v13; 
}

TSimPar ZHSSimPar(std::string PathRoot="Event_", std::string EvNum="4212") {
  // Input:
  // - RepoPath: directory of the event of interest minus the number and "_" at the end.
  // - EvNum: event number, std::string form.
  // Output:
  // TSimPar result.
  std::string EventDir = PathRoot + EvNum;
  std::string ZHSSummaryFile = EventDir + "/Event_" + EvNum + ".sry";
  std::string ZHSInputFile = EventDir + "/Event_" + EvNum + ".inp";
  std::cout << "Composed paths: " << ZHSSummaryFile << ", " << ZHSInputFile << std::endl;
  double HMax = -9999;
  double Ze = -9999;
  double Az = -9999;
  double AntX = -9999, AntY = -9999, AntZ = -9999;
  double Hg = -9999;

  // Grab Shower Ze, GroundAltitude (Hg) and AntXYZ (Antenna coordinates) from an event *.inp file.
  // Shower Az and Hmax will be taken from the event's *.sry file.
  WITH_LINES(
             ZHSInputFile.c_str(),
             ind,
             tokens,
             try {
               double ZeCandidate = ZenInp(tokens);
               if (ZeCandidate > -999) {Ze = ZeCandidate;}

               double HgCandidate = HgInp(tokens);
               if (HgCandidate > -999) {Hg = HgCandidate;}

               double AntXCandidate, AntYCandidate, AntZCandidate;
               int ret = AntXYZInp(tokens, AntXCandidate, AntYCandidate, AntZCandidate);
               if (ret > -999) {
                 AntX = AntXCandidate;
                 AntY = AntYCandidate;
                 AntZ = AntZCandidate;
                 // And we are done with this file:
                 break;
               }
             }
             catch (...) { std::cout << "<===" << std::endl << "BVV-EXCEPTION" << std::endl; break; }
             );

  // Grab Az and Hmax from ZHS *.sry file.
  WITH_LINES(
             ZHSSummaryFile.c_str(),
             ind,
             tokens,
             try {
               // *.sry doesn't provide precise enough output for Ze.
               // double ZeCandidate = ZenSry(tokens);
               // if (ZeCandidate > -999) {Ze = ZeCandidate;}

               double AzCandidate = AzSry(tokens);
               if (AzCandidate > -999) {Az = AzCandidate;}

               double HMaxCandidate = HMaxSry(tokens);
               if (HMaxCandidate > -999) {
                 HMax = HMaxCandidate;
                 // And we are done with this file:
                 break;
               }
             }
             catch (...) { std::cout << "<===" << std::endl << "BVV-EXCEPTION" << std::endl; break; }
             );
  return TSimPar(Ze, Az, HMax, AntX, AntY, AntZ, Hg);
}


// Compute direction of RF wave propagation simulated by ZHS.
void dir2bn(double& ux, double& uy, double& uz, std::string PathRoot="./Event_", std::string EvNumStr = "4212"){
  // Inputs:
  // PathRoot: name of the event directory with "_" and event number stripped from the end.
  // EvNumStr: ZHS event number in form of a string.
  // Outputs:
  // u[xyz]: Unit vector of the direction of RF wave propagation.

  TSimPar SimPar = ZHSSimPar(PathRoot, EvNumStr);
  atmosinit();
  // std::cout << "Altitude: " << AltFromDepth(SimPar.HMax) << std::endl;

  double Sx, Sy, Sz;
  double AiresEarthRad = 6370949;
  ShowerXYZ(Sx, Sy, Sz, SimPar.Ze, SimPar.Az, AiresEarthRad, SimPar.Hg, AltFromDepth(SimPar.HMax));
  UnitRA(ux, uy, uz, Sx, Sy, Sz, SimPar.AntX, SimPar.AntY, SimPar.AntZ);
}
