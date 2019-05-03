#ifndef  _ZHS_INPUT_H
#define _ZHS_INPUT_H 1

class TZHSSimPar {
  char buf[256];
public: 
  static double AiresEarthRad() { return 6370949; }
  double Ze;
  double Az;   // astro-ph/9911331v1: The azimuth angle phi is the angle between the horizontal projection of the shower axis and  the x-axis (0 <= phi < 360 deg).
  double HMax; // vertical, g/cm^2;
  double AntX; // Antenna X [m].
  double AntY; // Antenna Y [m].
  double AntZ; // Antenna Z [m].
  double Hg  ; // Ground level [m].

  TZHSSimPar(double Ze = -9999, double Az = -9999, double HMax = -9999, double AntX = -9999, double AntY = -9999, double AntZ = -9999,
          double Hg = -9999): Ze(Ze), Az(Az), HMax(HMax), AntX(AntX), AntY(AntY), AntZ(AntZ), Hg(Hg) {}
  std::string str() {
    sprintf(buf, "Ze: %9.4f, Az: %6.2f, HMax: %9.2f\nAntX: %6.1f, AntY: %6.1f, AntZ: %6.1f\nHg: %6.1f", Ze, Az, HMax, AntX, AntY, AntZ, Hg);
    return std::string(buf);
  }
};

TZHSSimPar dir2bn(double& ux, double& uy, double& uz, std::string RepoPath, std::string EvNumStr);

#endif /* !_HOT_LOOP_H */
