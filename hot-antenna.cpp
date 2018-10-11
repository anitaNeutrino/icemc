# include <iostream>
# include <stdlib.h>
# include "TCanvas.h"
# include "TSystem.h"
# include "game.h"
# include "TGraph.h"
# include "TF1.h"
# include "TH1F.h"
# include "TLegend.h"
# include "TView.h"
# include "TPolyLine3D.h"
# include "TPolyMarker3D.h"
#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_IMPLEMENTATION
#define NK_XLIB_IMPLEMENTATION
#include "nuklear.h"
#include "nuklear_xlib.h"
#include "buffer.hh"
#include "position.hh"
#include "vector.hh"
#include "Settings.h"
#include "anita.hh"
#include "balloon.hh"
#include "Constants.h"
#include "hot-antenna.h"

extern struct xlibstruct *gxlib;
struct nk_context *ctx = &(gxlib->ctx);

enum RunMode { m_init, m_reload, m_step, NRunModes };

// extern double tmp_vhz[2][Anita::NFREQ];
extern Settings* global_settings1;
extern Anita *global_anita1;
extern Balloon *global_bn1;
extern Vector direction2bn;
extern const double pi;
// bInteractive: can be assigned "false" once through game_init
// to instruct the module that it is in the batch mode.
// One time is enough since reloading should not happen in the batch mode.
bool bInteractive = true;

using namespace std;

#define INIT_VAR_ONCE(TYPE, NAME, INI)                                  \
  static void * p##NAME;                                                \
  if (mode == m_reload) {                                               \
    map<std::string, void*> &env = *penv;                               \
    std::map<std::string, void*>::iterator env_it = env.find(#NAME);    \
    if ( env_it == env.end() ) {                                        \
      p##NAME = new std::unique_ptr<TYPE>;                              \
      ((std::unique_ptr<TYPE> *) p##NAME)->reset(INI);                  \
      env[#NAME] = p##NAME;                                             \
    } else {                                                            \
      p##NAME = env_it->second;                                         \
    }                                                                   \
  }                                                                     \
  std::unique_ptr<TYPE> & NAME = *((std::unique_ptr<TYPE> *) p##NAME)

#define RESET_VAR_ONRELOAD(TYPE, NAME, INI)     \
  INIT_VAR_ONCE(TYPE, NAME, INI);               \
  if (mode == m_reload) {                       \
    NAME.reset(INI);                            \
  }

void property_int(struct nk_context* ctx, const char *name, int min, bvv::TBuffer <int> &buf, int max, int step, float inc_per_pixel){
  int val = buf; 
  nk_property_int(ctx, name, min, &val, max, step, inc_per_pixel);
  // cout << "val: " << val << endl;
  buf = val;
}

void property_double(struct nk_context* ctx, const char *name, double min, bvv::TBuffer <double> &buf, double max, double step, float inc_per_pixel){
  double val = buf; 
  nk_property_double(ctx, name, min, &val, max, step, inc_per_pixel);
  buf = val;
}

void checkbox_label(struct nk_context *ctx, const char *name, bvv::TBuffer <int> &buf) {
  int val = buf;
  nk_checkbox_label(ctx, name, &val);
  buf = val;
}

void PlotGain(std::map<std::string, void*> *penv, RunMode mode, struct nk_context *ctx) {

  if (mode == m_init) return;

  const double Dir2BalLen = 5;
  static int SelAntLayer = -1;
  static int SelAntFold = -1;
  const double RangeMax = 5;
  static Vector n_eplane, SelAnt_eplane;
  static Vector n_hplane, SelAnt_hplane;
  static Vector n_normal, SelAnt_normal;
  static Vector RotatedSel_eplane;
  static int LayerFromAntNum[48];
  static int PhiFromAntNum[48];
  static double SelAntX;
  static double SelAntY;
  static double SelAntZ;
  static int SelAntNum = 40;
  // int ReturnCode;

  static bvv::TBuffer <double> PolAngle(30);
  static bvv::TBuffer <int> SelLayer(0);
  static bvv::TBuffer <int> AntZoom(0);
  
  INIT_VAR_ONCE(TCanvas, cHotTest, new TCanvas());
  INIT_VAR_ONCE(TView, view, TView::CreateView(1));
  INIT_VAR_ONCE(TCanvas, cGain, new TCanvas());
  INIT_VAR_ONCE(bvv::TBuffer <int>, RotSpeed, new bvv::TBuffer <int> (0));
  RESET_VAR_ONRELOAD(vector <TPolyLine3D>, vAntNormals, new vector <TPolyLine3D> (48));
  RESET_VAR_ONRELOAD(vector <TPolyMarker3D>, vAntPos, new vector <TPolyMarker3D> (48));
  RESET_VAR_ONRELOAD(TPolyLine3D, lDir2Bal, new TPolyLine3D(2));
  RESET_VAR_ONRELOAD(TPolyLine3D, lPol, new TPolyLine3D(2));
  RESET_VAR_ONRELOAD(TPolyLine3D, lSelAnt_hplane, new TPolyLine3D(2));
  RESET_VAR_ONRELOAD(TPolyLine3D, lSelAnt_eplane, new TPolyLine3D(2));
  RESET_VAR_ONRELOAD(TGraph, gegain, new TGraph(Anita::NFREQ));
  RESET_VAR_ONRELOAD(TGraph, ghgain, new TGraph(Anita::NFREQ));

  cHotTest->cd();
  cHotTest->SetWindowPosition(200, 200);

  if (mode == m_reload) {
    int AntNum;
    // Map antenna number to (ilayer, ifold):
    AntNum = -1;
    for (int ilayer=0; ilayer < global_settings1->NLAYERS; ilayer++) { // loop over layers on the payload
      for (int ifold=0;ifold<global_anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
        AntNum++;
        LayerFromAntNum[AntNum] = ilayer;
        PhiFromAntNum[AntNum] = ifold;
      }
    }
    SelAntLayer = LayerFromAntNum[SelAntNum];
    SelAntFold = PhiFromAntNum[SelAntNum];
    SelAntX = global_anita1->antenna_positions[SelAntNum][0];
    SelAntY = global_anita1->antenna_positions[SelAntNum][1];
    SelAntZ = global_anita1->antenna_positions[SelAntNum][2];

    SelLayer.modified = true;
    // cout << "view: " << &view << endl;
    global_bn1->GetAntennaOrientation(global_settings1,  global_anita1,  SelAntLayer,  SelAntFold, SelAnt_eplane,  SelAnt_hplane,  SelAnt_normal);
    cout << "SelAnt_normal: " << SelAnt_normal << endl;
    
    view->SetRange(-RangeMax, -RangeMax, -RangeMax, +RangeMax, +RangeMax, +RangeMax);
    lDir2Bal->SetPoint(0, 0, 0, 0);
    lDir2Bal->SetPoint(1, SelAnt_normal[0] * Dir2BalLen, SelAnt_normal[1] * Dir2BalLen, SelAnt_normal[2] * Dir2BalLen);
    lDir2Bal->SetLineColor(kMagenta);
    lDir2Bal->Draw();
    AntNum = -1;
    for (int ilayer=0; ilayer < global_settings1->NLAYERS; ilayer++) { // loop over layers on the payload
      for (int ifold=0;ifold<global_anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
        AntNum++;
        // int antNum = global_anita1->GetRxTriggerNumbering(ilayer, ifold);
        // cout << "ilayer: " << ilayer << ", ifold: " << ifold << ", antNum: " << antNum << endl;

        global_bn1->GetAntennaOrientation(global_settings1,  global_anita1,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
        // vAntNormals->at(antNum) = TPolyLine3D(2);
        double AntX = global_anita1->antenna_positions[AntNum][0];
        double AntY = global_anita1->antenna_positions[AntNum][1];
        double AntZ = global_anita1->antenna_positions[AntNum][2];
        // double AntX = global_anita1->ANTENNA_POSITION_START[ilayer][ifold][0];
        // double AntY = global_anita1->ANTENNA_POSITION_START[ilayer][ifold][1];
        // double AntZ = global_anita1->ANTENNA_POSITION_START[ilayer][ifold][2];
        vAntPos->at(AntNum).SetPoint(0, AntX, AntY, AntZ);
        vAntNormals->at(AntNum).SetPoint(0, AntX, AntY, AntZ);
        vAntNormals->at(AntNum).SetPoint(1, AntX + n_normal[0], AntY + n_normal[1], AntZ + n_normal[2]);
        if (ilayer == SelAntLayer && ifold == SelAntFold) {
          lSelAnt_eplane->SetPoint(0, AntX, AntY, AntZ);
          lSelAnt_eplane->SetPoint(1, AntX + SelAnt_eplane[0], AntY + SelAnt_eplane[1], AntZ + SelAnt_eplane[2]);

          lSelAnt_hplane->SetPoint(0, AntX, AntY, AntZ);
          lSelAnt_hplane->SetPoint(1, AntX + SelAnt_hplane[0], AntY + SelAnt_hplane[1], AntZ + SelAnt_hplane[2]);

          lSelAnt_eplane->Draw();
          lSelAnt_hplane->Draw();
          cout << "AntXYZ: " << AntX << ", " << AntY << ", " << AntZ << endl;
          cout << "SelAnt_hplane[123]: " << SelAnt_hplane[0] << ", " << SelAnt_hplane[1] << ", " << SelAnt_hplane[2] << endl;
        }

        vAntPos->at(AntNum).SetMarkerSize(1);
        vAntPos->at(AntNum).SetMarkerColor(kBlue);
        vAntPos->at(AntNum).SetMarkerStyle(8);
        vAntNormals->at(AntNum).SetLineWidth(1);
        vAntNormals->at(AntNum).SetLineColor(kBlue);

        vAntNormals->at(AntNum).Draw();
        vAntPos->at(AntNum).Draw();
        // std::cout << "n_normal: " << n_normal << std::endl;
      }
    }

    cHotTest->Modified();
    cHotTest->Update();
  }

  cHotTest->cd();

  if (mode == m_step) {
    nk_layout_row_dynamic(ctx, 25, 1);
    property_int(ctx, "Rot Speed: ", 0, *RotSpeed, 50 /*max*/, 1 /*increment*/, 0.5 /*sensitivity*/);
    property_double(ctx, "Pol Angle: ", 0, PolAngle, 359 /*max*/, 1 /*increment*/, 0.5 /*sensitivity*/);
    checkbox_label(ctx, "Ant. Zoom: ", AntZoom);
    if (*AntZoom) {
      if (AntZoom) {
        const double RangeZoomMax = 1;
        view->SetRange(SelAntX - RangeZoomMax, SelAntY - RangeZoomMax, SelAntZ - RangeZoomMax, SelAntX + RangeZoomMax, SelAntY + RangeZoomMax, SelAntZ + RangeZoomMax);
      } else {
        view->SetRange(-RangeMax, -RangeMax, -RangeMax, +RangeMax, +RangeMax, +RangeMax);
      }
    }

    if (*SelLayer) { // This is either from the m_reload or previous step.
      int AntNum = -1;
      for (int ilayer=0; ilayer < global_settings1->NLAYERS; ilayer++) { // loop over layers on the payload
        for (int ifold=0;ifold<global_anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
          AntNum++;
          // int antNum = global_anita1->GetRxTriggerNumbering(ilayer, ifold);

          int NormalsColor;
          if (ilayer == SelLayer)
            NormalsColor = kRed;
          else
            NormalsColor = kBlue;
          vAntNormals->at(AntNum).SetLineColor(NormalsColor);
          vAntPos->at(AntNum).SetMarkerColor(NormalsColor);

          if (ilayer == SelAntLayer && ifold == SelAntFold){
            vAntNormals->at(AntNum).SetLineColor(kMagenta);
            vAntPos->at(AntNum).SetMarkerColor(kMagenta);
          }
          vAntNormals->at(AntNum).Draw();
        }
      }
      // I cannot wait until the end of m_step block because I rely on SelLayer.modified == true which expire after the first redrawing of the corresponding widget:
      cHotTest->Modified();
      cHotTest->Update();
    }
    // The following is shifted to the bottom since "modified" flag may have been inherited from the reload section and we
    // don't want to clear it by the call to the "property_int".
    property_int(ctx, "SelLayer: ", 0, SelLayer, global_settings1->NLAYERS - 1/*max*/, 1 /*increment*/, 0.5 /*sensitivity*/);
    int ReturnCode;
    if (*RotSpeed != 0) {
      view->SetView(view->GetLongitude() + 0.1 * *RotSpeed, view->GetLatitude() + 0, 0, ReturnCode);
    }
  }

  if (mode == m_reload || *PolAngle) {
    RotatedSel_eplane = SelAnt_eplane.Rotate(PolAngle * pi / 180, SelAnt_normal);
    lPol->SetPoint(0, SelAntX, SelAntY, SelAntZ);
    lPol->SetPoint(1, SelAntX + RotatedSel_eplane[0], SelAntY + RotatedSel_eplane[1], SelAntZ + RotatedSel_eplane[2]);
    lPol->SetLineColor(kGreen);
    lPol->Draw();
    // Computing gains:
    cGain->cd(); 
    double e_component_kvector, h_component_kvector, n_component_kvector;
    global_bn1->GetEcompHcompkvector(SelAnt_eplane,  SelAnt_hplane,  SelAnt_normal, /* -SelAnt_hplane */ -SelAnt_normal /* <-- direction2bn */, e_component_kvector,  h_component_kvector,  n_component_kvector);
    double e_component, h_component, n_component;
    global_bn1->GetEcompHcompEvector(global_settings1,  SelAnt_eplane,  SelAnt_hplane, /* SelAnt_eplane */ RotatedSel_eplane /* <-- polarization */,  e_component,  h_component,  n_component);
    double hitangle_e, hitangle_h;
    global_bn1->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);
    // double e_signal[global_anita1->NFREQ] = {1.0}, h_signal[global_anita1->NFREQ] = {1.0};

    for (int ifreq = 0; ifreq < Anita::NFREQ; ifreq++){
      double e_signal = 1.0, h_signal = 1.0;
      global_anita1->AntennaGain(global_settings1, hitangle_e, hitangle_h, e_component, h_component, ifreq, e_signal, h_signal);
      gegain->SetPoint(ifreq, ifreq, e_signal);
      ghgain->SetPoint(ifreq, ifreq, h_signal);
      // cout << "hitangles: " << hitangle_e << ", " << hitangle_h << endl;
      // cout << "eh_signal: " << e_signal << ", " << h_signal << endl;
    }
    gegain->SetLineColor(kBlue);
    gegain->SetLineWidth(2);
    ghgain->SetLineColor(kRed);
    ghgain->SetLineWidth(2);
    ghgain->Draw("AL");
    ghgain->GetYaxis()->SetRangeUser(0, 0.45);
    gegain->Draw("L");

  }

  if (mode == m_reload /* not needed: || *SelLayer */ || *AntZoom || *RotSpeed || *PolAngle) {
    cHotTest->Modified();
    cHotTest->Update();
    cGain->Modified();
    cGain->Update();
    cout << "Step updated canvases" << endl;
  }
}

  void PlotSomething(std::map<std::string, void*> *penv, RunMode mode, struct nk_context *ctx){
    map<std::string, void*> &env = *penv;
    if (mode == m_init) {
      env["cHotTest"] = new TCanvas();
      return;
    }

    static bvv::TBuffer <double> Deviation (1.0);

    if (mode == m_reload) {
      // BEGIN Introduce new persistent entries in "env":
      if ( env.find("gAmps0") == env.end() ) {
        // not found
        env["gAmps0"] = new std::unique_ptr<TGraph>;
      }
      if ( env.find("fHotTest") == env.end() ) {
        // not found
        env["fHotTest"] = new std::unique_ptr<TF1>;
      }
      // END   Introduce new persistent entries in "env".

      // BEGIN Create shortcuts to objects listed in "env":
      std::unique_ptr<TF1> & fHotTest = *((std::unique_ptr<TF1> *) env["fHotTest"]);
      std::unique_ptr<TGraph> & gAmps0 = *((std::unique_ptr<TGraph> *) env["gAmps0"]);
      // END   Create shortcuts to objects listed in "env".


      fHotTest.reset(new TF1("fHotTest","sin([0]*x)/x",0,10));
      fHotTest->SetParameter(0, Deviation);

      gAmps0.reset(new TGraph(Anita::NFREQ));

      for (int i = 0; i < Anita::NFREQ; i++) {
        gAmps0->SetPoint(i, i, 1.2345);
      }

      gAmps0->SetLineColor(kBlue);
      gAmps0->SetMarkerStyle(5);
      gAmps0->SetMarkerColor(kBlue);
      gAmps0->SetMarkerSize(2);
      gAmps0->SetLineWidth(6);
      gAmps0->Draw("AL");

    }

    if (mode == m_step) {
      // BEGIN Create shortcuts to objects listed in "env":
      std::unique_ptr<TF1> & fHotTest = *((std::unique_ptr<TF1> *) env["fHotTest"]);
      // END   Create shortcuts to objects listed in "env".
      nk_layout_row_dynamic(ctx, 25, 1);
      property_double(ctx, "STD: ", 0.0, Deviation, 15.0 /*max*/, 1.0 /*increment*/, 0.002 /*sensitivity*/);

      if (*Deviation) {
        fHotTest.reset(new TF1("fHotTest","sin([0]*x)/x",0,10));
        fHotTest->SetParameter(0, Deviation);
        fHotTest->Draw("L");
        fHotTest->SetLineColor(kRed);
        fHotTest->SetLineWidth(2);
        fHotTest->GetXaxis()->SetTitle("Time [ns]");
        fHotTest->GetYaxis()->SetTitle("E [V/m]");
      } 
    }

    if (mode == m_reload || *Deviation) {
      ((TCanvas *) env["cHotTest"])->Modified();
      ((TCanvas *) env["cHotTest"])->Update();
    }

    fflush(stdout);
  }

static std::map<std::string, void*> *game_init(bool bInteractive_arg)
{
  bInteractive = bInteractive_arg;

  std::map<std::string, void*> *env = new std::map<std::string, void*>;
  // PlotSomething(env, m_init, NULL);
  printf("Init Done\n");
  return env;
}

static void game_finalize(std::map<std::string, void*> *env)
{
  delete env;
}

static void game_reload(std::map<std::string, void*> *env)
{
  // PlotSomething(env, m_reload, ctx);
  PlotGain(env, m_reload, ctx);
  cout << "reloaded dl" << endl;
}

static void game_unload(std::map<std::string, void*> *env __attribute__ ((unused)))
{
  cout << "game unloaded" << endl;
}

static bool game_step_core(std::map<std::string, void*> *env) {
  // PlotSomething(env, m_step, ctx);
  PlotGain(env, m_step, ctx);
  return true;
}

static bool game_step(std::map<std::string, void*> *env)
{
  if (bInteractive) {
    gSystem->ProcessEvents();
 
    if (nk_begin(ctx, "Plot Controls", nk_rect(50, 50, 200, 200),
                 NK_WINDOW_BORDER|NK_WINDOW_MOVABLE|NK_WINDOW_SCALABLE|
                 NK_WINDOW_CLOSABLE|NK_WINDOW_MINIMIZABLE|NK_WINDOW_TITLE)) {
      game_step_core(env);
    }
    else
      cout << "nk_begin failed" << endl;
    nk_end(ctx);
  }
  else {
    game_step_core(env);
  }

  return true;
}

__attribute__((visibility("default"))) extern const struct game_api GAME_API = {
  (void *(*)(bool))  game_init, // game_init
  (void (*)(void *)) game_finalize, // game_finalize
  (void (*)(void *)) game_reload, // game_reload
  (void (*)(void *)) game_unload, // game_unload
  (bool (*)(void *)) game_step // game_step
};
