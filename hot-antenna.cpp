# include <iostream>
# include <stdlib.h>
# include "TCanvas.h"
# include "TSystem.h"
# include "game.h"
# include "TGraph.h"
# include "TF1.h"
# include "TH1F.h"
# include "TLegend.h"
#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_IMPLEMENTATION
#define NK_XLIB_IMPLEMENTATION
#include "nuklear.h"
#include "nuklear_xlib.h"
#include "buffer.hh"
#include "vector.hh"
#include "anita.hh"
#include "Settings.h"
#include "anita.hh"
#include "hot-antenna.h"

extern struct xlibstruct *gxlib;
struct nk_context *ctx = &(gxlib->ctx);

enum RunMode { m_init, m_reload, m_step, NRunModes };

// extern double tmp_vhz[2][Anita::NFREQ];
extern Settings* global_settings1;
extern Anita *global_anita1;

// bInteractive: can be assigned "false" once through game_init
// to instruct the module that it is in the batch mode.
// One time is enough since reloading should not happen in the batch mode.
bool bInteractive = true;

using namespace std;

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
  if (mode == m_reload) {
    for (int ilayer=0; ilayer < global_settings1->NLAYERS; ilayer++) { // loop over layers on the payload
      for (int ifold=0;ifold<global_anita1->NRX_PHI[ilayer];ifold++) { // ifold loops over phi
        cout << "ilayer: " << ilayer << ", ifold: " << ifold << ", antNum: " << global_anita1->GetRxTriggerNumbering(ilayer, ifold) << endl;
      }
    }
  }
  // bn1->GetAntennaOrientation(settings1,  anita1,  ilayer,  ifold, n_eplane,  n_hplane,  n_normal);
  // bn1->GetEcompHcompkvector(n_eplane,  n_hplane,  n_normal,  panel1->GetVec2bln(jpt), e_component_kvector,  h_component_kvector,  n_component_kvector);
  // bn1->GetEcompHcompEvector(settings1,  n_eplane,  n_hplane,  panel1->GetPol(jpt),  e_component,  h_component,  n_component);
  // bn1->GetHitAngles(e_component_kvector, h_component_kvector, n_component_kvector, hitangle_e, hitangle_h);

  // anita1->AntennaGain(settings1, hitangle_e, hitangle_h, e_component, h_component, k, tmp_vhz[0][k], tmp_vhz[1][k]);
}

void PlotSomething(std::map<std::string, void*> *penv, RunMode mode, struct nk_context *ctx){
  map<std::string, void*> &env = *penv;
  if (mode == m_init) {
    env["cHotTest"] = new TCanvas();
    return;
  }

  static bvv::TBuffer <double> Deviation(1.0);

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
  PlotSomething(env, m_init, NULL);
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
  PlotSomething(env, m_step, ctx);
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
