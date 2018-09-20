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
#include "hot-antenna.h"

extern struct xlibstruct *gxlib;
struct nk_context *ctx = &(gxlib->ctx);

enum RunMode { m_init, m_reload, m_step, NRunModes };

// extern double tmp_vhz[2][Anita::NFREQ];

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


void PlotSomething(std::map<std::string, void*> *env, RunMode mode, struct nk_context *ctx){
  if (mode == m_init) {
    (*env)["cHotTest"] = new TCanvas();
    (*env)["fHotTest"] = new std::unique_ptr<TF1>;
    return;
  }

  static bvv::TBuffer <double> Deviation(1.0);

  if (mode == m_reload) {
    ((std::unique_ptr<TF1> *) (*env)["fHotTest"])->reset(new TF1("fHotTest","sin([0]*x)/x",0,10));
    (*((std::unique_ptr<TF1> *) (*env)["fHotTest"]))->SetParameter(0, Deviation);

    if ( (*env).find("amps0") == (*env).end() ) {
      // not found
      (*env)["amps0"] = new std::unique_ptr<TGraph>;
    }
    ((std::unique_ptr<TGraph> *) (*env)["amps0"])->reset(new TGraph(Anita::NFREQ));
    for (int i = 0; i < Anita::NFREQ; i++) {
      // (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->SetPoint(i, i, tmp_vhz[0][i]);
      (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->SetPoint(i, i, 1.2345);
    }

    (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->SetLineColor(kBlue);
    (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->SetMarkerStyle(5);
    (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->SetMarkerColor(kBlue);
    (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->SetMarkerSize(2);
    (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->SetLineWidth(6);
    (*((std::unique_ptr<TGraph> *) (*env)["amps0"]))->Draw("AL");

  }

  if (mode == m_step) {
    nk_layout_row_dynamic(ctx, 25, 1);
    property_double(ctx, "STD: ", 0.0, Deviation, 15.0 /*max*/, 1.0 /*increment*/, 0.002 /*sensitivity*/);

    if (*Deviation) {
      ((std::unique_ptr<TF1> *) (*env)["fHotTest"])->reset(new TF1("fHotTest","sin([0]*x)/x",0,10));
      (*((std::unique_ptr<TF1> *) (*env)["fHotTest"]))->SetParameter(0, Deviation);
      (*((std::unique_ptr<TF1> *) (*env)["fHotTest"]))->Draw("L");
      (*((std::unique_ptr<TF1> *) (*env)["fHotTest"]))->SetLineColor(kBlue);
      (*((std::unique_ptr<TF1> *) (*env)["fHotTest"]))->SetLineWidth(2);
      (*((std::unique_ptr<TF1> *) (*env)["fHotTest"]))->GetXaxis()->SetTitle("Time [ns]");
      (*((std::unique_ptr<TF1> *) (*env)["fHotTest"]))->GetYaxis()->SetTitle("E [V/m]");
    } 
  }

  if (mode == m_reload || *Deviation) {
    ((TCanvas *) (*env)["cHotTest"])->Modified();
    ((TCanvas *) (*env)["cHotTest"])->Update();
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
  PlotSomething(env, m_reload, ctx);
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
