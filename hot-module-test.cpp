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
#include "hot-module-test.h"

extern struct xlibstruct *gxlib;
struct nk_context *ctx = &(gxlib->ctx);

enum RunMode { m_init, m_reload, m_step, NRunModes };

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


void PlotSomething(struct hot_test_state *state, RunMode mode, struct nk_context *ctx){
  if (mode == m_init) {
    state->cHotTest = new TCanvas();
  }

  static bvv::TBuffer <double> Deviation(1.0);

  if (mode == m_reload) {
    state->fHotTest.reset(new TF1("fHotTest","sin([0]*x)/x",0,10));
    state->fHotTest->SetParameter(0, Deviation);

  }

  if (mode == m_step) {
    nk_layout_row_dynamic(ctx, 25, 1);
    property_double(ctx, "STD: ", 0.0, Deviation, 15.0 /*max*/, 1.0 /*increment*/, 0.002 /*sensitivity*/);
    state->fHotTest->SetParameter(0, Deviation);
    state->fHotTest->Draw("L");
    state->fHotTest->SetLineColor(kRed);
    state->fHotTest->SetLineWidth(2);
    state->fHotTest->GetXaxis()->SetTitle("Time [ns]");
    state->fHotTest->GetYaxis()->SetTitle("E [V/m]");
    state->cHotTest->Modified();
    state->cHotTest->Update();
  }


  if (mode == m_reload || *Deviation) {
    state->fHotTest.reset(new TF1("fHotTest","sin([0]*x)/x",0,10));
    state->fHotTest->SetParameter(0, Deviation);
    state->cHotTest->Modified();
    state->cHotTest->Update();
  }

  fflush(stdout);
}

static struct hot_test_state *game_init(bool bInteractive_arg)
{
  bInteractive = bInteractive_arg;

  struct hot_test_state *state = new hot_test_state;
  PlotSomething(state, m_init, NULL);
  printf("Init Done\n");
  return state;
}

static void game_finalize(struct hot_test_state *state)
{
  delete state;
}

static void game_reload(struct hot_test_state *state)
{
  PlotSomething(state, m_reload, ctx);
  cout << "reloaded dl" << endl;
}

static void game_unload(struct hot_test_state *state __attribute__ ((unused)))
{
  cout << "game unloaded" << endl;
}

static bool game_step_core(struct hot_test_state *state) {
  PlotSomething(state, m_step, ctx);
  return true;
}

static bool game_step(struct hot_test_state *state)
{
  if (bInteractive) {
    gSystem->ProcessEvents();
 
    if (nk_begin(ctx, "Plot Controls", nk_rect(50, 50, 200, 200),
                 NK_WINDOW_BORDER|NK_WINDOW_MOVABLE|NK_WINDOW_SCALABLE|
                 NK_WINDOW_CLOSABLE|NK_WINDOW_MINIMIZABLE|NK_WINDOW_TITLE)) {
      game_step_core(state);
    }
    else
      cout << "nk_begin failed" << endl;
    nk_end(ctx);
  }
  else {
    game_step_core(state);
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
