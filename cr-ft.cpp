# include <iostream>
# include <stdlib.h>
# include "TCanvas.h"
# include "TSystem.h"
# include "game.h"
// # include "anita.hh"
# include "TGraph.h"
# include "TGaxis.h"
# include "TLine.h"
# include "TH1F.h"
# include "TLegend.h"
# include "FFTtools.h"
#define NK_INCLUDE_FIXED_TYPES
#define NK_INCLUDE_STANDARD_IO
#define NK_INCLUDE_STANDARD_VARARGS
#define NK_INCLUDE_DEFAULT_ALLOCATOR
#define NK_IMPLEMENTATION
#define NK_XLIB_IMPLEMENTATION
#include "nuklear.h"
#include "nuklear_xlib.h"
#include "buffer.hh"
#include "cr-ft.h"
#include "Tools.h"

extern struct xlibstruct *gxlib;
struct nk_context *ctx = &(gxlib->ctx);

enum RunMode { m_init, m_reload, m_step, NRunModes };
enum ComplexNumForm { cf_polar, cf_rect, NComplexNumForm };
namespace unit {
  const double ns = 1e-9;
}

// bInteractive: can be assigned "false" once through game_init
// to instruct the module that it is in the batch mode.
// One time is enough since reloading should not happen in the batch mode.
bool bInteractive = true;

using namespace std;

class TGraph1: public TGraph {
public:
  using TGraph::TGraph;
  ~TGraph1() { cout<<"Destructing TGraph1 " << fName << endl; } 
};

class TPad1: public TPad {
public:
  using TPad::TPad;
  ~TPad1() { cout<<"Destructing TPad1!!!" << endl; } 
};



extern const double pi;
extern int ZhsTimeN;
extern double ZhsTimeStart; 
extern double ZhsTimeDelta;
extern vector<double> ZhsTimeArr;
extern vector<double> ZhsTimeE;
extern vector<double> ZhsAlpha;
extern double NrFT[ANITA_TIME_SAMPLES];
bool FWHM(long int n, double *x, double *y, double &xmin, double &xmax, int &ind_maxval, double &xmaxval, double &ymaxval, double threshold_rel = 0.5);


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

void PlotIFT(struct cr_ft_state *state, RunMode mode, struct nk_context *ctx){
  // We shouldn't be here if in batch mode:
  if (!bInteractive) return;

  if (mode == m_init) {
    state->cZhsIFft = new TCanvas();
    return;
  }

  static bvv::TBuffer <int> Phi(0);
  static bvv::TBuffer <int> FreqLoCut(0);

  if (mode == m_step) {
    nk_layout_row_dynamic(ctx, 25, 2);
    // property_int(ctx, "Phase: ", -180 /*min*/, Phi, +180 /*max*/, 1 /*increment*/, 1.0 /*sensitivity*/);
    property_int(ctx, "First non-zero freq: ", 0 /*min*/, FreqLoCut, +2000 /*max*/, 1 /*increment*/, 1.0 /*sensitivity*/);
  }

  if (mode == m_reload || *state->vis_nbins || *Phi || *state->BinShift || *FreqLoCut) {
    bool bHiResIFT = false;
    bool bConstPhi = false;
    state->cZhsIFft->Clear();
    if (bHiResIFT) {
      if (bConstPhi) {
        for (int i = 0; i < state->vis_nbins / 2 + 1; i++) {
          double re = state->ZhsFft[i].re;
          double im = state->ZhsFft[i].im;
          double rho = TMath::Sqrt(re * re + im * im);
          double new_re = rho * TMath::Cos(Phi * TMath::Pi() / 180.);
          double new_im = rho * TMath::Sin(Phi * TMath::Pi() / 180.);
          state->ZhsFft[i].re = new_re;
          state->ZhsFft[i].im = new_im;
        }
      }
      if (mode == m_reload || *FreqLoCut) {
        int FreqLoCutBin = FreqLoCut * 1e6 / state->dfreq + 1;
        cout << "The actual FreqLoCut: " << FreqLoCutBin * state->dfreq / 1e6 << endl;
        cout << "FreqLoCut: " << FreqLoCut << endl;
        for (int i = 0; i <= (int) FreqLoCutBin; i++) {
          state->ZhsFft[i].re = 0;
          state->ZhsFft[i].im = 0;
        }
      }
      state->ZhsIFft.reset(FFTtools::doInvFFT(state->vis_nbins, state->ZhsFft.get()));
      state->grIFft.reset(new TGraph(state->vis_nbins));
      for (int i = 0; i < state->vis_nbins; i++) {
        state->grIFft->SetPoint(i, (i + state->vis_xmin_bin) * ZhsTimeDelta + ZhsTimeStart, state->ZhsIFft[i]);
      }
    }
    else { // bHiResIFT == false
      // Filling state->AnitaFT using state->ZhsFft:
      for (int i = 0; i < ANITA_FT_SAMPLES; i++) {
        double target_freq = ANITA_FREQ_HIGH / ANITA_FT_BINS * i;
        double target_freq_src_units = target_freq / state->dfreq;
        int int_target_freq_src_units = target_freq_src_units + 0.5;
        state->AnitaFT[i].re = state->ZhsFft[int_target_freq_src_units].re;
        state->AnitaFT[i].im = state->ZhsFft[int_target_freq_src_units].im;
      }

      // BEGIN Filling state->AnitaNRInp from state->AnitaFT and computing IFT with realft in place of state->AnitaNRInp:
      for (int i = 0; i < ANITA_TIME_SAMPLES / 2; i++) {
        state->AnitaNRInp[i * 2] = state->AnitaFT[i].re;
        state->AnitaNRInp[i * 2 + 1] = -state->AnitaFT[i].im;
      }

      state->AnitaNRInp[1] = state->AnitaFT[ANITA_TIME_SAMPLES / 2].re;

      Tools::realft(state->AnitaNRInp, -1, ANITA_TIME_SAMPLES);
      // END  Filling state->AnitaNRInp from state->AnitaFT and computing IFT with realft in place of state->AnitaNRInp.

      state->grINRFft.reset(new TGraph(ANITA_TIME_SAMPLES));
      // Filling state->grINRFft with state->AnitaNRInp[] containing IFT after a call to realft:
      for (int i = 0; i < ANITA_TIME_SAMPLES; i++){
        // By my design, ANITA_TIME_SAMPLES / 2 + 1 corresponds to ANITA_FREQ_HIGH.
        // That means ANITA_TIME_SAMPLES / 2 * (1 / T) = ANITA_FREQ_HIGH
        // Therefore T = ANITA_TIME_SAMPLES / 2 / ANITA_FREQ_HIGH
        // And the time step is T / (ANITA_TIME_SAMPLES - 1)
        double TimeDelta = ANITA_TIME_SAMPLES / 2 / ANITA_FREQ_HIGH / (ANITA_TIME_SAMPLES - 1); // Should be valid for any even number of time samples, not just ANITA_TIME_SAMPLES.
        state->grINRFft->SetPoint(i, i * TimeDelta / unit::ns + state->vis_xmin_bin * ZhsTimeDelta + ZhsTimeStart, state->AnitaNRInp[i] * 2.0 / ANITA_TIME_SAMPLES * ZhsTimeDelta * unit::ns / TimeDelta);
      }

      state->ZhsIFft.reset(FFTtools::doInvFFT(ANITA_TIME_SAMPLES, state->AnitaFT));
      // GDB: plot1d_opt *(state->ZhsIFft.get())@512 'with lines'

      state->grIFft.reset(new TGraph(ANITA_TIME_SAMPLES));
      for (int i = 0; i < ANITA_TIME_SAMPLES; i++){
        // By my design, ANITA_TIME_SAMPLES / 2 + 1 corresponds to ANITA_FREQ_HIGH.
        // That means ANITA_TIME_SAMPLES / 2 * (1 / T) = ANITA_FREQ_HIGH
        // Therefore T = ANITA_TIME_SAMPLES / 2 / ANITA_FREQ_HIGH
        // And the time step is T / (ANITA_TIME_SAMPLES - 1)
        double TimeDelta = ANITA_TIME_SAMPLES / 2 / ANITA_FREQ_HIGH / (ANITA_TIME_SAMPLES - 1); // Should be valid for any even number of time samples, not just ANITA_TIME_SAMPLES.
        // For your records:
        // ZhsTimeDelta * unit::ns / TimeDelta = 0.778477.
        // TimeDelta = 3.85368e-10 s.
        state->grIFft->SetPoint(i, i * TimeDelta / unit::ns + state->vis_xmin_bin * ZhsTimeDelta + ZhsTimeStart, state->ZhsIFft[i] * ZhsTimeDelta * unit::ns / TimeDelta);
      }
    }
    // GDB: plot1d_opt *state->grIFft->GetY()@512 'with lines'

    state->grZhsTimeERec.reset(new TGraph(state->vis_nbins));

    state->grIFft->SetLineColor(kBlue);
    state->grIFft->SetMarkerStyle(5);
    state->grIFft->SetMarkerColor(kRed);
    state->grIFft->SetMarkerSize(2);
    state->grIFft->SetLineWidth(6);
    state->grIFft->Draw("AL");

    state->grINRFft->SetLineColor(kGreen);
    state->grINRFft->SetMarkerStyle(5);
    state->grINRFft->SetMarkerColor(kGreen);
    state->grINRFft->SetMarkerSize(2);
    state->grINRFft->SetLineWidth(2);
    state->grINRFft->Draw("L");

    // Reference waveform shifted such that peaks are split between the first and the last bins:
    for (int i = 0; i < state->vis_nbins; i++) {
      state->grZhsTimeERec->SetPoint(i, ((i + state->BinShift) % state->vis_nbins + state->vis_xmin_bin) * ZhsTimeDelta + ZhsTimeStart, ZhsTimeE[state->vis_xmin_bin + i]);
    }
    state->grZhsTimeERec->SetMarkerSize(2);
    state->grZhsTimeERec->SetMarkerStyle(4);
    state->grZhsTimeERec->SetMarkerColor(kRed);
    state->grZhsTimeERec->Draw("L");
    state->grZhsTimeERec->SetLineColor(kRed);
    // state->grZhsTimeERec->Draw("L");
    state->cZhsIFft->Modified(); state->cZhsIFft->Update(); 
  }
}

void PlotFT(struct cr_ft_state *state, RunMode mode, struct nk_context *ctx){
  if (mode == m_init) {
    if (!bInteractive) return;

    state->cZhsFft   = new TCanvas();
    state->panel_ft_phi = NULL;
    state->panel_ft_rho = NULL;
    return;
  }
  static bvv::TBuffer <int> ShiftPeakToZero(1);
  static bvv::TBuffer <int> cf(cf_polar);
  if (ShiftPeakToZero == 1) {
    state->BinShift = state->vis_nbins - (state->ind_maxval - state->vis_xmin_bin) -1; 
    // printf("state->ind_maxval: %d, vis_xmin_bin: %d\n", state->ind_maxval, state->vis_xmin_bin);
    // printf("ShiftPeakToZero: %d, by %d\n", int(ShiftPeakToZero), int(state->BinShift));
  } else {
    property_int(ctx, "BinShift: ", 0 /*min*/, state->BinShift, +1000 /*max*/, 1 /*increment*/, 1.0 /*sensitivity*/);
  }
  if (mode == m_step) {
    if (bInteractive) {
      nk_layout_row_dynamic(ctx, 30, 2);
      if (nk_option_label(ctx, "Polar", cf == cf_polar )) { cf = cf_polar; }
      if (nk_option_label(ctx, "Rectg", cf == cf_rect))   { cf = cf_rect; }
      nk_layout_row_dynamic(ctx, 30, 2);
      checkbox_label(ctx, "To Zero", ShiftPeakToZero);
    }
  }


  if (mode == m_reload || *cf || *state->BinShift || *state->vis_nbins || *ShiftPeakToZero) {
    // printf("Before state->cZhsFft->Clear()\n");
    // I want to start from fresh: clean canvas cZhsFft, no panels that can or cannot go into cZhsFft.
    // Panels will be _deleted_ (object destruction) as a result of TCanvas.Clear(). Graphs need to be deleted manually.

    state->ZhsFftInp.reset( new double[state->vis_nbins]);
    for (int i = 0; i < state->vis_nbins; i++) {
      state->ZhsFftInp[(i + state->BinShift) % state->vis_nbins] = ZhsTimeE[int(state->vis_xmin_bin) + i];
    }
    state->ZhsFft.reset(FFTtools::doFFT(state->vis_nbins, state->ZhsFftInp.get()));
    cout << "vis_nbins: " << state->vis_nbins << endl;

    if (bInteractive) {
      state->cZhsFft->Clear(); // No panels any more if there were any.
      if (cf == cf_polar) {
        state->grFftRho.reset( new TGraph(state->vis_nbins / 2 + 1));
        state->grFftPhi.reset( new TGraph(state->vis_nbins / 2 + 1));
      } else {
        state->grFftRe.reset(  new TGraph(state->vis_nbins / 2 + 1));
        state->grFftIm.reset(  new TGraph(state->vis_nbins / 2 + 1));
      }
    }
    else { // Non-interactive mode:
      state->FftRho = new double[state->vis_nbins / 2 + 1];
      state->FftPhi = new double[state->vis_nbins / 2 + 1];
    }

    state->dfreq = 1.0 / ((state->vis_xmax_bin - state->vis_xmin_bin) * ZhsTimeDelta * unit::ns);
    int FreqLoCutBin = 200 * 1e6 / state->dfreq + 1;

    // for (int i = 0; i <= (int) FreqLoCutBin; i++) {
    //   cout << i << " re: " << state->ZhsFft[i].re << " im: " << state->ZhsFft[i].im << endl;
    // }

    for (int i = 0; i <= (int) FreqLoCutBin; i++) {
      state->ZhsFft[i].re = 0;
      state->ZhsFft[i].im = 0;
    }
    int FreqHiCutBin = 1200 * 1e6 / state->dfreq;
    for (int i = FreqHiCutBin; i < state->vis_nbins / 2 + 1; i++) {
      state->ZhsFft[i].re = 0;
      state->ZhsFft[i].im = 0;
    }
    for (int i = 0; i < state->vis_nbins / 2 + 1; i++){
      double re = state->ZhsFft[i].re * /* --> */ ZhsTimeDelta  * unit::ns /* <-- to make discrete ft comparable to analytical one */;
      double im = state->ZhsFft[i].im * /* --> */ ZhsTimeDelta  * unit::ns /* <-- to make discrete ft comparable to analytical one */;
      if (bInteractive) {
        if (cf == cf_polar) {
          state->grFftRho->SetPoint(i, i * state->dfreq, TMath::Sqrt(im * im + re * re));
          // state->grFftPhi->SetPoint(i, i * state->dfreq, TMath::ATan2(re, -im) * 180.0 / TMath::Pi());
          double phi = TMath::ATan2(-im, re) * 180.0 / TMath::Pi();
          if (phi < 0) phi = phi + 360;
          state->grFftPhi->SetPoint(i, i * state->dfreq, phi);
          // cout << i << " Rho: " << TMath::Sqrt(im * im + re * re) << " phi: " << phi << endl;
        }
        else {
          state->grFftRe->SetPoint(i, i * state->dfreq, re);
          state->grFftIm->SetPoint(i, i * state->dfreq, im);
        }
      }
      else { // Non-interactive:
        state->FftRho[i] = TMath::Sqrt(im * im + re * re);
        // IceMC's "realft" apparently uses the opposite sign of the power of the exponent:
        double phi = TMath::ATan2(-im, re) * 180.0 / TMath::Pi();
        if (phi < 0) phi = phi + 360;
        state->FftPhi[i] = phi;
      }
    }

    // The rest of this function serves visualization only:
    if (!bInteractive) {
      return;
    } else {
      // At this point, we have all plots that are chosen for
      // visualization and even those that are not needed for the
      // visualization (to simplify logic).  Lets visualize those that are
      // chosen.
      printf("*cf, *state->vis_nbins, state->vis_nbins: %d, %d, %d\n", *cf, *state->vis_nbins, int(state->vis_nbins));
      state->cZhsFft->cd();
      if (cf == cf_rect) {
        state->grFftRe->Draw("AL");
        state->grFftIm->Draw("L");

        state->grFftIm->SetLineWidth(2);
        state->grFftIm->SetLineColor(kMagenta);

        state->grFftRe->SetLineWidth(2);
        state->grFftRe->SetLineColor(kBlue);
      } else {
        // Memory is conserved because of the previous call to state->cZhsFft->Clear():
        state->panel_ft_rho = new TPad1("ft_rho","",0,0,1,1);
        state->panel_ft_rho->cd();
        state->grFftRho->SetLineWidth(2);
        state->grFftRho->SetLineColor(kBlue);
        state->grFftRho->Draw("AL");
        state->panel_ft_rho->Modified();
        state->cZhsFft->cd();
        state->panel_ft_rho->Draw();

        // Memory is conserved because of the previous call to state->cZhsFft->Clear().
        state->panel_ft_phi = new TPad1("ft_phi","",0,0,1,1);
        state->panel_ft_phi->SetFillStyle(4000);
        state->panel_ft_phi->SetFrameFillStyle(0);
        state->panel_ft_phi->cd();
        state->grFftPhi->SetLineWidth(2);
        state->grFftPhi->SetLineColor(kMagenta);
        state->grFftPhi->Draw("ALY+");
        state->panel_ft_phi->Modified();
        state->cZhsFft->cd();
        state->panel_ft_phi->Draw();
      }
      state->cZhsFft->Modified(); state->cZhsFft->Update(); 
    }
  }
}

void PlotWaveform(struct cr_ft_state *state, RunMode mode, struct nk_context *ctx){
  if (mode == m_init) {
    bool fwhm_res = FWHM(ZhsTimeN, ZhsTimeArr.data(), ZhsTimeE.data(), state->fwhm_xmin, state->fwhm_xmax, state->ind_maxval, state->fwhm_xmaxval, state->fwhm_ymaxval);
    if (!fwhm_res) printf("FWHM was not successful\n");
    // printf("PlotWaveform: state->ind_maxval = %d\n", state->ind_maxval);
    if (!bInteractive) return;

    state->cZhsEAndAlpha = new TCanvas();
    state->px1 = new TPad("px1","",0,0,1,1);
    state->px1->Draw();
    state->px1->cd();
    state->grZhsTimeE = new TGraph(ZhsTimeN, ZhsTimeArr.data(), ZhsTimeE.data());
    state->grZhsTimeE->Draw("AL");

    // This update is needed here to get coordinates for the lines marking FWHM:
    gPad->Update(); 
    state->line_min = new TLine(state->fwhm_xmin, gPad->GetUymin(), state->fwhm_xmin, gPad->GetUymax());
    state->line_max = new TLine(state->fwhm_xmax, gPad->GetUymin(), state->fwhm_xmax, gPad->GetUymax());
    state->line_min->Draw();
    state->line_max->Draw();

    // Transparent panel on the top of px1:
    state->px2 = new TPad("px2","",0,0,1,1);
    state->px2->SetFillStyle(4000);
    state->px2->SetFrameFillStyle(0);
    state->px2->Draw();
    state->px2->cd();

    state->grZhsAlpha = new TGraph(ZhsTimeN, ZhsTimeArr.data(), ZhsAlpha.data());
  
    state->grZhsAlpha->Draw("ALY+");
    state->legend = new TLegend(0.1,0.8,0.3,0.9);
    state->legend->AddEntry(state->grZhsTimeE, "E_{proj}","l");
    state->legend->AddEntry(state->grZhsAlpha, "pol. angle","l");
    state->legend->Draw();
    return;
  }

  if (mode == m_reload && bInteractive) {
    TGaxis::SetExponentOffset(0.02, -0.04, "x");

    state->grZhsTimeE->SetLineColor(kRed);
    state->grZhsTimeE->SetLineWidth(2);
    state->grZhsTimeE->GetXaxis()->SetTitle("Time [ns]");
    state->grZhsTimeE->GetYaxis()->SetTitle("E [V/m]");

    state->grZhsAlpha->SetLineColor(kBlue);
    state->grZhsAlpha->SetLineWidth(2);
    state->grZhsAlpha->GetHistogram()->SetMaximum(180);
    state->grZhsAlpha->GetHistogram()->SetMinimum(0);
    state->legend->SetX1NDC(0.12);
    state->legend->SetX2NDC(0.3);
    state->legend->SetY1NDC(0.75);
  }

  static bvv::TBuffer <int> ZhsTimePlotHalfWidth(413);
  static bvv::TBuffer <double> ZhsTimePlotY2NDC(0.85);
  if (mode == m_step && bInteractive) {
    nk_layout_row_dynamic(ctx, 25, 1);
    property_int(ctx, "Plt wid: ", 0, ZhsTimePlotHalfWidth, 3000 /*max*/, 10 /*increment*/, 1.0 /*sensitivity*/);
    //    property_double(ctx, "Y2NDC: ", 0.0, ZhsTimePlotY2NDC, 1.0 /*max*/, 0.02 /*increment*/, 0.002 /*sensitivity*/);
  }


  if (mode == m_reload || *ZhsTimePlotHalfWidth || *ZhsTimePlotY2NDC) {
    // vis_xmin, vis_xmax: which range to visualize, should encompass (fwhm_xmin, fwhm_xmax):
    state->vis_xmin = state->fwhm_xmaxval - ZhsTimePlotHalfWidth*(state->fwhm_xmaxval - state->fwhm_xmin);
    state->vis_xmax = state->fwhm_xmaxval + ZhsTimePlotHalfWidth*(state->fwhm_xmax - state->fwhm_xmaxval); 

    state->vis_xmin_bin = (state->vis_xmin - ZhsTimeStart) / ZhsTimeDelta;
    state->vis_xmax_bin = (state->vis_xmax - ZhsTimeStart) / ZhsTimeDelta;
    // state->vis_nbins = state->vis_xmax_bin - state->vis_xmin_bin + 1;
    state->vis_nbins.set_val(state->vis_xmax_bin - state->vis_xmin_bin + 1);

    if (bInteractive) {
      state->grZhsTimeE->GetXaxis()->SetLimits(state->vis_xmin, state->vis_xmax);
      state->grZhsAlpha->GetXaxis()->SetLimits(state->vis_xmin, state->vis_xmax);

      state->legend->SetY2NDC(ZhsTimePlotY2NDC);

      state->px1->Modified();
      state->px2->Modified();
      state->cZhsEAndAlpha->Modified();
      state->cZhsEAndAlpha->Update();
    }
  }

  if (mode == m_step && bInteractive) {
    // To detect "unmodified" state of vis_nbins.
    // Without this statement, vis_nbins will stay in "modified" state.
    state->vis_nbins.modified = ZhsTimePlotHalfWidth.modified;
  }

  fflush(stdout);
  
}

static struct cr_ft_state *game_init(bool bInteractive_arg)
{
  bInteractive = bInteractive_arg;

  // struct cr_ft_state *state = (cr_ft_state *) malloc(sizeof(*state));
  struct cr_ft_state *state = new cr_ft_state;
  PlotWaveform(state, m_init, NULL);

  PlotFT(state, m_init, NULL);

  if (bInteractive) PlotIFT(state, m_init, NULL);

  printf("Init Done\n");

  return state;
}

static void game_finalize(struct cr_ft_state *state)
{
  delete state;
}

static void game_reload(struct cr_ft_state *state)
{
  PlotWaveform(state, m_reload, ctx);
  PlotFT(state, m_reload, ctx);
  if (bInteractive) PlotIFT(state, m_reload, ctx);
  
  cout << "reloaded dl" << endl;
  
}

static void game_unload(struct cr_ft_state *state __attribute__ ((unused)))
{
  cout << "game unloaded" << endl;
}

static bool game_step_core(struct cr_ft_state *state) {
  PlotWaveform(state, m_step, ctx);
  PlotFT(state, m_step, ctx);
  if (bInteractive) PlotIFT(state, m_step, ctx);
  return true;
}

static bool game_step(struct cr_ft_state *state)
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

  // Try something along these lines as NEXT STEP (causes segfault without changing):
  // *gxlib = xlib;
    
  return true;
}

__attribute__((visibility("default"))) extern const struct game_api GAME_API = {
  (void *(*)(bool))  game_init, // game_init
  (void (*)(void *)) game_finalize, // game_finalize
  (void (*)(void *)) game_reload, // game_reload
  (void (*)(void *)) game_unload, // game_unload
  (bool (*)(void *)) game_step // game_step
};