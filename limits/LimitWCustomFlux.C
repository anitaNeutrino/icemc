void LimitWCustomFlux()
{
//=========Macro generated from canvas: cConst_2/A Simple Graph Example
//=========  (Wed Sep  9 14:49:38 2020) by ROOT version 6.13/08
   TCanvas *cConst_2 = new TCanvas("cConst_2", "A Simple Graph Example",1966,104,800,800);
   cConst_2->Range(17.39948,-19.75,21.40006,-12.25);
   cConst_2->SetFillColor(0);
   cConst_2->SetBorderMode(0);
   cConst_2->SetBorderSize(2);
   cConst_2->SetLogx();
   cConst_2->SetLogy();
   cConst_2->SetLeftMargin(0.15);
   cConst_2->SetFrameBorderMode(0);
   cConst_2->SetFrameBorderMode(0);
   
   Double_t Graph0_fx1[79] = {
   2.710816e+14,
   3.906609e+14,
   5.810319e+14,
   8.550667e+14,
   1.325562e+15,
   1.930634e+15,
   2.506686e+15,
   3.538344e+15,
   5.207149e+15,
   7.989139e+15,
   1.264736e+16,
   2.023485e+16,
   3.204055e+16,
   4.715199e+16,
   6.795165e+16,
   9.691701e+16,
   1.471635e+17,
   2.077304e+17,
   3.057032e+17,
   4.269726e+17,
   6.348922e+17,
   1.139462e+18,
   2.293507e+18,
   4.030881e+18,
   5.931984e+18,
   8.822669e+18,
   1.325562e+19,
   1.992049e+19,
   2.782275e+19,
   3.927354e+19,
   5.601444e+19,
   7.348522e+19,
   1.00531e+20,
   1.291517e+20,
   1.694338e+20,
   2.154269e+20,
   2.654606e+20,
   3.271147e+20,
   3.906609e+20,
   4.813933e+20,
   5.810319e+20,
   7.234358e+20,
   5.152286e+19,
   4.739145e+19,
   4.314197e+19,
   3.688926e+19,
   3.088873e+19,
   2.506686e+19,
   2.165708e+19,
   1.757519e+19,
   1.426264e+19,
   1.194263e+19,
   9.007411e+18,
   6.449113e+18,
   4.865192e+18,
   3.631617e+18,
   2.546244e+18,
   1.823056e+18,
   1.200605e+18,
   7.82348e+17,
   4.641945e+17,
   2.993643e+17,
   2.055417e+17,
   1.411562e+17,
   8.729714e+16,
   6.185859e+16,
   3.519654e+16,
   2.023485e+16,
   1.418731e+16,
   9.055241e+15,
   6.153185e+15,
   3.727348e+15,
   2.305685e+15,
   1.456465e+15,
   1.010416e+15,
   7.387543e+14,
   5.457579e+14,
   3.747141e+14,
   2.768216e+14};
   Double_t Graph0_fy1[79] = {
   6.085552e-15,
   6.549377e-15,
   6.235912e-15,
   5.516962e-15,
   4.762116e-15,
   3.818563e-15,
   3.255367e-15,
   2.708944e-15,
   1.898017e-15,
   1.250835e-15,
   7.202778e-16,
   4.047622e-16,
   2.138947e-16,
   1.309785e-16,
   1.076465e-16,
   9.63829e-17,
   1e-16,
   1.063163e-16,
   1.130577e-16,
   1.089683e-16,
   9.756634e-17,
   7.181249e-17,
   3.842377e-17,
   2.212585e-17,
   1.354565e-17,
   8.396533e-18,
   5.333349e-18,
   3.305978e-18,
   2.024417e-18,
   1.194538e-18,
   6.878599e-19,
   4.317178e-19,
   2.425493e-19,
   1.379749e-19,
   7.473086e-20,
   4.199523e-20,
   2.388911e-20,
   1.309785e-20,
   7.269424e-21,
   3.235937e-21,
   1.668784e-21,
   6.732866e-22,
   6.732866e-22,
   8.818607e-22,
   1.338135e-21,
   2.532212e-21,
   4.970499e-21,
   1e-20,
   1.736601e-20,
   3.408788e-20,
   6.217274e-20,
   1.002998e-19,
   2.172201e-19,
   4.703271e-19,
   8.578276e-19,
   1.564588e-18,
   2.684108e-18,
   4.072865e-18,
   6.256047e-18,
   7.99466e-18,
   1.009253e-17,
   1.099765e-17,
   1.099765e-17,
   1.034189e-17,
   1.009253e-17,
   1.059986e-17,
   1.4054e-17,
   2.212585e-17,
   2.969614e-17,
   4.343102e-17,
   5.55009e-17,
   7.269424e-17,
   1e-16,
   1.172735e-16,
   1.325867e-16,
   1.480471e-16,
   1.554891e-16,
   1.480471e-16,
   1.409613e-16};
   TGraph *graph = new TGraph(79,Graph0_fx1,Graph0_fy1);
   graph->SetName("Graph0");
   graph->SetTitle(";E (eV);E dN/dE dA d#Omega dt (cm^{-2} sr ^{-1} s^{-1} )");
   graph->SetFillColor(15);
   graph->SetLineColor(0);
   
   TH1F *Graph_Graph01 = new TH1F("Graph_Graph01","",100,9.99e+17,1e+21);
   Graph_Graph01->SetMinimum(1e-19);
   Graph_Graph01->SetMaximum(1e-13);
   Graph_Graph01->SetDirectory(0);
   Graph_Graph01->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph01->SetLineColor(ci);
   Graph_Graph01->GetXaxis()->SetTitle("E (eV)");
   Graph_Graph01->GetXaxis()->SetLabelFont(42);
   Graph_Graph01->GetXaxis()->SetLabelOffset(0.006);
   Graph_Graph01->GetXaxis()->SetTitleOffset(1.3);
   Graph_Graph01->GetXaxis()->SetTitleFont(42);
   Graph_Graph01->GetYaxis()->SetTitle("E dN/dE dA d#Omega dt (cm^{-2} sr ^{-1} s^{-1} )");
   Graph_Graph01->GetYaxis()->SetLabelFont(42);
   Graph_Graph01->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph01->GetYaxis()->SetTitleOffset(1.8);
   Graph_Graph01->GetYaxis()->SetTitleFont(42);
   Graph_Graph01->GetZaxis()->SetLabelFont(42);
   Graph_Graph01->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph01->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph01);
   
   graph->Draw("af");
   
   TH1F *Graph_copy = new TH1F("Graph_copy","",100,9.99e+17,1e+21);
   Graph_copy->SetMinimum(1e-19);
   Graph_copy->SetMaximum(1e-13);
   Graph_copy->SetDirectory(0);
   Graph_copy->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_copy->SetLineColor(ci);
   Graph_copy->GetXaxis()->SetTitle("E (eV)");
   Graph_copy->GetXaxis()->SetLabelFont(42);
   Graph_copy->GetXaxis()->SetLabelOffset(0.006);
   Graph_copy->GetXaxis()->SetTitleOffset(1.3);
   Graph_copy->GetXaxis()->SetTitleFont(42);
   Graph_copy->GetYaxis()->SetTitle("E dN/dE dA d#Omega dt (cm^{-2} sr ^{-1} s^{-1} )");
   Graph_copy->GetYaxis()->SetLabelFont(42);
   Graph_copy->GetYaxis()->SetLabelOffset(0.007);
   Graph_copy->GetYaxis()->SetTitleOffset(1.8);
   Graph_copy->GetYaxis()->SetTitleFont(42);
   Graph_copy->GetZaxis()->SetLabelFont(42);
   Graph_copy->GetZaxis()->SetLabelSize(0.035);
   Graph_copy->GetZaxis()->SetTitleSize(0.035);
   Graph_copy->GetZaxis()->SetTitleFont(42);
   Graph_copy->Draw("sameaxis");
   
   Double_t Graph1_fx2[32] = {
   1.021175e+15,
   1.871113e+15,
   3.538344e+15,
   6.62064e+15,
   9.151663e+15,
   1.360818e+16,
   1.961101e+16,
   3.237427e+16,
   5.514422e+16,
   9.19814e+16,
   1.325562e+17,
   1.930634e+17,
   2.696497e+17,
   4.094492e+17,
   6.483358e+17,
   9.844644e+17,
   1.479449e+18,
   2.002627e+18,
   3.137618e+18,
   5.126253e+18,
   7.946941e+18,
   1.194263e+19,
   1.871113e+19,
   2.641192e+19,
   3.612435e+19,
   5.152286e+19,
   8.07235e+19,
   1.05901e+20,
   1.463862e+20,
   1.86123e+20,
   2.269865e+20,
   2.826181e+20};
   Double_t Graph1_fy2[32] = {
   7.634841e-17,
   6.590221e-17,
   5.61953e-17,
   4.850651e-17,
   4.506091e-17,
   4.506091e-17,
   4.506091e-17,
   4.396428e-17,
   4.451435e-17,
   4.343102e-17,
   4.451435e-17,
   4.732602e-17,
   4.970499e-17,
   4.617427e-17,
   4.186007e-17,
   3.196686e-17,
   2.239752e-17,
   1.569278e-17,
   9.492924e-18,
   4.894406e-18,
   2.618786e-18,
   1.526863e-18,
   7.495489e-19,
   4.534193e-19,
   2.396073e-19,
   1.147889e-19,
   4.576147e-20,
   2.061579e-20,
   7.269424e-21,
   3.081059e-21,
   1.440456e-21,
   6.817105e-22};
   graph = new TGraph(32,Graph1_fx2,Graph1_fy2);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetLineStyle(9);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","Graph",100,9.190571e+14,3.108798e+20);
   Graph_Graph12->SetMinimum(6.135395e-22);
   Graph_Graph12->SetMaximum(8.398318e-17);
   Graph_Graph12->SetDirectory(0);
   Graph_Graph12->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph12->SetLineColor(ci);
   Graph_Graph12->GetXaxis()->SetLabelFont(42);
   Graph_Graph12->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleFont(42);
   Graph_Graph12->GetYaxis()->SetLabelFont(42);
   Graph_Graph12->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleOffset(0);
   Graph_Graph12->GetYaxis()->SetTitleFont(42);
   Graph_Graph12->GetZaxis()->SetLabelFont(42);
   Graph_Graph12->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph12);
   
   graph->Draw("l");
   
   Double_t Graph2_fx3[11] = {
   1e+16,
   3.162278e+16,
   1e+17,
   3.162278e+17,
   1e+18,
   3.162278e+18,
   1e+19,
   3.162278e+19,
   1e+20,
   3.162278e+20,
   1e+21};
   Double_t Graph2_fy3[11] = {
   1.360912e-15,
   7.909882e-16,
   1.234099e-15,
   1.675714e-15,
   6.4217e-16,
   1.764371e-16,
   5.992941e-17,
   3.185884e-17,
   1.382038e-17,
   4.135043e-18,
   7.773561e-19};
   graph = new TGraph(11,Graph2_fx3,Graph2_fy3);
   graph->SetName("Graph2");
   graph->SetTitle(";log_{10} #left(#frac{E_{#nu}}{eV}#right);EF(E), cm^{-2} s^{-1} str^{-1}");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#ff00ff");
   graph->SetLineColor(ci);
   graph->SetLineStyle(9);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph23 = new TH1F("Graph_Graph23","",100,9e+15,1.099999e+21);
   Graph_Graph23->SetMinimum(6.996205e-19);
   Graph_Graph23->SetMaximum(1.843208e-15);
   Graph_Graph23->SetDirectory(0);
   Graph_Graph23->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph23->SetLineColor(ci);
   Graph_Graph23->GetXaxis()->SetTitle("log_{10} #left(#frac{E_{#nu}}{eV}#right)");
   Graph_Graph23->GetXaxis()->SetLabelFont(42);
   Graph_Graph23->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph23->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph23->GetXaxis()->SetTitleFont(42);
   Graph_Graph23->GetYaxis()->SetTitle("EF(E), cm^{-2} s^{-1} str^{-1}");
   Graph_Graph23->GetYaxis()->SetLabelFont(42);
   Graph_Graph23->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph23->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph23->GetYaxis()->SetTitleOffset(0);
   Graph_Graph23->GetYaxis()->SetTitleFont(42);
   Graph_Graph23->GetZaxis()->SetLabelFont(42);
   Graph_Graph23->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph23->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph23->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph23);
   
   graph->Draw("l");
   
   Double_t Graph3_fx4[52] = {
   8.912509e+16,
   1.122018e+17,
   1.412538e+17,
   1.778279e+17,
   2.238721e+17,
   2.818383e+17,
   3.548134e+17,
   4.466836e+17,
   5.623413e+17,
   7.079458e+17,
   8.912509e+17,
   1.122018e+18,
   1.412538e+18,
   1.778279e+18,
   2.238721e+18,
   2.818383e+18,
   3.548134e+18,
   4.466836e+18,
   5.623413e+18,
   7.079458e+18,
   8.912509e+18,
   1.122018e+19,
   1.412538e+19,
   1.778279e+19,
   2.238721e+19,
   2.818383e+19,
   3.548134e+19,
   4.466836e+19,
   5.623413e+19,
   7.079458e+19,
   8.912509e+19,
   1.122018e+20,
   1.412538e+20,
   1.778279e+20,
   2.238721e+20,
   2.818383e+20,
   3.548134e+20,
   4.466836e+20,
   5.623413e+20,
   7.079458e+20,
   8.912509e+20,
   1.122018e+21,
   1.412538e+21,
   1.778279e+21,
   2.238721e+21,
   2.818383e+21,
   3.548134e+21,
   4.466836e+21,
   5.623413e+21,
   7.079458e+21,
   8.912509e+21,
   1.122018e+22};
   Double_t Graph3_fy4[52] = {
   4.59198e-17,
   5.248075e-17,
   6.053409e-17,
   6.839116e-17,
   7.638358e-17,
   8.433348e-17,
   9.036495e-17,
   9.484185e-17,
   9.549926e-17,
   9.289664e-17,
   8.749838e-17,
   7.961594e-17,
   7.063176e-17,
   6.095369e-17,
   5.140437e-17,
   4.295364e-17,
   3.580964e-17,
   2.951209e-17,
   2.42661e-17,
   1.981527e-17,
   1.648162e-17,
   1.342765e-17,
   1.111732e-17,
   9.120108e-18,
   7.430191e-18,
   5.861382e-18,
   4.830588e-18,
   3.775722e-18,
   2.951209e-18,
   2.213095e-18,
   1.721869e-18,
   1.294196e-18,
   9.977001e-19,
   7.603263e-19,
   5.754399e-19,
   4.265795e-19,
   3.206269e-19,
   2.208005e-19,
   1.770109e-19,
   1.282331e-19,
   7.726806e-20,
   5.794287e-20,
   4.246196e-20,
   2.844461e-20,
   1.399587e-20,
   1.34586e-20,
   5.432503e-21,
   4.74242e-21,
   2.654606e-21,
   7.726806e-22,
   8.810489e-22,
   7.638358e-22};
   graph = new TGraph(52,Graph3_fx4,Graph3_fy4);
   graph->SetName("Graph3");
   graph->SetTitle(";log_{10} #left(#frac{E_{#nu}}{eV}#right);EF(E), cm^{-2} s^{-1} str^{-1}");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#00cc00");
   graph->SetLineColor(ci);
   graph->SetLineStyle(9);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph34 = new TH1F("Graph_Graph34","",100,8.021258e+16,1.234219e+22);
   Graph_Graph34->SetMinimum(6.874522e-22);
   Graph_Graph34->SetMaximum(1.050491e-16);
   Graph_Graph34->SetDirectory(0);
   Graph_Graph34->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph34->SetLineColor(ci);
   Graph_Graph34->GetXaxis()->SetTitle("log_{10} #left(#frac{E_{#nu}}{eV}#right)");
   Graph_Graph34->GetXaxis()->SetLabelFont(42);
   Graph_Graph34->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph34->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph34->GetXaxis()->SetTitleFont(42);
   Graph_Graph34->GetYaxis()->SetTitle("EF(E), cm^{-2} s^{-1} str^{-1}");
   Graph_Graph34->GetYaxis()->SetLabelFont(42);
   Graph_Graph34->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph34->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph34->GetYaxis()->SetTitleOffset(0);
   Graph_Graph34->GetYaxis()->SetTitleFont(42);
   Graph_Graph34->GetZaxis()->SetLabelFont(42);
   Graph_Graph34->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph34->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph34->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph34);
   
   graph->Draw("l");
   
   Double_t Graph4_fx5[31] = {
   3.819835e+15,
   6.183073e+15,
   1.040372e+16,
   1.752584e+16,
   2.826956e+16,
   4.460515e+16,
   6.240037e+16,
   8.240834e+16,
   1.088257e+17,
   1.40649e+17,
   1.857262e+17,
   2.452772e+17,
   3.010526e+17,
   4.614069e+17,
   7.603215e+17,
   1.266472e+18,
   2.055144e+18,
   3.291979e+18,
   5.599494e+18,
   9.536065e+18,
   1.523234e+19,
   2.111381e+19,
   2.729611e+19,
   3.368756e+19,
   3.968064e+19,
   4.566709e+19,
   5.256539e+19,
   6.051548e+19,
   6.805568e+19,
   7.656307e+19,
   8.413952e+19};
   Double_t Graph4_fy5[31] = {
   1.233136e-19,
   1.892336e-19,
   2.553814e-19,
   3.003688e-19,
   3.526557e-19,
   4.82996e-19,
   7.036039e-19,
   1.035779e-18,
   1.551814e-18,
   2.366169e-18,
   3.607879e-18,
   5.311178e-18,
   6.999286e-18,
   1.092055e-17,
   1.593044e-17,
   2.000996e-17,
   2.206393e-17,
   2.206393e-17,
   1.989488e-17,
   1.537349e-17,
   1.028017e-17,
   6.656502e-18,
   4.329982e-18,
   2.799711e-18,
   1.929836e-18,
   1.328128e-18,
   8.666731e-19,
   5.369284e-19,
   3.534943e-19,
   2.071843e-19,
   1.29567e-19};
   graph = new TGraph(31,Graph4_fx5,Graph4_fy5);
   graph->SetName("Graph4");
   graph->SetTitle("takami_cosmogenic.csv");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#ffcc00");
   graph->SetLineColor(ci);
   graph->SetLineStyle(9);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph45 = new TH1F("Graph_Graph45","takami_cosmogenic.csv",100,3.437852e+15,9.255309e+19);
   Graph_Graph45->SetMinimum(1.109823e-19);
   Graph_Graph45->SetMaximum(2.425799e-17);
   Graph_Graph45->SetDirectory(0);
   Graph_Graph45->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph45->SetLineColor(ci);
   Graph_Graph45->GetXaxis()->SetLabelFont(42);
   Graph_Graph45->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph45->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph45->GetXaxis()->SetTitleFont(42);
   Graph_Graph45->GetYaxis()->SetLabelFont(42);
   Graph_Graph45->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph45->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph45->GetYaxis()->SetTitleOffset(0);
   Graph_Graph45->GetYaxis()->SetTitleFont(42);
   Graph_Graph45->GetZaxis()->SetLabelFont(42);
   Graph_Graph45->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph45->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph45->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph45);
   
   graph->Draw("l");
   
   Double_t Graph5_fx6[7] = {
   1e+18,
   3.162278e+18,
   1e+19,
   3.162278e+19,
   1e+20,
   3.162278e+20,
   1e+21};
   Double_t Graph5_fy6[7] = {
   9.513554e-14,
   6.401941e-15,
   1.697179e-16,
   2.612106e-17,
   4.846184e-18,
   1.571316e-18,
   6.750369e-19};
   graph = new TGraph(7,Graph5_fx6,Graph5_fy6);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#990000");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph56 = new TH1F("Graph_Graph56","Graph",100,9e+17,1.0999e+21);
   Graph_Graph56->SetMinimum(6.075332e-19);
   Graph_Graph56->SetMaximum(1.04649e-13);
   Graph_Graph56->SetDirectory(0);
   Graph_Graph56->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph56->SetLineColor(ci);
   Graph_Graph56->GetXaxis()->SetLabelFont(42);
   Graph_Graph56->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph56->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph56->GetXaxis()->SetTitleFont(42);
   Graph_Graph56->GetYaxis()->SetLabelFont(42);
   Graph_Graph56->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph56->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph56->GetYaxis()->SetTitleOffset(0);
   Graph_Graph56->GetYaxis()->SetTitleFont(42);
   Graph_Graph56->GetZaxis()->SetLabelFont(42);
   Graph_Graph56->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph56->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph56->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph56);
   
   graph->Draw("l");
   
   Double_t Graph6_fx7[7] = {
   1e+18,
   3.162278e+18,
   1e+19,
   3.162278e+19,
   1e+20,
   3.162278e+20,
   1e+21};
   Double_t Graph6_fy7[7] = {
   3.0372e-14,
   1.148099e-15,
   6.003807e-17,
   7.816794e-18,
   1.427859e-18,
   4.83258e-19,
   2.005803e-19};
   graph = new TGraph(7,Graph6_fx7,Graph6_fy7);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#00ffff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph67 = new TH1F("Graph_Graph67","Graph",100,9e+17,1.0999e+21);
   Graph_Graph67->SetMinimum(1.805222e-19);
   Graph_Graph67->SetMaximum(3.340917e-14);
   Graph_Graph67->SetDirectory(0);
   Graph_Graph67->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph67->SetLineColor(ci);
   Graph_Graph67->GetXaxis()->SetLabelFont(42);
   Graph_Graph67->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph67->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph67->GetXaxis()->SetTitleFont(42);
   Graph_Graph67->GetYaxis()->SetLabelFont(42);
   Graph_Graph67->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph67->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph67->GetYaxis()->SetTitleOffset(0);
   Graph_Graph67->GetYaxis()->SetTitleFont(42);
   Graph_Graph67->GetZaxis()->SetLabelFont(42);
   Graph_Graph67->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph67->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph67->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph67);
   
   graph->Draw("l");
   
   Double_t Graph7_fx8[7] = {
   1e+18,
   3.162278e+18,
   1e+19,
   3.162278e+19,
   1e+20,
   3.162278e+20,
   1e+21};
   Double_t Graph7_fy8[7] = {
   3.87842e-14,
   1.689472e-15,
   1.222469e-16,
   1.685489e-17,
   2.525102e-18,
   7.614465e-19,
   3.117629e-19};
   graph = new TGraph(7,Graph7_fx8,Graph7_fy8);
   graph->SetName("Graph7");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph78 = new TH1F("Graph_Graph78","Graph",100,9e+17,1.0999e+21);
   Graph_Graph78->SetMinimum(2.805866e-19);
   Graph_Graph78->SetMaximum(4.266258e-14);
   Graph_Graph78->SetDirectory(0);
   Graph_Graph78->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph78->SetLineColor(ci);
   Graph_Graph78->GetXaxis()->SetLabelFont(42);
   Graph_Graph78->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph78->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph78->GetXaxis()->SetTitleFont(42);
   Graph_Graph78->GetYaxis()->SetLabelFont(42);
   Graph_Graph78->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph78->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph78->GetYaxis()->SetTitleOffset(0);
   Graph_Graph78->GetYaxis()->SetTitleFont(42);
   Graph_Graph78->GetZaxis()->SetLabelFont(42);
   Graph_Graph78->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph78->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph78->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph78);
   
   graph->Draw("l");
   
   Double_t Graph8_fx9[30] = {
   5.905432e+16,
   6.83567e+16,
   8.345808e+16,
   1.005433e+17,
   1.227403e+17,
   1.539021e+17,
   1.781379e+17,
   2.175631e+17,
   2.656704e+17,
   3.376832e+17,
   4.015765e+17,
   4.904127e+17,
   5.754165e+17,
   7.614568e+17,
   9.681733e+17,
   1.18264e+18,
   1.425558e+18,
   1.788432e+18,
   2.463951e+18,
   3.091649e+18,
   3.984641e+18,
   5.273791e+18,
   7.562943e+18,
   1.084529e+19,
   1.707487e+19,
   2.384366e+19,
   3.329571e+19,
   4.775395e+19,
   6.84877e+19,
   9.821559e+19};
   Double_t Graph8_fy9[30] = {
   5.574434e-15,
   3.898582e-15,
   2.560025e-15,
   1.680985e-15,
   1.036244e-15,
   6.524425e-16,
   4.467882e-16,
   3.472245e-16,
   2.480467e-16,
   1.772187e-16,
   1.4364e-16,
   1.070265e-16,
   8.316621e-17,
   6.602495e-17,
   5.582842e-17,
   4.720082e-17,
   4.162161e-17,
   3.445895e-17,
   3.039823e-17,
   2.737897e-17,
   2.572255e-17,
   2.221574e-17,
   1.960018e-17,
   1.693218e-17,
   1.373569e-17,
   1.290784e-17,
   1.212989e-17,
   1.139975e-17,
   1.049028e-17,
   9.255216e-18};
   graph = new TGraph(30,Graph8_fx9,Graph8_fy9);
   graph->SetName("Graph8");
   graph->SetTitle("auger_2017.dat");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph89 = new TH1F("Graph_Graph89","auger_2017.dat",100,5.314889e+16,1.080312e+20);
   Graph_Graph89->SetMinimum(8.329694e-18);
   Graph_Graph89->SetMaximum(6.130952e-15);
   Graph_Graph89->SetDirectory(0);
   Graph_Graph89->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph89->SetLineColor(ci);
   Graph_Graph89->GetXaxis()->SetLabelFont(42);
   Graph_Graph89->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph89->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph89->GetXaxis()->SetTitleFont(42);
   Graph_Graph89->GetYaxis()->SetLabelFont(42);
   Graph_Graph89->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph89->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph89->GetYaxis()->SetTitleOffset(0);
   Graph_Graph89->GetYaxis()->SetTitleFont(42);
   Graph_Graph89->GetZaxis()->SetLabelFont(42);
   Graph_Graph89->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph89->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph89->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph89);
   
   graph->Draw("l");
   
   Double_t Graph9_fx10[23] = {
   1.080692e+16,
   1.488644e+16,
   1.996685e+16,
   2.979384e+16,
   4.328505e+16,
   5.503137e+16,
   7.478389e+16,
   1.057696e+17,
   1.709227e+17,
   2.653892e+17,
   3.241247e+17,
   6.067345e+17,
   8.24578e+17,
   1.351149e+18,
   2.273946e+18,
   4.036429e+18,
   5.271214e+18,
   8.299361e+18,
   1.554268e+19,
   2.414568e+19,
   3.238739e+19,
   4.645498e+19,
   7.82143e+19};
   Double_t Graph9_fy10[23] = {
   2.710045e-15,
   2.197542e-15,
   1.858469e-15,
   1.356719e-15,
   9.903518e-16,
   8.028677e-16,
   5.859424e-16,
   4.014947e-16,
   2.325468e-16,
   1.434588e-16,
   1.1149e-16,
   6.597925e-17,
   5.022398e-17,
   3.825274e-17,
   2.91373e-17,
   2.040428e-17,
   1.725457e-17,
   1.429513e-17,
   1.066522e-17,
   8.651468e-18,
   7.472299e-18,
   7.324627e-18,
   6.887098e-18};
   graph = new TGraph(23,Graph9_fx10,Graph9_fy10);
   graph->SetName("Graph9");
   graph->SetTitle("icecube_2018.dat");
   graph->SetFillStyle(1000);

   ci = TColor::GetColor("#0000ff");
   graph->SetLineColor(ci);
   graph->SetLineWidth(3);
   
   TH1F *Graph_Graph910 = new TH1F("Graph_Graph910","icecube_2018.dat",100,9.726227e+15,8.603465e+19);
   Graph_Graph910->SetMinimum(6.198388e-18);
   Graph_Graph910->SetMaximum(2.980361e-15);
   Graph_Graph910->SetDirectory(0);
   Graph_Graph910->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph910->SetLineColor(ci);
   Graph_Graph910->GetXaxis()->SetLabelFont(42);
   Graph_Graph910->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph910->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph910->GetXaxis()->SetTitleFont(42);
   Graph_Graph910->GetYaxis()->SetLabelFont(42);
   Graph_Graph910->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph910->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph910->GetYaxis()->SetTitleOffset(0);
   Graph_Graph910->GetYaxis()->SetTitleFont(42);
   Graph_Graph910->GetZaxis()->SetLabelFont(42);
   Graph_Graph910->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph910->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph910->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph910);
   
   graph->Draw("l");
   
   TLegend *leg = new TLegend(0.5,0.7,0.89,0.89,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph5","ANITA-III","l");

   ci = TColor::GetColor("#990000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph7","ANITA IV","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph6","ANITA I-IV","l");

   ci = TColor::GetColor("#00ffff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph8","Auger 2017","l");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph9","IceCube 2018","l");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   leg = new TLegend(0.16,0.12,0.43,0.254,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph3","Simulation '20","l");

   ci = TColor::GetColor("#00cc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(9);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","KKSS '02 ","l");

   ci = TColor::GetColor("#ff00ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(9);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph4","Takami et al '09","l");

   ci = TColor::GetColor("#ffcc00");
   entry->SetLineColor(ci);
   entry->SetLineStyle(9);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Ahlers '12, E_{min}=10^{18.5} eV","l");
   entry->SetLineColor(43);
   entry->SetLineStyle(9);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph0","GZK, Kotera '10","f");
   entry->SetFillColor(15);
   entry->SetFillStyle(1001);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   cConst_2->Modified();
   cConst_2->cd();
   cConst_2->SetSelected(cConst_2);
}
