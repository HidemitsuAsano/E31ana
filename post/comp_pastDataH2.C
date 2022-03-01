void comp_pastDataH2()
{
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  //gStyle->SetErrorX(0);
  // ==================================== //
  // reference data
  // NPB90(1975)349 @ LBL, bubble chamber
  // d(sigma)/d(Omega)
  // ==================================== //
  // P1001.0 MeV,   REK- P --> PI- SIGMA+,    SQRT(S)1.794 GeV
  // Plot: p7071_d40x1y1
  double p7071_d40x1y1_xval[] = { -0.975, -0.925, -0.875, -0.825, -0.775, -0.7, -0.625, -0.575, -0.525, 
				  -0.475, -0.35, -0.15, 0.0, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 
				  0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 
				  0.875, 0.925, 0.975 };
  double p7071_d40x1y1_xerrminus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				       0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				       0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				       0.025, 0.025, 0.025 };
  double p7071_d40x1y1_xerrplus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				      0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				      0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				      0.025, 0.025, 0.025 };
  double p7071_d40x1y1_yval[] = { 0.215, 0.129, 0.145, 0.236, 0.172, 0.063, 0.086, 0.095, 0.057, 
				  0.083, 0.02, 0.014, 0.061, 0.058, 0.117, 0.109, 0.157, 0.137, 0.341, 
				  0.311, 0.28, 0.324, 0.407, 0.233, 0.381, 0.375, 0.301, 0.257, 0.257, 
				  0.127, 0.297, 0.141 };
  double p7071_d40x1y1_yerrminus[] = { 0.057, 0.039, 0.046, 0.053, 0.05, 0.02, 0.032, 0.034, 0.025, 
				       0.031, 0.008, 0.006, 0.019, 0.026, 0.037, 0.036, 0.044, 0.04, 0.065, 
				       0.062, 0.058, 0.064, 0.072, 0.053, 0.071, 0.071, 0.064, 0.061, 0.062, 
				       0.045, 0.074, 0.053 };
  double p7071_d40x1y1_yerrplus[] = { 0.057, 0.039, 0.046, 0.053, 0.05, 0.02, 0.032, 0.034, 0.025, 
				      0.031, 0.008, 0.006, 0.019, 0.026, 0.037, 0.036, 0.044, 0.04, 0.065, 
				      0.062, 0.058, 0.064, 0.072, 0.053, 0.071, 0.071, 0.064, 0.061, 0.062, 
				      0.045, 0.074, 0.053 };
  double p7071_d40x1y1_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0 };
  double p7071_d40x1y1_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0 };
  int p7071_d40x1y1_numpoints = 32;
  for(int i=0; i<p7071_d40x1y1_numpoints; i++){
    p7071_d40x1y1_yval[i] *= 1000;
    p7071_d40x1y1_yerrminus[i] *= 1000;
    p7071_d40x1y1_yerrplus[i] *= 1000;
  }
  TGraphAsymmErrors *p7071_d40x1y1 = new TGraphAsymmErrors(p7071_d40x1y1_numpoints, p7071_d40x1y1_xval, p7071_d40x1y1_yval, p7071_d40x1y1_xerrminus, p7071_d40x1y1_xerrplus, p7071_d40x1y1_yerrminus, p7071_d40x1y1_yerrplus);
  p7071_d40x1y1->SetName("/HepData/7071/d40x1y1");
  p7071_d40x1y1->SetTitle("/HepData/7071/d40x1y1");
  p7071_d40x1y1->SetLineColor(2);
  p7071_d40x1y1->SetMarkerColor(2);
  p7071_d40x1y1->SetMarkerStyle(24);
  //TCanvas *ctest1 = new TCanvas("ctest1","ctest1");
  //p7071_d40x1y1->Print();
  //p7071_d40x1y1->Draw("AP*");

  //P1001.0 MeV,   REK- P --> PI- SIGMA+,    SQRT(S)1.794 GeV
  // Plot: p7071_d60x1y1
  double p7071_d60x1y1_xval[] = { -0.975, -0.925, -0.875, -0.825, -0.775, -0.725, -0.675, -0.625, -0.575, 
				  -0.525, -0.475, -0.425, -0.375, -0.325, -0.275, -0.225, -0.175, -0.125, -0.075, 
				  -0.025, 0.025, 0.075, 0.125, 0.175, 0.225, 0.3, 0.375, 0.425, 0.5, 
				  0.625, 0.75, 0.85, 0.925, 0.975 };
  double p7071_d60x1y1_xerrminus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				       0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				       0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				       0.025, 0.025, 0.025, 0.025, 0.025 };
  double p7071_d60x1y1_xerrplus[] = { 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				      0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				      0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 0.025, 
				      0.025, 0.025, 0.025, 0.025, 0.025 };
  double p7071_d60x1y1_yval[] = { 0.048, 0.135, 0.309, 0.324, 0.28, 0.393, 0.266, 0.357, 0.164, 
				  0.202, 0.217, 0.143, 0.077, 0.081, 0.068, 0.078, 0.067, 0.068, 0.08, 
				  0.05, 0.08, 0.089, 0.108, 0.059, 0.069, 0.05, 0.061, 0.092, 0.047, 
				  0.018, 0.055, 0.076, 0.187, 0.204 };
  double p7071_d60x1y1_yerrminus[] = { 0.021, 0.036, 0.055, 0.057, 0.052, 0.063, 0.052, 0.059, 0.04, 
				       0.045, 0.046, 0.037, 0.027, 0.029, 0.026, 0.028, 0.025, 0.026, 0.028, 
				       0.022, 0.028, 0.03, 0.033, 0.024, 0.026, 0.016, 0.025, 0.031, 0.016, 
				       0.008, 0.017, 0.021, 0.048, 0.053 };
  double p7071_d60x1y1_yerrplus[] = { 0.021, 0.036, 0.055, 0.057, 0.052, 0.063, 0.052, 0.059, 0.04, 
				      0.045, 0.046, 0.037, 0.027, 0.029, 0.026, 0.028, 0.025, 0.026, 0.028, 
				      0.022, 0.028, 0.03, 0.033, 0.024, 0.026, 0.016, 0.025, 0.031, 0.016, 
				      0.008, 0.017, 0.021, 0.048, 0.053 };
  double p7071_d60x1y1_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0 };
  double p7071_d60x1y1_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0 };
  int p7071_d60x1y1_numpoints = 34;
  for(int i=0; i<p7071_d60x1y1_numpoints; i++){
    p7071_d60x1y1_yval[i] *= 1000;
    p7071_d60x1y1_yerrminus[i] *= 1000;
    p7071_d60x1y1_yerrplus[i] *= 1000;
  }
  TGraphAsymmErrors *p7071_d60x1y1 = new TGraphAsymmErrors(p7071_d60x1y1_numpoints, p7071_d60x1y1_xval, p7071_d60x1y1_yval, p7071_d60x1y1_xerrminus, p7071_d60x1y1_xerrplus, p7071_d60x1y1_yerrminus, p7071_d60x1y1_yerrplus);
  p7071_d60x1y1->SetName("/HepData/7071/d60x1y1");
  p7071_d60x1y1->SetTitle("/HepData/7071/d60x1y1");
  p7071_d60x1y1->SetLineColor(2);
  p7071_d60x1y1->SetMarkerColor(2);
  p7071_d60x1y1->SetMarkerStyle(24);
  //TCanvas *ctest1 = new TCanvas("ctest1","ctest1");
  //p7071_d60x1y1.Draw("AP");


  // ==================================== //
  // reference data
  // NPB105(1976)189 @ CERN, bubble chamber
  // d(sigma)/d(cos)
  // ==================================== //
  // P1005.0 MeV,   REK- P --> PI- SIGMA+
  double p7006_d21x1y2_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				  -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				  0.95 };
  double p7006_d21x1y2_xerrminus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				       0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				       0.05 };
  double p7006_d21x1y2_xerrplus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05 };
  double p7006_d21x1y2_yval[] = { 1.11, 1.41, 0.94, 0.88, 0.62, 0.21, 0.21, 0.05, 0.11, 
				  0.38, 0.38, 0.79, 0.88, 1.11, 1.14, 2.12, 1.94, 1.52, 2.5, 
				  1.21 };
  double p7006_d21x1y2_yerrminus[] = { 0.24, 0.27, 0.22, 0.21, 0.18, 0.1, 0.1, 0.07, 0.08, 
				       0.15, 0.14, 0.21, 0.21, 0.24, 0.25, 0.34, 0.33, 0.31, 0.41, 
				       0.31 };
  double p7006_d21x1y2_yerrplus[] = { 0.24, 0.27, 0.22, 0.21, 0.18, 0.1, 0.1, 0.07, 0.08, 
				      0.15, 0.14, 0.21, 0.21, 0.24, 0.25, 0.34, 0.33, 0.31, 0.41, 
				      0.31 };
  double p7006_d21x1y2_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0 };
  double p7006_d21x1y2_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  int p7006_d21x1y2_numpoints = 20;
  for(int i=0; i<p7006_d21x1y2_numpoints; i++){
    p7006_d21x1y2_yval[i] *= 1000;
    p7006_d21x1y2_yerrminus[i] *= 1000;
    p7006_d21x1y2_yerrplus[i] *= 1000;
    // d(sigma)/d(cos) -> d(sigma)/d(Omega) = *(2./NBIN)/(2pi*2./NBIN) = *(1./2pi)
    p7006_d21x1y2_yval[i] *= 1./(2.*3.1415);
    p7006_d21x1y2_yerrminus[i] *= 1./(2.*3.1415);
    p7006_d21x1y2_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *p7006_d21x1y2 = new TGraphAsymmErrors(p7006_d21x1y2_numpoints, p7006_d21x1y2_xval, p7006_d21x1y2_yval, p7006_d21x1y2_xerrminus, p7006_d21x1y2_xerrplus, p7006_d21x1y2_yerrminus, p7006_d21x1y2_yerrplus);
  p7006_d21x1y2->SetName("/HepData/7006/d21x1y2");
  p7006_d21x1y2->SetTitle("/HepData/7006/d21x1y2");
  p7006_d21x1y2->SetLineColor(4);
  p7006_d21x1y2->SetMarkerColor(4);
  p7006_d21x1y2->SetMarkerStyle(25);
  //p7006_d21x1y2.Draw("AP");


  // P1005.0 MeV,   REK- P --> PI+ SIGMA-
  double p7006_d20x1y2_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				  -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				  0.95 };
  double p7006_d20x1y2_xerrminus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				       0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				       0.05 };
  double p7006_d20x1y2_xerrplus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05 };
  double p7006_d20x1y2_yval[] = { 1.01, 1.54, 1.9, 1.36, 1.2, 0.78, 0.47, 0.47, 0.41, 
				  0.23, 0.45, 0.31, 0.42, 0.38, 0.32, 0.18, 0.12, 0.14, 0.59, 
				  1.52 };
  double p7006_d20x1y2_yerrminus[] = { 0.15, 0.19, 0.21, 0.17, 0.16, 0.13, 0.1, 0.1, 0.1, 
				       0.07, 0.1, 0.08, 0.1, 0.09, 0.09, 0.06, 0.05, 0.06, 0.12, 
				       0.2 };
  double p7006_d20x1y2_yerrplus[] = { 0.15, 0.19, 0.21, 0.17, 0.16, 0.13, 0.1, 0.1, 0.1, 
				      0.07, 0.1, 0.08, 0.1, 0.09, 0.09, 0.06, 0.05, 0.06, 0.12, 
				      0.2 };
  double p7006_d20x1y2_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0 };
  double p7006_d20x1y2_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  int p7006_d20x1y2_numpoints = 20;
  for(int i=0; i<p7006_d20x1y2_numpoints; i++){
    p7006_d20x1y2_yval[i] *= 1000;
    p7006_d20x1y2_yerrminus[i] *= 1000;
    p7006_d20x1y2_yerrplus[i] *= 1000;
    // d(sigma)/d(cos) -> d(sigma)/d(Omega) = *(2./NBIN)/(2pi*2./NBIN) = *(1./2pi)
    p7006_d20x1y2_yval[i] *= 1./(2.*3.1415);
    p7006_d20x1y2_yerrminus[i] *= 1./(2.*3.1415);
    p7006_d20x1y2_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *p7006_d20x1y2 =  new TGraphAsymmErrors(p7006_d20x1y2_numpoints, p7006_d20x1y2_xval, p7006_d20x1y2_yval, p7006_d20x1y2_xerrminus, p7006_d20x1y2_xerrplus, p7006_d20x1y2_yerrminus, p7006_d20x1y2_yerrplus);
  p7006_d20x1y2->SetName("/HepData/7006/d20x1y2");
  p7006_d20x1y2->SetTitle("/HepData/7006/d20x1y2");
  p7006_d20x1y2->SetLineColor(4);
  p7006_d20x1y2->SetMarkerColor(4);
  p7006_d20x1y2->SetMarkerStyle(25);
  //p7006_d20x1y2.Draw("AP");

  // ==================================== //
  // reference data
  // NPB8(1968)233 @ CERN, bubble chamber
  // 2pi*d(sigma)/d(Omega)
  // ==================================== //
  // P991.0 MeV,   PI- SIGMA+
  double past_data1_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				 -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				 0.95 };
  double past_data1_xerrminus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05 };
  double past_data1_xerrplus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05 };
  double past_data1_yval[] = { 0.61, 0.81, 0.86, 0.65, 0.39, 0.26, 0.21, 0.10, 0.10, 0.48,
			       0.79, 0.66, 0.70, 2.03, 1.30, 2.08, 1.69, 1.60, 1.03, 0.75 };
  double past_data1_yerrminus[] = { 0.18, 0.22, 0.22, 0.19, 0.15, 0.12, 0.11, 0.07, 0.07, 0.16,
				    0.21, 0.20, 0.20, 0.36, 0.29, 0.38, 0.34, 0.34, 0.27, 0.27 };
  double past_data1_yerrplus[] = { 0.18, 0.22, 0.22, 0.19, 0.15, 0.12, 0.11, 0.07, 0.07, 0.16,
				   0.21, 0.20, 0.20, 0.36, 0.29, 0.38, 0.34, 0.34, 0.27, 0.27 };
  double past_data1_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  double past_data1_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0 };
  int past_data1_numpoints = 20;
  for(int i=0; i<past_data1_numpoints; i++){
    past_data1_yval[i] *= 1000;
    past_data1_yerrminus[i] *= 1000;
    past_data1_yerrplus[i] *= 1000;
    // 2pi*d(sigma)/d(Omega) -> d(sigma)/d(Omega)
    past_data1_yval[i] *= 1./(2.*3.1415);
    past_data1_yerrminus[i] *= 1./(2.*3.1415);
    past_data1_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *past_data1 = new TGraphAsymmErrors(past_data1_numpoints, past_data1_xval, past_data1_yval, past_data1_xerrminus, past_data1_xerrplus, past_data1_yerrminus, past_data1_yerrplus);
  past_data1->SetName("past_data1");
  past_data1->SetTitle("past_data1");
  past_data1->SetLineColor(6);
  past_data1->SetMarkerColor(6);
  past_data1->SetMarkerStyle(26);
 
  //past_data1.Draw("AP");

  // P1022.0 MeV,   PI- SIGMA+
  double past_data2_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				 -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				 0.95 };
  double past_data2_xerrminus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05 };
  double past_data2_xerrplus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05 };
  double past_data2_yval[] = { 0.38, 1.26, 0.69, 0.48, 0.51, 0.07, 0.0, 0.09, 0.08, 0.19,
			       0.70, 0.69, 0.78, 1.34, 0.90, 1.08, 1.94, 1.08, 1.70, 0.43 };
  double past_data2_yerrminus[] = { 0.20, 0.35, 0.25, 0.20, 0.21, 0.07, 0.06, 0.09, 0.08, 0.13,
				    0.26, 0.25, 0.27, 0.37, 0.30, 0.33, 0.47, 0.34, 0.47, 0.25 };
  double past_data2_yerrplus[] = { 0.20, 0.35, 0.25, 0.20, 0.21, 0.07, 0.06, 0.09, 0.08, 0.13,
				   0.26, 0.25, 0.27, 0.37, 0.30, 0.33, 0.47, 0.34, 0.47, 0.25 };
  double past_data2_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  double past_data2_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0 };
  int past_data2_numpoints = 20;
  for(int i=0; i<past_data2_numpoints; i++){
    past_data2_yval[i] *= 1000;
    past_data2_yerrminus[i] *= 1000;
    past_data2_yerrplus[i] *= 1000;
    // 2pi*d(sigma)/d(Omega) -> d(sigma)/d(Omega)
    past_data2_yval[i] *= 1./(2.*3.1415);
    past_data2_yerrminus[i] *= 1./(2.*3.1415);
    past_data2_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *past_data2 = new TGraphAsymmErrors(past_data2_numpoints, past_data2_xval, past_data2_yval, past_data2_xerrminus, past_data2_xerrplus, past_data2_yerrminus, past_data2_yerrplus);
  past_data2->SetName("past_data2");
  past_data2->SetTitle("past_data2");
  past_data2->SetLineColor(8);
  past_data2->SetMarkerColor(8);
  past_data2->SetMarkerStyle(27);
  //past_data2.Draw("AP");


  // P991.0 MeV,   PI+ SIGMA-
  double past_data3_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				 -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				 0.95 };
  double past_data3_xerrminus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05 };
  double past_data3_xerrplus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05 };
  double past_data3_yval[] = { 0.83, 1.51, 2.24, 1.97, 1.24, 0.75, 0.22, 0.37, 0.22, 0.41,
			       0.36, 0.67, 0.27, 0.69, 0.28, 0.14, 0.00, 0.15, 0.42, 1.56 };
  double past_data3_yerrminus[] = { 0.20, 0.29, 0.35, 0.33, 0.25, 0.19, 0.10, 0.13, 0.10, 0.14,
				    0.13, 0.18, 0.11, 0.18, 0.12, 0.08, 0.04, 0.09, 0.15, 0.32 };
  double past_data3_yerrplus[] = { 0.20, 0.29, 0.35, 0.33, 0.25, 0.19, 0.10, 0.13, 0.10, 0.14,
				   0.13, 0.18, 0.11, 0.18, 0.12, 0.08, 0.04, 0.09, 0.15, 0.32 };
  double past_data3_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  double past_data3_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0 };
  int past_data3_numpoints = 20;
  for(int i=0; i<past_data3_numpoints; i++){
    past_data3_yval[i] *= 1000;
    past_data3_yerrminus[i] *= 1000;
    past_data3_yerrplus[i] *= 1000;
    // 2pi*d(sigma)/d(Omega) -> d(sigma)/d(Omega)
    past_data3_yval[i] *= 1./(2.*3.1415);
    past_data3_yerrminus[i] *= 1./(2.*3.1415);
    past_data3_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *past_data3 = new TGraphAsymmErrors(past_data3_numpoints, past_data3_xval, past_data3_yval, past_data3_xerrminus, past_data3_xerrplus, past_data3_yerrminus, past_data3_yerrplus);
  past_data3->SetName("past_data3");
  past_data3->SetTitle("past_data3");
  past_data3->SetLineColor(6);
  past_data3->SetMarkerColor(6);
  past_data3->SetMarkerStyle(26);
  //past_data3.Draw("AP");

  // P1022.0 MeV,   PI+ SIGMA-
  double past_data4_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				 -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				 0.95 };
  double past_data4_xerrminus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05 };
  double past_data4_xerrplus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				     0.05 };
  double past_data4_yval[] = { 0.76, 1.83, 1.39, 1.89, 0.86, 0.65, 0.35, 0.14, 0.14, 0.37,
			       0.70, 0.51, 0.52, 0.43, 0.37, 0.16, 0.00, 0.08, 0.82, 1.33 };
  double past_data4_yerrminus[] = { 0.24, 0.39, 0.34, 0.40, 0.26, 0.22, 0.16, 0.10, 0.10, 0.17,
				    0.23, 0.20, 0.20, 0.18, 0.17, 0.12, 0.06, 0.08, 0.27, 0.31 };
  double past_data4_yerrplus[] = { 0.24, 0.39, 0.34, 0.40, 0.26, 0.22, 0.16, 0.10, 0.10, 0.17,
				    0.23, 0.20, 0.20, 0.18, 0.17, 0.12, 0.06, 0.08, 0.27, 0.31 };
  double past_data4_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  double past_data4_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0 };
  int past_data4_numpoints = 20;
  for(int i=0; i<past_data4_numpoints; i++){
    past_data4_yval[i] *= 1000;
    past_data4_yerrminus[i] *= 1000;
    past_data4_yerrplus[i] *= 1000;
    // 2pi*d(sigma)/d(Omega) -> d(sigma)/d(Omega)
    past_data4_yval[i] *= 1./(2.*3.1415);
    past_data4_yerrminus[i] *= 1./(2.*3.1415);
    past_data4_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *past_data4 = new TGraphAsymmErrors(past_data4_numpoints, past_data4_xval, past_data4_yval, past_data4_xerrminus, past_data4_xerrplus, past_data4_yerrminus, past_data4_yerrplus);
  past_data4->SetName("past_data4");
  past_data4->SetTitle("past_data4");
  past_data4->SetLineColor(8);
  past_data4->SetMarkerColor(8);
  past_data4->SetMarkerStyle(27);
  //past_data4.Draw("AP");


  // ==================================== //
  // reference data
  // NPB193(1981)21./RL-80-073 @ CERN, bubble chamber
  // d(sigma)/d(cos)
  // ==================================== //
  // P990.0 MeV,   REK- P --> PI- SIGMA+,    SQRT(S)1.789 GeV
  double p7266_d12x1y8_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				  -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				  0.99 };
  double p7266_d12x1y8_xerrminus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				       0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				       0.05 };
  double p7266_d12x1y8_xerrplus[] = { 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
				      0.05 };
  double p7266_d12x1y8_yval[] = { 0.88, 0.98, 0.94, 0.65, 0.58, 0.41, 0.17, 0.07, 0.2, 
				  0.27, 0.5, 0.87, 0.85, 1.06, 1.86, 1.93, 1.77, 1.53, 1.24, 
				  0.8 };
  double p7266_d12x1y8_yerrminus[] = { 0.17, 0.18, 0.18, 0.15, 0.14, 0.12, 0.08, 0.05, 0.08, 
				       0.09, 0.13, 0.17, 0.17, 0.18, 0.23, 0.24, 0.23, 0.22, 0.19, 
				       0.16 };
  double p7266_d12x1y8_yerrplus[] = { 0.17, 0.18, 0.18, 0.15, 0.14, 0.12, 0.08, 0.05, 0.08, 
				      0.09, 0.13, 0.17, 0.17, 0.18, 0.23, 0.24, 0.23, 0.22, 0.19, 
				      0.16 };
  double p7266_d12x1y8_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0 };
  double p7266_d12x1y8_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  int p7266_d12x1y8_numpoints = 20;
  for(int i=0; i<p7266_d12x1y8_numpoints; i++){
    p7266_d12x1y8_yval[i] *= 1000.;
    p7266_d12x1y8_yerrminus[i] *= 1000.;
    p7266_d12x1y8_yerrplus[i] *= 1000.;
    // d(sigma)/d(cos) -> d(sigma)/d(Omega) = *(2/NBIN)/(2pi*2/NBIN) = *(1./2pi)
    p7266_d12x1y8_yval[i] *= 1./(2.*3.1415);
    p7266_d12x1y8_yerrminus[i] *= 1./(2.*3.1415);
    p7266_d12x1y8_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *p7266_d12x1y8 = new TGraphAsymmErrors(p7266_d12x1y8_numpoints, p7266_d12x1y8_xval, p7266_d12x1y8_yval, p7266_d12x1y8_xerrminus, p7266_d12x1y8_xerrplus, p7266_d12x1y8_yerrminus, p7266_d12x1y8_yerrplus);
  p7266_d12x1y8->SetName("/HepData/7266/d12x1y8");
  p7266_d12x1y8->SetTitle("/HepData/7266/d12x1y8");
  p7266_d12x1y8->SetLineColor(3);
  p7266_d12x1y8->SetMarkerColor(3);
  p7266_d12x1y8->SetMarkerStyle(28);
  //p7266_d12x1y8.Draw("AP");


  // P990.0 MeV,   REK- P --> PI+ SIGMA-,    SQRT(S)1.789 GeV
  double p7266_d17x1y8_xval[] = { -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, 
				  -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
				  0.95 };
  double p7266_d17x1y8_xerrminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  double p7266_d17x1y8_xerrplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				      0.0 };
  double p7266_d17x1y8_yval[] = { 0.77, 1.51, 1.73, 1.84, 1.21, 0.77, 0.5, 0.46, 0.37, 
				  0.48, 0.31, 0.47, 0.37, 0.44, 0.17, 0.17, 0.09, 0.23, 0.54, 
				  1.61 };
  double p7266_d17x1y8_yerrminus[] = { 0.16, 0.22, 0.23, 0.23, 0.2, 0.16, 0.13, 0.12, 0.11, 
				       0.13, 0.1, 0.12, 0.11, 0.12, 0.08, 0.07, 0.05, 0.09, 0.13, 
				       0.22 };
  double p7266_d17x1y8_yerrplus[] = { 0.16, 0.22, 0.23, 0.23, 0.2, 0.16, 0.13, 0.12, 0.11, 
				      0.13, 0.1, 0.12, 0.11, 0.12, 0.08, 0.07, 0.05, 0.09, 0.13, 
				      0.22 };
  double p7266_d17x1y8_ystatminus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					0.0 };
  double p7266_d17x1y8_ystatplus[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
				       0.0 };
  int p7266_d17x1y8_numpoints = 20;
  for(int i=0; i<p7266_d17x1y8_numpoints; i++){
    p7266_d17x1y8_yval[i] *= 1000;
    p7266_d17x1y8_yerrminus[i] *= 1000;
    p7266_d17x1y8_yerrplus[i] *= 1000;
    // d(sigma)/d(cos) -> d(sigma)/d(Omega) = *(2/NBIN)/(2pi*2/NBIN) = *(1./2pi)
    p7266_d17x1y8_yval[i] *= 1./(2.*3.1415);
    p7266_d17x1y8_yerrminus[i] *= 1./(2.*3.1415);
    p7266_d17x1y8_yerrplus[i] *= 1./(2.*3.1415);
  }
  TGraphAsymmErrors *p7266_d17x1y8 = new TGraphAsymmErrors(p7266_d17x1y8_numpoints, p7266_d17x1y8_xval, p7266_d17x1y8_yval, p7266_d17x1y8_xerrminus, p7266_d17x1y8_xerrplus, p7266_d17x1y8_yerrminus, p7266_d17x1y8_yerrplus);
  p7266_d17x1y8->SetName("/HepData/7266/d17x1y8");
  p7266_d17x1y8->SetTitle("/HepData/7266/d17x1y8");
  p7266_d17x1y8->SetLineColor(3);
  p7266_d17x1y8->SetMarkerColor(3);
  p7266_d17x1y8->SetMarkerStyle(28);
  //p7266_d17x1y8.Draw("AP");


  //TCanvas *w1;
  //w1 = new TCanvas("w1", "", 600, 300);
  //w1->Divide(2,1);
  //for( int x=0; x<2; x++ ){
  //  w1->cd(x+1);
  //  his = (TH1F*)Corr_Cospi[x];
  //  his->SetStats(0);
  //  his->SetTitle(0);
  //  his->Draw();
  //  his->SetXTitle("cos#theta_{#pi}^{CM}");
  //  his->SetYTitle("d#sigma/d#Omega [#mub/sr]");
  //  his->SetMinimum(0);
  //  his->SetMaximum(600);
  //  his->SetMarkerStyle(20);
  TFile *fcs = TFile::Open("CSsigma_H2.root","READ");
  TFile *fK0 = TFile::Open("CS_K0contami_H2.root","READ");
  // TGraphErrors *gCS_Sp = (TGraphErrors*) fcs->Get("Graph_from_CS_Sp");
//  TGraphErrors *gCS_Sm = (TGraphErrors*) fcs->Get("Graph_from_CS_Sm");
  //TGraphAsymmErrors *gCS_Sp = (TGraphAsymmErrors*) fcs->Get("Graph_from_Cospicm_pi_Sp");
  //TGraphAsymmErrors *gCS_Sm = (TGraphAsymmErrors*) fcs->Get("Graph_from_Cospicm_pi_Sm");
  TGraphAsymmErrors *gCS_Sp = (TGraphAsymmErrors*) fcs->Get("gCS_Sp_K0sub");
  TGraphAsymmErrors *gCS_Sm = (TGraphAsymmErrors*) fcs->Get("gCS_Sm_K0sub");
  TGraphAsymmErrors *gCS_Spsys = (TGraphAsymmErrors*) fcs->Get("gCS_Sp_syssum_K0sub");
  TGraphAsymmErrors *gCS_Smsys = (TGraphAsymmErrors*) fcs->Get("gCS_Sm_syssum_K0sub");
  TGraphErrors *gCSK0_Sp = (TGraphErrors*)fK0->Get("gCSK0_Sp");
  TGraphErrors *gCSK0_Sm = (TGraphErrors*)fK0->Get("gCSK0_Sm");
  TCanvas *cCS_Sp = new TCanvas("cCS_Sp","cCS_Sp",1000,800);
  TLegend *leg = new TLegend(0.3, 0.6, 0.6, 0.9);
  leg->AddEntry(gCS_Sp, "This analysis", "LP");
  //  if(x==0){ // K- P --> PI- SIGMA+
  cCS_Sp->cd();
  gCS_Sp->GetXaxis()->SetRangeUser(0.5,1);
  gCS_Sp->GetYaxis()->SetRangeUser(0,600);
  gCS_Sp->SetMarkerStyle(20);
  //gCS_Sp->SetFillColor(0);
  gCS_Sp->SetMarkerColor(1);
  gCS_Sp->SetLineColor(1);
  gCS_Sp->GetXaxis()->SetTitle("Cos#theta_{CM} miss-#pi^{-}");
  gCS_Sp->GetYaxis()->SetTitle("d#sigma/d#Omega [#mub/sr]");
  gCS_Sp->GetXaxis()->CenterTitle();
  gCS_Sp->GetYaxis()->CenterTitle();
  gCS_Sp->SetTitle("#Sigma^{+}#pi^{-}");
  gCS_Sp->Draw("AP");
  gCS_Spsys->SetLineColor(1);
  gCS_Spsys->Draw("5");
  p7071_d40x1y1->Draw("PZ same");//p7071_d40x1y1->Print();
  p7006_d21x1y2->Draw(" PZ same");
  past_data1->Draw(" PZ same");
  past_data2->Draw(" PZ same");
  p7266_d12x1y8->Draw(" PZ same");
      leg->AddEntry(p7266_d12x1y8, "NPB193(1981)21   0.99 GeV/c", "LP");
      leg->AddEntry(past_data1,    "NPB8(1968)233,   0.991 GeV/c", "LP");
      leg->AddEntry(p7071_d40x1y1, "NPB90(1975)349,  1.001 GeV/c", "LP");
      leg->AddEntry(p7006_d21x1y2, "NPB105(1976)189, 1.005 GeV/c", "LP");
      leg->AddEntry(past_data2,    "NPB8(1968)233,   1.022 GeV/c", "LP");
  //  }
   //   leg->Draw();
   // else if(x==1){ // K- P --> PI+ SIGMA-
  TCanvas *cCS_Sm = new TCanvas("cCS_Sm","cCS_Sm",1000,800);
  cCS_Sm->cd();
  gCS_Sm->GetXaxis()->SetRangeUser(0.5,1);
  gCS_Sm->GetYaxis()->SetRangeUser(0,300);
  gCS_Sm->SetMarkerStyle(20);
  gCS_Sm->GetXaxis()->SetTitle("Cos#theta_{CM} miss-#pi^{+}");
  gCS_Sm->GetYaxis()->SetTitle("d#sigma/d#Omega [#mub/sr]");
  gCS_Sm->GetXaxis()->CenterTitle();
  gCS_Sm->GetYaxis()->CenterTitle();
  gCS_Sm->SetTitle("#Sigma^{-}#pi^{+}");
  gCS_Sm->Draw("AP");
  gCS_Smsys->SetLineColor(1);
  gCS_Smsys->Draw("5");
      p7071_d60x1y1->Draw("same PZ");
      p7006_d20x1y2->Draw("same PZ");
      past_data3->Draw("same PZ");
      past_data4->Draw("same PZ");
      p7266_d17x1y8->Draw("same PZ");
      TLegend *leg2 = new TLegend(0.3, 0.6, 0.6, 0.9);
      leg2->AddEntry(gCS_Sm, "This Analysis ", "LP");
      leg2->AddEntry(p7266_d17x1y8, "NPB193(1981)21   0.99 GeV/c", "LP");
      leg2->AddEntry(past_data3,    "NPB8(1968)233,   0.991 GeV/c", "LP");
      leg2->AddEntry(p7071_d60x1y1, "NPB90(1975)349,  1.001 GeV/c", "LP");
      leg2->AddEntry(p7006_d20x1y2, "NPB105(1976)189, 1.005 GeV/c", "LP");
      leg2->AddEntry(past_data4,    "NPB8(1968)233,   1.022 GeV/c", "LP");
    //}
    leg2->Draw();
    //his->Draw("same");
 // }
 // w1->Print("tmp.pdf");
  return;
  // CS obtained by MM & IM fits [bg_piS.C] 20191113
  // mean CS from 2/4/6 MeVee cut
  double val[2][6] = {{ 0, 485.9, 487.9, 328.6, 273.4, 165.8 },
		      { 0, 0, 0, 76.9, 104.7, 200.7 }};

  // mean CS error from 2/4/6 MeVee cut
  double sta[2][6] = {{ 0, 95.9, 69.0, 39.5, 26.0, 13.7 },
		      { 0, 0, 0, 13.7, 12.1, 13.2 }};

  // syst - each value (max-cent etc.)
  double sys[2][6] = {{ 0, 51.4, 44.0, 17.9, 29.3, 8.5 },
		      { 0, 0, 0, 8.2, 13.8, 15.3 }};
 

  TCanvas *w2;
  w2 = new TCanvas("w2", "", 600, 300);
  w2->Divide(2,1);
  for( int x=0; x<2; x++ ){
    w2->cd(x+1);
    his = (TH1F*)Corr_Cospi[x]->Clone("his");
    TH1* hise = (TH1F*)his->Clone("hise");
    hise->Reset("ICES");
    for( int i=14; i<20; i++ ){
      double v = his->GetBinContent(i+1);
      double e = his->GetBinError(i+1);
      std::cerr<<i<<" "<<v<<" "<<e<<std::endl;
      his->SetBinContent(i+1, val[x][i-14]);
      his->SetBinError(i+1, sta[x][i-14]);
      hise->SetBinContent(i+1, val[x][i-14]);
      hise->SetBinError(i+1, sys[x][i-14]);
    }
    his->SetStats(0);
    his->SetTitle(0);
    his->Draw();
    his->SetXTitle("cos#theta_{#pi}^{CM}");
    his->SetYTitle("d#sigma/d#Omega [#mub/sr]");
    his->SetMinimum(0);
    his->SetMaximum(600);
    his->SetMarkerStyle(20);
    TLegend *leg = new TLegend(0.3, 0.6, 0.6, 0.9);
    leg->AddEntry(his, "DATA", "LP");
    if(x==0){ // K- P --> PI- SIGMA+
      p7071_d40x1y1.Draw("same PZ");
      p7006_d21x1y2.Draw("same PZ");
      past_data1.Draw("same PZ");
      past_data2.Draw("same PZ");
      p7266_d12x1y8.Draw("same PZ");
      leg->AddEntry(&p7266_d12x1y8, "NPB193(1981)21   0.99 GeV/c", "LP");
      leg->AddEntry(&past_data1,    "NPB8(1968)233,   0.991 GeV/c", "LP");
      leg->AddEntry(&p7071_d40x1y1, "NPB90(1975)349,  1.001 GeV/c", "LP");
      leg->AddEntry(&p7006_d21x1y2, "NPB105(1976)189, 1.005 GeV/c", "LP");
      leg->AddEntry(&past_data2,    "NPB8(1968)233,   1.022 GeV/c", "LP");
    }
    else if(x==1){ // K- P --> PI+ SIGMA-
      p7071_d60x1y1.Draw("same PZ");
      p7006_d20x1y2.Draw("same PZ");
      past_data3.Draw("same PZ");
      past_data4.Draw("same PZ");
      p7266_d17x1y8.Draw("same PZ");
      leg->AddEntry(&p7266_d17x1y8, "NPB193(1981)21   0.99 GeV/c", "LP");
      leg->AddEntry(&past_data3,    "NPB8(1968)233,   0.991 GeV/c", "LP");
      leg->AddEntry(&p7071_d60x1y1, "NPB90(1975)349,  1.001 GeV/c", "LP");
      leg->AddEntry(&p7006_d20x1y2, "NPB105(1976)189, 1.005 GeV/c", "LP");
      leg->AddEntry(&past_data4,    "NPB8(1968)233,   1.022 GeV/c", "LP");
    }
    leg->Draw();
    his->Draw("same");
    hise->SetFillStyle(0);
    hise->Draw("same,e2");

  }
  w2->Print("tmp.pdf)");



}
