
double Legendre(int n, float x);
void PlotAngularKbarN()
{
  //K-n -> K-n elastic
  const double CSKm = 20000.0;//ub
  const double CSKmerr = 1.0;
  double LegendreC_Kmn[7]    = {1.0,  1.09,  2.53,  1.13,  1.08,  0.29,  0.28} ;
  double LegendreC_Kmnerr[7] = {0,     0.2,   0.2,  0.2,    0.2,   0.1,  0.1} ;
  //K-p -> K0n charge ex.
  const double CSK0  = 6454.;//ub
  const double CSK0err = 0.042;
  double LegendreC_K0n[7] = { 0.44,  -0.096,  0.506,  -0.373,  0.527,  -0.204,  -0.007};
  double LegendreC_K0nerr[7] = { 0.015,0.032,0.044,0.048,0.053,0.054,0.059 };

  const int nparamKm = 7;
  const int nparamK0 = 7;

  std::vector<double> coscmKmn, yKmn, yKmnehi, yKmnelow;
  std::vector<double> coscmK0n, yK0n, yK0nehi, yK0nelow;
  double normKm = CSKm/7.7461;//integral of K-
  double normK0 = CSK0/2.05008;//integral of K0

  for(double coscm=-1.;coscm<1;coscm += 0.001){
    double S_Kmn=0;
    double S_Kmnup=0;
    double S_Kmndown=0;
    for(int il=0;il<nparamKm;il++){
      S_Kmn += LegendreC_Kmn[il]*Legendre(il,coscm);
      S_Kmnup += (LegendreC_Kmn[il]+LegendreC_Kmnerr[il]  ) *Legendre(il,coscm);
      S_Kmndown += (LegendreC_Kmn[il]-LegendreC_Kmnerr[il]  ) *Legendre(il,coscm);
    }
    coscmKmn.push_back(coscm);
    yKmn.push_back(S_Kmn*normKm);
    if(S_Kmnup > S_Kmndown){
      yKmnehi.push_back((S_Kmnup-S_Kmn)*normKm);
      yKmnelow.push_back((S_Kmn-S_Kmndown)*normKm);
    }else{
      yKmnehi.push_back((S_Kmndown-S_Kmn)*normKm);
      yKmnelow.push_back((S_Kmn-S_Kmnup)*normKm);
    }
    double S_K0n=0;
    double S_K0nup=0;
    double S_K0ndown=0;
    for(int il=0;il<nparamK0;il++){
      S_K0n += LegendreC_K0n[il]*Legendre(il,coscm);
      S_K0nup += (LegendreC_K0n[il]+LegendreC_K0nerr[il]  ) *Legendre(il,coscm);
      S_K0ndown += (LegendreC_K0n[il]-LegendreC_K0nerr[il]  ) *Legendre(il,coscm);
    }
    coscmK0n.push_back(coscm);
    yK0n.push_back(S_K0n*normK0);
    if(S_K0nup > S_K0ndown){
      yK0nehi.push_back((S_K0nup-S_K0n)*normK0);
      yK0nelow.push_back((S_K0n-S_K0ndown)*normK0);
    }else{
      yK0nehi.push_back((S_K0ndown-S_K0n)*normK0);
      yK0nelow.push_back((S_K0n-S_K0nup)*normK0);
    }
  }

  std::vector<double> coscmsumn, ysumn, ysumnehi, ysumnelow;//n angle (K0n + K-n sum) 
  const int ndata = coscmKmn.size();
  for(int i=0;i<ndata;i++){
    double S_nsum=0;
    double S_nsumup=0;
    double S_nsumdown=0;
    coscmsumn.push_back(-1.0*coscmKmn.at(i));
    ysumn.push_back(yKmn.at(i)+yK0n.at(i));
    ysumnehi.push_back(yKmnehi.at(i)+yK0nehi.at(i));
    ysumnelow.push_back(yKmnelow.at(i)+yK0nelow.at(i));
  }


  //K-n -> K-n elastic
  TCanvas *cKm = new TCanvas("cKm","cKm");
  TGraphAsymmErrors *grKmn = new TGraphAsymmErrors(coscmKmn.size(),&coscmKmn[0],&yKmn[0],0,0,&yKmnelow[0],&yKmnehi[0]);
  grKmn->SetTitle("K^{-}n #rightarrow K^{-}n");
  grKmn->SetLineColor(4);
  grKmn->SetFillColor(4);
  grKmn->SetFillStyle(3001);
  grKmn->GetXaxis()->SetTitle("cosCM K^{-}");
  grKmn->GetXaxis()->CenterTitle();
  grKmn->Draw("a4");
  std::cout << "K- " << grKmn->Integral(0,coscmKmn.size())  << std::endl;


  //K-p ->K0n 
  TCanvas *cK0 = new TCanvas("cK0","cK0");
  TGraphAsymmErrors *grK0n = new TGraphAsymmErrors(coscmK0n.size(),&coscmK0n[0],&yK0n[0],0,0,&yK0nelow[0],&yK0nehi[0]);
  grK0n->SetTitle("K^{-}p #rightarrow #bar{K}^{0}n");
  grK0n->SetLineColor(4);
  grK0n->SetFillColor(4);
  grK0n->SetFillStyle(3001);
  grK0n->GetXaxis()->SetTitle("cosCM #bar{K}^{0}");
  grK0n->GetXaxis()->CenterTitle();
  grK0n->Draw("a4");
  std::cout << "K0 " << grK0n->Integral(0,coscmK0n.size())  << std::endl;

  TCanvas *csum = new TCanvas("csum","csum");
  TGraphAsymmErrors *grsumn = new TGraphAsymmErrors(coscmsumn.size(),&coscmsumn[0],&ysumn[0],0,0,&ysumnelow[0],&ysumnehi[0]);
  grsumn->SetTitle("K^{-}p #rightarrow #bar{K}^{0}n  &  K^{-}n #rightarrow K^{-}n  sum");
  grsumn->SetLineColor(2);
  grsumn->SetFillColor(2);
  grsumn->SetFillStyle(3001);
  grsumn->GetXaxis()->SetTitle("cosCM n");
  grsumn->GetXaxis()->CenterTitle();
  grsumn->Draw("a4");
  grsumn->SetName("grsumn");
  grK0n->SetName("grK0n");
  grKmn->SetName("grKmn");
  std::cout << "n " << grsumn->Integral(0,coscmsumn.size())  << std::endl;
  TFile *feleangle = new TFile("feleangle.root","RECREATE");
  grsumn->Write();
  grK0n->Write();
  grKmn->Write();


}


double Legendre(int n, float x)
{
  float f = 0.0;

  if (n==0){
    f = 1.0;
  } else if(n==1){
    f = x;
  } else if(n==2){
    f = 0.5*(3.0*x*x-1.0);
  } else if(n==3){
    f = 0.5*(5.0*x*x*x-3.0*x);
  } else if (n==4){
    f = 0.125*(35.0*x*x*x*x-30.0*x*x + 3.0);
  } else if (n==5){
    f = 0.125*(63.0*x*x*x*x*x-70.0*x*x*x+15.0*x);
  } else if (n==6){
    f = 0.0625*(231.0*x*x*x*x*x*x - 315.0*x*x*x*x + 105.0*x*x - 5.0);
  } else if (n==7){
    f = 0.0625*(429.0*x*x*x*x*x*x*x - 693.0*x*x*x*x*x + 315.0*x*x*x - 35.0*x);
  } else if (n==8){
    f = 0.0078125*(6435.0*x*x*x*x*x*x*x*x - 12012.0*x*x*x*x*x*x + 6930.0*x*x*x*x - 1260.0*x*x + 35);
  } else {
    // G4cout << "n > 8 is not supported!!!" << G4endl;
  }
  return f;
}
