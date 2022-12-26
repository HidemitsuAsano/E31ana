
double Legendre(int n, float x);
void PlotAngularKbarN()
{
  //K-n -> K-n elastic
  const double CSKm = 20000.0;//ub
  const double CSKmerr = 1000.0;
  const double CSKm800 = 9560; 
  const double CSKm800err =  960; 

  //963 MeV
  double LegendreC_Kmn[7]    = {1.0,  1.09,  2.53,  1.13,  1.08,  0.29,  0.28} ;
  double LegendreC_Kmnerr[7] = {0,     0.17,   0.2,  0.2,    0.2,   0.11,  0.12} ;
  //806 MeV
  double LegendreC_Kmn800[7]    = {1.0,  1.07,  2.15,  0.77,  0.90,  -0.04,  0.31} ;
  double LegendreC_Kmnerr800[7] = {0,     0.23,   0.28,  0.23,    0.21,   0.18,  0.20} ;
  //K-p -> K0n charge ex.
  const double CSK0  = 7840.;//ub
  const double CSK0err = 380.;
  const double CSK0800 = 4790.;
  const double CSK0800err = 270.;  

  double LegendreC_K0n[7] = { 0.44,  -0.096,  0.506,  -0.373,  0.527,  -0.204,  -0.007};
  double LegendreC_K0nerr[7] = { 0.015,0.032,0.044,0.048,0.053,0.054,0.059 };
  //K-p -> K0n charge ex. //862 MeV
  double LegendreC_K0n800[7] = { 0.215,  -0.114, 0.194,  -0.315,  0.216,  -0.078,  -0.032};
  double LegendreC_K0nerr800[7] = { 0.010,0.020,0.028,0.030,0.032,0.034,0.036 };

  const int nparamKm = 7;
  const int nparamK0 = 7;

  std::vector<double> coscmKmn, yKmn, yKmnehi, yKmnelow;
  std::vector<double> coscmKmn800, yKmn800, yKmnehi800, yKmnelow800;
  std::vector<double> coscmK0n, yK0n, yK0nehi, yK0nelow;
  std::vector<double> coscmK0n800, yK0n800, yK0nehi800, yK0nelow800;
  double normKm = CSKm/7.7461;//integral of K-
  double normK0 = CSK0/2.05008;//integral of K0
  double normKm800 = CSKm800/6.69481;//integral of K-
  double normK0800 = CSK0800/1.36187;//integral of K0
  
  TH2D *h2Kmn = new TH2D("h2Kmn","h2Kmn",2000,-1,1,4000,0,20000);
  TH2D *h2Kmn800 = new TH2D("h2Kmn800","h2Kmn800",2000,-1,1,4000,0,20000);
  TH2D *h2K0n = new TH2D("h2K0n","h2K0n",2000,-1,1,4000,0,20000);
  TH2D *h2K0n800 = new TH2D("h2K0n800","h2K0n800",2000,-1,1,4000,0,20000);
  for(double coscm=-1.;coscm<1;coscm += 0.001){
    double S_Kmn=0;
    //compute default val.
    for(int il=0;il<nparamKm;il++){
      S_Kmn += LegendreC_Kmn[il]*Legendre(il,coscm);
    }
   
    for(int imc=0;imc<10000;imc++){ 
      double val  = 0;
      for(int il=0;il<nparamKm;il++){
        val += (LegendreC_Kmn[il]+gRandom->Gaus(0,LegendreC_Kmnerr[il]))*Legendre(il,coscm);
      }
      h2Kmn->Fill(-1.*coscm,val*normKm);
    }
    coscmKmn.push_back(-1.*coscm);
    yKmn.push_back(S_Kmn*normKm);

    double S_Kmn800=0;
    for(int imc=0;imc<10000;imc++){ 
      double val  = 0;
      for(int il=0;il<nparamKm;il++){
        val += (LegendreC_Kmn800[il]+gRandom->Gaus(0,LegendreC_Kmnerr800[il]))*Legendre(il,coscm);
      }
      h2Kmn800->Fill(-1.*coscm,val*normKm800);
    }
    coscmKmn800.push_back(-1.0*coscm);
    yKmn800.push_back(S_Kmn800*normKm800);
    double S_K0n=0;
    for(int il=0;il<nparamKm;il++){
      S_K0n += LegendreC_K0n[il]*Legendre(il,coscm);
    }
    for(int imc=0;imc<10000;imc++){ 
      double val  = 0;
      for(int il=0;il<nparamK0;il++){
        val += (LegendreC_K0n[il]+gRandom->Gaus(0,LegendreC_K0nerr[il]))*Legendre(il,coscm);
      }
      h2K0n->Fill(-1.*coscm,val*normK0);
    }
    coscmK0n.push_back(-1.0*coscm);
    yK0n.push_back(S_K0n*normK0);
    
    double S_K0n800=0;
    for(int imc=0;imc<10000;imc++){ 
      double val  = 0;
      for(int il=0;il<nparamK0;il++){
        val += (LegendreC_K0n800[il]+gRandom->Gaus(0,LegendreC_K0nerr800[il]))*Legendre(il,coscm);
      }
      h2K0n800->Fill(-1.*coscm,val*normK0800);
    }
    coscmK0n800.push_back(-1.0*coscm);
    yK0n800.push_back(S_K0n800*normK0800);
  }
   
  auto *ch2Kmn = new TCanvas("ch2Kmn","ch2Kmn");
  h2Kmn->Draw("colz");
  auto *ch2K0n = new TCanvas("ch2K0n","ch2K0n");
  h2K0n->Draw("colz");
  
  //
  //TGraphAsymmErrors *gKmn = new TGraphAsymmErrors();
  //gKmn->SetName("gKmn");
  //TGraphAsymmErrors *gK0n = new TGraphAsymmErrors();
  //gK0n->SetName("gK0n");
 
  const int ndata = coscmKmn.size();
  std::cout << "ndata" << ndata << std::endl;
  TH1D* pyKmn[ndata];
  for(int ip = 0;ip<ndata;ip++){
    pyKmn[ip] = (TH1D*)h2Kmn->ProjectionY(Form("pyKmn%d",ip),ip+1,ip+2);
    double rms = pyKmn[ip]->GetRMS();
    yKmnelow.push_back( rms);
    yKmnehi.push_back( rms);
  }

  TH1D* pyK0n[ndata];
  for(int ip = 0;ip<ndata;ip++){
    pyK0n[ip] = (TH1D*)h2K0n->ProjectionY(Form("pyK0n%d",ip),ip+1,ip+2);
    double rms = pyK0n[ip]->GetRMS();
    yK0nelow.push_back( rms);
    yK0nehi.push_back( rms);
  }


  //K-n -> K-n elastic
  TCanvas *cKm = new TCanvas("cKm","cKm");
  TGraphAsymmErrors *grKmn = new TGraphAsymmErrors(coscmKmn.size(),&coscmKmn[0],&yKmn[0],0,0,&yKmnelow[0],&yKmnehi[0]);
  grKmn->SetTitle("K^{-}n #rightarrow K^{-}n");
  grKmn->SetLineColor(4);
  grKmn->SetFillColor(4);
  grKmn->SetFillStyle(3001);
  //grKmn->GetXaxis()->SetTitle("cosCM K^{-}");
  grKmn->GetXaxis()->SetTitle("cosCM n");
  grKmn->GetXaxis()->CenterTitle();
  grKmn->GetYaxis()->SetTitle("A.U.");
  grKmn->GetYaxis()->SetMaxDigits(3);
  grKmn->GetYaxis()->CenterTitle();
  grKmn->Draw("ap4");
  //grKmn->Print();
  std::cout << "K- " << grKmn->Integral(0,coscmKmn.size())  << std::endl;
 

  TCanvas *cK0 = new TCanvas("cK0","cK0");
  TGraphAsymmErrors *grK0n = new TGraphAsymmErrors(coscmK0n.size(),&coscmK0n[0],&yK0n[0],0,0,&yK0nelow[0],&yK0nehi[0]);
  grK0n->SetTitle("K^{-}p #rightarrow #bar{K}^{0}n");
  grK0n->SetLineColor(4);
  grK0n->SetFillColor(4);
  grK0n->SetFillStyle(3001);
  //grK0n->GetXaxis()->SetTitle("cosCM #bar{K}^{0}");
  grK0n->GetXaxis()->SetTitle("cosCM n");
  grK0n->GetXaxis()->CenterTitle();
  grK0n->GetYaxis()->SetTitle("A.U.");
  grK0n->GetYaxis()->SetMaxDigits(3);
  grK0n->GetYaxis()->CenterTitle();
  grK0n->Draw("ap4");
  std::cout << "K0 " << grK0n->Integral(0,coscmK0n.size())  << std::endl;

  std::vector<double> coscmsumn, ysumn, ysumnehi, ysumnelow;//n angle (K0n + K-n sum) 
  //std::vector<double> coscmsumn800, ysumn800, ysumnehi800, ysumnelow800;//n angle (K0n + K-n sum) 
  for(int i=0;i<ndata;i++){
    double S_nsum=0;
    double S_nsumup=0;
    double S_nsumdown=0;
    coscmsumn.push_back(coscmKmn.at(i));
    ysumn.push_back(yKmn.at(i)+yK0n.at(i));
    ysumnehi.push_back(yKmnehi.at(i)+yK0nehi.at(i));
    ysumnelow.push_back(yKmnelow.at(i)+yK0nelow.at(i));
    //double S_nsum800=0;
    //double S_nsumup800=0;
    //double S_nsumdown800=0;
    //coscmsumn800.push_back(coscmKmn800.at(i));
    //ysumn800.push_back(yKmn800.at(i)+yK0n800.at(i));
    //ysumnehi800.push_back(yKmnehi800.at(i)+yK0nehi800.at(i));
    //ysumnelow800.push_back(yKmnelow800.at(i)+yK0nelow800.at(i));
  }


  //K-n -> K-n elastic
  //TCanvas *cKm800 = new TCanvas("cKm800","cKm800");
  //TGraphAsymmErrors *grKmn800 = new TGraphAsymmErrors(coscmKmn800.size(),&coscmKmn800[0],&yKmn800[0],0,0,&yKmnelow800[0],&yKmnehi800[0]);
  //grKmn800->SetTitle("K^{-}n #rightarrow K^{-}n");
  //grKmn800->SetLineColor(3);
  //grKmn800->SetFillColor(3);
  //grKmn800->SetFillStyle(3001);
  //grKmn800->GetXaxis()->SetTitle("cosCM K^{-}");
  //grKmn800->GetXaxis()->SetTitle("cosCM n");
  //grKmn800->GetXaxis()->CenterTitle();
  //grKmn800->GetYaxis()->SetTitle("A.U.");
  //grKmn800->GetYaxis()->SetMaxDigits(3);
  //grKmn800->GetYaxis()->CenterTitle();
  //grKmn800->Draw("a4");
  //std::cout << "K-800: " << grKmn800->Integral(0,coscmKmn800.size())  << std::endl;
  //cKm->cd();
  //grKmn800->Draw("4");

  //K-p ->K0n 
  //TGraphAsymmErrors *grK0n800 = new TGraphAsymmErrors(coscmK0n800.size(),&coscmK0n800[0],&yK0n800[0],0,0,&yK0nelow800[0],&yK0nehi800[0]);
  //grK0n800->SetTitle("K^{-}p #rightarrow #bar{K}^{0}n");
  //grK0n800->SetLineColor(3);
  //grK0n800->SetFillColor(3);
  //grK0n800->SetFillStyle(3001);
  //grK0n800->GetXaxis()->SetTitle("cosCM n");
  //grK0n800->GetXaxis()->CenterTitle();
  //grK0n800->Scale(1.8);
  //grK0n800->Draw("c");
  //std::cout << "K0800 " << grK0n800->Integral(0,coscmK0n800.size())  << std::endl;
  
  TCanvas *csum = new TCanvas("csum","csum");
  TGraphAsymmErrors *grsumn = new TGraphAsymmErrors(coscmsumn.size(),&coscmsumn[0],&ysumn[0],0,0,&ysumnelow[0],&ysumnehi[0]);
  //TGraphAsymmErrors *grsumn800 = new TGraphAsymmErrors(coscmsumn800.size(),&coscmsumn800[0],&ysumn800[0],0,0,&ysumnelow800[0],&ysumnehi800[0]);
  grsumn->SetTitle("K^{-}p #rightarrow #bar{K}^{0}n  &  K^{-}n #rightarrow K^{-}n  sum");
  grsumn->SetLineColor(2);
  grsumn->SetFillColor(2);
  grsumn->SetFillStyle(3001);
  grsumn->GetXaxis()->SetTitle("cosCM n");
  grsumn->GetXaxis()->CenterTitle();
  grsumn->GetYaxis()->SetTitle("A.U.");
  grsumn->GetYaxis()->SetMaxDigits(3);
  grsumn->GetYaxis()->CenterTitle();
  grsumn->Draw("ac");

  //grsumn800->SetLineColor(3);
  //grsumn800->SetFillColor(3);
  //grsumn800->SetFillColor(3);
  //grsumn800->SetFillStyle(3002);
  //grsumn800->Scale(1.9);
  //grsumn800->Draw("c");

 
  grsumn->SetName("grsumn");
  //grsumn->SetName("grsumn800");
  grK0n->SetName("grK0n");
  grKmn->SetName("grKmn");
  //grKmn800->SetName("grKmn");
  std::cout << "n " << grsumn->Integral(0,coscmsumn.size())  << std::endl;
  TFile *feleangle = new TFile("feleangle.root","RECREATE");
  grsumn->Write();
  //grsumn800->Write();
  grK0n->Write();
  grKmn->Write();
  //grKmn800->Write();
  //grK0n800->Write();

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
