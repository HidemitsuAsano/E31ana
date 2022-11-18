
double Legendre(int n, float x);
void PlotAngularKbarN()
{

  //K-n -> K-n elastic
  double LegendreC_Kmn[7]    = {1.0,  1.09,  2.53,  1.13,  1.08,  0.29,  0.28} ;
  double LegendreC_Kmnerr[7] = {0,     0.2,   0.2,  0.2,    0.2,   0.1,  0.1} ;
  //K-p -> K0n charge ex.
  double LegendreC_K0n[7] = { 0.44,  -0.096,  0.506,  -0.373,  0.527,  -0.204,  -0.007};
  double LegendreC_K0nerr[7] = { 0.015,0.032,0.044,0.048,0.053,0.054,0.059 };


   const int nparamKm = 7;
   const int nparamK0 = 7;
   
   std::vector<double> coscmKmn, yKmn, yKmnehi, yKmnelow;
   std::vector<double> coscmK0n, yK0n, yK0nehi, yK0nelow;
   for(double coscm=-1;coscm<1;coscm += 0.01){
     double S_Kmn=0;
     double S_Kmnup=0;
     double S_Kmndown=0;
     for(int il=0;il<nparamKm;il++){
       S_Kmn += LegendreC_Kmn[il]*Legendre(il,coscm);
       S_Kmnup += (LegendreC_Kmn[il]+LegendreC_Kmnerr[il]  ) *Legendre(il,coscm);
       S_Kmndown += (LegendreC_Kmn[il]-LegendreC_Kmnerr[il]  ) *Legendre(il,coscm);
     }
     coscmKmn.push_back(coscm);
     yKmn.push_back(S_Kmn);
     if(S_Kmnup > S_Kmndown){
       yKmnehi.push_back(S_Kmnup-S_Kmn);
       yKmnelow.push_back(S_Kmn-S_Kmndown);
     }else{
       yKmnehi.push_back(S_Kmndown-S_Kmn);
       yKmnelow.push_back(S_Kmn-S_Kmnup);
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
     yK0n.push_back(S_K0n);
     if(S_K0nup > S_K0ndown){
       yK0nehi.push_back(S_K0nup-S_K0n);
       yK0nelow.push_back(S_K0n-S_K0ndown);
     }else{
       yK0nehi.push_back(S_K0ndown-S_K0n);
       yK0nelow.push_back(S_K0n-S_K0nup);
     }
   }


   TGraphAsymmErrors *grKmn = new TGraphAsymmErrors(coscmKmn.size(),&coscmKmn[0],&yKmn[0],0,0,&yKmnelow[0],&yKmnehi[0]);
   grKmn->SetLineColor(4);
   grKmn->SetFillColor(4);
   grKmn->SetFillStyle(4);
   grKmn->Draw("ap5");
   
   TCanvas *cK0 = new TCanvas("cK0","cK0");
   TGraphAsymmErrors *grK0n = new TGraphAsymmErrors(coscmK0n.size(),&coscmK0n[0],&yK0n[0],0,0,&yK0nelow[0],&yK0nehi[0]);
   grK0n->SetLineColor(4);
   grK0n->SetFillColor(4);
   grK0n->SetFillStyle(4);
   grK0n->Draw("ap5");
   
  
}


double Legendre(int n, float x)
{
  float f = 0.0;

  if (n==0){
    f = 1.0;
  } else
  if (n==1){
    f = x;
  } else
  if (n==2){
    f = 0.5*(3.0*x*x-1.0);
  } else
  if (n==3){
    f = 0.5*(5.0*x*x*x-3.0*x);
  } else
  if (n==4){
    f = 0.125*(35.0*x*x*x*x-30.0*x*x + 3.0);
  } else
  if (n==5){
    f = 0.125*(63.0*x*x*x*x*x-70.0*x*x*x+15.0*x);
  } else
  if (n==6){
    f = 0.0625*(231.0*x*x*x*x*x*x - 315.0*x*x*x*x + 105.0*x*x - 5.0);
  } else
  if (n==7){
    f = 0.0625*(429.0*x*x*x*x*x*x*x - 693.0*x*x*x*x*x + 315.0*x*x*x - 35.0*x);
  } else
  if (n==8){
    f = 0.0078125*(6435.0*x*x*x*x*x*x*x*x - 12012.0*x*x*x*x*x*x + 6930.0*x*x*x*x - 1260.0*x*x + 35);
  }

 else {
   // G4cout << "n > 8 is not supported!!!" << G4endl;
  }
  return f;
}
