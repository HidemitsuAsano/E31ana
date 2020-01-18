//Missing mass
Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(x[0]<1.116){
    return par[0]*exp(-0.5*pow((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2.0)); 
  }else if(1.116<=x[0] && x[0]<1.5){
    return par[4]*exp(-0.5*pow(((x[0]-par[5])/par[6]),2.0)); 
  }else{
    return 0.;
  }
}

Double_t param_MMnmiss[7]={
  1.64736,
  0.906359,
  0.235915,
  -0.929384,
  2.58101,
  1.48951,
  0.265225};

//IM(npip+)
Double_t func_IMnpip(Double_t *x,Double_t *par)
{
  if(x[0]<1.08){
    return 0.;
  } else if(1.08 <= x[0] && x[0]<1.10){
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
  } else if(1.10 <= x[0] && x[0]<1.25){
    return exp(par[3]+par[4]*x[0]); 
  } else if(1.25 <= x[0] && x[0]<2.0){
    return par[5]*exp(-0.5*pow(((x[0]-par[6])/par[7]),2.0)); 
  } else {
    return 0.;
  }
}

Double_t param_IMnpip[8]={
     1.10885,
     1.09604,
     0.0172468,
     0.815124,
    -0.715797,
     1.23375,
     1.56849,
     0.43139
};

Double_t func_IMnpim(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 0.;
  } else if(1.07 <= x[0] && x[0]<1.10){
    return par[0]+par[1]*x[0]+pow(par[2]*x[0],2.0)+paw(par[3]*x[0],3.0);
  } else if(1.10 <= x[0] && x[0]<2.00){
    return par[4]+par[5]*x[0]+pow(par[6]*x[0],2.0)+paw(par[7]*x[0],3.0);
  } else {
    return 0.;
  }
}

Double_t param_IMnpim[8]={
  77240.2,
   -214642.,
   198787.,
   -61356.4,
   8.92664,
   -6.14661,
   -3.5443,
   2.62163,
};
