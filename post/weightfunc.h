//Missing mass
Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(x[0]<1.116){
    return par[0]*exp(-0.5*pow((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2.0)); 
  }else if(1.116<=x[0] && x[0]<1.5){
    return par[4]*exp(-0.5*pow(((x[0]-par[5])/par[6]),2.0)); 
  }else{
    return 1.;
  }
}


Double_t param_MMnmiss[7]={
  2.54287,
  0.919059,
  0.151276,
  -0.874727,
  2.83954,
  1.48577,
  0.258338};

Double_t param_MMnmiss_corr[7]={
  1.07364,
  2.01487,
  45.9014,
  22.4603,
  1.06126,
  1.01389,
  0.528322};




//IM(npip+)
Double_t func_IMnpip(Double_t *x,Double_t *par)
{
  if(x[0]<1.08){
    return 1.;
  } else if(1.08 <= x[0] && x[0]<1.10){
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
  } else if(1.10 <= x[0] && x[0]<1.25){
    return exp(par[3]+par[4]*x[0]); 
  } else if(1.25 <= x[0] && x[0]<2.0){
    return par[5]*exp(-0.5*pow(((x[0]-par[6])/par[7]),2.0)); 
  } else {
    return 1.;
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

Double_t param_IMnpip_wK0[8]={
  4.30739,
  0.993153,
  0.0563368,
  -0.328542,
  0.25817,
  1.51837,
  1.37419,
  0.118939
};



Double_t func_IMnpim(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 1.;
  } else if(1.07 <= x[0] && x[0]<1.10){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  } else if(1.10 <= x[0] && x[0]<1.70){
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0);
  } else {
    return 1.;
  }
}

Double_t param_IMnpim[8]={
  183924.,
  -506893.,
  465649.,
  -142581.,
  23.0502,
  -50.0064,
  37.8527,
  -9.57951
};

Double_t param_IMnpim_wK0[8]={
  670470.,
  -1.84802e+06,
  1.69783e+06,
  -519928.,
  21.839,
  -47.0231,
  35.9704,
  -9.37302
};


Double_t func_MMom(Double_t *x,Double_t *par)
{
  if(x[0]<0.4){
    return 1.0;
  }else if(0.4<= x[0] && x[0]<1.5){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return 1.0;
  }
}


Double_t param_MMom[6]={
  -5.86273,
    49.697,
    -133.313,
  166.548,
    -97.7351,
    21.737
};

Double_t param_MMom_corr[6]={
  -5.48873,
  38.0513,
  -83.0645,
  86.0181,
  -42.9564,
  8.37741
};

Double_t param_MMom_wK0[6]={
   -0.162034,
    6.15245,
    15.2708,
    -71.5352,
    82.0394,
    -28.1011
};

Double_t param_MMom_wK0_corr[6]={
  -1.0349,
  25.1577,
  -73.4263,
  89.5101,
  -49.5942,
  10.3317
};


Double_t func_nmom(Double_t *x,Double_t *par)
{
  if(x[0]<0.14){
    return 1.0;
  }else if(0.14<= x[0] && x[0]<1.0){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_nmom[9]={
  6.52201,
  -102.617,
  755.024,
  -2881.34,
  6339.79,
  -8366.75,
  6529.5,
  -2766.3,
  486.944
};

Double_t param_nmom_wK0[9]={
 -0.972747,
   57.1471,
   -639.926,
   3615.68,
   -11347.4,
   20614.5,
   -21550.4,
   12025.9,
   -2774.2
};


Double_t func_pipmom(Double_t *x,Double_t *par)
{
  if(x[0]<0.08){
    return 1.0;
  }else if(0.08<= x[0] && x[0]<0.70){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_pipmom[9]={
  6.69768,
  -186.48,
  2330.16,
  -14279.8,
  45669.6,
  -71999.8,
  35235.8,
  34255.4,
  -33836.8
};


Double_t param_pipmom_wK0[9]={
  -4.38653,
  227.618,
  -3861.76,
  33849.9,
  -169632.,
  503872.,
  -875238.,
  819170.,
  -317802.
};

Double_t param_pipmom_wK0_corr[9]={
 -2.94713,
 108.708,
 -1308.73,
 9142.9,
 -40229.2,
 111437.,
 -185455.,
 167671.,
 -62943.4
};



Double_t func_pimmom(Double_t *x,Double_t *par)
{
  if(x[0]<0.06){
    return 1.0;
  }else if(0.06<= x[0] && x[0]<0.73){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_pimmom[9]={
  8.1307,    
 -232.745,    
  2867.69,    
 -17454.4,    
  58065.8,    
  -107834,    
   104767,    
 -41082.4,    
  -365.787   
};


Double_t param_pimmom_wK0[9]={
  -3.38479,
    142.391,
    -2272.15,
    20489.,
    -107791.,
    335405.,
    -606102.,
    585953.,
    -233120.
};




Double_t func_Mompippim(Double_t *x,Double_t *par)
{
  if(x[0]<0.98){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  }else{
    return 1.0;
  }
}

Double_t param_Mompippim[4]={
 1.01184,
 1.44482,
  -5.32594,
 3.39873
};

Double_t func_IMpippim(Double_t *x,Double_t *par)
{
  if(x[0]<0.28){
    return 1.0;
  }else if(0.28 <= x[0] && x[0] < 0.97){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0);
  }else{
    return 1.0;
  }
}

Double_t param_IMpippim[7]={
  22.6532,
  -259.532,
  1236.75,
  -2977.95,
  3771.67,
  -2340.11,
  542.993
};


Double_t func_IMnpipi(Double_t *x,Double_t *par)
{
  if(x[0]<1.22){
    return 1.;
  }else if(1.22<=x[0] && x[0]<2.00) {
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return 1.;
  }
}

Double_t func_IMnpipi_wK0(Double_t *x,Double_t *par)
{
  if(x[0]<1.22){
    return 1.;
  }else if(1.22<=x[0] && x[0]<2.00) {
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0);
  }else{
    return 1.;
  }
}

Double_t param_IMnpipi[6]={
 3655.85,
 -11877.6,
 15385.6,
 -9924.94,
 3186.74,
 -407.218
};


Double_t param_IMnpipi_wK0[5]={
  -3481.04,
  8597.94,
  -7929.73,
  3236.3,
  -492.935
};

