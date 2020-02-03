//Missing mass

/*
Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(x[0]<1.116){
    return par[0]*exp(-0.5*pow((x[0]-par[1])/(par[2]+(x[0]<par[1])*par[3]*(x[0]-par[1])),2.0)); 
  }else if(1.116<=x[0] && x[0]<1.5){
    return par[4]*exp(-0.5*pow(((x[0]-par[5])/par[6]),2.0)); 
  }else{
    return 1.;
  }
}*/

/*
Double_t param_MMnmiss[7]={
  2.54287,
  0.919059,
  0.151276,
  -0.874727,
  2.83954,
  1.48577,
  0.258338};
*/
Double_t param_MMnmiss_corr[7]={
  1.07364,
  2.01487,
  45.9014,
  22.4603,
  1.06126,
  1.01389,
  0.528322};


Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(x[0]<1.12) {
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  } else if(1.12<=x[0] && x[0]<1.5) {
    return par[6]+par[7]*x[0]+par[8]*pow(x[0],2.0)+par[9]*pow(x[0],3.0)
    +par[10]*pow(x[0],4.0)+par[11]*pow(x[0],5.0);
  } else {
    return 1.;
  }
}


Double_t param_MMnmiss[12]={
   1.269308,
  -1.832203, 
  30.448986, 
 -90.510740, 
 105.535941, 
 -42.719303, 
 -121.197566,
 123.372754, 
  66.063951, 
 -40.664855, 
 -64.274182, 
  33.698206 
};



Double_t func_MMnmiss_wK0(Double_t *x,Double_t *par)
{
  if(x[0]<1.05){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0);
  }else if(1.05<=x[0] && x[0]<1.50) {
    return par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)
    +par[9]*pow(x[0],4.0);
  }else{ 
    return 1.;
  }
}


Double_t param_MMnmiss_wK0[10]={
   0.961862,
  -2.193081, 
  13.915484, 
 -21.849738, 
  10.338688, 
 403.554988, 
 -684.177521, 
 -19.923610, 
 523.680608, 
 -220.485964 
};



/*

Double_t param_MMnmiss_wK0[10]={
   1.269308,
  -1.832203, 
  30.448986, 
 -90.510740, 
 105.535941, 
 -42.719303, 
 -121.197566,
 123.372754, 
  66.063951, 
 -40.664855, 
 -64.274182, 
  33.698206 
};
*/







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
  14.564993,
 -250.889223, 
 1970.987705,
 -8167.852054, 
 19812.064637,
 -29213.144032, 
 25803.695052,
 -12551.544466,
 2583.399144 
};

Double_t param_nmom_wK0[9]={
  5.637132,
  -87.411220,
  627.422705, 
  -2207.331058, 
  4284.083705, 
  -4738.747351, 
  2849.782330, 
  -774.223646, 
  41.481221 
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

/*
Double_t param_pipmom[9]={
   7.483433,
 -215.010358, 
 2775.303785,
 -17909.009610, 
 62443.095579, 
 -117666.604158, 
 108053.041757, 
 -28828.746802,
 -10899.659054 
};*/

Double_t param_pipmom[9]={
   7.452202,
 -213.554201, 
 2750.025390, 
 -17696.704073, 
 61465.830531, 
 -115080.561086,
 104140.974163, 
 -25694.488752, 
 -11928.126127 
};


Double_t param_pipmom_wK0[9]={
 /*
  4.338103,
  -108.036329, 
  1223.867445, 
  -7072.068393, 
  24309.791523, 
  -52495.222080, 
  70038.843418 ,
  -52560.063945 ,
  16873.372632 */
  4.607702,
  -114.813459, 
  1286.104696, 
  -7324.184562, 
  24692.647823, 
  -52136.361698,
  67997.988007, 
  -49975.656940,
  15753.832819 
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
   9.095308,
 -269.228884, 
 3432.977598, 
 -21943.903612, 
 78417.721706, 
 -162868.224400,
 192696.478111, 
 -117678.128409,
 27613.785804 
};


Double_t param_pimmom_wK0[9]={
    2.267339,
  -49.837594, 
  691.129233, 
  -4264.817683, 
  13304.697807, 
  -21446.999941,
  15303.410719, 
  -478.665434, 
  -3348.043095 
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


/*
Double_t func_cosn(Double_t *x,Double_t *p)
{
  if(-0.92<= x[0] && x[0]<0.30){
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    //+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
    return (p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3)+p[4]*TMath::Power(x[0],4)+p[5]*TMath::Power(x[0],5)+p[6]*TMath::Power(x[0],6)
    +p[7]*TMath::Power(x[0],7)+p[8]*TMath::Power(x[0],8)) ; 
  }else{
    return 1.0;
  }
}

Double_t param_cosn[9]={
    1.015500, 
    1.006929, 
   -2.654014, 
  -28.385825, 
   10.852024, 
  384.093426, 
  945.706089, 
  919.509377, 
  323.637467 
};


Double_t func_cosn_wK0(Double_t *x,Double_t *p)
{
  if(-0.92<= x[0] && x[0]<0.20){
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    //+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
    return (p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3)+p[4]*TMath::Power(x[0],4)+p[5]*TMath::Power(x[0],5)+p[6]*TMath::Power(x[0],6)
    +p[7]*TMath::Power(x[0],7)+p[8]*TMath::Power(x[0],8)) ; 
  }else{
    return 1.0;
  }
}

Double_t param_cosn_wK0[9]={
   0.977648,
   0.490488, 
  -7.280459, 
  -4.424170, 
 321.226554, 
 1457.236988, 
 2663.463669, 
 2241.231879, 
 718.822831 
};

*/

Double_t func_cosn(Double_t *x,Double_t *p)
{
  if(x[0]<-0.90){
    return TMath::Exp(p[0]+p[1]*x[0]);
  }else{
    return (p[2]+p[3]*x[0]+p[4]*TMath::Power(x[0],2)+p[5]*TMath::Power(x[0],3)) ; 
  }
}

Double_t param_cosn[6]={
  -8.454380,
  -9.339147, 
  0.948693, 
  0.943733, 
  2.495556, 
  1.473135 
};


Double_t func_cosn_wK0(Double_t *x,Double_t *p)
{
  if(x[0]<-0.90){
    return TMath::Exp(p[0]+p[1]*x[0]);
  }else{
    return (p[2]+p[3]*x[0]+p[4]*TMath::Power(x[0],2)+p[5]*TMath::Power(x[0],3)) ; 
  }
}

Double_t param_cosn_wK0[6]={
  -8.052033,
  -8.903550, 
  0.981118, 
  0.614340, 
  1.344470, 
  0.661165 
};



Double_t func_cospip(Double_t *x,Double_t *p)
{
  if(x[0]<-0.92){
    return 1.0;
  }else if(-0.92<=x[0] && x[0]<-0.75){
    return TMath::Exp(p[0]+p[1]*x[0]);
  }else{
    return (p[2]+p[3]*x[0]+p[4]*TMath::Power(x[0],2)+p[5]*TMath::Power(x[0],3)) ;
  }
}

Double_t param_cospip[6]={
  -9.694997,
 -12.635173,
   0.967719,
  -0.043963,
   0.362578,
   0.364421 
};


Double_t func_cospip_wK0(Double_t *x,Double_t *par)
{
  if(0.5<x[0]){
    return 1.0;
  }else{
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }
}

Double_t param_cospip_wK0[6]={
  0.940005,
  0.225100, 
  0.179505, 
  -0.567310, 
  4.950085, 
  6.509011 
};




Double_t func_cospim(Double_t *x,Double_t *p)
{
  if(x[0]<-0.92){
    return 1.0;
  }else if(-0.92<= x[0] && x[0]<0.50){
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    //+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
    return (p[0]+p[1]*x[0]+p[2]*TMath::Power(x[0],2)+p[3]*TMath::Power(x[0],3)+p[4]*TMath::Power(x[0],4)+p[5]*TMath::Power(x[0],5)+p[6]*TMath::Power(x[0],6)
    +p[7]*TMath::Power(x[0],7)+p[8]*TMath::Power(x[0],8)) ; 
  }else{
    return 1.0;
  }
}

Double_t param_cospim[9]={
   1.021878,
   0.190030, 
   0.854445, 
   4.520616, 
   0.026460, 
 -25.521797, 
   1.003050, 
  99.555661, 
  84.895447 
};


Double_t func_cospim_wK0(Double_t *x,Double_t *par)
{
  if(x[0]<-0.92){
    return 1.0;
  }else if(-0.92<= x[0] && x[0]<0.55){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_cospim_wK0[9]={
  0.979292,
  0.189653,
  1.102991,
  1.069020,
  -16.429364,
  -29.377114,
  92.255499,
  247.149189,
  147.666808 
};



Double_t func_phinpip(Double_t *x,Double_t *par)
{
  if(x[0]<-0.4){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  }else{
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0);
  }
}


Double_t param_phinpip[8]={
   0.527852,
  -0.782251, 
  -0.384799, 
  -0.058251, 
   1.079872, 
  -0.006227, 
  -0.042964, 
   0.011449 
};


Double_t param_phinpip_wK0[8]={
  0.444693, 
 -0.967532, 
 -0.516172, 
 -0.088949, 
  1.175380, 
 -0.354065, 
  0.144912, 
 -0.014048 
};



Double_t func_phinpim(Double_t *x,Double_t *par)
{
  if(x[0]<0.5){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return par[6]+par[7]*x[0]+par[8]*pow(x[0],2.0)+par[9]*pow(x[0],3.0);
  }
}


Double_t param_phinpim[10]={
  1.359675,
  0.019439,
 -0.843374,
 -0.751765,
 -0.252200,
 -0.029723,
  1.100313,
  0.350875,
 -0.350422,
  0.067273 
};



Double_t func_phinpim_wK0(Double_t *x,Double_t *par)
{
  if(x[0]<0.3){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return par[6]+par[7]*x[0]+par[8]*pow(x[0],2.0)+par[9]*pow(x[0],3.0);
  }
}


Double_t param_phinpim_wK0[10]={
  1.263597,
  -0.245141, 
  -1.660550, 
  -1.572056, 
  -0.567819, 
  -0.071471, 
  1.492275, 
  -0.542323, 
  0.152484, 
  -0.015358 
};
