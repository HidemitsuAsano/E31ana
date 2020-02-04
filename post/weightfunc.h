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

Double_t param_MMnmiss_corr[7]={
  1.07364,
  2.01487,
  45.9014,
  22.4603,
  1.06126,
  1.01389,
  0.528322};
*/

Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
  if(0<x[0] && x[0]<1.12) {
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


Double_t func_MMnmiss_corr(Double_t *x,Double_t *par)
{
  if(0.0 <x[0] && x[0]<0.86){
    return 1.0;
  }else if(0.86<=x[0] && x[0]<1.40){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0);
  }else{
    return 1.0;
  }
}


Double_t param_MMnmiss_corr[7]={
  -5262.058359,
  28413.402661, 
  -63537.801498, 
  75330.303498, 
  -49941.176533, 
  17554.106342, 
  -2555.797844 
};

Double_t func_MMnmiss_corr2(Double_t *x,Double_t *par)
{
  if(0.0 <x[0] && x[0]<0.78){
    return 1.0;
  }else if(0.78<=x[0] && x[0]<1.50){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0);
  }else{
    return 1.0;
  }
}

Double_t param_MMnmiss_corr2[7]={
  -856.534093,
  4911.977326,
  -11598.489427,
  14447.996381,
  -10014.411373,
  3663.111015,
  -552.664798 
};



Double_t func_MMnmiss_mul(Double_t *x,Double_t *par)
{
  return func_MMnmiss(&x[0],&par[0])*func_MMnmiss_corr(&x[0],&par[12])*func_MMnmiss_corr2(&x[0],&par[19]) ;
}



Double_t func_MMnmiss_wK0(Double_t *x,Double_t *par)
{
  if(0<x[0] && x[0]<1.05){
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

Double_t param_MMnmiss_wK0_corr[10]={
   0.729318,
   3.499398, 
 -13.489569, 
  19.241882, 
  -8.890493, 
   0.071010, 
   0.652360, 
   0.673780, 
   0.213987,
  -0.602705 
};

Double_t param_MMnmiss_wK0_corr2[10]={
   0.950016,
  -2.030146, 
  13.543845, 
 -21.590511, 
  10.298567, 
 406.365779, 
 -687.941032, 
 -20.959304, 
 526.633366, 
 -221.406781 
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
 /*
 14.564993,
 -250.889223, 
 1970.987705,
 -8167.852054, 
 19812.064637,
 -29213.144032, 
 25803.695052,
 -12551.544466,
 2583.399144 
 */
  7.902928,
  -135.560771, 
  1057.258135, 
  -4308.463166, 
  10200.236327, 
  -14569.209201,
  12366.972812, 
  -5731.878085, 
  1113.876156 
};


Double_t param_nmom_corr[9]={
   -0.708220,
   30.930350, 
  -254.858919, 
  1167.488566, 
  -3142.463775, 
  5085.774947, 
  -4863.326637, 
  2530.254126, 
  -551.985785 
};




Double_t param_nmom_wK0[9]={
  /*
  5.637132,
  -87.411220,
  627.422705, 
  -2207.331058, 
  4284.083705, 
  -4738.747351, 
  2849.782330, 
  -774.223646, 
  41.481221 
  */
  5.435864,
  -89.784427 ,
  676.759745, 
  -2569.240720,
  5578.457146, 
  -7254.775564,
  5577.731920, 
  -2329.123836,
  405.226024 
};


Double_t param_nmom_wK0_corr[9]={
  -1.204856,
  41.161437, 
 -338.246799, 
 1551.636399, 
 -4235.483964, 
 6984.450592, 
 -6802.010060, 
 3597.041146, 
 -796.191342 
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
   7.864272,
 -231.438415, 
 3003.086228, 
 -19581.795118, 
 69983.446486, 
 -139003.416918, 
 144597.079010, 
 -63170.795188, 
 2591.178812 
};

/*
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
};*/





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
  16873.372632 
  4.607702,
  -114.813459, 
  1286.104696, 
  -7324.184562, 
  24692.647823, 
  -52136.361698,
  67997.988007, 
  -49975.656940,
  15753.832819 */
      4.761786,
  -117.597953, 
  1301.086943, 
  -7324.008728, 
  24403.268865, 
  -50960.716944,
  65873.916717, 
  -48120.574334, 
  15118.954446 
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
   8.785302,
 -259.702985, 
 3263.641163, 
 -20431.247986, 
 71345.355180, 
 -144557.884564, 
 166351.254056, 
 -98116.074470, 
 21832.308887 
};


Double_t param_pimmom_wK0[9]={
    2.128888,
  -42.982985, 
  580.745659, 
  -3411.868426, 
  9651.958685, 
  -12300.654097,
  1983.825305, 
  9964.374290, 
  -6753.731277 
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
  /*
  -8.454380,
  -9.339147, 
  0.948693, 
  0.943733, 
  2.495556, 
  1.473135 */
 -9.353473,
-10.342388,
  0.906928,
  0.915011,
  2.610918,
  1.560637 
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
  /*
  -8.052033,
  -8.903550, 
  0.981118, 
  0.614340, 
  1.344470, 
  0.661165 */
  -8.473595,
 -9.376222, 
  0.879538, 
  0.198493, 
  0.877567, 
  0.512796 
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
  -9.462239,
 -12.379514, 
   0.957114, 
  -0.011929, 
   0.526920, 
   0.446129 
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
  0.933200,
  0.183558,
  0.141299,
  -0.585394,
  5.228950,
  6.826591 
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
  0.997606,
  0.111185,
  0.938944,
  5.478411,
  3.586188,
  -24.128151,
  -13.398514,
  74.738372,
  73.161008 
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
   0.964987,
   0.129265,
   1.293863,
   2.040254,
 -17.122606,
 -36.781413,
  87.956922,
 259.287208,
 159.063798 
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
  0.515040, 
 -0.819747, 
 -0.414733, 
 -0.064466, 
  1.088241, 
  0.011982, 
 -0.061337, 
  0.015112 
};


Double_t param_phinpip_wK0[8]={
   0.483800, 
  -0.931034, 
  -0.507744, 
  -0.088314, 
   1.188477, 
  -0.352878, 
   0.137703, 
  -0.012871 
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
   1.369371,
  -0.021330,
  -0.911575,
  -0.797019,
  -0.266403,
  -0.031387,
   1.047935,
   0.418273,
  -0.373949,
   0.068988 
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
    1.262950,
   -0.350961,
   -1.892044,
   -1.756252,
   -0.631022,
   -0.079336,
    1.513497,
   -0.572610,
    0.166750,
   -0.017844 
};
