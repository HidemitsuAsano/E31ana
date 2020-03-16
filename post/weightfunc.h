//Missing mass


Double_t func_MMnmiss(Double_t *x,Double_t *par)
{
   if(0.0<x[0] && x[0]<1.06){
     return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0);
   }else if(1.06<=x[0] && x[0]<1.115 ){
     return par[10]+par[11]*x[0]+par[12]*pow(x[0],2.0);
   }else if(1.115<x[0] && x[0]<1.5){
     return par[13]+par[14]*x[0]+par[15]*pow(x[0],2.0)+par[16]*pow(x[0],3.0)+par[17]*pow(x[0],4.0);
   }else{
     return 0;
   }
}

Double_t param_MMnmiss[18]={
  0.00229985,
  -0.078266,
  7.28136,
  -56.1799,
  196.571,
  -305.011,
  108.226,
  245.874,
  -280.618,
  85.8073,
  595.196,
  -1120.27,
  529.193,
  -46.2783,
  81.7302,
  1.25405,
  -54.0075,
  20.3614
};

Double_t param_MMnmiss_corr[18]={
  1.96335,
  -3.24185,
  15.6824,
  -18.8532,
  -8.62532,
  13.4245,
  14.0861,
  -5.64434,
  -18.9119,
  11.1562,
  3.37067,
  -3.2522,
  0.9139,
  40.0071,
  -57.8304,
  -3.90708,
  36.6173,
  -13.2461
};

Double_t param_MMnmiss_corr2[18]={
  1.06569,
  2.1173,
  -10.0731,
  13.8839,
  -0.979701,
  -9.12568,
  -1.73836,
  5.7668,
  3.40999,
  -3.45749,
  3.64316,
  -5.69929,
  2.924,
  -27.2928,
  44.084,
  1.0866,
  -28.7472,
  11.3445
};


Double_t param_MMnmiss_corr3[18]={
  1.06531,
  -0.403542,
  0.722876,
  1.93327,
  -5.97501,
  0.618494,
  7.67952,
  0.739576,
  -10.7638,
  5.37883,
  6.41249,
  -9.71456,
  4.3228,
  -18.9226,
  32.9773,
  -0.248451,
  -22.1338,
  9.10923
};

//wK0 v3 use func_MMnmiss
Double_t param_MMnmiss_wK0_corr3[18]={
  1.72962,
  -9.61454,
  44.1084,
  -83.7342,
  52.1106,
  23.1782,
  -23.9721,
  -13.0442,
  8.30979,
  1.85185,
  176.179,
  -324.257,
  149.945,
  50.3157,
  -81.7262,
  -0.401614,
  56.2598,
  -23.1369
};

Double_t func_MMnmiss_wK0(Double_t *x,Double_t *par)
{
  if(0<x[0] && x[0]<1.05){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0);
  }else if(1.05<=x[0] && x[0]<1.40) {
    return par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)
    +par[9]*pow(x[0],4.0);
  }else{ 
    return 0;
  }
}


Double_t param_MMnmiss_wK0[10]={
-0.008548,
0.670184,
-1.7269,
3.40415,
-0.300939,
-20.7379,
11.2535,
16.9602,
6.34851,
-12.3974
};



Double_t param_MMnmiss_wK0_corr[10]={
  0.553894,
  9.68177,
  -27.7101,
  32.0468,
  -13.4982,
  2.30503,
  -0.504196,
  -0.874001,
  -0.277716,
  0.361689
};

Double_t param_MMnmiss_wK0_corr2[10]={
  1.15529,
  1.32269,
  -4.70185,
  5.39287,
  -2.18035,
  21.5779,
  -35.5298,
  0.344495,
  26.1662,
  -11.4934
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


Double_t func_IMnpip(Double_t *x,Double_t *par)
{
  if(1.06<=x[0] && x[0]<1.25){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  }else if(1.25<=x[0] && x[0]<1.92) {
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)
    +par[8]*pow(x[0],4.0);
  }else{ 
    return 1.;
  }
};

Double_t param_IMnpip[9]={
  -319.022,
  740.821,
  -567.762,
  143.975,
  -215.28,
  535.924,
  -488.017,
  193.642,
  -28.3287
};

Double_t param_IMnpip_wK0[9]={
-1185.59,
3042.75,
-2603.28,
743.055,
-1226.74,
3267.49,
-3233.84,
1411.57,
-229.532
};




Double_t func_IMnpim(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<1.11) {
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
  } else if(1.11 <= x[0] && x[0]<=1.7) {
    return par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0)+par[6]*pow(x[0],3.0)+par[7]*pow(x[0],4.0);
  } else {
    return 1.;
  }
}

Double_t param_IMnpim[8]={
1.74,
1.114,
0.02519,
126.126,
-316.818,
301.132,
-127.767,
20.3552
};



Double_t func_IMnpim_wK0(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<1.11) {
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
  } else if(1.11 <= x[0] && x[0]<=1.7) {
    return par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0)+par[6]*pow(x[0],3.0)+par[7]*pow(x[0],4.0)+par[8]*pow(x[0],5.0);
  } else {
    return 1.;
  }
}

Double_t param_IMnpim_wK0[9]={
1.3699,
1.09471,
0.0108092,
3201.6,
-11318,
15913.4,
-11116.7,
3857.44,
-531.86
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
  66.7953,
  -1067.75,
  7583.42,
  -30074,
  72041.4,
  -106532,
  95090.7,
  -46952.3,
  9844.54
};


Double_t param_nmom_v305[9]={
  6.84201,
  -108.686,
  862.596,
  -3690.65,
  9276.35,
  -14119.6,
  12787.2,
  -6332.86,
  1319.26
};



Double_t param_nmom_wK0[9]={
  52.4252,
  -829.667,
  5870.17,
  -23144.3,
  54993.3,
  -80529.5,
  71096.2,
  -34693.1,
  7184.83
};


Double_t param_nmom_wK0_v305[9]={
5.94715,
-92.4088,
757.239,
-3350.8,
8684.88,
-13606.8,
12678,
-6462.93,
1387.31
};




Double_t func_pipmom(Double_t *x,Double_t *par)
{
  if(x[0]<0.06){
    return 1.0;
  }else if(0.06<= x[0] && x[0]<0.70){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0);
  }else{
    return 0.0;
  }
}


Double_t param_pipmom[9]={
0.934566,
-60.6745,
1217.44,
-10210.8,
45380.8,
-115656,
169541,
-132863,
43110.8
};



Double_t param_pipmom_wK0[9]={
  -11.4993,
    339.057,
    -3744.95,
    22299.8,
    -76098.1,
    152471,
    -177179,
    110449,
    -28540.1
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
  if(x[0]<1.0){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return 1.0;
  }
}

Double_t param_Mompippim[6]={
  1.10873,
  -3.93937,
  34.8285,
  -103.117,
  117.397,
  -46.2839
};


Double_t param_Mompippim_corr[6]={
 1.07468,
 -1.43303,
 9.66428,
 -29.7112,
 41.4303,
 -20.8171
};


Double_t param_Mompippim_wK0[6]={
  1.19725,
  -16.9661,
  159.908,
  -581.328,
  916.705,
  -508.948
};


Double_t param_Mompippim_wK0_corr[6]={
 0.785575,
 -1.67717,
 19.0277,
 -68.8766,
 122.602,
 -77.1012
};





Double_t func_IMpippim(Double_t *x,Double_t *par)
{
  if(x[0]<0.28){
    return 1.0;
  }else if(0.28 <= x[0] && x[0] <= 0.97){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)
    +par[6]*pow(x[0],6.0);
  }else{
    return 1.0;
  }
}


Double_t param_IMpippim[7]={
5.87624,
-58.0674,
178.204,
-75.3271,
-380.574,
531.672,
-201.818
};



Double_t func_IMnpipi(Double_t *x,Double_t *par)
{
  if(x[0]<1.19) {
    return 1.;
  } else if(1.19 <= x[0] && x[0]<1.35) {
    //return exp(par[0]+par[1]*x[0]);
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  } else if(1.35 <= x[0] && x[0]<2.0) {
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)
          +par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0)  
          +par[11]*pow(x[0],7.0);   
  } else {
    return 1.;
  }
}

Double_t param_IMnpipi[12]={
1335.63,
-3224.24,
2584.02,
-687.339,
979.52,
649.894,
-10007,
18988.2,
-16787,
7957.51,
-1963.1,
198.693
};

Double_t func_IMnpipi_wK0(Double_t *x,Double_t *par)
{
  if(x[0]<1.4){
    return 1.;
  }else if(1.4<=x[0] && x[0]<=2.00) {
  //  return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0);
  }else{
    return 1.;
  }
}


Double_t param_IMnpipi_wK0[6]={
  -2689.01,
  7308.9
  -7875.29,
  4208.6,
  -1116.09,
  117.533
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


Double_t func_q(Double_t *x,Double_t *par)
{
  if(0<x[0] && x[0]<1.5){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0);
  }else{
    return 1.0;
  }
}


Double_t param_q[8]={
  6.65152,
  8.14869,
  -143.676,
  425.004,
  -636.562,
  541.661,
  -248.891,
  47.9206
};

Double_t param_q_corr[8]={
  0.785026,
  2.57318,
  -22.4933,
  107.921,
  -278.341,
  375.986,
  -244.352,
  61.9188,
};

Double_t param_q_corr2[8]={
0.569318,
0.261346,
-0.802152,
8.46911,
-2.09522,
-8.4265,
4.1212,
0.273668
};

Double_t param_q_corr3[8]={
};

Double_t param_q_corr4[8]={
};

Double_t param_q_corr5[8]={
};

Double_t param_q_corr6[8]={
};

Double_t param_q_corr7[8]={
};




Double_t param_q_wK0[8]={
  3.35328,
  -0.101967,
  -59.486,
  207.349,
  -340.556,
  307.91,
  -147.653,
  29.3524
};

Double_t param_q_wK0_corr[8]={
  0.698491,
  6.046,
  -56.4227,
  247.988,
  -578.943,
  735.773,
  -471.459,
  119.356
};

Double_t param_q_wK0_corr2[8]={
0.505648,
-0.3244,
11.7009,
-58.2898,
155.612,
-197.424,
117.117,
-26.2648
};

Double_t param_q_wK0_corr3[8]={
};

Double_t param_q_wK0_corr4[8]={
};


Double_t param_q_wK0_corr5[8]={
};



Double_t param_q_wK0_corr6[8]={
};
