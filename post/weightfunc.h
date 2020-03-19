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
1.08141,
-1.6342,
10.8374,
-12.96,
-5.34959,
10.0688,
9.65441,
-5.64761,
-15.2146,
10.1763,
4.5168,
-5.72916,
2.2231,
31.997,
-45.6284,
-3.37016,
28.8531,
-10.3279
};

Double_t param_MMnmiss_corr2[18]={
  1.42416,
  -0.339543,
  0.502144,
  -0.170362,
  -0.927739,
  -0.562897,
  0.406446,
  1.00229,
  0.595765,
  -0.929904,
  1.1192,
  0.199145,
  -0.324193,
  1.1662,
  0.0793446,
  -0.247439,
  -0.139973,
  0.118117
};


Double_t param_MMnmiss_corr3[18]={
};

//wK0 v3 use func_MMnmiss
Double_t param_MMnmiss_wK0[18]={
  0.00726697,
  0.0673846,
  4.28727,
  -27.1111,
  94.5407,
  -147.735,
  39.534,
  166.387,
  -183.989,
  56.1687,
  131.436,
  -242.709,
  113.887,
  218.859,
  -394.94,
  27.3593,
  278.09,
  -126.597
};


Double_t param_MMnmiss_wK0_corr[18]={
0.724401,
0.440247,
0.458313,
-0.0358145,
-0.200083,
-0.176673,
-0.175163,
-0.199547,
-0.100104,
0.309052,
1.49351,
0.302198,
-0.737416,
-0.243336,
0.889151,
0.906723,
0.223969,
-0.833955
};


Double_t param_MMnmiss_wK0_corr2[18]={
1.24054,
4.09863,
-32.3319,
98.914,
-109.221,
-32.3324,
110.896,
35.7169,
-130.335,
54.2821,
-27.1712,
52.4349,
-24.5099,
101.228,
-164.906,
-1.37626,
113.435,
-46.5863
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


Double_t param_IMnpip_corr[9]={
-118.684,
301.742,
-252.959,
70.5453,
-153.475,
424.237,
-431.773,
193.042,
-32.0311
};


Double_t param_IMnpip_wK0[9]={
-1779.86,
4566.16,
-3903.96,
1113.01,
-1656.31,
4487.98,
-4524.16,
2013.4,
-334.067
};

Double_t param_IMnpip_wK0_corr[9]={
1576.47,
-3973.79,
3342.09,
-937.169,
-92.396,
299.593,
-352.617,
181.006,
-34.3025
};

Double_t func_IMnpim(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<1.11) {
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
  } else if(1.11 <= x[0] && x[0]<=2.0) {
    return par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0)+par[6]*pow(x[0],3.0)+par[7]*pow(x[0],4.0);
  } else {
    return 1.;
  }
}

Double_t param_IMnpim[8]={
1.74,
1.114,
0.02519,
131.607,
-333.194,
319.371,
-136.741,
22.0011
};

Double_t param_IMnpim_corr[8]={
1.092,
1.114,
0.06987,
5.30638,
-13.0001,
16.1339,
-9.15064,
1.86824
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
2.12266,
1.11762,
0.0272097,
5356.81,
-19303.5,
27699.9,
-19776.9,
7024.27,
-992.862
};


Double_t param_IMnpim_wK0_corr[9]={
1.21442,
1.0981,
0.0546432,
7269.69,
-27461.6,
41352.5,
-31016.6,
11585.2,
-1723.58
};


Double_t func_IMnpip_wK0_corr(Double_t *x,Double_t *par)
{
   if(1<x[0] && x[0]<1.14){
     return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
   }else if(1.14<=x[0] && x[0]<=1.30){
    //pol4
     return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)
           +par[8]*pow(x[0],4.0);
   }else if(1.30<x[0] && x[0]<1.71){
     return par[9]+par[10]*x[0]+par[11]*pow(x[0],2.0)+par[12]*pow(x[0],3.0);
   }else{
     return 1;
   }
}


Double_t param_IMnpip_wK0_corr2[13]={
  14578.1,
  -39072.2,
  34907.8,
  -10395.1,
  21331.6,
  -69746.0,
  85473.5,
  -46528.6,
  9492.76,
  51.1006,
  -102.339,
  69.2752,
  -15.5523
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

Double_t param_nmom_v317[9]={
3.7885,
-60.5717,
534.321,
-2519.07,
6965.54,
-11616,
11477,
-6180.95,
1397.09
};


Double_t param_nmom_wK0[9]={
48.5158,
-756.394,
5307.45,
-20791,
49109.1,
-71502.8,
62779.5,
-30473.6,
6279.58
};


Double_t param_nmom_wK0_v313[9]={
6.48565,
-102.943,
841.758,
-3712.86,
9596.8,
-14996.7,
13935.9,
-7084.48,
1516.45
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
  }else if(0.28 <= x[0] && x[0] <= 1.0){
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

Double_t func_IMpippim_corr(Double_t *x,Double_t *par)
{
  if(x[0]<=0.29){
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0));
  }else if(0.29<x[0] && x[0]<0.32){
    return par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0);
  }else if(0.32<=x[0] && x[0]<0.70){
   //pol5
    return par[6]+par[7]*x[0]+par[8]*pow(x[0],2.0)+par[9]*pow(x[0],3.0)+par[10]*pow(x[0],4.0)+par[11]*pow(x[0],5.0);
  }else{
    return 1.;
  }
};

Double_t param_IMpippim_corr[12]={
1.2413,
0.292752,
0.00867818,
-3.55155,
37.8065,
-74.3704,
-54.0954,
578.806,
-2380.41,
4783.04,
-4693.25,
1800.26
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
0.823398,
-0.394939,
-1.19818,
19.1463,
-55.1249,
82.4005,
-59.0306,
15.753
};

Double_t param_q_corr4[8]={
1.05505,
1.20589,
-10.8729,
41.848,
-81.0985,
80.9196,
-40.2394,
7.96846
};

Double_t param_q_corr5[8]={
1.02577,
0.553953,
-5.17268,
21.2187,
-42.8392,
44.8756,
-23.7955,
5.05474
};

Double_t param_q_corr6[8]={
};

Double_t param_q_corr7[8]={
};




Double_t param_q_wK0[8]={
  2.60581,
  -0.608159,
  -42.6983,
  150.097,
  -245.271,
  220.252,
  -105.109,
  20.8622
};

Double_t param_q_wK0_corr[8]={
  0.654005,
  7.41498,
  -66.5307,
  280.971,
  -632.764,
  780.543,
  -489.537,
  122.115
};

Double_t param_q_wK0_corr2[8]={
0.508119,
-0.16851,
9.77129,
-50.5494,
140.106,
-181.39,
108.972,
-24.6617
};

Double_t param_q_wK0_corr3[8]={
0.825218,
-0.0970836,
3.89259,
-20.8223,
56.3774,
-73.6802,
45.6857,
-10.7515
};

Double_t param_q_wK0_corr4[8]={
0.796584,
0.287783,
3.61358,
-12.8859,
14.6609,
-0.811114,
-7.85312,
3.38527
};


Double_t param_q_wK0_corr5[8]={
1.0909,
1.83473,
-13.6843,
50.2759,
-100.122,
104.935,
-54.7527,
11.2442
};



Double_t param_q_wK0_corr6[8]={
};
