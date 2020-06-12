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
2.05093,
-1.48646,
-0.236963,
0.541371,
0.464892,
0.0802679,
-0.230139,
-0.313202,
-0.140528,
0.268801,
0.935134,
0.290191,
-0.233375,
0.84013,
0.222917,
-0.0346273,
-0.0691156,
0.0113051
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


Double_t param_MMnmiss_wK0_corr3[18]={
  1.0105,
  -0.27061,
  2.4707,
  -5.53168,
  1.18957,
  5.16548,
  1.86557,
  -4.1672,
  -5.50124,
  4.73453,
  1.60959,
  -0.788262,
  0.140453,
  1.05666,
  0.21581,
  -0.185559,
  -0.217815,
  0.0904736
};

Double_t param_MMnmiss_wK0_corr4[18]={
1.1373,
-0.188608,
-0.085393,
0.0286818,
0.00295522,
-0.0227861,
0.00822664,
0.0428339,
0.0292489,
-0.0375888,
0.822598,
0.261248,
-0.170677,
0.889593,
0.171011,
-0.109043,
-0.105096,
0.0665909
};

Double_t param_MMnmiss_wK0_corr5[18]={
0.940773,
0.0324743,
0.00701859,
-0.0202257,
-0.062378,
-0.0484597,
0.0196156,
0.0664951,
0.0398941,
-0.0581025,
0.830253,
0.262358,
-0.177713,
0.773121,
0.233155,
-0.0155865,
-0.071451,
-0.00862406
};

Double_t param_MMnmiss_wK0_corr6[18]={
0.950869,
-0.032816,
0.0173145,
0.0254298,
-0.0411465,
-0.0547954,
0.00805618,
0.0636743,
0.0445691,
-0.0556762,
0.786981,
0.2664,
-0.128765,
0.825653,
0.199314,
-0.0604462,
-0.0811634,
0.0412611
};

Double_t param_MMnmiss_wK0_corr7[18]={
0.939284,
-0.00957458,
0.00774752,
0.0146914,
-0.0424951,
-0.048764,
0.0158822,
0.0681307,
0.0428932,
-0.0631986,
0.788746,
0.265698,
-0.131069,
0.811152,
0.208241,
-0.0481434,
-0.077585,
0.0298788
};

Double_t param_MMnmiss_wK0_corr8[18]={
1.07506,
-0.18663,
0.554224,
-0.370442,
-0.908513,
-0.391632,
0.553045,
1.03094,
0.534216,
-0.985947,
0.747851,
0.230939,
-0.073633,
0.722746,
0.228487,
-0.0010042,
-0.0531887,
0.00724505
};


Double_t param_MMnmiss_wK0_corr9[18]={
2.10406,
-4.63889,
14.2547,
-18.0977,
-7.52286,
16.7752,
15.8671,
-8.61103,
-23.1818,
14.0537,
2.03065,
-1.70971,
0.665223,
0.935659,
0.279347,
-0.122071,
-0.213727,
0.0984271
};


Double_t func_MMnmiss_wK0(Double_t *x,Double_t *par)
{
   if(0.8<x[0] && x[0]<1.06){
     return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0);
   }else if(1.06<=x[0] && x[0]<1.115 ){
     return par[10]+par[11]*x[0]+par[12]*pow(x[0],2.0);
   }else if(1.115<x[0] && x[0]<1.5){
     return par[13]+par[14]*x[0]+par[15]*pow(x[0],2.0)+par[16]*pow(x[0],3.0)+par[17]*pow(x[0],4.0);
   }else{
     return 1.0;
   }
}


Double_t func_MMnmiss_mod(Double_t *x,Double_t *par)
{
   //connection using woods-saxon
   if(0.0<x[0] && x[0]<1.5){
  //   return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
  //     +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
  //     +par[7]*pow(x[0],7.0))*(1./(1.0+exp((x[0]-1.10)/par[8])))
     return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
       +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)
       +par[7]*pow(x[0],7.0)+par[8]*pow(x[0],8.0)+par[9]*pow(x[0],9.0))*(1./(1.0+exp((x[0]-1.08)/par[10])))
       +(1.0-1./(1.0+exp((x[0]-1.08)/par[10])))*(par[11]*exp(-0.5*pow(((x[0]-par[12])/par[13]),2.0)))*(1./(1.0+exp((x[0]-1.14)/par[14])))
       +(1.0-1./(1.0+exp((x[0]-1.14)/par[14])))*
        (par[15]+par[16]*x[0]+par[17]*pow(x[0],2.0)+par[18]*pow(x[0],3.0)+par[19]*pow(x[0],4.0));
   }else{
     return 1.0;
   }
}

/*
Double_t param_MMnmiss_mod[18]={
-3.43372,
58.0474,
-386.213,
1363.79,
-2711.64,
3046.57,
-1793.34,
428.104,
0.015,
42.9194,
0.40951,
0.303197,
0.01,
722.62,
-2158.24,
2433.18,
-1222.05,
230.02
};
*/
/*
Double_t param_MMnmiss_mod[18]={
-3.40395,
57.6393,
-385.855,
1363.85,
-2712.14,
3046.57,
-1792.4,
427.548,
0.015,
314.216,
-0.108428,
0.398183,
0.01,
723.056,
-2159.13,
2433.32,
-1221.57,
229.823
};
*/

/*
Double_t param_MMnmiss_mod[18]={
-3.90353,
63.7699,
-414.916,
1422.4,
-2746.62,
3001.35,
-1721.86,
401.488,
0.015,
1063.05,
-0.333118,
0.419609,
0.01,
739.969,
-2185.25,
2435.47,
-1209.43,
225.167
};
*/

/*
Double_t param_MMnmiss_mod[20]={
0.00193703,
0.0108121,
5.89338,
-38.6911,
94.4113,
46.6664,
-609.268,
1098.4,
-824.566,
228.764,
2.3678,
1.11385,
0.0665073,
771.754,
-2322.27,
2632.98,
-1328.19,
250.877
};
*/
/*
Double_t param_MMnmiss_mod[20]={
0.00130506,
0.101503,
3.63065,
-18.7114,
15.189,
191.464,
-689.597,
995.684,
-664.774,
168.681,
0.01,
2.57533,
9.33939,
24.3035,
0.01,
260.528,
-802.777,
942.813,
-494.31,
96.904
};
*/
/*
Double_t param_MMnmiss_mod[20]={
0.00131172,
0.276019,
3.7462,
-18.7678,
15.4473,
191.425,
-690.015,
995.172,
-664.941,
169.308,
0.01,
3.52661,
-15.7252,
18.5641,
0.01,
260.33,
-802.781,
942.869,
-494.271,
96.8947
};
*/

Double_t param_MMnmiss_mod[20]={
0.000281332,
0.490481,
2.6641,
-15.438,
13.4367,
188.441,
-689.654,
998.37,
-662.618,
166.036,
0.01,
3.59474,
0.712467,
0.432273,
0.01,
260.568,
-803.045,
942.772,
-494.219,
96.9351
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
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0);
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

Double_t param_IMnpip_corr2[9]={
-34.5925,
91.374,
-77.7511,
21.9431,
-1.6005,
4.90642,
-2.29113,
-0.301563,
0.257101
};




Double_t func_IMnpipmul_s(Double_t *x,Double_t *par)
{
  if(1.06<=x[0] && x[0]<1.92){
    return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))*(1./(1.0+exp((x[0]-1.25)/par[11])))+
           (par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))*(1.0-1./(1.0+exp((x[0]-1.25)/par[11])));
  }else{
    return 1.;
  }
}

Double_t param_IMnpip_s[12]={
-327.854,
748.394,
-561.008,
138.25,
-3277.16,
12041.2,
-18359.5,
14891.1,
-6782.67,
1645.87,
-166.272,
0.01
};



/*
Double_t func_IMnpip_corr(Double_t *x,Double_t *par)
{
  if(1.06<=x[0] && x[0]<=2.00){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0);
  }else{
    return 1.0;
  } 
}


Double_t param_IMnpip_corr3[5]={
3.23373,
-2.89104,
-0.474793,
1.93943,
-0.685505
};
*/

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

Double_t func_IMnpim_corr(Double_t *x,Double_t *par)
{  
  if(1.00 <= x[0] && x[0]<1.11) {
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
    //return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
    return 1.0;
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



Double_t param_IMnpim_corr2[8]={
1.092,
1.114,
0.06987,
2.54129,
-3.00676,
2.34562,
-0.944198,
0.151111
};



Double_t param_IMnpim_corr3[8]={
1.092,
1.114,
0.06987,
49.5531,
-142.585,
157.349,
-77.2074,
14.1624
};


Double_t func_IMnpim_mod(Double_t *x,Double_t *par)
{
  if(x[0]<1.07){
    return 0.0;
  }else if(1.07 <= x[0]){
    return (par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.0)))
          *(1./(1.0+exp((x[0]-1.10)/par[3])))
          +(1.0 - 1./(1.0+exp((x[0]-1.10)/par[3])))
          //*(exp(par[4]+par[5]*x[0]));
          *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)
           +par[9]*pow(x[0],5.0));//+par[10]*pow(x[0],6.0)+par[11]*pow(x[0],7.0));
  }else{
   return 0;
  }
}


Double_t param_IMnpim_mod[10]={
1.93919,
1.10969,
0.0209215,
0.02,
111.098,
-210.672,
89.7735,
55.8835,
-53.8587,
11.4735
};




Double_t func_IMnpim_wK0(Double_t *x,Double_t *par)
{
  if(1.00 <= x[0] && x[0]<1.11) {
    //return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
    return par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)); 
  } else if(1.11 <= x[0] && x[0]<=1.7) {
    return par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0)+par[6]*pow(x[0],3.0)+par[7]*pow(x[0],4.0)+par[8]*pow(x[0],5.0);
  } else {
    return 0.;
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

Double_t param_IMnpim_wK0_corr2[9]={
1.06217,
1.10574,
0.084444,
10.5054,
-98.0943,
256.286,
-283.087,
142.609,
-27.0793
};

Double_t param_IMnpim_wK0_corr3[9]={
1.08479,
1.08658,
0.0844952,
881.291,
-3380.98,
5184.96,
-3966.96,
1513.51,
-230.271
};

Double_t param_IMnpim_wK0_corr4[9]={
1.03365,
1.0955,
0.0846515,
1513.43,
-5817.75,
8923.28,
-6820.81,
2597.97,
-394.424
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
     return 0;
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


Double_t func_IMnpip_wK0_corr2(Double_t *x,Double_t *par)
{
  if(1.0<=x[0] && x[0]<1.15){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0);
  }else if(1.15<=x[0] && x[0]<1.25){
    return par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0);
  }else if(1.25<=x[0] && x[0]<1.71){
    return par[8]+par[9]*x[0]+par[10]*pow(x[0],2.0)+par[11]*pow(x[0],3.0);
  }else{
    return 1.0;
  }
}


Double_t param_IMnpip_wK0_corr3[12]={
  -8143.0,
  21923.1,
  -19668.4,
  5880.91,
  3652.29,
  -9101.91,
  7561.11,
  -2093.17,
  38.4039,
  -76.8394,
  52.4779,
  -11.9415
};


Double_t param_IMnpip_wK0_corr4[12]={
  -1322.95,
  3543.57,
  -3160.53,
  939.427,
  93.8916,
  -228.603,
  188.06,
  -51.6619,
  9.33782,
  -16.5109,
  10.9703,
  -2.42934
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

Double_t param_nmom_v341[9]={
0.891033,
1.77151,
-5.04277,
-34.8959,
286.365,
-822.087,
1149.45,
-789.999,
213.998
};

Double_t param_nmom_v345[9]={
3.8321,
-64.0065,
589.792,
-2913.44,
8461.09,
-14844.7,
15435.0,
-8734.75,
2068.74
};


Double_t param_nmom_v353[9]={
30.2372,
-313.925,
1412.82,
-3246.26,
3356.36,
456.405,
-4583.45,
4091.76,
-1203.77
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


Double_t param_nmom_wK0_v326[9]={
3.90496,
-62.1428,
548.478,
-2595,
7202.31,
-12045.8,
11911.7,
-6403.85,
1441.19
};


Double_t param_nmom_wK0_v333[9]={
1.57273,
-11.0729,
101.821,
-498.737,
1431.25,
-2475.1,
2528.24,
-1401.52,
324.551
};


Double_t param_nmom_wK0_v339[9]={
2.87928,
-42.1032,
386.902,
-1909.21,
5549.52,
-9765.02,
10202.5,
-5804.67,
1381.19
};


Double_t param_nmom_wK0_v341[9]={
-2.36278,
77.1795,
-699.762,
3376.34,
-9495.56,
16009.6,
-15949.4,
8660.18,
-1975.94
};

Double_t param_nmom_wK0_v345[9]={
4.29847,
-73.5664,
667.905,
-3252.01,
9319.22,
-16160.5,
16644.7,
-9354.26,
2206.7
};

Double_t param_nmom_wK0_v349[9]={
2.19649,
-14.6197,
97.3219,
-324.903,
518.34,
-188.501,
-515.487,
664.399,
-237.424
};


Double_t func_nmom_mod(Double_t *x,Double_t *par)
{
  if(x[0]<0.140){
    return 1.0;
  }else if(0.14<= x[0] && x[0]<1.0){
 //   return exp(par[0]+par[1]*x[0])*((1./(1.0+exp((x[0]-0.2)/par[2]))))
 //   + (par[3]+par[4]*x[0]+par[5]*pow(x[0],2.0)+par[6]*pow(x[0],3.0)+par[7]*pow(x[0],4.0)+par[8]*pow(x[0],5.0)
 //   +par[9]*pow(x[0],6.0))*(1.-1./(1.0+exp((x[0]-0.2)/par[2])));
    return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0))*((1./(1.0+exp((x[0]-0.2)/par[4]))))
    + (par[5]+par[6]*x[0]+par[7]*pow(x[0],2.0)+par[8]*pow(x[0],3.0)+par[9]*pow(x[0],4.0)+par[10]*pow(x[0],5.0)
    +par[11]*pow(x[0],6.0))*(1.-1./(1.0+exp((x[0]-0.2)/par[4])));
  }else{
    return 1.0;
  }
}

Double_t param_nmom_mod[12]={
126.564,
-1837.3,
9351.47,
-16281.6,
0.01,
11.039,
-61.0051,
118.523,
-49.6603,
-125.447,
165.959,
-59.1096
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

Double_t param_IMpippim_corr2[12]={
  0.994132,
  0.290115,
  0.0188121,
  -1.07393,
  13.0808,
  -21.1081,
  1.04552,
  3.0547,
  -28.1953,
  82.9608,
  -99.6095,
  42.8865
};


Double_t func_IMpippim_mod(Double_t *x,Double_t *par)
{
  if(0.27<x[0]){
  return (par[0]*exp(-0.5*pow(((x[0]-par[1])/par[2]),2.0)))
  //      *(1./(1.0+exp((x[0]-0.32)/par[3])))            
        *(1./(1.0+exp((x[0]-0.31)/par[3])))            
  //    +(1.0 - 1./(1.0+exp((x[0]-0.32)/par[3])))
      +(1.0 - 1./(1.0+exp((x[0]-0.31)/par[3])))
       *(par[4]+par[5]*x[0]+par[6]*pow(x[0],2.0)+par[7]*pow(x[0],3.0)+par[8]*pow(x[0],4.0)+par[9]*pow(x[0],5.0)+par[10]*pow(x[0],6.0))
       *(1./(1.0+exp((x[0]-0.73)/par[11])))
      +(1.0 - 1./(1.0+exp((x[0]-0.73)/par[11])))
       *(par[12]*exp(-0.5*pow(((x[0]-par[13])/par[14]),2.0)))
       ;
   }else{
    return 0.0;
   }
};

/*
Double_t param_IMpippim_mod[15]={
0.395041,
0.296922,
0.0118245,
0.02,
5.07512,
-23.1285,
10.27,
84.9313,
1.30107,
-231.161,
140.749,
0.03,
1.78472,
0.753667,
0.0668154
};
*/
/*
Double_t param_IMpippim_mod[15]={
0.468509,
0.297271,
0.0130168,
0.02,
-16.0994,
157.248,
-429.852,
-190.526,
2851.48,
-4601.61,
2316.36,
0.03,
1.62175,
0.791874,
-0.0483207
};
*/

Double_t param_IMpippim_mod[15]={
0.431559,
0.297098,
0.0124366,
0.02,
2.1313,
-5.02899,
-6.19126,
17.1507,
50.3553,
26.175,
-154.199,
0.03,
4.39784,
0.671838,
0.088268
};


Double_t func_IMpippim_wK0(Double_t *x,Double_t *par)
{
  return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)+par[4]*pow(x[0],4.0);
};

Double_t param_IMpippim_wK0[5]={
5.81027e+06,
-4.6716e+07,
1.40843e+08,
-1.88708e+08,
9.48076e+07
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
};




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
1.02276,
-0.107561,
1.2506,
-8.54088,
28.0589,
-44.0838,
31.7414,
-8.53317
};

Double_t param_q_corr7[8]={
1.07683,
-1.27815,
11.458,
-48.4124,
106.792,
-127.263,
76.7718,
-18.3157
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
1.02266,
1.35377,
-9.40087,
35.7952,
-73.1855,
78.2075,
-41.4799,
8.62952
};



Double_t param_q_wK0_corr7[8]={
1.07843,
0.211107,
-1.0433,
2.31208,
-2.71374,
1.10978,
0.403476,
-0.297473
};

/*
Double_t param_q_mul[8]={
3.8072,
-5.52768,
9.87903,
-26.5909,
75.0323,
-123.557,
99.1959,
-29.4453
};*/
/*
Double_t param_q_mul[8]={
4.21641,
-15.0937,
102.561,
-448.393,
1090.17,
-1450.83,
983.821,
-262.107
};*/

Double_t func_q_mod(Double_t *x,Double_t *par)
{
  if(0<x[0]) {  //&& x[0]<1.17){
    return (par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0));
    //*(1./(1.0+exp((x[0]-1.17)/par[8]))) 
    //+(1.0 - 1./(1.0+exp((x[0]-1.17)/par[8])))
    //*par[9];
    //*(par[9]+par[10]*x[0]+par[11]*pow(x[0],2.0)+par[12]*pow(x[0],3.0)+par[13]*pow(x[0],4.0)+par[14]*pow(x[0],5.0));
  }
}

/*
Double_t param_q_mod[10]={
4.22238,
-11.1141,
46.9233,
-117.049,
91.6876,
148.49,
-314.236,
159.295,
0.01,
-3.5824
};*/

/*
Double_t param_q_mod[8]={
4.10422,
-4.40469,
-0.437715,
15.9143,
-20.1113,
-11.6589,
36.7471,
-16.5354
};*/

/*
Double_t param_q_mod[8]={
4.25462,
-4.71068,
0.126449,
14.8645,
-20.1902,
-9.06257,
32.105,
-14.488
};
*/

Double_t param_q_mod[8]={
4.37725,
-5.26915,
1.8429,
12.0264,
-19.8566,
-5.04731,
27.5322,
-12.9654
};

Double_t func_q_wK0(Double_t *x,Double_t *par)
{
  if(0<x[0] && x[0]<1.5){
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2.0)+par[3]*pow(x[0],3.0)
    +par[4]*pow(x[0],4.0)+par[5]*pow(x[0],5.0)+par[6]*pow(x[0],6.0)+par[7]*pow(x[0],7.0) 
    +par[8]*pow(x[0],8.0);
  }else{
    return 1.0;
  }
}

Double_t param_q_wK0_corr8[9]={
1.19774,
-3.94891,
40.229,
-190.51,
488.979,
-719.832,
606.061,
-270.856,
49.7548
};

Double_t param_q_wK0_corr9[9]={
0.979371,
1.2872,
-11.4152,
39.8568,
-49.3429,
-19.2135,
96.0009,
-77.2793,
19.9878
};

Double_t param_q_wK0_corr10[9]={
0.994817,
1.46609,
-13.9497,
64.9922,
-169.423,
254.779,
-220.419,
101.976,
-19.5526
};
